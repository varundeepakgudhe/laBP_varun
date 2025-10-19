/**
 * World.cpp
 *  InvCoal
 *  Version:  labp_v17
 *  Edited by Rafael Guerrero 2008-2016
 *
 * This is the main class in our simulation
 * holds the vector of populations, the vectors of carriers per population, all parameter values
 * Methods for this class carry out the main events of the simulation (migration, recombination, coalescence)
 *
 *
 */

//#define DEBUGGER
// Includes from STL:
	#include <iostream>
		using std::cout;
		using std::endl;
	#include <vector>
		using std::vector;
	#include <map>
		using std::map;
	#include <math.h>
	
// Includes from our files:
	#include "chromosome.h"
	#include "ran_mk.h"
	#include "world.h"
	#include "argnode.h"
	#include <algorithm>    // for sort()


void World::rebuildClustersAndCarriers() {
    // Gather all chromosomes
    std::vector< std::shared_ptr<Chromosome> > allChr;
    allChr.reserve(totalNCarriers());
    for (auto &vec : *worldData->carriers)
        for (auto &chr : vec)
            allChr.push_back(chr);

    // Rebuild cluster map
    cluster.clear();
    makecluster(worldData->popSize.size());
    worldData->nClust = cluster.size();

    // Reset carriers buckets
    worldData->carriers->assign(worldData->nClust, {});

    // Reassign carriers to new clusters
    for (auto &chr : allChr) {
        auto it = cluster.find(chr->getContext());
        if (it == cluster.end()) {
            std::cerr << "[rebuildClustersAndCarriers] ERROR: unknown context after relabel: "
                      << "pop=" << chr->getContext().pop
                      << " inv=" << chr->getContext().inversion << "\n";
            continue; // or assert(false);
        }
        int cid = it->second;
        worldData->carriers->at(cid).push_back(chr);
    }
}



World::World(shared_ptr<Parameters::ParameterData> p){	
	// Constructor for World
	//

	DBG("Population:: Constructing Population(size, shared_ptr<Chromosome>)...")

	unsigned int nPops = p->popSizeVec.size();
	
	makecluster(nPops); // builds a map of all contexts
	
	unsigned int nClust=cluster.size();
	

	// Wipe clean WorldData:
	//
	worldData.reset(new WorldData);
	worldData->carriers.reset( new vector < vector< shared_ptr<Chromosome> > > );
	worldData->generation = 0;
	worldData->nPops = nPops;
	worldData->nClust = nClust;
    worldData->kingman=p->kingman;
	worldData->popSize=p->popSizeVec;
	worldData->drift=p->drift;
    worldData->snpSites = p->neut_site;
	worldData->invRange = p->invRange;
	worldData->phi=p->phi;
	worldData->freq=p->initialFreqs;

    worldData->freq.resize(nClust);
    for( cluster_t::iterator iter = cluster.begin(); iter != cluster.end(); ++iter )
    {    worldData->freq[(*iter).second]= 1 - p->initialFreqs[(*iter).first.pop];
        if((*iter).first.inversion==1)
            worldData->freq[(*iter).second]= p->initialFreqs[(*iter).first.pop];
    }
    
    worldData->current_epoch = 0;

    struct Event { double time; int type; double value; vector<unsigned int> popSizes; };
    vector<Event> events;

    // 1) Demography events → type=2
    if (p->demography.at(0) == 1) {
        size_t np     = p->popSizeVec.size();
        size_t stride = 1 + np;
        size_t count  = (p->demography.size() - 1) / stride;
        for (size_t e = 0; e < count; ++e) {
            size_t idx = 1 + e * stride;
            double t = p->demography[idx];
            vector<unsigned int> sizes;
            for (size_t k = 0; k < np; ++k) {
                sizes.push_back(static_cast<unsigned int>(p->demography[idx + 1 + k]));
            }
            events.push_back({ t, 2, 0.0, sizes });
        }
    }

    // --- speciation parsing (new 1 A B T F R format; backward compatible) ---
	if (!p->speciation.empty() && p->speciation.at(0) == 1) {
		// Prefer new stride = 5 (A,B,T,F,R). If only 4 fields remain, treat like legacy (A,B,T,F) with no R.
		size_t i = 1;
		while (i < p->speciation.size()) {
			if (i + 4 < p->speciation.size()) {
				// New format with R
				unsigned int A = (unsigned int)p->speciation[i];
				unsigned int B = (unsigned int)p->speciation[i+1];
				double       T =               p->speciation[i+2];
				double       F =               p->speciation[i+3];
				unsigned int R = (unsigned int)p->speciation[i+4];
				events.push_back({ T, 1, F, { A, B, R } });
				i += 5;
			} else if (i + 3 < p->speciation.size()) {
				// Legacy format (no R)
				unsigned int A = (unsigned int)p->speciation[i];
				unsigned int B = (unsigned int)p->speciation[i+1];
				double       T =               p->speciation[i+2];
				double       F =               p->speciation[i+3];
				events.push_back({ T, 1, F, { A, B } });
				i += 4;
			} else {
				break; // trailing junk; ignore
			}
		}
	}




    // sort all events by their time
    sort(events.begin(), events.end(),
         [](auto &a, auto &b){
             if (a.time != b.time) return a.time < b.time;
             return a.type  > b.type;    // demography (2) comes before speciation (1)
         });

    // populate worldData vectors
    for (auto &ev : events) {
        worldData->epoch_breaks.push_back(ev.time);
        worldData->epochType.push_back(ev.type);
        worldData->epochValues.push_back(ev.value);
        worldData->epochPopSizes.push_back(ev.popSizes);
    }

    worldData->epochs_over = worldData->epoch_breaks.empty();

    // ---------------- Inversion-origin initialization ----------------
    inv_age = p->inv_age;
    originPop = p->originPop;
    // NOTE: origin_start_time is *re-inferred* here. If your pipeline guarantees
    // originPop is a present-day leaf (deme end_time == 0), then origin_start_time = 0.
    // If you later permit non-zero origin_start_time, pass it explicitly from Python.
    origin_start_time       = 0.0;
    originPopInitialized    = false;
    inversionOriginHandled  = false;

    Context originCtx(0,0);
	worldData->originCtx=originCtx;
	// Copy the population sizes and resize vectors to the number of contexts:
	//
	worldData->carriers->resize(nClust);
	
	DBG("n of contexts "<<cluster.size())
	for( cluster_t::iterator iter = cluster.begin(); iter != cluster.end(); ++iter ) 
	{
		DBG((*iter).first.inversion<<" inversion, "<<(*iter).first.pop<<" pop is: "<< (*iter).second<<" freq = "<<worldData->freq.at((*iter).second));
	}
	
	
	vector<int> initNCarriers;
	initNCarriers.resize(nClust);
	shared_ptr< ARGNode > nulle ;
	
    for(int i = 0; i < p->initChr.size(); ++i){
        assert(!p->initChr.at(i)->isEmpty());
        shared_ptr<Chromosome> chr= p->initChr.at(i)->copy_values(worldData->invRange);
        shared_ptr < ARGNode > newNode ;newNode.reset( new ARGNode(worldData->argNodeVec.size(), chr) );// Make a terminal ARG node
        chr->setDescendant(newNode);
        worldData->carriers->at(cluster[chr->getContext()]).push_back( chr );					// Adds this chromosome to the carriers array
        worldData->argNodeVec.push_back(newNode);							// Add the new node to the vector of ARG nodes
        
        initialARGnodes.push_back(newNode); //added X-13, saved for alternative output
    }
    
    
    worldData->generation=1.0;

	// debug: what epochs do we have?
	std::cerr << ">>> Epoch breaks: ";
	for (double t : worldData->epoch_breaks) std::cerr << t << " ";
	std::cerr << "\n>>> Epoch types:  ";
	for (unsigned ty : worldData->epochType)   std::cerr << ty << " ";
	std::cerr << "\n";

} // End constructor for World



World::~World(){
	DBG("Population:: Destructing...")
}



bool World::simulationFinished(){
    //
    // Determine if only one carrier remains in the world
    //
    if(totalNCarriers()==1) return true;

//C.Cheng, change end condition
//	if( worldData->LeftTotalSegSize <= 0 || totalNCarriers()==1) return true;
//	if (worldData->recomb_breakpoints.size() == 0) return true;

    return false;
}

bool World::sitesCoalesced(){
    //
    // Determine if all sites of interest have coalesced
    //
    vector<int> carrierCounter(worldData->snpSites.size(),0);

	for( cluster_t::iterator iter = cluster.begin(); iter != cluster.end(); ++iter ) {
        for(int j=0; j< worldData->carriers->at((*iter).second).size(); ++j){
            for(int s=0; s< worldData->snpSites.size(); ++s){
                double site= worldData->snpSites[s];
                if(worldData->carriers->at((*iter).second)[j]->carriesSite(site)) ++carrierCounter[s];
            }
        }
    }
    bool out = true;
    for (int i=0; i< carrierCounter.size();++i){
        if(carrierCounter[i]>1) out= false;
    }
    return out;
}

void World::forceAllCoal(){
 
	vector<shared_ptr<Chromosome> > sons;
	for( cluster_t::iterator iter = cluster.begin(); iter != cluster.end(); ++iter ) {
			unsigned long n=ctxtNcarriers((*iter).first);
			
			for(int j=0; j<n;++j){
				shared_ptr<Chromosome> chr=worldData->carriers->at((*iter).second).at(j);
                chr->setStdSegs(chr->get_All_SegVec());
				chr->setContext(worldData->originCtx);
				sons.push_back(chr);
            }
			worldData->carriers->at((*iter).second).clear();
	}
    
	if(sons.size()>0){
		shared_ptr<Chromosome> cain=sons.at(0);
		
		shared_ptr < ARGNode > cainNode ;
		cainNode.reset(  new ARGNode(worldData->argNodeVec.size(), cain, worldData->generation) );
		
		for(int i=1; i<sons.size(); ++i){
			shared_ptr<Chromosome> abel=sons.at(i);
			shared_ptr < ARGNode > newNode ;
			newNode.reset(  new ARGNode(worldData->argNodeVec.size(), cain, abel, worldData->generation) );
			worldData->argNodeVec.push_back(newNode);		// Add a pointer to this node to the vector of all nodes in	the simulation
			cain->merge( abel, newNode); // merging includes setting decendant to new ARGNode.
		}
		
		worldData->carriers->at(cluster[worldData->originCtx]).push_back(cain);
	}

}


unsigned long World::ctxtNcarriers(const Context ctxt){
	unsigned long total = worldData->carriers->at(cluster[ctxt]).size();
	return total;
}


unsigned long World::popNcarriers(unsigned long pop){
// Returns the number of carriers in all contexts of population pop	
	unsigned long total = 0;
	for( cluster_t::iterator iter = cluster.begin(); iter != cluster.end(); ++iter ) {
		if (pop==(*iter).first.pop)
			total = worldData->carriers->at((*iter).second).size() +total;
	}
	return total;
	
}

unsigned long World::totalNCarriers(){
//
// Returns the total number of carriers in the world
//
	unsigned long total = 0;
	for(unsigned long i = 0; i < worldData->nClust; i++)
	{	total = worldData->carriers->at(i).size() +total;
		}
	return total;
}



double World::nGenerations(){
	return worldData->generation;
}

void World::Generation_pp (){
	worldData->generation++;
}

void World::makecluster (unsigned long nPops){
	
	int iter=0;
	for (int i=0; i<nPops; ++i){			//loop through populations
		for (int j=0; j < 2; ++j){			//loop through inversion states (I,S)	//loop through alleles in siteA

				Context c(i, j); 
				cluster[c]=iter;
				++iter;
			
			}
		}

}

vector< shared_ptr < ARGNode > >& World::getARGVec(){
	return worldData->argNodeVec;
}

double World::getFreqI(unsigned long pop){
	double total = 0;
	for( cluster_t::iterator iter = cluster.begin(); iter != cluster.end(); ++iter ) {
		if (1 ==(*iter).first.inversion && pop==(*iter).first.pop)
			total+= worldData->freq.at((*iter).second);
	}
	return total;
	
}

vector<double> World::getFreqbyPop(Context c){
	vector<double> freqs (worldData->nPops,0);
	Context geno=c;
	for (int i=0; i<worldData->nPops;++i){
		geno.pop=i;
		freqs.at(i)=worldData->freq.at(cluster[geno]);
	}
	return freqs;
}

vector<double> World::getSizebyPop(Context c){
    vector<double> sizes (worldData->nPops,0);
    Context geno=c;
    for (int i=0; i<worldData->nPops;++i){
        geno.pop=i;
        sizes.at(i)=worldData->freq.at(cluster[geno]) * worldData->popSize.at(geno.pop);
    }
    return sizes;
}


Context World::getRecombContext(Context c, bool hetero){
	// this function gives the context of a chromosome that will recombine with a carrier
	// given that the recombination is homokaryotypic or heterokaryotipic
	//  the function randomly selects one of the possible genotypes within standard or inverted contexts.
	// 
	
	//if(hetero)cout<<hetero<<'\n';
	double f= getFreqI(c.pop);
	int inv=1;
	if ((c.inversion==1 && hetero) || (c.inversion==0 && !hetero)){ 
		f=1-f;
		inv=0;
	}
	double draw = randreal(0, f);
	double total = 0;
	for( cluster_t::iterator iter = cluster.begin(); iter != cluster.end(); ++iter ) {
		if (inv ==(*iter).first.inversion && c.pop==(*iter).first.pop){
			total+= worldData->freq.at((*iter).second);
			if(draw<total) {
					return (*iter).first;
			}
		}
	}
	cout<<"Error in getRecombContext. hetero= "<<hetero<<", f= "<<f<<", inv= "<<inv<<", total= "<<total<<"\n";
	return c;
}

void World::updateContexts(Context originCtx){
	// this function will transfer all carriers on inverted contexts to standard context. It's meant to be done 
	// at the origin of the inversion
	
	vector<shared_ptr<Chromosome> > sons;
	for( cluster_t::iterator iter = cluster.begin(); iter != cluster.end(); ++iter ) {
		if (1 ==(*iter).first.inversion && ctxtNcarriers((*iter).first)>0 ){
			unsigned long n=ctxtNcarriers((*iter).first);
			
			for(int j=0; j<n;++j){
				shared_ptr<Chromosome> chr=worldData->carriers->at((*iter).second).at(j);
				chr->setStdSegs(chr->get_All_SegVec());
				chr->clearInvSegs();
				chr->setContext(originCtx);
				sons.push_back(chr);
			}
			worldData->carriers->at((*iter).second).clear();
		}
		
	}
	if(sons.size()>0){
		shared_ptr<Chromosome> cain=sons.at(0);
		
		//UPDATE NOV/10 RG included these two lines to make sure the node itself changes context when switched to Standard (even without coalescence)
		shared_ptr < ARGNode > cainNode ;
		cainNode.reset(  new ARGNode(worldData->argNodeVec.size(), cain, worldData->generation) );
		
		for(int i=1; i<sons.size(); ++i){
			shared_ptr<Chromosome> abel=sons.at(i);
			shared_ptr < ARGNode > newNode ;
			newNode.reset(  new ARGNode(worldData->argNodeVec.size(), cain, abel, worldData->generation) );
			worldData->argNodeVec.push_back(newNode);		// Add a pointer to this node to the vector of all nodes in	the simulation	
			cain->merge( abel, newNode); // merging includes setting decendant to new ARGNode.
		}
		
		worldData->carriers->at(cluster[originCtx]).push_back(cain);
	}
}

void World::updateFreqs(vector<double> newF){
	worldData->freq=newF;
}
	
Context World::randCtx(unsigned long pop){
	double draw = randreal(0, 1);
	Context c(0,0);
	double total = 0;
	for( cluster_t::iterator iter = cluster.begin(); iter != cluster.end(); ++iter ) {
		if (pop==(*iter).first.pop){
			total+= worldData->freq.at((*iter).second);
			if(draw<total) {
				return (*iter).first;
			}
		}
	}
	cout<<"Error in randCtx\n";
	return c;
	
}

//----------------------------------------------------------
// Inversion-origin handling
//   t = 0: present
//   t = origin_start_time: move all inverted carriers into originPop
//   origin_start_time < t < inv_age: inversion exists only in originPop
//   t = inv_age: collapse all inverted lineages in originPop and convert to standard
//   t > inv_age: inverted chromosomes forbidden
//----------------------------------------------------------
void World::handleInversionOrigin() {
    double t = worldData->generation;

	if (inv_age < 0) return;
	
    // 1) Before origin_start_time → nothing special
    if (t < origin_start_time)
        return;

	// 2) First time reaching/clearing origin_start_time → move ALL inversions to originPop
	if (!originPopInitialized && t >= origin_start_time) {
		std::cerr << "[InversionOrigin] t=" << t
				<< " >= origin_start_time (" << origin_start_time
				<< "): moving all inverted carriers into originPop=" << originPop << "\n";

		for (auto &kv : cluster) {
            const Context &ctx = kv.first;
            if (ctx.inversion == 1 && ctx.pop != originPop) {
                auto &vec = worldData->carriers->at(kv.second);
                for (auto &chr : vec) chr->setContext(Context(originPop, 1));
            }
        }
		// Re-bucket carriers according to the new contexts
		rebuildClustersAndCarriers();

		// ensures collapse converts to STANDARD in originPop (even if another event also happens at this t)
		worldData->originCtx = Context(originPop, 0);

		originPopInitialized = true;
		return;
	}


    // 3) Between origin_start_time and inv_age → enforce that inversions remain in originPop
    if (originPopInitialized && !inversionOriginHandled && t > origin_start_time && t < inv_age) {
        bool moved_any = false;
        for (auto &kv : cluster) {
            const Context &ctx = kv.first;
            if (ctx.inversion == 1 && ctx.pop != originPop) {
                auto &vec = worldData->carriers->at(kv.second);
                if (!vec.empty()) {
                    std::cerr << "[Warning] Inverted carriers outside originPop at t=" << t
                              << " → moving to originPop=" << originPop << "\n";
                    for (auto &chr : vec) {
                        chr->setContext(Context(originPop, 1));
                    }
                    moved_any = true;
                }
            }
        }
        if (moved_any) rebuildClustersAndCarriers();
        return;
    }

    // 4) At or beyond inv_age (first time) → collapse all inverted lineages and convert to standard
    if (!inversionOriginHandled && t >= inv_age) {
		std::cerr << "[InversionOrigin] t=" << t
				<< " >= inv_age (" << inv_age
				<< "): collapsing and converting to standard in originPop=" << originPop << "\n";
		// make sure collapse target is the standard context of originPop
		worldData->originCtx = Context(originPop, 0);
		// Note: zero out all inverted (and non-origin) frequencies, so the simulator cannot re-create inversions after inv_age.
		freqStepToLoss();  // this calls updateContexts(originCtx) internally
		inversionOriginHandled = true;
		return;
	}

    // 5) After inv_age → no inversions allowed
    if (inversionOriginHandled && t > inv_age) {
        for (auto &kv : cluster) {
            const Context &ctx = kv.first;
            if (ctx.inversion == 1 && !worldData->carriers->at(kv.second).empty()) {
                std::cerr << "[Error] Inverted chromosomes detected at t=" << t
                          << " (> inv_age=" << inv_age << "). Aborting.\n";
                exit(1);
            }
        }
    }
}

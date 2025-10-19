/** simulate.cpp
 *  InvCoal
 *  Version:  labp_v17
 *  transfered from World.cpp by RG 03.09
 *  This file has the workhorse function of the simulation, simulateGeneration(), that runs in a loop and determines which events happen, controling the demographic epochs. 
 *
 */

//#define DEBUGGER
//#define NDEBUG
//#define CURRENTDEBUGGER
//#define COALDBGR
//#define RCDBGR

// Includes from STL:


#include <iostream>
using std::cout;
using std::endl;
#include <vector>
using std::vector;
#include <map>
using std::map;
#include <math.h>
using std::pair;
using std::make_pair;
#include <algorithm>
using std::max;
using std::min;


#include <boost/math/special_functions/binomial.hpp>
using boost::math::binomial_coefficient;

// Includes from our files:
#include "chromosome.h"
#include "ran_mk.h"
#include "world.h"
#include "argnode.h"
#include "sitenode.h"


//	Implimentation notes:

//		Homokaryotypic recombination is done with interference: there's a maximum of one recombination event per chromosome.
//		Autosomal inheritance assumed
//		Recombination assumes the partner chromosome is *not* a carrier.

unsigned short World::simulateGeneration(vector < vector <double> > & mig_prob){
    //
    //
    DBG(">> Simulating Generation "<<worldData->generation);
    if (inv_age >= 0) handleInversionOrigin();
    
    vector <vector<double> > migmap;
    vector<double> mRate;
    vector<double> cRate;
    vector<double> rRate;
    vector<double> hRate;
    vector<double> dRate;
    
    double totalM=0;
    double totalC=0;
    double totalR=0;
    double totalH=0;
    double totalD=0;
    
    for( cluster_t::iterator i = cluster.begin(); i != cluster.end(); ++i ) {
        Context cxt =i->first;
        unsigned int k=ctxtNcarriers(cxt);
        
        if(k>0){
            
            unsigned long clustID=i->second;
            int pop= i->first.pop;															// note that it goes through contexts, not populations
            double  freqI = getFreqI(pop);														//
            double clustSize= worldData->freq.at(clustID) * worldData->popSize.at(pop);
            if (cxt.inversion==0) freqI= 1- freqI;
            
            //UPDATE DEC10, RG: Migration rates adjusted to account for nonconservative migration
            //typical of scenarios like local adaptation
            //originally, the migration rate per individual was m
            // corrected, backward migration rate S1->S2 = m*q1/q2 where qi is the size of the context i
            
            vector<double> mfreq = getSizebyPop(cxt);
            vector<double> migs;
            double b_mig_total=0;
            for(int l=0; l<worldData->nPops;++l){
                if(l!=pop){
                    double q1q2 = 0;
                    if(mfreq.at(pop) > 0) q1q2 = mfreq.at(l)/mfreq.at(pop);
                    b_mig_total += mig_prob.at(pop).at(l) * q1q2;
                    migs.push_back(b_mig_total);
                }
                else migs.push_back(0);
            }
            migmap.push_back(migs);
            
            mRate.push_back(k * b_mig_total);	// rate is number of carriers in context, times probability of migration in population
            totalM += mRate[clustID];		//
            
            if(k > 1){
                unsigned int a=2;
                if(clustSize==0) {std::cerr<<"Error in SimulateGeneration(): Carriers in seemly empty context\n";}
                double ka= binomial_coefficient<double>(k,a)/clustSize;
                cRate.push_back( ka );
                totalC += ka;		//
            }
            
            // --- Recombination rates per carrier with correct partner karyotype ---
            for (int carrierID = 0; carrierID < k; ++carrierID) {
                // Is THIS carrier inverted?
                bool invCarrier = worldData->carriers->at(clustID).at(carrierID)->getInv() == 1;

                // Partner probabilities in this population:
                //   partner_same: prob partner has same karyotype as this carrier
                //   partner_opp : prob partner has opposite karyotype
                double partner_same = invCarrier ? freqI : (1.0 - freqI);
                double partner_opp  = invCarrier ? (1.0 - freqI) : freqI;

                // After collapse (t >= inv_age), no heterokaryotypes anywhere.
                if (inversionOriginHandled) {
                    partner_opp  = 0.0;
                    partner_same = 1.0;
                }

                // Homokaryotypic recombination rate (uses same-karyotype partner)
                double length = partner_same *
                    worldData->carriers->at(clustID).at(carrierID)->getHomoLength(worldData->invRange);
                rRate.push_back(length);
                totalR += length;
                CDBG("                homo = " << length)

                // Heterokaryotypic recombination rate (uses opposite-karyotype partner)
                double hetero = partner_opp *
                    worldData->carriers->at(clustID).at(carrierID)->getHeteroLength(worldData->invRange);
                hRate.push_back(hetero);
                totalH += hetero;
                CDBG("                hetero = " << hetero)

                // Double-cut / gene-flux rate (also needs opposite-karyotype partner)
                double doubrec = partner_opp * worldData->phi;
                dRate.push_back(doubrec);
                totalD += doubrec;
                CDBG("                doubrec = " << doubrec)
            }
            
        }
        else{vector<double> dum (1,0); migmap.push_back(dum);mRate.push_back(0);}
    }
    
    ///////////////////
    
    int nEvents=0;
    double Rate= totalM+ totalC+ totalR+ totalH + totalD;
    double waiting_t= randexp(Rate);
    
    CDBG("total Rate ="<<Rate)
    
    if(worldData->kingman){ //simulation is running on Kingman approx
        if(timeCheck(waiting_t, Rate)) {// waiting time exceeds time to next epoch, need to update freqs/sizes and skip the events
            updateToNextEpoch();
        }
        else{ //waiting time looks good, run Kingman events
            worldData->generation += waiting_t;
            double event= randreal(0,1);
            CDBG("event = "<<event<<"\n")
            if (event < totalM/Rate) nEvents+= migrateEvent(migmap, mRate, totalM);
            else if (event < (totalM+totalC)/Rate) nEvents+=coalesceEvent(cRate, totalC);
            else if (event < (totalM+totalC+totalR)/Rate)			nEvents+= recombineEvent(rRate, totalR, false, false);
            else if (event < (totalM+totalC+totalR+totalH)/Rate)	nEvents+= recombineEvent(hRate, totalH, true, false);
            else													nEvents+= recombineEvent(dRate, totalD, true, true);
        }
    }
    else { // Simulation is running gen-by-gen
        if(timeCheck(1, Rate)) {
            updateToNextEpoch(); // at boundary of epochs, updating freqs/sizes
        }
        else{
            nEvents += migratePoisson(mig_prob);
            nEvents += recombine_all();
            nEvents += coalescePoisson();
            worldData->generation ++;
        }
    }
    
    if( sitesCoalesced()) forceAllCoal();
    
    return nEvents;
}

bool World::timeCheck(double waiting, double Rate){
    double t = waiting + worldData->generation;
    if(!worldData->epochs_over){ //haven't reached end of epochs (or we're simulating equilibrium to infinity, inv_age==0)
        if( t >= worldData->epoch_breaks.at(worldData->current_epoch)) {//waiting time exceeds epoch change
            return true;
        }
        if( Rate ==0) {
            return true; // one carrier remains per context, need to move to next epoch (done in updateToNextEpoch()
        }
    }
    else if(Rate ==0) std::cerr<<"Error in World::timeCheck(), zero rate will lead to infinite waiting time\n";
    return false;
}

void World::updateToNextEpoch(){
    int out = worldData->current_epoch;
    std::cerr << "*** updateToNextEpoch() called for epoch " << out
              << " (time=" << worldData->epoch_breaks[out]
              << ", type=" << worldData->epochType[out] << ")\n";


    int outof_epoch = worldData->current_epoch;
    int into_epoch = worldData->current_epoch + 1;
    worldData->generation = worldData->epoch_breaks[outof_epoch];
    
    if(into_epoch <= worldData->epoch_breaks.size()){
        switch (worldData->epochType[outof_epoch]) {
            // case 0://Inversion age
            //     freqStepToLoss();
            //     break;
            case 1://Speciation event
                DBG("Going into speciation()\t")
                speciation();
                break;
            case 2://Demography change
                DBG("Going into demoChange()\t")
                demoChange();
                break;
            default:
                std::cerr<< "Error in World::updateToNextEpoch(). Epoch type not recognized\n"; exit(1);
                break;
        }
    }
    else{ std::cerr<< "Error in World::updateToNextEpoch(). Current epoch > epoch vector size\n"; exit(1);}

    if(into_epoch == worldData->epoch_breaks.size()){//reached last epoch
        worldData->epochs_over = true;
    }

    worldData->current_epoch = into_epoch;
    DBG("In updateToNextEpoch() from: "<<outof_epoch<<" to: "<<into_epoch<<"...epochType= "<< worldData->epochType[outof_epoch])
}

void World::demoChange() {
    const auto &sizes = worldData->epochPopSizes[worldData->current_epoch];
    if (sizes.empty()) return;

    for (size_t i = 0; i < sizes.size() && i < worldData->popSize.size(); ++i) {
        if (sizes[i] != 0) {
            worldData->popSize[i] = sizes[i];
        }
    }

    std::cerr << "[demoChange] generation " << worldData->generation
              << ": resulting popSizes = ";
    for (auto s : worldData->popSize) std::cerr << s << " ";
    std::cerr << "\n";
}




void World::speciation() {
    // Epoch payload: [A, B, (optional R)]
    const auto& payload = worldData->epochPopSizes[worldData->current_epoch];
    if (payload.size() < 2) {
        std::cerr << "[Speciation] ERROR: event payload too short.\n";
        return;
    }
    unsigned int A = payload[0]; // sink
    unsigned int B = payload[1]; // source
    double newFreq = worldData->epochValues[worldData->current_epoch];

    bool hasR = (payload.size() >= 3);
    unsigned int R = hasR ? payload[2] : 0;

    // Snapshot old cluster & freqs BEFORE we change anything
    const auto oldCluster = cluster;
    const auto oldFreq    = worldData->freq;

    // Debug (before we erase B)
    size_t beforeA = popNcarriers(A);
    size_t beforeB = popNcarriers(B);
    std::cerr << "[Speciation] BEFORE merge: pop "
              << A << "=" << beforeA << " carriers, pop "
              << B << "=" << beforeB << " carriers\n";

    if (hasR) {
        std::cerr << "[Speciation] Using result size R=" << R << " from demes.\n";
        worldData->popSize[A] = R;
    } else {
        std::cerr << "[Speciation] No result size provided; using A += B fallback.\n";
        worldData->popSize[A] += worldData->popSize[B];
    }
    worldData->popSize.erase(worldData->popSize.begin() + B);

    unsigned int newA = (A > B) ? (A - 1) : A;

    // Remap originPop across the A <- B merge and index shift removing B
    if (inv_age >= 0) {
        unsigned int oldOrigin = originPop;
        if (oldOrigin == A || oldOrigin == B) {
            originPop = newA;
        } else if (oldOrigin > B) {
            originPop = oldOrigin - 1;
        }
        if (originPop != oldOrigin) {
            std::cerr << "[Speciation] originPop remapped: "
                    << oldOrigin << " -> " << originPop << "\n";
        }
        // keep the collapse target in sync
        worldData->originCtx = Context(originPop, 0);
    }


    std::vector< std::shared_ptr<Chromosome> > allChr;
    allChr.reserve(totalNCarriers());
    for (auto &vec : *worldData->carriers)
        for (auto &chr : vec)
            allChr.push_back(chr);

    for (auto &chr : allChr) {
        int p = chr->getContext().pop;
        if      (p == (int)B) chr->setPopulation(newA);
        else if (p  > (int)B) chr->setPopulation(p - 1);
    }

    // --- DEBUG: clusters BEFORE rebuild (map still reflects old layout) ---
    std::cerr << "Clusters BEFORE rebuild:\n";
    for (auto &kv : cluster) {
        std::cerr << "  pop=" << kv.first.pop
                  << ", inv=" << kv.first.inversion
                  << " -> cid=" << kv.second << "\n";
    }

    rebuildClustersAndCarriers();

    // --- DEBUG: clusters AFTER rebuild ---
    std::cerr << "Clusters AFTER rebuild:\n";
    for (auto &kv : cluster) {
        std::cerr << "  pop=" << kv.first.pop
                  << ", inv=" << kv.first.inversion
                  << " -> cid=" << kv.second << "\n";
    }


    std::vector<double> newFreqVec(worldData->nClust, 0.0);
    for (auto &kv : cluster) {
        const Context &ctxNew = kv.first;
        int cidNew = kv.second;

        if (ctxNew.pop == (int)newA) {
            newFreqVec[cidNew] = ctxNew.inversion ? newFreq : (1.0 - newFreq);
        } else {
            int oldPop = ctxNew.pop;
            if (oldPop >= (int)B) oldPop += 1;  
            Context ctxOld(oldPop, ctxNew.inversion);
            auto itOld = oldCluster.find(ctxOld);
            if (itOld != oldCluster.end()) {
                newFreqVec[cidNew] = oldFreq[itOld->second];
            }
        }
    }
    worldData->freq.swap(newFreqVec);

    worldData->nPops = worldData->popSize.size();

    size_t afterA = popNcarriers(newA);
    size_t afterB = (B < worldData->nPops) ? popNcarriers(B) : 0;

    std::cerr << "[Speciation] AFTER  merge: pop "
              << newA << "=" << afterA << " carriers, pop "
              << B   << "=" << afterB << " carriers\n";

    std::cerr << "[Speciation] pop " << B << "→" << newA << "; popSizes = ";
    for (auto s : worldData->popSize) std::cerr << s << " ";
    std::cerr << " ; nClust=" << worldData->nClust << "\n";
}




void World::freqStepToLoss(){
    
    for( cluster_t::iterator iter = cluster.begin(); iter != cluster.end(); ++iter ){
        if(worldData->originCtx ==(*iter).first){
            worldData->freq[(*iter).second] = 1.0;
        }
        else{
            worldData->freq[(*iter).second] = 0.0;
        }
    }
    updateContexts(worldData->originCtx);
}

unsigned short World::migrateEvent(vector < vector< double> >& mig_prob, vector<double>& rate, double total ){
    //
    // Move one carrier from its current population to a new one.
    // The function figures out in which population is the migration event, who migrates and where to.
    //
    // Each population has a rate of migration equal to ri=(Ci * Mi) ("number of carriers in that context" times "probability of migration for one carrier in that context")
    // The function calculates these rates, and divides them by the sum of them through all contexts (scales to 1).
    //	A cumulative limit is calculated, to be used then with a random [0,1).
    //
    //			0				r1			r1+r2		r1+r2+r3	  1
    //			|---------------]-------------]-------------]---------]
    //					^				^				^			^
    //			[migration in 1]	[in 2]			[in 3]		[in 4]
    //
    //
    //	This assignment is done with a map< KEY_type double, VALUE_type int> . the map makes KEY/VALUE pairs. We use the lower_bound(x) function, that
    //	returns an iterator to the first element in the map whose key does not compare less than x, i.e. it is either equal or greater.
    //
    //
    //  Modified DEC-10 by RG
    //	DBG("Migration happened in gen "<< nGenerations())
    
    // Initialize things:
    //
    
    map<double, int> whichClust;
    
    double add=0;
    for( cluster_t::iterator i = cluster.begin(); i != cluster.end(); ++i ) {			// this loop scales the rates to one, and constructs the map
        if(ctxtNcarriers(i->first)!=0){													// used to determine in which context migration happened
            whichClust[(rate[i->second] + add) / total] = i->second;
            add += rate[i->second];
        }
    }
    
    unsigned long c = whichClust.lower_bound(randreal(0,1))->second;				// get context where mig happened, using a random number on the map
    unsigned long who= randint(0, worldData->carriers->at(c).size()-1);			// random number that will be the index of the migrant carrier
    
    shared_ptr<Chromosome> chrom = worldData->carriers->at(c).at(who);	// Get the migrant carrier
    
    
    vector<shared_ptr<Chromosome> >::iterator pos = worldData->carriers->at(c).begin() + who;
    worldData->carriers->at(c).erase(pos);
    
    //THIS FIX ASSUMES TWO POPS ONLY!!! MUST CHANGE
    //int whereto=0; if( chrom->getContext().pop==0) whereto=1;
    // Draw a destination using the row-stochastic migration matrix for the source pop
    int from = chrom->getContext().pop;
    double u = randreal(0,1);
    double cum = 0.0;
    int whereto = from;
    for (int j = 0; j < worldData->nPops; ++j) {
        cum += mig_prob.at(from).at(j);
        if (u <= cum) { whereto = j; break; }
    }
    // If we drew "stay in same pop", undo the erase and treat as no-op migration
    if (whereto == from) {
        worldData->carriers->at(c).insert(worldData->carriers->at(c).begin() + who, chrom);
        return 0;
    }
    
    chrom->setPopulation(whereto);													// change the context of this carrier to its new population
    
    int newC= cluster[chrom->getContext()];
    worldData->carriers->at(newC).push_back( chrom );								// add the carrier to the vector of its new population
    
    shared_ptr < ARGNode > newNode ;
    newNode.reset(new ARGNode(worldData->argNodeVec.size(), chrom, worldData->generation));
    worldData->argNodeVec.push_back(newNode);
    chrom->setDescendant(newNode);
    
    
    return 1;
}	//****** end of migrateEvent

unsigned short World::coalesceEvent(vector<double>& rate, double total){
    //
    //	This function executes one simple coalescent event.
    //
    // Written by Rafael Guerrero, VI-09
    //
    
    DBG("Simple coalescent occurred in gen "<< nGenerations())
    map<double, int> whichClust;
    
    int id=0;
    double add=0;
    for( cluster_t::iterator i = cluster.begin(); i != cluster.end(); ++i ) {	// this loop builds the map of cumulative fractions of coalescence (binomial coeff/ total)
        COALDBG("carriers in "<< i->second<<" "<<ctxtNcarriers(i->first))
        if(ctxtNcarriers(i->first) > 1){
            whichClust[(rate.at(id) + add) / total] = i->second;
            COALDBG("rate "<<(rate.at(id) + add))
            COALDBG(total)
            COALDBG("whichClust[] "<<whichClust[(rate.at(id) + add) / total] )
            add += rate.at(id);
            ++id;
            
        }
    }
    
    if(id==0) // no clusters have more than 1 carrier.
    {
        COALDBG("no coalescence possible, all clusters have < 2 carriers")
        return 0;
    }
    
    else{ // carry out the rest of the function.
        int c = whichClust.lower_bound(randreal(0,1))->second;					// randomly decide in which context is coalescence happening
        //
        vector<int> carr_idx;
        for (int u=0; u< worldData->carriers->at(c).size(); ++u) carr_idx.push_back(u);
        random_shuffle(carr_idx.begin(), carr_idx.end());
        
        shared_ptr<Chromosome> chrom = worldData->carriers->at(c).at(carr_idx.at(0));
        shared_ptr<Chromosome> chrom2 = worldData->carriers->at(c).at(carr_idx.at(1));
        
        shared_ptr < ARGNode > newNode ;
        newNode.reset(  new ARGNode(worldData->argNodeVec.size(), chrom, chrom2, worldData->generation) );
        worldData->argNodeVec.push_back(newNode);		// Add a pointer to this node to the vector of all nodes in	the simulation
        
        
        chrom->merge( chrom2 , newNode); // merging includes setting decendant to new ARGNode.
        
        worldData->carriers->at(c).erase(remove(worldData->carriers->at(c).begin(), worldData->carriers->at(c).end(), chrom2), worldData->carriers->at(c).end());
        DBG(totalNCarriers()<<" in "<<worldData->generation<<"\n")
        
        
        //###CCD change imported from sexCoal
        /*
         vector<double>::size_type outit = 0;
         while ( outit < worldData->recomb_breakpoints.size() )
         {
         int myCopy = 0;
         for( cluster_t::iterator i = cluster.begin(); i != cluster.end(); ++i )
         {
         unsigned int clustID=i->second;
         int nLocalCarriers  = worldData->carriers->at(clustID).size();
         
         for(int carrierID = 0; carrierID < nLocalCarriers; ++carrierID)
         {
         if( siteQuery( worldData->recomb_breakpoints.at( outit ), worldData->carriers->at( clustID ).at( carrierID )->getSegmentVector()  ) )
         myCopy ++;
         }
         }
         if( myCopy == 1 )
         {
         std::pair <double, shared_ptr < ARGNode > > myFinishedSegment;
         myFinishedSegment = std::make_pair ( worldData->recomb_breakpoints.at( outit ), newNode );
         worldData->outputSegments.push_back( myFinishedSegment );
         worldData->recomb_breakpoints.erase( worldData->recomb_breakpoints.begin() + outit );
         }else
         {
         ++outit;
         }
         //if myCopy == 1, remove breakpoint from list, and put it into outlist with the argnode
         }
         //cout << endl;
         */
        //##############################
        
        return 1;
    }
}//************************************end of coalesceEvent()
unsigned short World::recombineEvent(vector<double>& rate, double total, bool hetero, bool gflux ){
    
    /*
     // this function runs one recombination event between a carrier and a non-carrier
     // It determines in which context, and who recombines.
     // It adds the resulting new chromosome to the vector of carriers.
     //	It also checks for "empty" carriers after recombination.
     //	In this case, the carriers are not included in the Carriers vector.
     //
     //
     // Given that recombination happened, each carrier has a probability of being the recombinant (Ri) equal to
     //	the total length of its segments divided by the total length of segments in the world.
     //
     //	A cumulative limit is calculated, to be used then with a random [0,1).
     //
     // For example, in a case with four carriers:
     //
     //			0				R1			R1+R2		R1+R2+R3	  1
     //			|---------------]-------------]-------------]---------]
     //					^				^				^			^
     //			[recombine 1]		[ 2]			[ 3]		[ 4]
     //
     //
     //	This assignment is done with a map< KEY_type double, VALUE_type pair< unsigned long, int> > . the map makes KEY/VALUE pairs.
     //	In the example, we would have a map with the pairs < R1 , 1 > , < R1+R2, 2 >, < R1+R2+R3, 3 >, < 1== R1+R2+R3+R4, 4 >
     //	We use the lower_bound(x) function, that returns an iterator to the first element in the map whose key does not
     //	compare less than x, i.e. it is either equal or greater.
     //	In this particular map, the VALUE will be a pair <unsigned long, int>, the complete "address" of the carrier.
     //	This pair is:
     //	- the cluster[context] ( number of the context in the map of all contexts )
     //	- the postion of the carrier in the vector of carriers in that context.
     
     */
    
    RCDBG( "Recombination (hetero= "<<(int)hetero<<", double ="<<(int) gflux <<") happened in gen "<<nGenerations())
    
    map<double, pair<unsigned long,int> > whichCarrier;
    
    int id=0;
    double add=0;
    for( cluster_t::iterator i = cluster.begin(); i != cluster.end(); ++i ) {
        unsigned long clustID=i->second;
        unsigned long nLocalCarriers = worldData->carriers->at(clustID).size();
        
        for(int carrierID = 0; carrierID < nLocalCarriers; ++carrierID){
            if (rate[id]!=0) whichCarrier[(rate[id] + add) / total] = make_pair(clustID, carrierID) ;
            add += rate[id];
            ++id;
        }
    }
    
    pair<unsigned long,int> who = whichCarrier.lower_bound(randreal(0,1))->second;		// randomly decide which carrier will recombine, using the map.
    // "who" is a pair
    // who->first is the number of the context
    // who->second is the number of the carrier in that context
    
    shared_ptr<Chromosome> chrom = worldData->carriers->at(who.first).at(who.second);			// Get the recombining carrier
    
    // create new chromosome
    shared_ptr<Chromosome> chrom2 = recomb_Wrap(chrom, hetero, gflux);
    
    vector<shared_ptr<Chromosome> >::iterator pos = worldData->carriers->at(who.first).begin() + who.second;	// an iterator to the position of chrom. needed only for the function erase()
    worldData->carriers->at(who.first).erase(pos); //erase chrom from old context
    
    if(!chrom->isEmpty()){
        worldData->carriers->at(cluster[chrom->getContext()]).push_back(chrom);	 //add chrom to new context
    }
    
    if(!chrom2->isEmpty())
    {
        worldData->carriers->at(cluster[chrom2->getContext()]).push_back(chrom2);	//	add chrom2 only if it holds segments
    }
    
    RCDBG(totalNCarriers()<<" in "<<worldData->generation<<"\n")
    return 1;
    
}
shared_ptr<Chromosome> World::recomb_Wrap(shared_ptr<Chromosome> chrom, bool hetero, bool gflux){
    
    // Make a new ARG node
    //
    shared_ptr < ARGNode > newNode ;
    newNode.reset(new ARGNode(worldData->argNodeVec.size(), chrom, worldData->generation));
    
    
    // Make two recombinant chromosomes	
    shared_ptr<Chromosome> chrom2;
    int h= static_cast<int> (hetero); if (chrom->getInv()==h) {h=0;} else {h=1;}
    Context other_ctx (chrom->getPopulation(), h);				
    
    if(gflux){
        double mid= (worldData->invRange.R - worldData->invRange.L)/2 + worldData->invRange.L;
        double bp1= randreal(worldData->invRange.L, mid);
        double bp2= randreal(mid, worldData->invRange.R);
        /*	cout<<"bp1 "<<bp1<<", bp2 "<<bp2<<'\n'; 
         cout<<"before"<<'\n';
         cout<<"InvSegs ";for(int i=0;i<chrom->get_Inv_Segs().size();++i)cout<<chrom->get_Inv_Segs().at(i).L<<" "; cout<<'\n';
         cout<<"StdSegs ";for(int i=0;i<chrom->get_Std_Segs().size();++i)cout<<chrom->get_Std_Segs().at(i).L<<" ";cout<<'\n';
         */
        chrom2 = chrom->doubrecombine(bp1, bp2, newNode, other_ctx, worldData->invRange);
        
        /*		cout<<"after "<<chrom->getInv()<<" "<<chrom->get_Inv_Segs().size()<<", chrom2 ->"<<chrom2->getInv()<<" "<<chrom2->get_Inv_Segs().size()<<'\n';
         cout<<"InvSegs ";for(int i=0;i<chrom->get_Inv_Segs().size();++i)cout<<chrom->get_Inv_Segs().at(i).L<<" "; 
         cout<<"InvSegs2 ";for(int i=0;i<chrom2->get_Inv_Segs().size();++i)cout<<chrom2->get_Inv_Segs().at(i).L<<" ";cout<<'\n';
         cout<<"StdSegs ";for(int i=0;i<chrom->get_Std_Segs().size();++i)cout<<chrom->get_Std_Segs().at(i).L<<" ";
         cout<<"StdSegs2 ";for(int i=0;i<chrom2->get_Std_Segs().size();++i)cout<<chrom2->get_Std_Segs().at(i).L<<" ";cout<<'\n';
         */
    }
    else{
        double breakpoint= chrom->calcBreakpoint(worldData->invRange, hetero);
        chrom2 = chrom->recombine(breakpoint, newNode, other_ctx, worldData->invRange);	
        
    }
    
    // add ARGnode to the vector of nodes
    if(!chrom2->isEmpty() && !chrom->isEmpty()) {
        worldData->argNodeVec.push_back(newNode);
    }
    
    return chrom2;
}



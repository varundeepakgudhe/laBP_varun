/*
 *  InvCoal
 *  Version:  labp_v17
 *
 *  This is the main file for the simulation.
 *  Note: Parameter values are now updated via YAML files (other.yaml and demes.yaml)
 *        using the Python wrapper.
 */

// Includes from STL:
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <ctime>
using std::cout;
using std::endl;
#include <fstream>
using std::ifstream;
using std::ofstream;
using std::ostream;
#include <sstream>
using std::istringstream;
using std::stringstream;
#include <string>
using std::string;
#include <random>
#include <chrono>
#include <vector>
#include <map>
#include <iterator>
#include <algorithm>
#include <numeric>
using namespace std;

#include <memory>
using std::shared_ptr;
using std::unique_ptr;
using std::weak_ptr;

// Includes for our files:
#include "parameters.h"
#include "world.h"
#include "ran_mk.h"
#include "chromosome.h"
#include "sitenode.h"
#include "snptree.h"
#include "tajima.h"

// Global declarations:
std::random_device rd;
auto seed = rd();

std::mt19937_64 gen(seed); // Random generator declared globally.

vector<vector<double>> buildMigMatrix(Parameters &p)
{
    // A simple stepping-stone migration matrix.
    unsigned int nPops = static_cast<int>(p.paramData->popSizeVec.size());
    vector<vector<double>> mig_prob;
    mig_prob.resize(nPops);
    double m = p.paramData->migRate.at(0) / (2 * p.paramData->totalPopSize);

    for (unsigned int i = 0; i < nPops; i++)
        mig_prob.at(i).resize(nPops);
    for (unsigned int i = 0; i < nPops; i++)
    {
        for (unsigned int j = 0; j < nPops; j++)
        {
            if (nPops == 1)
                mig_prob.at(i).at(j) = 1;
            else if (i == j)
                mig_prob.at(i).at(j) = 1 - m;
            else if (i == j - 1)
            {
                mig_prob.at(i).at(j) = m / 2;
                if (i == 0)
                    mig_prob.at(i).at(j) = m;
            }
            else if (i == j + 1)
            {
                mig_prob.at(i).at(j) = m / 2;
                if (i == nPops - 1)
                    mig_prob.at(i).at(j) = m;
            }
            else
                mig_prob.at(i).at(j) = 0;
        }
    }

    vector<map<double, int>> mig;
    mig.resize(nPops);
    for (unsigned int i = 0; i < nPops; i++)
    {
        map<double, int> migrate;
        migrate[mig_prob.at(i).at(i)] = i;
        double cumulative = mig_prob.at(i).at(i);
        for (unsigned int j = 0; j < nPops; j++)
        {
            if (i != j && mig_prob.at(i).at(j) != 0)
            {
                cumulative += mig_prob.at(i).at(j);
                migrate[cumulative] = j;
            }
        }
        mig.at(i) = migrate;
    }
    return mig_prob;
}

int main(int argc, const char *argv[])
{
    std::cerr << "Random Seed: " << seed << '\n';

    auto input_seed = std::stoi(argv[1]);
    if (input_seed != -1)
    {
        std::cerr << "Input Seed: " << input_seed << '\n';
        seed = input_seed;
        gen.seed(input_seed);
    }
    // Timer start
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    unsigned int totalEvents = 0;

    stringstream infile;
    stringstream ms_ss;
    stringstream sstat_ss;

    infile << "inLABP.pars";
    ms_ss << "tests/outLABP_" << seed << ".sites";
    sstat_ss << "tests/outLABP_" << seed << ".stats";

    std::vector<std::string> param_vec(argv + 2, argv + argc);
    Parameters params(infile.str().c_str(), param_vec);

    unsigned int nRuns = params.paramData->nRuns;
    unsigned nSites = params.paramData->n_SNPs;

    vector<vector<double>> mig_prob = buildMigMatrix(params);
    vector<double> outTime(params.paramData->n_SNPs, 0);

    std::ofstream msout;
    if (params.paramData->msOutput)
    {
        msout.open(ms_ss.str().c_str());
    }

    std::ofstream stout;
    if (params.paramData->msOutput)
    {
        stout.open(sstat_ss.str().c_str());
    }
    stout << "ET totL S piTotal pi1 pi2 fst dxy tajD\n";

    double LDsum = 0;
    int ticker = (int)(nRuns / 10);

    if (ticker > 0)
        std::cerr << "Progress ticker (1 tick = 10% of runs): ";
    for (int timer = 0; timer < nRuns; ++timer)
    {
        if (ticker > 0 && timer % ticker == 0)
            std::cerr << "+";

        params.setPhi();
        params.setSNPs();
        params.setCarriers();
        unsigned nCarriers = params.paramData->initChr.size();

        World *world = new World(params.getpData());

        while (!world->simulationFinished())
        {
            totalEvents += world->simulateGeneration(mig_prob);
        }

        vector<shared_ptr<ARGNode>> allNodes = world->getARGVec();
        vector<double> tempLD(params.paramData->n_SNPs, 0);

        vector<SNPtree> trees;
        unsigned novar = 0;
        vector<double> varPos;
        double lengthLastSite = 0;

        for (int k = 0; k < nSites; ++k)
        {
            unsigned pos = k;
            if (!params.paramData->fixedS)
            {
                pos = 0;
            }
            SiteNode geneTree = SiteNode(params.paramData->neut_site[pos], allNodes.at(allNodes.size() - 1));
            std::cout << "Value:" << params.paramData->neut_site[pos];
            SNPtree tmp(geneTree, nCarriers, params.paramData->theta);
            unsigned mutantCount = accumulate(tmp.SNPvalues.begin(), tmp.SNPvalues.end(), 0);
            if (mutantCount > 0)
            {
                trees.push_back(tmp);
                varPos.push_back(k);
            }
            else
            {
                ++novar;
            }
            lengthLastSite = tmp.totalLength / (double)params.paramData->totalPopSize;
            double totalmrca = geneTree.getTime() / (double)params.paramData->totalPopSize;
            outTime.at(k) += totalmrca;
            tempLD[k] = totalmrca;
        }

        double pi = 0, piP1 = 0, piP2 = 0, dxy = 0, fst = 0;
        double varSites = nSites - novar;

        for (int i = 0; i < varSites; ++i)
        {
            vector<unsigned> snps = trees.at(i).SNPvalues;
            vector<vector<unsigned>> split = split_by_population(snps, params.getSamplePerPop());
            double totalPi = heterozygosity(snps);
            pi += totalPi;
            vector<double> pipop = pi_by_pop(split);
            piP1 += pipop[0];
            piP2 += pipop[1];
            fst += fst_nei(totalPi, pipop);
            dxy += calcdxy(split);
        }
        pi = pi / nSites;
        piP1 = piP1 / nSites;
        piP2 = piP2 / nSites;
        fst = fst / varSites;
        dxy = dxy / nSites;

        vector<double> posout = params.paramData->neut_site;
        if (!params.paramData->fixedS)
            posout = varPos;

        msout << "// \nsegsites: " << varSites << "\npositions: ";
        for (int k = 0; k < varSites; ++k)
        {
            msout << " " << posout.at(k) << " ";
        }
        msout << '\n';
        for (int i = 0; i < nCarriers; ++i)
        {
            for (int k = 0; k < varSites; ++k)
            {
                msout << trees[k].SNPvalues[i];
            }
            msout << '\n';
        }

        double tempProd = tempLD[0] * tempLD[nSites - 1];
        LDsum += tempProd;

        stout << tempLD[0] << " " << lengthLastSite << " " << varSites << " " << pi << " " << piP1 << " " << piP2 << " " << fst << " " << dxy << " " << tajd(nCarriers, varSites, pi) << '\n';

        delete world;
    }

    std::cerr.precision(6);
    std::cerr << "\nMean LD (E[T1,n]) = " << LDsum / (double)nRuns << '\n';

    if (!params.paramData->fixedS)
    {
        std::cerr << "Mean TMRCA = " << outTime.at(0) / nRuns << '\n';
    }
    else
    {
        std::cerr << "Mean TMRCA per site (E[T1]...E[Tn])\n";
        for (int k = 0; k < nSites; ++k)
        {
            std::cerr << " " << outTime.at(k) / nRuns;
        }
        std::cerr << "\n";
    }

    msout.close();
    stout.close();

    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cerr << "Elapsed time: " << elapsed_seconds.count() << "s\n";

    return 0;
}

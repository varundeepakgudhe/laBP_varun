/**
 * Parameters.cpp
 *  InvCoal
 *  Version:  labp_v17
 *
 *  This file handles the input parameters.
 */

// Includes from STL:
#include <string>
using std::string;
#include <fstream>
using std::ifstream;
#include <iostream>
using std::cout;
using std::endl;
#include <sstream>
using std::istringstream;
#include <vector>
using std::vector;
#include <iterator>
#include <algorithm>
using std::sort;
#include <cassert>
#include <memory>
using std::shared_ptr;
#include <cmath>
#include <numeric>

// Includes for our files:
#include "typedefs.h"
#include "parameters.h"
#include "chromosome.h"
#include "ran_mk.h"

//======================================================================

Parameters::Parameters(const char *insstring, const std::vector<std::string> &param_vec)
{
    DBG("Parameters:: Constructing...")

    paramData.reset(new ParameterData);

    int temp;
    double doubtemp;
    std::istringstream iss;

    // Read in different params (order must match the list passed from the Python wrapper)

    // nruns (param_vec[0])
    paramData->nRuns = std::stoi(param_vec[0]);
    std::cerr << "Number of runs: " << paramData->nRuns << '\n';

    // kingman coalescent flag (param_vec[1])
    paramData->kingman = (param_vec[1] == "1");
    std::cerr << "Assuming Kingman? " << paramData->kingman << '\n';

    // drift simulation flag (param_vec[2])
    paramData->drift = (param_vec[2] == "1");
    std::cerr << "Simulating drift? " << paramData->drift << '\n';
    if (paramData->drift)
    {
        std::cerr << "Sorry, drift not implemented in this version yet\n";
        exit(1);
    }

    // msOutput flag (param_vec[3])
    paramData->msOutput = (param_vec[3] == "1");
    std::cerr << "Output each set in ms format? " << paramData->msOutput << '\n';

    // Population Sizes (param_vec[4])
    iss.str(param_vec[4]);
    paramData->popSizeVec = vector<unsigned int>(std::istream_iterator<int>(iss), std::istream_iterator<int>());
    std::cerr << "Population Sizes: ";
    for (size_t n = 0; n < paramData->popSizeVec.size(); n++)
    {
        std::cerr << n << " = " << paramData->popSizeVec.at(n) << " | ";
    }
    std::cerr << '\n';
    paramData->totalPopSize = std::accumulate(paramData->popSizeVec.begin(), paramData->popSizeVec.end(), 0);

    unsigned int pops = paramData->popSizeVec.size();

    // Inversion initial frequencies (param_vec[5])
    iss.clear();
    iss.str(param_vec[5]);
    paramData->initialFreqs = vector<double>(std::istream_iterator<double>(iss), std::istream_iterator<double>());
    if (paramData->initialFreqs.size() != pops)
    {
        std::cerr << "Inconsistent number of sample sizes found (should be equal to number of pops)\n";
        exit(1);
    }
    std::cerr << "Frequencies of inverted chromosomes: ";
    for (size_t n = 0; n < paramData->initialFreqs.size(); n++)
    {
        std::cerr << n << " = " << paramData->initialFreqs.at(n) << " | ";
    }
    std::cerr << '\n';

    // Speciation (param_vec[6])
    iss.clear();
    iss.str(param_vec[6]);
    paramData->speciation = vector<double>(std::istream_iterator<double>(iss), std::istream_iterator<double>());
    if (paramData->speciation.at(0) == 1)
    {
        if (paramData->speciation.size() < 3)
        {
            std::cerr << "Error; Simulating speciation but no time or invFreq parameters found\n";
            exit(1);
        }
        std::cerr << "Speciation time: " << paramData->speciation.at(1) << "\n";
    }
    else
    {
        std::cerr << "No speciation.\n";
    }

    // Demography (param_vec[7])
    iss.clear();
    iss.str(param_vec[7]);
    paramData->demography = vector<double>(std::istream_iterator<double>(iss), std::istream_iterator<double>());
    if (paramData->demography.at(0) == 1)
    {
        if (paramData->demography.size() < 3)
        {
            std::cerr << "Error; Simulating demography but no time parameter found\n";
            exit(1);
        }
        std::cerr << "Demography change time: " << paramData->demography.at(1) << "\n";
        std::cerr << "Reminder: this coeff multiplies all pop sizes in the preceding epoch\n";
    }
    else
    {
        std::cerr << "No demographic event.\n";
    }

    // Age of inversion (param_vec[8])
    paramData->inv_age = std::stoi(param_vec[8]);
    std::cerr << "Age of inversion = " << paramData->inv_age << '\n';

    // Migration rate (param_vec[9])
    iss.clear();
    iss.str(param_vec[9]);
    paramData->migRate = vector<double>(std::istream_iterator<double>(iss), std::istream_iterator<double>());
    std::cerr << "Migration rate (pop0<-->1) = " << paramData->migRate[0] << '\n';

    // Bases per Morgan (param_vec[10])
    paramData->BasesPerMorgan = std::stod(param_vec[10]);
    std::cerr << "Bases per Morgan in homokaryotypic recombination " << paramData->BasesPerMorgan << '\n';

    // Random phi flag (param_vec[11])
    paramData->randPhi = (param_vec[11] == "1");
    std::cerr << "Random phi values? " << paramData->randPhi << '\n';

    // Phi range (param_vec[12])
    iss.clear();
    iss.str(param_vec[12]);
    paramData->phi_range = vector<double>(std::istream_iterator<double>(iss), std::istream_iterator<double>());
    if (paramData->randPhi)
    {
        std::cerr << "Random gene flux (phi) range = " << paramData->phi_range[0] << " - " << paramData->phi_range[1] << '\n';
    }
    else
    {
        std::cerr << "Gene flux (phi) = " << paramData->phi_range[0] << '\n';
    }

    // Inversion range (param_vec[13])
    vector<double> invtemp;
    iss.clear();
    iss.str(param_vec[13]);
    invtemp = vector<double>(std::istream_iterator<double>(iss), std::istream_iterator<double>());
    paramData->invRange.L = invtemp[0] / paramData->BasesPerMorgan;
    paramData->invRange.R = invtemp[1] / paramData->BasesPerMorgan;
    std::cerr << "Inversion from: " << invtemp[0] << " to " << invtemp[1] << " (" << paramData->invRange.L << " - " << paramData->invRange.R << " recUnits)\n";

    // Fixed S and Theta (param_vec[14])
    vector<double> stemp;
    iss.clear();
    iss.str(param_vec[14]);
    stemp = vector<double>(std::istream_iterator<double>(iss), std::istream_iterator<double>());
    paramData->fixedS = static_cast<bool>(stemp[0]);
    paramData->n_SNPs = static_cast<int>(stemp[1]);
    if (paramData->fixedS)
    {
        paramData->theta = 0;
        std::cerr << "Number of markers to simulate: " << paramData->n_SNPs << '\n';
    }
    else
    {
        paramData->theta = stemp[2];
        std::cerr << "Number of bases (non-recombining) to simulate: " << paramData->n_SNPs << ", with mutation rate: " << paramData->theta << '\n';
    }

    // Random SNP positions flag (param_vec[15])
    paramData->randSNP = (param_vec[15] == "1");
    std::cerr << "Markers in random locations? " << paramData->randSNP << '\n';

    // SNP positions (param_vec[16])
    iss.clear();
    iss.str(param_vec[16]);
    while (iss >> doubtemp)
    {
        paramData->snpPositions.push_back(doubtemp / paramData->BasesPerMorgan);
    }
    std::cerr << paramData->snpPositions.size() << " Site positions read.\n";

    if (paramData->randSNP)
    {
        paramData->snpRange.push_back(paramData->snpPositions[0]);
        paramData->snpRange.push_back(paramData->snpPositions[paramData->snpPositions.size() - 1]);
        if (paramData->snpPositions.size() != 2)
        {
            std::cerr << "Warning: Asking for random SNPs, but number of SNPs read !=2. Will only take first and last SNP positions as the range for random SNP locations\n";
        }
        if (!paramData->fixedS)
        {
            std::cerr << "Warning: Number of segregating sites (S) is not fixed. SNPRange will be ignored. Random SNPs will be simulated in non-recombining region, using theta parameter.\n";
        }
        std::cerr << "Window of SNP locations (in recUnits): " << paramData->snpRange[0] << " - " << paramData->snpRange[paramData->snpRange.size() - 1] << '\n';
    }

    if (!paramData->fixedS)
    {
        paramData->snpPositions.resize(1);
        paramData->randSNP = false;
    }

    // Random sample flag (param_vec[17])
    bool randomSample = (param_vec[17] == "1");
    std::cerr << "Random sample of carriers? " << randomSample << '\n';

    paramData->nCarriers.resize(pops);
    if (randomSample)
    {
        vector<int> tempRead;
        iss.clear();
        iss.str(param_vec[18]);
        tempRead = vector<int>(std::istream_iterator<int>(iss), std::istream_iterator<int>());
        if (tempRead.size() != pops)
        {
            std::cerr << "Inconsistent number of sample sizes found (should be equal to number of pops)\n";
            exit(1);
        }
        for (unsigned int p = 0; p < pops; ++p)
        {
            int invcount = randbinom(tempRead[p], paramData->initialFreqs.at(p));
            paramData->nCarriers.at(p).push_back(tempRead[p] - invcount);
            paramData->nCarriers.at(p).push_back(invcount);
        }
    }
    else
    {
        int totalSample = 0;
        for (unsigned int p = 0; p < pops; ++p)
        {
            std::cerr << "Sample in Pop " << p << ": ";
            iss.clear();
            iss.str(param_vec[19]);
            while (iss >> temp)
            {
                paramData->nCarriers.at(p).push_back(temp);
                totalSample += temp;
                std::cerr << temp << " ";
            }
            std::cerr << '\n';
        }
        if (totalSample == 2)
            std::cerr << "Sample size = 2 || Warning: Informative sites length function won't work.\n";
    }
}

Parameters::~Parameters(){
    DBG("Parameters:: Destructing...")
}

vector<unsigned int> Parameters::getPopulationSizes()
{
    return paramData->popSizeVec;
}

vector<unsigned> Parameters::getSamplePerPop()
{
    vector<unsigned> popVal;
    for (unsigned i = 0; i < paramData->nCarriers.size(); ++i)
    {
        unsigned int ncarr = std::accumulate(paramData->nCarriers.at(i).begin(), paramData->nCarriers.at(i).end(), 0);
        for (unsigned j = 0; j < ncarr; ++j)
            popVal.push_back(i);
    }
    return popVal;
}

vector<shared_ptr<Chromosome>> Parameters::getChromVec()
{
    return paramData->initChr;
}

shared_ptr<Parameters::ParameterData> Parameters::getpData()
{
    return paramData;
}

void Parameters::setCarriers()
{
    vector<Segment> initSegVec;
    for (size_t j = 0; j < paramData->neut_site.size(); ++j)
    {
        Segment initSeg;
        initSeg.L = paramData->neut_site[j];
        initSeg.R = paramData->neut_site[j];
        initSegVec.push_back(initSeg);
    }
    paramData->initChr.clear();
    for (unsigned int p = 0; p < paramData->popSizeVec.size(); ++p)
    {
        for (int inv = 0; inv < 2; ++inv)
        {
            Context initCtxt(p, inv);
            shared_ptr<ARGNode> nulle;
            for (int n = 0; n < paramData->nCarriers.at(p).at(inv); ++n)
            {
                shared_ptr<Chromosome> chr(new Chromosome(initCtxt, initSegVec, paramData->invRange, nulle));
                paramData->initChr.push_back(chr);
            }
        }
    }
}

void Parameters::setPhi()
{
    if (paramData->randPhi)
    {
        double phiExponent = randreal(paramData->phi_range[0], paramData->phi_range[1]);
        paramData->phi = pow(10, phiExponent);
    }
    else
    {
        paramData->phi = paramData->phi_range[0];
    }
}

void Parameters::setSNPs()
{
    paramData->neut_site.clear();
    if (paramData->randSNP)
    {
        for (int p = 0; p < paramData->n_SNPs; ++p)
        {
            paramData->neut_site.push_back(randreal(paramData->snpRange[0], paramData->snpRange[1]));
        }
        sort(paramData->neut_site.begin(), paramData->neut_site.end());
    }
    else
    {
        paramData->neut_site = paramData->snpPositions;
    }
}

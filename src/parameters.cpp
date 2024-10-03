/**
 * Parameters.cpp
 *  InvCoal
 *  Version:  labp_v17
 *
 *
 */

// Includes from STL:
#include <string>
using std::string;
#include <fstream>
using std::ifstream;
#include <iostream>
using std::cout;
using std::endl;
#include <string>
using std::string;
#include <sstream>
using std::istringstream;
#include <vector>
using std::vector;
#include <algorithm> // std::sort
using std::sort;
#include <assert.h>
#include <memory>
using std::shared_ptr;
#include <math.h>

// Includes from our files:
#include "typedefs.h"
#include "parameters.h"
#include "chromosome.h"
#include "ran_mk.h"

//======================================================================

Parameters::Parameters(const char *insstring, const std::vector<std::string> &param_vec)
{ // Constructor.  Reads parameter values from the file parameters.txt

    DBG("Parameters:: Constructing...")

    paramData.reset(new ParameterData);

    int temp;
    double doubtemp;
    std::istringstream iss;

    // Read in different params
    //  number of runs
    paramData->nRuns = std::stoi(param_vec[0]);
    std::cerr << "Number of runs: " << paramData->nRuns << '\n';

    // Should we simulate kingman coalescent? (i.e. exponential waiting times). If false, does generation by generation
    paramData->kingman = (param_vec[1] == "1"); // this can also be assigned to true
    std::cerr << "Assuming Kingman? " << paramData->kingman << '\n';

    // Simulate drift trajectory?
    paramData->drift = (param_vec[2] == "1");
    std::cerr << "Simulating drift? " << paramData->drift << '\n';

    if (paramData->drift)
    {
        std::cerr << "Sorry, drift not implemented in this version yet\n";
        exit(1);
    }

    // Output ms-formatted haplotypes?
    paramData->msOutput = (param_vec[3] == "1");
    std::cerr << "Output each set in ms format? " << paramData->msOutput << '\n';

    // Population Sizes
    iss.str(param_vec[4]);
    paramData->popSizeVec = std::vector<unsigned int>(std::istream_iterator<int>(iss), std::istream_iterator<int>());

    std::cerr << "Population Sizes: ";
    for (int n = 0; n < (int)paramData->popSizeVec.size(); n++)
    {
        std::cerr << (n) << " = " << paramData->popSizeVec.at(n) << " | ";
    }
    std::cerr << '\n';
    paramData->totalPopSize = std::accumulate(paramData->popSizeVec.begin(), paramData->popSizeVec.end(), 0);

    unsigned int pops = paramData->popSizeVec.size();

    // Inversion initial frequencies
    iss.clear();
    iss.str(param_vec[5]);
    paramData->initialFreqs = std::vector<double>(std::istream_iterator<double>(iss), std::istream_iterator<double>());

    if (paramData->initialFreqs.size() != pops)
    {
        std::cerr << "Inconsistent number of sample sizes found (should be equal to number of pops)\n";
        exit(1);
    }

    std::cerr << "Frequencies of inverted chromosomes: ";
    for (int n = 0; n < (int)paramData->initialFreqs.size(); n++)
    {
        std::cerr << (n) << " = " << paramData->initialFreqs.at(n) << " | ";
    }
    std::cerr << '\n';

    // Simulate speciation?
    iss.clear();
    iss.str(param_vec[6]);
    paramData->speciation = std::vector<double>(std::istream_iterator<double>(iss), std::istream_iterator<double>());

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

    // Simulate demography?
    iss.clear();
    iss.str(param_vec[7]);
    paramData->demography = std::vector<double>(std::istream_iterator<double>(iss), std::istream_iterator<double>());

    if (paramData->demography.at(0) == 1)
    {
        if (paramData->speciation.size() < 3)
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

    // Age of inversion in number of generations (EX: 1500)
    paramData->inv_age = std::stoi(param_vec[8]);

    // change the below ones also to paramData->inv_age instead of param_vec[8]
    std::cerr << "Age of inversion = " << paramData->inv_age << '\n';

    if (paramData->inv_age > 0 && paramData->inv_age < std::max(paramData->demography.at(1), paramData->speciation.at(1)))
    {
        std::cerr << "Error: Current implementations only allows for inversions that predate speciation and demographic changes\n";
        exit(1);
    }

    // Migration rate (4Nm)
    // why is migRate a vector?
    iss.clear();
    iss.str(param_vec[9]);
    paramData->migRate = std::vector<double>(std::istream_iterator<double>(iss), std::istream_iterator<double>());

    std::cerr << "Migration rate (pop0<-->1) = " << paramData->migRate[0] << '\n';

    // bases per Morgan in the hh genotype.
    paramData->BasesPerMorgan = std::stod(param_vec[10]);
    std::cerr << "Bases per Morgan in homokaryotypic recombination " << paramData->BasesPerMorgan << '\n';

    // Random phi?
    paramData->randPhi = (param_vec[11] == "1");
    std::cerr << "Random phi values? (if 1, range below. if 0, single value read below) " << paramData->randPhi << '\n';

    // why is phi_range a  vector?
    iss.clear();
    iss.str(param_vec[12]);
    paramData->phi_range = std::vector<double>(std::istream_iterator<double>(iss), std::istream_iterator<double>());

    if (paramData->randPhi)
    {
        std::cerr << "Random gene flux (phi) range = " << paramData->phi_range[0] << " - " << paramData->phi_range[1] << '\n';
    }
    else
    {
        std::cerr << "Gene flux (phi) = " << paramData->phi_range[0] << '\n';
    }

    // the range of the inversion
    vector<double> invtemp;

    iss.clear();
    iss.str(param_vec[13]);
    invtemp = std::vector<double>(std::istream_iterator<double>(iss), std::istream_iterator<double>());

    paramData->invRange.L = invtemp[0] / paramData->BasesPerMorgan;
    paramData->invRange.R = invtemp[1] / paramData->BasesPerMorgan;
    std::cerr << "Inversion from: " << invtemp[0] << " to " << invtemp[1] << " (" << paramData->invRange.L << " - " << paramData->invRange.R << " recUnits)\n";

    // Fixed S? Value of S or sequence length (bases), Theta
    vector<double> stemp;

    iss.clear();
    iss.str(param_vec[14]);
    stemp = std::vector<double>(std::istream_iterator<double>(iss), std::istream_iterator<double>());

    paramData->fixedS = (bool)stemp[0];
    paramData->n_SNPs = (int)stemp[1];
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

    // random positions of SNPs?
    paramData->randSNP = (param_vec[15] == "1");
    std::cerr << "Markers in random locations? " << paramData->randSNP << '\n';

    // range (positions of begin-end) where the SNPs are -- in bases
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

    // should the sample be random?
    bool randomSample = false;

    randomSample = (param_vec[17] == "1");
    std::cerr << "Random sample of carriers? " << randomSample << '\n';

    paramData->nCarriers.resize(pops);

    if (randomSample)
    {
        vector<int> tempRead;

        iss.clear();
        iss.str(param_vec[18]);
        tempRead = std::vector<int>(std::istream_iterator<int>(iss), std::istream_iterator<int>());

        if (tempRead.size() != pops)
        {
            std::cerr << "Inconsistent number of sample sizes found (should be equal to number of pops)\n";
            exit(1);
        }

        for (int p = 0; p < pops; ++p)
        {
            int invcount = randbinom(tempRead[p], paramData->initialFreqs.at(p));
            paramData->nCarriers.at(p).push_back(tempRead[p] - invcount);
            paramData->nCarriers.at(p).push_back(invcount);
        }
    }
    else
    { // number of S, I carriers in each pop
        int totalSample = 0;
        for (int p = 0; p < pops; ++p) // why for loop, its gonna execute only once right?
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

} // end constructor

Parameters::~Parameters(){// Destructor
                          DBG("Parameters:: Destructing...")}

vector<unsigned int> Parameters::getPopulationSizes()
{ // Returns vector of population sizes (over populations)
    return paramData->popSizeVec;
}

vector<unsigned> Parameters::getSamplePerPop()
{
    vector<unsigned> popVal;
    for (unsigned i = 0; i < paramData->nCarriers.size(); ++i)
    {
        unsigned int ncarr = accumulate(paramData->nCarriers.at(i).begin(), paramData->nCarriers.at(i).end(), 0);
        for (unsigned j = 0; j < ncarr; ++j)
            popVal.push_back(i);
    }
    return popVal;
}

vector<shared_ptr<Chromosome>> Parameters::getChromVec()
{ // Returns the prototype chromosome for the initial carriers
    return paramData->initChr;
}

shared_ptr<Parameters::ParameterData> Parameters::getpData()
{
    return paramData;
}

void Parameters::setCarriers()
{

    vector<Segment> initSegVec;
    for (int j = 0; j < paramData->neut_site.size(); ++j)
    {
        Segment initSeg;
        initSeg.L = paramData->neut_site[j];
        initSeg.R = paramData->neut_site[j];
        initSegVec.push_back(initSeg);
    }
    paramData->initChr.clear();

    for (int p = 0; p < paramData->popSizeVec.size(); ++p)
    {
        for (int inv = 0; inv < 2; ++inv)
        {

            Context initCtxt(p, inv);
            shared_ptr<ARGNode> nulle;

            for (int n = 0; n < paramData->nCarriers.at(p).at(inv); ++n)
            {
                // Make a prototype chromosome:
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
    // std::cerr.precision(12);
    // std::cerr<<"\nSetting SNPs at: ";for(int p=0;p< paramData->n_SNPs;++p){std::cerr<<paramData->neut_site[p]<<" ";}
    // std::cerr<<'\n';
}

void Parameters::setTimeStamp(double t)
{
    paramData->timestamp = t;
}

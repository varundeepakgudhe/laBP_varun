/**
 *  parameters.h
 *  InvCoal
 *  Version:  labp_v17
 *
 *  This header defines the Parameters class.
 */

 #ifndef PARAMETERS_H_
 #define PARAMETERS_H_
 
 #include <memory>
 using std::shared_ptr;
 using std::unique_ptr;
 #include <vector>
 #include <numeric>
 
 #include "typedefs.h"
 
 // Forward declaration:
 class Chromosome;
 
 class Parameters
 {
 public:
     struct ParameterData;
     shared_ptr<ParameterData> paramData;
 
     Parameters();
     Parameters(const char *insstring, const std::vector<std::string> &param_vec);
     ~Parameters();
 
     std::vector<unsigned int> getPopulationSizes();
     std::vector<unsigned int> getSamplePerPop();
     std::vector<shared_ptr<Chromosome>> getChromVec();
     shared_ptr<ParameterData> getpData();
     void setPhi();
     void setCarriers();
     void setSNPs();
     void setTimeStamp(double t);
 };
 
 struct Parameters::ParameterData
 {
     unsigned int nRuns;
     std::vector<unsigned int> popSizeVec;
     std::vector<double> migRate;
     std::vector<double> initialFreqs;
     std::vector<double> speciation;
     std::vector<double> demography;
     unsigned int inv_age;
     std::vector<double> phi_range;
     Segment invRange;
     double phi;
     double theta;
     double BasesPerMorgan;
     bool kingman;
     bool drift;
     bool msOutput;
     unsigned int n_SNPs;
     bool randSNP;
     bool randPhi;
     bool fixedS;
     std::vector<double> snpRange;
     std::vector<double> snpPositions;
     std::vector<std::vector<int>> nCarriers;
     std::vector<shared_ptr<Chromosome>> initChr;
     std::vector<double> neut_site;
     double timestamp;
     unsigned long totalPopSize;
 };
 
 #endif /* PARAMETERS_H_ */
 
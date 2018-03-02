//
//   @file    fastnlo_utils.cxx
//
//
//   @author D. Britzger, K. Rabbertz
//
//
// Get preprocessor defines from autotools setup
#include "amconfig.h"

// Require fastNLO
#ifdef HAVE_FASTNLO

// System includes
#include <cfloat>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <algorithm> //c++98
#include <utility> //c++11
#include <sys/stat.h> //c++98

// NNLO bridge includes
#include "nnlo_common.h"
#include "nnlo_utils.h"

// fastNLO includes
#include "fastnlo_utils.h"
//#include "fastnlotk/fastNLOCreate.h"



//-----------------------------------------------------------------------------
//
// --- --- --- fnloUtils --- --- --
//
//-----------------------------------------------------------------------------
extern "C" void   evolvepdf_(const double& x, const double& Q, double* xf);
extern "C" double alphaspdf_(const double& Q);
extern "C" bool getunitphase_();



//
// --- fastNLO specific functions
//

// _____________________________________________________________________ //
fastNLO::GeneratorConstants  fnloUtils::GetNNLOJET_GenConsts() {
   //!< Get NNLOJET specific generator constants for this interface
   fastNLO::GeneratorConstants gc;
   gc.Name = "NNLOJET_rev_4585"; //!< Name and version of generator
   gc.References.push_back("Single jet inclusive production (adapt in case of other process):"); //!< References for generator
   gc.References.push_back("Currie, J.; Glover, E. W. N. & Pires, J., Phys. Rev. Lett., 2017, 118, 072002 [arXiv:1611.01460]"); //!< References for generator
   gc.References.push_back("Currie, J.; Glover, E. W. N.; Gehrmann-De Ridder, A.; Gehrmann, T.; Huss, A. & Pires, J., In Proceedings 23rd Cracow Epiphany Conference, 2017, [arXiv:1704.00923]"); //!< References for generator
   gc.References.push_back("Currie, J.; Gehrmann-De Ridder, A.; Gehrmann, T.; Glover, E. W. N.; Huss, A. & Pires, J., In Proceedings 52nd Rencontres de Moriond QCD, 2017, [arXiv:1705.08205]"); //!< References for generator
   gc.UnitsOfCoefficients = 15; //!< fb
   return gc;
}



// _____________________________________________________________________ //
fastNLO::ProcessConstants fnloUtils::GetNNLOJET_ProcConstsPP() {
   //!< Get default process constants for pp processes

   // --- external input from NNLOJET
   std::vector<std::vector<int> > partLiCos = fnloUtils::sp.GetPartonCombinations();
   if ( partLiCos.empty() ) {
      std::cerr<<"[fnloUtils::GetNNLOJET_ProcConstsPP] Error. No parton combinations found. "
               << "Please read config file into subprocess-class: fnloUtils::sp. Aborting!"<<std::endl;
      exit(5);
   }

   fastNLO::ProcessConstants pc;
   pc.LeadingOrder     = nnlo::GetPowLO(); //!< Power in alpha_s of LO process
   pc.NPDF             = 2;                //!< No. of PDFs involved
   int nprocs = partLiCos.size();          // TODO Why only LO?
   pc.NSubProcessesLO  = nprocs;           //!< No. of LO   subprocesses
   pc.NSubProcessesNLO = nprocs;           //!< No. of NLO  subprocesses
   pc.NSubProcessesNNLO= nprocs;           //!< No. of NNLO subprocesses
   pc.IPDFdef1         = 3;                //!< Flag 1 to define PDF linear combinations of partonic subprocesses
   pc.IPDFdef2         = 0;                //!< Flag 2 to define PDF linear combinations of partonic subprocesses
   pc.IPDFdef3LO       = nprocs;           //!< Flag 3 to define PDF LCs at   LO
   pc.IPDFdef3NLO      = nprocs;           //!< Flag 3 to define PDF LCs at  NLO
   pc.IPDFdef3NNLO     = nprocs;           //!< Flag 3 to define PDF LCs at NNLO
   pc.NPDFDim          = 2;                //!< Internal storage mode for PDF LCs, 2: full-matrix ==> we do not need 'asymmetric processes'
   //!< PDF Linear combinations in two different formats (used only if IPDFdef2==0)
   pc.PDFCoeffLO       = VectorToPairs(partLiCos); //!< PDF Linear combinations
   pc.PDFCoeffNLO      = pc.PDFCoeffLO;
   pc.PDFCoeffNNLO     = pc.PDFCoeffLO;
   pc.PDFLiCoInLO      = partLiCos;                //!< PDF Linear combinations
   pc.PDFLiCoInNLO     = pc.PDFLiCoInLO;
   pc.PDFLiCoInNNLO    = pc.PDFLiCoInLO;
   pc.AsymmetricProcesses.clear();                 //!< Specify processes that need to be exchanged in half-matrix notation, when xmin>xmax (only if NPDFDim==1)
   pc.Name             = nnlo::GetProcessName();   //!< Name of process as given in NNLOJET runcard
   pc.References.push_back("Please check and modify NNLOJET references according to selected process."); // references for processes
   return pc;
}



// _____________________________________________________________________ //
fastNLO::ProcessConstants fnloUtils::GetNNLOJET_ProcConstsDIS() {
   //!< Get default process constants for DIS processes

   // --- external input from NNLOJET
   std::vector<std::vector<int> > partLiCos = fnloUtils::sp.GetPartonCombinations();
   if ( partLiCos.empty() ) {
      std::cerr<<"[fnloUtils::GetNNLOJET_ProcConstsDIS] Error. No parton combinations found. "
               << "Please read config file into subprocess-class: fnloUtils::sp. Aborting!"<<std::endl;
      exit(5);
   }

   fastNLO::ProcessConstants pc;
   pc.LeadingOrder     = nnlo::GetPowLO(); //!< Power in alpha_s of LO process
   pc.NPDF             = 1;                //!< No. of PDFs involved
   int nprocs = partLiCos.size();          // TODO Why only LO?
   pc.NSubProcessesLO  = nprocs;           //!< No. of LO   subprocesses
   pc.NSubProcessesNLO = nprocs;           //!< No. of NLO  subprocesses
   pc.NSubProcessesNNLO= nprocs;           //!< No. of NNLO subprocesses
   pc.IPDFdef1         = 2;                //!< Flag 1 to define PDF linear combinations of partonic subprocesses
   pc.IPDFdef2         = 0;                //!< Flag 2 to define PDF linear combinations of partonic subprocesses
   pc.IPDFdef3LO       = nprocs;           //!< Flag 3 to define PDF LCs at   LO
   pc.IPDFdef3NLO      = nprocs;           //!< Flag 3 to define PDF LCs at  NLO
   pc.IPDFdef3NNLO     = nprocs;           //!< Flag 3 to define PDF LCs at NNLO
   pc.NPDFDim          = 0;                //!< Internal storage mode for PDF LCs, 2: full-matrix ==> we do not need 'asymmetric processes'
   pc.IPDFdef1         = 2; // Flag 1 to define PDF linear combinations
   pc.IPDFdef2         = 0; // Flag 2 to define PDF linear combinations
   pc.IPDFdef3LO       = nprocs; // No. of  subprocesses
   pc.IPDFdef3NLO      = nprocs;
   pc.IPDFdef3NNLO     = nprocs;
   //!< PDF Linear combinations in two different formats (used only if IPDFdef2==0)
   pc.PDFCoeffLO       = VectorToPairs(partLiCos); //!< PDF Linear combinations
   pc.PDFCoeffNLO      = pc.PDFCoeffLO;
   pc.PDFCoeffNNLO     = pc.PDFCoeffLO;
   pc.PDFLiCoInLO      = partLiCos; //!< PDF Linear combinations
   pc.PDFLiCoInNLO     = pc.PDFLiCoInLO;
   pc.PDFLiCoInNNLO    = pc.PDFLiCoInLO;
   pc.AsymmetricProcesses.clear();                 //!< Specify processes that need to be exchanged in half-matrix notation, when xmin>xmax (only if NPDFDim==1)
   pc.Name             = nnlo::GetProcessName();   //!< Name of process as given in NNLOJET runcard
   pc.References.push_back("Please check and modify NNLOJET references according to selected process."); // references for processes
   return pc;
}



// _____________________________________________________________________ //
fastNLO::ScenarioConstants fnloUtils::GetNNLOJET_ScenConsts() {

   //!< Get default process constants for pp processes
   fastNLO::ScenarioConstants sc;

   sc.ScenarioName = "fnlXYZ";     // No white space here!
   sc.ScenarioDescription.clear();
   sc.ScenarioDescription.push_back("Add here full description of observable"); //< Description of the scenario
   sc.PublicationUnits = 12;       // LHC Rivet analyses are usually working in units of [pb]
   sc.DifferentialDimension = 1;   // NNLOJET always provides 1d histograms

   // Labels (symbol and unit) for the measurement dimensions (from outer to inner "loop"),
   // e.g. "|y|" and "p_T_[GeV]".
   // This may also help to define the observables to be calculated in an automatized way!
   sc.DimensionLabels.clear();
   sc.DimensionLabels.push_back("varName_[unit]"); // Could eventually be extracted from NNLOJET runcard (grid name)

   // Specify:
   //  0 : the cross section is NOT differential,
   //      i.e. there are two bin borders (but NO division (normalization) by bin width);
   //  1 : the cross section is point-wise differential, i.e. only one point is given;
   //  2 : the cross section is bin-wise differential,   i.e. there are two bin borders and division by bin width
   sc.DimensionIsDifferential.clear();
   sc.DimensionIsDifferential.push_back(2);

   sc.CalculateBinSize = true; //< Calculate bin width from lower and upper bin boundaries
   sc.BinSizeFactor = 1; //< Possibility to provide additional normalization factor, e.g. of 2 for bins in |y|
   //sc.BinSize = ; //< If 'CalculateBinSize' is 'false' provide table with bin widths for normalization.

   // Could eventually be extracted from NNLOJET runcard (SCALES)
   sc.ScaleDescriptionScale1 = "scale1"; //< "<pT_1,2>_[GeV]" # This defines the scale to be used (Note: The 1st scale should always be in units of [GeV]!)
   sc.ScaleDescriptionScale2 = "scale2"; //< "pT_max_[GeV]"   # Specify 2nd scale name and unit (ONLY for flexible-scale tables)

   // Binning is read from NNLOJET in warmup run and later defined via warmup files
   //
   // Observable binning Use either 'SingleDifferentialBinning' or
   //    'DoubleDifferentialBinning' or 'TripleDifferentialBinning' in
   //    accordance with 'DifferentialDimension' above
   // sc.SingleDifferentialBinning;
   // sc.DoubleDifferentialBinning; //< Observable binning
   // sc.TripleDifferentialBinning; //< Observable binning

   sc.CenterOfMassEnergy = inputpar_.roots; //< Center-of-mass energy in GeV from NNLOJET.
   sc.PDF1 = 2212; //< PDF of 1st hadron (following PDG convention: proton 2212).
   sc.PDF2 = 2212; //< PDF of 2nd hadron (following PDG convention: proton 2212).

   sc.OutputFilename = "test";//< Filename of fastNLO output table
   sc.OutputPrecision = 8; //< Number of decimal digits to store in output table (def.=8).
   //   sc.OutputCompression = true; // Output compression not set here, since availability of zlib unknown

   // Create table fully flexible in mu_f
   //    larger size, and requires scale independent weights during creation
   //    or table with fixed number of mu_f scale factors, def.=false.
   sc.FlexibleScaleTable = false;

   // Factorization scale variations
   //   only needed for fixed-scale tables,
   //   List of scale factors must include factor '1',
   //   Scale factors will be ordered according to fastNLO convention:
   //   (1, min, ... , max). Defaults: {0.5, 1.0, 2.0}
   sc.ScaleVariationFactors = {1.0};
   sc.ReadBinningFromSteering = true; // Specify if binning is read from fScenConst or from warmup
   sc.IgnoreWarmupBinningCheck = true; // Do not check once more the binning between warmup and steering
   sc.ApplyPDFReweighting  = true; // Apply reweighting of pdfs for an optimized interpolation, def.=true.

   // For warmup run! Set limits for scale nodes to bin borders, if possible
   sc.CheckScaleLimitsAgainstBins = true;

   // Choose default fastNLO interpolation settings
   // To be checked and eventually modified for each observable via additional steering file!
   sc.X_Kernel = "Lagrange";
   sc.X_DistanceMeasure = "sqrtlog10"; //"3rdrtlog10"
   //   sc.X_NNodes = 30; // equivalent to APPLgrid setting; too large tables
   sc.X_NNodes = 20;
   sc.X_NNodeCounting = "NodesPerBin"; // Corresponds to v14 & v21 standard
   if ( nnlo::IsDIS() ) {
      sc.X_Kernel = "Catmull";
      sc.X_DistanceMeasure = "log10"; // we like to have many nodes at LOW-x !
      sc.X_NNodes = 18;
   }
   sc.Mu1_Kernel = "Lagrange";
   sc.Mu1_DistanceMeasure = "loglog025";//"loglog025";
   sc.Mu1_NNodes = 7;
   if ( nnlo::IsDIS() ) {
      sc.Mu1_NNodes = 7; // This is Q2 for DIS
      sc.Mu1_DistanceMeasure = "loglog025";//"loglog025";
   }
   // Do this for a constant scale like M_Z
   // sc.Mu1_Kernel = "OneNode";
   // sc.Mu1_DistanceMeasure = "log10";//"loglog025";
   // sc.Mu1_NNodes = 1;
   sc.Mu2_Kernel = "Lagrange"; //"Lagrange";//<   Lagrange     # Scale2 not used for fixed-scale tables
   sc.Mu2_DistanceMeasure = "loglog025";//"loglog025";//< "loglog025"
   sc.Mu2_NNodes = 7;

   return sc;
}



// _____________________________________________________________________ //
void fnloUtils::SetNNLOJETDefaultBinning(fastNLO::ScenarioConstants& sc, const int& nbins, const double& lo) {
   //!< set binning as used by NNLOJET
   sc.SingleDifferentialBinning.resize(nbins+1);
   for (int i=0; i<=nbins; i++) sc.SingleDifferentialBinning[i] = (&lo)[i];
   std::cout<<"DefaultBinning: "<<lo;
   for (int i=1; i<=nbins; i++) std::cout<<"\t"<<sc.SingleDifferentialBinning[i];
   std::cout<<std::endl;
   return ;
   // double binwidth = (hi-lo)/nbins;
   // sc.SingleDifferentialBinning[0] = lo;
   // std::cout<<"DefaultBinning: "<<lo;
   // for (int i=1; i<=nbins; i++) {
   //    sc.SingleDifferentialBinning[i] = sc.SingleDifferentialBinning[i-1]+binwidth;
   //    std::cout<<"\t"<<sc.SingleDifferentialBinning[i];
   // }
   // std::cout<<std::endl;
}



// // _____________________________________________________________________ //
// std::string fnloUtils::GetProcessName() {
//    //!< get name of the process from NNLOJET
//    std::string process(process_.sproc);
//    process.erase(std::remove(process.begin(), process.end(), ' '),process.end());
//    return process;
// }



// _____________________________________________________________________ //
std::string fnloUtils::GetContribName(int nloops, int ne) {
   //!<
   //!< Get name of contributions
   //!<

 /*
 * LO = #looops=0  #extra-emissions=0
 * V = #looops=1  #extra-emissions=0
 * R = #looops=0  #extra-emissions=1
 * VV = #looops=2  #extra-emissions=0
 * RV = #looops=1  #extra-emissions=1
 * RR = #looops=0  #extra-emissions=2
 */

   if ( nloops==-1 && ne==-1 ) {
      std::cout<<"Works only for having a CURRENT event !!!"<<std::endl;
      nloops = ndimcurrent_.nv;
      static const int& nord   = nnlo::GetOrder();
      // int njet     = inppar_.njets;
      // if ( nnlo::IsDIS() ) njet -= 1;
      ne = nord - nloops;
   }
   if ( nloops == 0 && ne == 0) return "LO";
   else if ( nloops == 1 && ne == 0) return "V";
   else if ( nloops == 0 && ne == 1) return "R";
   else if ( nloops == 2 && ne == 0) return "VV";
   else if ( nloops == 1 && ne == 1) return "RV";
   else if ( nloops == 0 && ne == 2) return "RR";
   else return "Contrib Name NOT IDENTIFIED: "+std::to_string(nnlo::GetOrder());

}



// // _____________________________________________________________________ //
// int fnloUtils::GetOrder() {
//    static const int& nord = ndimcurrent_.norder ;
//    static const int& nv   = ndimcurrent_.nv ;
//    return nord + nv;
// }



// _____________________________________________________________________ //
std::string fnloUtils::GetOrderName( std::set<std::pair<int,int> >* nlne ) {
   //!< get order of calculation as 'string'

   if ( nlne == NULL) {
      int nord   = nnlo::GetOrder();
      //static const int& nord   = ndimcurrent_.norder;
      if ( nord == 0 ) return "LO";
      else if ( nord == 1 ) return "NLO";
      else if ( nord == 2 ) return "NNLO";
   }
   else if ( fnloUtils::fIDs.size()==1 ) {
      // only one channel was calculated
      // take:  'chZYZ'
      char buf[5];
      sprintf(buf, "%03d", *fnloUtils::fIDs.begin() );
      return "ch"+std::string(buf);
   }
   else {
      std::set<std::string> names;
      for ( auto i : *nlne ) names.insert(fnloUtils::GetContribName(i.first,i.second));
      if ( names.count("R") && names.count("V") ) {
         names.erase("R");
         names.erase("V");
         names.insert("NLO");
      }
      if ( names.count("RR") && names.count("VV") && names.count("RV") ) {
         names.erase("RR");
         names.erase("RV");
         names.erase("VV");
         names.insert("NNLO");
      }
      std::string ret;
      for ( auto i : names )
         if ( ret.empty() ) ret = i;
         else ret += "+"+i;
      return ret;
   }
   return ""; //?
}



// _____________________________________________________________________ //
int fnloUtils::GetSeed() {
   //!< get seed value

   // new
   return iseed_hook; // nnlo_common.

   // old
   /*
   std::vector<std::string> tok;
   char clog[101];
   strcpy(clog,log_.slogname);
   char *p = strtok(clog, ".");
   while (p) {
      tok.push_back(p);
      p = strtok(NULL, ".");
   }
   if ( tok.size()<=2 ) return -1;
   std::string sN = tok[tok.size()-2];
   sN.erase(0,1);
   int ret = atoi( sN.c_str());
   return ret;
   */
}



// _____________________________________________________________________ //
std::string fnloUtils::GetConfigDir() {
   //!< Get directory where config files are stored
   // CONFIGDIR is set by NNLOJET at 'make' time! Later changes per environment variable not possible!
   std::string configdir(getenv("CONFIGDIR")); // Read from NNLOJET setting
   std::cout << "[fnloUtils::GetConfigDir()] CONFIGDIR of NNLOJET found to be: " << configdir << std::endl;
   if ( configdir == "" ) {
      std::cout << "[fnloUtils::GetConfigDir()] WARNING! CONFIGDIR of NNLOJET not set, "
                << "using default $PWD/process instead." << std::endl;
      configdir = std::string(getenv("PWD"))+"/process";
   }
   // --- test config directory
   struct stat myStat;
   std::cout << "[fnloUtils::GetConfigDir()] Status of CONFIGDIR found to be: " << stat(configdir.c_str(), &myStat) << std::endl;
   if ( !((stat(configdir.c_str(), &myStat) == 0) && (((myStat.st_mode) & S_IFMT) == S_IFDIR)) ) {
      std::cout << "[fnloUtils::GetConfigDir()] WARNING! Could not find the config directory " <<
         configdir << std::endl;
      // try default in case CONFIGDIR was set differently above
      configdir = std::string(getenv("PWD"))+"/process";
      std::cout << "[fnloUtils::GetConfigDir()] Trying default setting for CONFIGDIR of NNLOJET: " << configdir << std::endl;
      if ( !((stat(configdir.c_str(), &myStat) == 0) && (((myStat.st_mode) & S_IFMT) == S_IFDIR)) ) {
         std::cout << "[fnloUtils::GetConfigDir()] WARNING! Could not find the config directory " <<
            configdir << std::endl;
         // another and next to last guess (Grid site)
         configdir = "/cvmfs/etp.kit.edu/nnlo/src/driver/process";
         std::cout << "[fnloUtils::GetConfigDir()] Trying grid site setting for CONFIGDIR: " << configdir << std::endl;
         if ( !((stat(configdir.c_str(), &myStat) == 0) && (((myStat.st_mode) & S_IFMT) == S_IFDIR)) ) {
            std::cout << "[fnloUtils::GetConfigDir()] WARNING! Could not find the config directory " <<
               configdir << std::endl;
            // really last guess (local user install); procedure to be rethought ...
            configdir = std::string(getenv("HOME"))+"/local/src/NNLOJET-rev3678/driver/process";
            std::cout << "[fnloUtils::GetConfigDir()] Trying local install setting for CONFIGDIR: " << configdir << std::endl;
            if ( !((stat(configdir.c_str(), &myStat) == 0) && (((myStat.st_mode) & S_IFMT) == S_IFDIR)) ) {
               std::cerr << "[fnloUtils::GetConfigDir()] ERROR! Really could not find a config directory " <<
                  configdir << std::endl;
               exit(2);
            }
         }
      }
   }
   std::cout<<"[fnloUtils::GetConfigDir()] Using config directory: "<<configdir<<std::endl;
   return configdir;
}



// _____________________________________________________________________ //
std::string fnloUtils::GetConfigFile() {
   //!< get name of the process
   std::string pname = nnlo::GetProcessName();
   std::cout << "[fnloUtils::GetConfigFile()] The process name is reported to be: " << pname << std::endl;
   if ( pname == "1jet" || pname == "2jet" ) {
     pname = pname.substr(1,std::string::npos);
     std::cout << "[fnloUtils::GetConfigFile()] For subprocess configuration shortening process name to: " << pname << std::endl;
   }
   std::string fname = pname + "-" + fnloUtils::GetOrderName() + ".config";
   fname = fnloUtils::GetConfigDir() + "/" + pname + "/" + fname;
   std::cout << "[fnloUtils::GetConfigFile()] Using config file: " << fname << std::endl;
   return fname;
}



// _____________________________________________________________________ //
std::vector<std::string> fnloUtils::GetConfigFiles() {
   //!< get name of the process
   std::string pname = nnlo::GetProcessName();
   std::cout << "[fnloUtils::GetConfigFiles()] The process name is reported to be: " << pname << std::endl;
   if ( pname == "1jet" || pname == "2jet" ) {
     pname = pname.substr(1,std::string::npos);
     std::cout << "[fnloUtils::GetConfigFiles()]  For subprocess configuration shortening process name to: " << pname << std::endl;
   }
   std::vector<std::string> fnames;
   std::string pref = fnloUtils::GetConfigDir()+"/"+pname+"/"+pname;
   if ( nnlo::GetOrder() == 0 ) {
      fnames.push_back( pref+"-LO.config");
   }
   else if ( nnlo::GetOrder() == 1 ){
      fnames.push_back( pref+"-R.config");
      fnames.push_back( pref+"-V.config");
   }
   else if ( nnlo::GetOrder() == 2 ){
      fnames.push_back( pref+"-RR.config");
      fnames.push_back( pref+"-RV.config");
      fnames.push_back( pref+"-VV.config");
   }
   std::cout<<"[fnloUtils::GetConfigFiles()]  Found order "<<nnlo::GetOrder()<<". Using config files: "<<std::endl;
   for ( auto f : fnames ) std::cout<<"\t"<<f<<std::endl;
   return fnames;
}



// _____________________________________________________________________ //
void fnloUtils::RememberContribution() {
   //!< fill the variable fNlNe
   //!< and remember the calculated contributions

   static const int& iproc = currentprocess_.iproc;
   int nloops = channel_.iloops[iproc-1];
   int nreal  = ndimcurrent_.norder;

   // static const int& nloops = ndimcurrent_.nv;
   // int nord                 = nnlo::GetOrder();//ndimcurrent_.norder;
   // // int njet                 = inppar_.njets;
   // // if ( nnlo::IsDIS() ) njet -= 1;
   // // int ne = nord - nloops;
   // //std::cout<<"nl="<<nloops<<"\tne="<<ne<<"\t\tnord="<<nord<<"\tIsDIS="<<nnlo::IsDIS()<<"\tnjet="<<njet<<std::endl;
   fnloUtils::fNlNe.insert({ nloops,nreal });
   fnloUtils::fIDs.insert(iproc);

}



// _____________________________________________________________________ //
double fnloUtils::CalculateReferenceWeight(double wgt, double x1, double x2) {
   //!< calculate the reference weight
   //!< and compare it later with the ref-weight from NNLOJET
   //!<
   //!< Set value: fnloUtils::_myfillwgt
   //!< Need access to evolvepdf_()
   //!< needs either strong_.as or alphaspdf_()
   //!<

   double PDF = 0;
   double muf2  = muf2_hook(1);
   const std::vector<std::vector<int > >&  pc =  fnloUtils::sp.GetInitialValues();
   if ( nnlo::IsDIS() ) {
      double xfx[13] = {0};
      evolvepdf_(x1,sqrt(muf2),xfx);
      for ( unsigned int ip = 0 ; ip < pc[currentprocess_.iproc].size() ; ip+=2 ) {
            PDF += xfx[pc[currentprocess_.iproc][ip]+6];
      }
   }
   else {
      double xfx1[13] = {0};
      double xfx2[13] = {0};
      evolvepdf_(x1,sqrt(muf2),xfx1);
      evolvepdf_(x2,sqrt(muf2),xfx2);
      for ( unsigned int ip = 0 ; ip < pc[currentprocess_.iproc].size() ; ip+=2 ) {
         PDF += xfx1[pc[currentprocess_.iproc][ip]+6]*xfx2[pc[currentprocess_.iproc][ip+1]+6];
      }
   }

   static const double& as = strong_.as; //alphaspdf_(sqrt(mur2));
   fnloUtils::_myfillwgt = wgt * PDF * pow( as/(2*M_PI) , nnlo::GetNPow() ) ;
   return fnloUtils::_myfillwgt;
}



// _____________________________________________________________________ //
std::vector<double> fnloUtils::SclIndepWgt(const double* wt, double muf, double mur) {
   //! Calculate scale independent weights for fastNLO

   // Do not use debug mode for mass production!
   static const bool Debug = false;
   static bool First = true;
   if ( First && Debug ) {
      First = false;
      std::cout << "[fnloUtils::SclIndepWgt]: Running in DEBUG mode. Do not use for mass production!" << std::endl;
   }

   std::vector<double> ret{0,0,0,0,0,0};  // return value
   std::vector<double> sret{0,0,0,0,0,0}; // simple sum
   std::vector<double> oret{0,0,0,0,0,0}; // ordered sum
   std::vector<double> kret{0,0,0,0,0,0}; // Kahan sum

   if ( nnlo::GetOrder()==0 ) {
      //! LO does not have scale dependent terms
      sret[0] = wt[0];
      oret[0] = wt[0];
      kret[0] = wt[0];
      ret[0]  = kret[0];
      return ret;
   }

   double wamax = DBL_MAX;
   const double EPS   = 1.e-3;  // aim for permille precision
   const double SMALL = 1.e-12; // relative weight limit to consider some deviations > EPS negligible
   static const unsigned int nSclVNNLO = 6;
   const unsigned int nSclV = nSclVNNLO;
   // Simplify in case of NLO, to be checked!!
   // unsigned int nSclV = nnlo::GetOrder()==1  ? 3 : 6; // only NLO

   if ( nSclV==3 ) {
      std::cerr << "[fnloUtils::SclIndepWgt]: Do not use yet NLO with nSclV = " << nSclV << "! Aborted!" << std::endl;
      exit(1);
      // -- not working, to be checked!
      // static const double SclInvNLO[3][3] = {{3, -2, 0},{0, 1, -1},{-1, 0, 1}}; // i+1
      // for ( unsigned int j = 0 ; j<nSclV ; j++ ){
      //    for ( unsigned int i = 0 ; i<nSclV ; i++ ){
      //       ret[j] += wt[i+1]*SclInvNLO[j][i];
      //    }
      //    //if ( fabs(ret[j]/wt[1]) < EPS ) ret[j] = 0;
      // }
      // -- not working, to be checked!
      // const double SclInvNLO[3][3] = {{5.5,-2,-2.5},{-1.5,1,0.5},{-0.5,0,0.5}}; // i+3
      // for ( unsigned int j = 0 ; j<nSclV ; j++ ){
      //         for ( unsigned int i = 0 ; i<nSclV ; i++ ){
      //            ret[j] += wt[i+3]*SclInvNLO[j][i];
      //         }
      // }
   } else {
      // This decomposition matrix is based on i.a. fixed scale choices of
      // exp^1, exp^(3/2), exp^(5/2);
      static const double SclInvNNLO[nSclVNNLO][nSclVNNLO] = {
         {    13.,  4.,    -9.,    -9.,     1., 1.},
         {-14./3., -2.,     2., 11./2.,     0., -5./6.},
         {-14./3., -2., 11./2.,     2., -5./6., 0.},
         {  1./3.,  0.,     0., -1./2.,     0., 1./6.},
         {  1./3.,  0., -1./2.,     0.,  1./6., 0.},
         {     1.,  1.,    -1.,    -1.,     0., 0.} };
      double sum[nSclVNNLO] = {0., 0., 0., 0., 0., 0.};
      int nRet = ( nnlo::GetOrder()==1 ) ? 4 : nSclV;

      // static unsigned int nprnt = 0;
      // static unsigned int npmod = 1;
      // if ( nprnt < 10 && nprnt % npmod == 0 ) {
      //    nprnt++;
      //    for ( int i = 0 ; i<=nRet ; i++ ){
      //       std::cout << std::setprecision(16) << "AAA i = " << i << ", wt[i] = " << wt[i] << std::endl;
      //    }
      // } else if ( nprnt == 10 ) {
      //    nprnt = 0;
      //    npmod *= 10;
      // }

      wamax = std::abs(wt[0]);
      for ( int j = 0 ; j<nRet ; j++ ){
         wamax = std::max(std::abs(wt[j+1]),wamax);
         double wasum = 0.;
         for ( unsigned int i = 0 ; i<nSclV ; i++ ){ // DB
            //         for ( unsigned int i = 0 ; i<nRet ; i++ ){ // KR TEST
            sum[i]   = wt[i+1]*SclInvNNLO[j][i];
            //            sret[j] += sum[i];
            wasum   += std::abs(sum[i]);
         }
         //         oret[j]  = fnloUtils::OrderedSum(sum,nSclVNNLO);
         kret[j]  = fnloUtils::KahanSum(sum,nSclVNNLO);
         // std::cout << "RD1 j = " << j << ", normal sum[j]  = " << sret[j] << std::endl;
         // std::cout << "RD2 j = " << j << ", ordered sum[j] = " << oret[j] << ", o/s = " << oret[j]/sret[j] << std::endl;
         // std::cout << "RD3 j = " << j << ", Kahan sum[j]   = " << kret[j] << ", k/s = " << kret[j]/sret[j] <<std::endl;
         // Numerically Kahan summing is the most precise
         ret[j] = kret[j];

         if ( Debug ) {
            // Set small weights to zero? Not necessary. Does not occur anyway for DBL_MIN ...
            static unsigned int nzwgt = 0;
            static unsigned int nzmss = 0;
            static unsigned int nzmod = 1;
            if ( ret[j] != 0. && std::abs(ret[j])/wasum < DBL_MIN ) {
               nzwgt++;
               if ( nzmss < 10 && nzwgt % nzmod == 0 ) {
                  nzmss++;
                  std::cout << "[fnloUtils::SclIndepWgt]: " << nzwgt << "th occurrence of tiny scale weight = " << ret[j] << std::endl;
               } else if ( nzmss == 10 ) {
                  nzmss = 0;
                  nzmod *= 10;
               }
            }
         }
      }

      // if ( nnlo::GetOrder()==1 ) {
      //    // ret[3]=0;
      //    ret[4]=0;// must be exact zero for fastNLO
      //    ret[5]=0;// must be exact zero for fastNLO
      // }

   }

   // Try to reconstruct original weight wt[0] as cross check
   if ( Debug ) {
      double r = log(mur*mur);
      double f = log(muf*muf);
      double sum[6] = { ret[0], r*ret[1], f*ret[2], r*r*ret[3], f*f*ret[4], r*f*ret[5] };
      // double sMyWgt0 = 0.;
      // for ( unsigned int i=0; i<nSclVNNLO; i++) {
      //    sMyWgt0 += sum[i];
      // }
      // double oMyWgt0 = OrderedSum(sum,nSclVNNLO);
      double kMyWgt0 = KahanSum(sum,nSclVNNLO);
      // Numerically Kahan summing is the most precise
      double MyWgt0  = kMyWgt0;

      //double MyWgtNLO   = ret[0] + r*ret[1] + f*ret[2];
      //double MyRefWgt = ret[0] + 2*r*ret[1] + 2*f*ret[2] + 4*r*r*ret[3] + 4*f*f*ret[4] + 4*r*f*ret[5];

      // Count misreconstructed weights for first scale
      static unsigned int nmwgt = 0;
      static unsigned int nmmss = 0;
      static unsigned int nmmod = 1;
      if  ( std::abs(MyWgt0/wt[0]-1) > EPS || !std::isfinite(MyWgt0)) {
         if ( std::abs(MyWgt0)/wamax > SMALL || std::abs(wt[0])/wamax > SMALL ) {
            nmwgt++;
            if ( nmmss < 10 && nmwgt % nmmod == 0 ) {
               nmmss++;
               std::cout << "[fnloUtils::SclIndepWgt]: " << nmwgt << "th occurrence of large misreconstructed weight!" << std::endl;
               std::cout<<"[fnloUtils::SclIndepWgt]: LARGE discrepant weight! Original weight = " << wt[0] << ", reconstructed weight = " << MyWgt0 << ", deviation (per mille) = " << 1000*(std::abs(MyWgt0/wt[0])-1) << ", weight max = " << wamax << std::endl;
               std::cout<<"[fnloUtils::SclIndepWgt]: Original weight/max = " << std::abs(wt[0])/wamax << ", reconstructed weight/max = " << std::abs(MyWgt0)/wamax << std::endl;
            } else if ( nmmss == 10 ) {
               nmmss = 0;
               nmmod *= 10;
            }
         }
         std::cerr<<"[fnloUtils::SclIndepWgt]. Warning!"<<std::endl;
         std::cerr<<"[fnloUtils::SclIndepWgt]. Discrepant weight! Original weight = " << wt[0] << ", reconstructed weight = " << MyWgt0 << ", weight max = " << wamax << std::endl;
         std::cerr<<"[fnloUtils::SclIndepWgt]. Original weight/max = " << std::abs(wt[0])/wamax << ", reconstructed weight/max = " << std::abs(MyWgt0)/wamax << std::endl;
         std::cout<<"  muf = " << muf << ", mur = " << mur << std::endl;
         std::cout<<"  wt [] =";
         for ( unsigned int i = 0 ; i<7 ; i++ )  std::cout<<"\t"<<wt[i];
         std::cout<<"\t\t"<<std::endl;

         std::cout<<"  ln(mur^2) []";
         for ( unsigned int i = 0 ; i<7 ; i++ )  std::cout<<"\t"<<log(mur2_hook(i+1));
         std::cout<<"\t\t"<<std::endl;

         std::cout<<"  ln(muf^2) []";
         for ( unsigned int i = 0 ; i<7 ; i++ )  std::cout<<"\t"<<log(muf2_hook(i+1));
         std::cout<<"\t\t"<<std::endl;

         std::cout<<"  ret[]";
         for ( unsigned int j = 0 ; j<6 ; j++ )std::cout<<"\t"<<ret[j];
         std::cout<<std::endl;
         std::cout<<" oret[]";
         for ( unsigned int j = 0 ; j<6 ; j++ )std::cout<<"\t"<<oret[j];
         std::cout<<std::endl;
         std::cout<<" kret[]";
         for ( unsigned int j = 0 ; j<6 ; j++ )std::cout<<"\t"<<kret[j];
         std::cout<<std::endl;

         //         exit(3);
      }
   }

   return ret;
}


// Kahan summing to improve numerical precision
double fnloUtils::KahanSum( double in[], unsigned int insize ) {
   double ksum = in[0];
   double c = 0.;
   for ( unsigned int i = 1; i<insize ; i++ ){
      double y = in[i] - c;
      double t = ksum + y;
      c        = (t  - ksum)  - y;
      ksum     = t;
   }
   return ksum;
}


// Ordered summing (in descending order here) to improve numerical precision
double fnloUtils::OrderedSum( double in[], unsigned int insize ) {
   double osum = 0.;
   std::vector<double> vin;
   for ( unsigned int i = 0 ; i<insize ; i++ ){
      vin.push_back(in[i]);
   }
   auto func=[](double a, double b) { return std::abs(a) > std::abs(b); };
   std::sort(std::begin(vin), std::end(vin), func);
   for ( unsigned int i = 0 ; i<insize ; i++ ){
      osum += vin[i];
   }
   return osum;
}



// _____________________________________________________________________ //
void fnloUtils::InitFastNLO( const int& id, const double* wt ) {
   // static const int& iproc = currentprocess_.iproc;
   // int nloops = channel_.iloops[iproc-1]; /// loop order
   // int ipars  = channel_.ipars[iproc-1];  /// number of additional particles + 1
   // std::cout<<"        nloops          "<<nloops<< std::endl;
   // std::cout<<"        ipars           "<<ipars<< std::endl;
   // std::cout<<"        ndim.norder     "<<ndimcurrent_.norder<< std::endl;
   // std::cout<<"        inppar.njets    "<<inppar_.njets<< std::endl;
   // std::cout<<"        order_.ieorder  "<<order_.ieorder<< std::endl;
   // std::cout<<"        GetNPos         "<<nnlo::GetNPow()<< std::endl;
   // std::cout<<"        GetOrder()      "<<nnlo::GetOrder()<< std::endl;
   // std::cout<< std::endl;

   // --- welcome
   std::cout << "[fnloUtils::InitFastNLO] fastNLO initialisation!"  << std::endl;
   if ( fnloUtils::ftable[id] != NULL )  return;
   std::cout << "[fnloUtils::InitFastNLO] NNLOJET provides process: " << nnlo::GetProcessName() << std::endl;

   // --- get settings of histogram [id] from NNLOJET
   const std::string& gridname = fnloUtils::fInits[id].gridname;
   const int& nbins = fnloUtils::fInits[id].nbins;
   const double& lo = *(fnloUtils::fInits[id].lo);
   //const double& hi = *(fnloUtils::fInits[id].hi);
   fnloUtils::_myfillwgt  = 0;

   // --- read config file with  parton combinations
   if ( fnloUtils::sp.GetConfigFile().empty() )
     fnloUtils::sp.ReadConfigFiles(fnloUtils::GetConfigFiles());

   // --- get NNLOJET specific 'GeneratorConstants'
   fastNLO::GeneratorConstants gc = fnloUtils::GetNNLOJET_GenConsts();

   // --- build process specific 'ProcessConstants' (DIS or pp?)
   fastNLO::ProcessConstants pc = nnlo::IsDIS() ? //fnloUtils::GetProcess() == "DIS" ) ?
     fnloUtils::GetNNLOJET_ProcConstsDIS() :
     fnloUtils::GetNNLOJET_ProcConstsPP() ;

   // --- build scenario specific 'ScenarioConstants'
   fastNLO::ScenarioConstants sc = fnloUtils::GetNNLOJET_ScenConsts();
   if ( gridname.find("iff") != std::string::npos ) { // TODO: Fix this very ugly hack
     sc.X_NNodes = 50;
     sc.X_DistanceMeasure = "sqrtlog10"; // we like to have increasingly many nodes at high-x
   }
   if ( nnlo::IsDIS() && gridname.find("q2") == std::string::npos ) {
     std::cout<<"[fnloUtils::InitFastNLO] IsDIS, but not differential as function of q2, using more q2 nodes, less x-nodes"<<std::endl;
     sc.Mu1_NNodes = 14; // more q2 nodes
     sc.Mu2_NNodes = 6;
     if ( gridname.find("iff") != std::string::npos ) sc.X_NNodes = 42;
   }
   // TODO: Fix hack to select between flexible- and fixed-scale table
   std::cout<<"[fnloUtils::InitFastNLO] NNLOJET scale choices: muf2(1) = " << muf2_hook(1) << "\tmur2(1) = " << mur2_hook(1) << std::endl;
   std::cout<<"[fnloUtils::InitFastNLO] NNLOJET scale choices: muf2(2) = " << muf2_hook(2) << "\tmur2(2) = " << mur2_hook(2) << std::endl;
   if ( muf2_hook(2) == mur2_hook(2) && fabs(sqrt(mur2_hook(2))/exp(1)-1)<1.e-6) {
      std::cout<<"[fnloUtils::InitFastNLO] Flexible-scale table chosen."<< std::endl;
      sc.FlexibleScaleTable = true;
   } else if ( muf2_hook(1) == mur2_hook(1) ) {
      std::cout<<"[fnloUtils::InitFastNLO] Fixed-scale table chosen."<< std::endl;
      sc.FlexibleScaleTable = false;
   } else {
      std::cout<<"[fnloUtils::InitFastNLO] Unrecognized scale settings. Please choose either a fixed-scale or a flexible-scale table."<< std::endl;
      std::cout<<" Exiting."<< std::endl;
      exit(4);
   }
   fnloUtils::SetNNLOJETDefaultBinning(sc, nbins, lo);
   // if ( nnlo::IsDIS() ){
   //    sc.CheckScaleLimitsAgainstBins = sc.FlexibleScaleTable;
   // }


   // TODO: The scale setting interface is extremely intransparent ==> IMPROVE
   // --- adapt for muf-variations for fixed-scale tables
   if ( sc.FlexibleScaleTable ) {
      static const double sclfac[1] = {1.0};
   } else {
      static const double sclfac[maxscl] = {1,0.5,2.,0.505, 202,0.005,200}; // ugly convention: muf*0.5, muf*2, mur+muf * 0.5, mur+muf* 2, mur*0.5 mur*2
   }
   std::cout<<"[fnloUtils::InitFastNLO] wt  []:";
   for ( int is = 0 ; is<maxscl ; is++ ) printf("\t%5.4e",wt[is]);
   std::cout<< std::endl;
   std::cout<<"[fnloUtils::InitFastNLO] muf2[]:";
   for ( int is = 0 ; is<maxscl ; is++ ) printf("\t%5.4e",muf2_hook(is+1));
   std::cout<< std::endl;
   std::cout<<"[fnloUtils::InitFastNLO] mur2[]:";
   for ( int is = 0 ; is<maxscl ; is++ ) printf("\t%5.4e",mur2_hook(is+1));
   std::cout<< std::endl;

   //if ( !getunitphase_() ) {
   // KR: Use C++ macro constants like DBL_MIN instead
   //      for ( int is = 1 ; is<maxscl && fabs(wt[is])> TINY; is++ )  {
      for ( int is = 1 ; is<maxscl && fabs(wt[is]) > DBL_MIN; is++ )  {
   //         std::cout<<"[fnloUtils::InitFastNLO] wt["<<is<<"]="<<wt[is]<< std::endl;
         sc.ScaleVariationFactors.push_back(sclfac[is]);
      }
      std::cout<<"[fnloUtils::InitFastNLO] Number of scale factors found: "<<sc.ScaleVariationFactors.size()<<"\t Please respect the convention how to specify those in the run-card."<< std::endl;
      std::cout<<"[fnloUtils::InitFastNLO] Scale variation factors: "<<sc.ScaleVariationFactors<< std::endl;

      if ( sc.FlexibleScaleTable ) {
         // double q2 = muf2_hook(1);
         // double p2 = mur2_hook(1);
         // std::cout<<"mur: ";
         // for ( int is = 0 ; is<maxscl ; is++ ) std::cout<<"\t"<<sqrt(mur2_hook(1+is)) ;
         // std::cout<<std::endl;
         // std::cout<<"muf: ";
         // for ( int is = 0 ; is<maxscl ; is++ ) std::cout<<"\t"<<sqrt(muf2_hook(1+is)) ;
         // std::cout<<std::endl;

         // double muf2[7];
         // double mur2[7];
         // for ( int is = 0 ; is<maxscl ; is++ ) {
         //    muf2[is] = muf2_hook(is+1);
         //    mur2[is] = mur2_hook(is+1);
         // }
         // std::cout<<fabs( mur2[1] / q2 - 1. )<<std::endl;
         // std::cout<<fabs( muf2[2] / p2 - 1. )<<std::endl;
         // std::cout<<fabs( muf2[3] / p2 - 1. )<<std::endl;
         // std::cout<<fabs( mur2[4] / (p2*p2) - 1. )<<std::endl;
         // std::cout<<fabs( muf2[5] / (q2*q2) - 1. )<<std::endl;
         // std::cout<<fabs( muf2[6] / (q2*q2) - 1. )<<std::endl;
         // const double EPS = 1.e-10;
         // if ( fabs( muf2[1] / q2 - 1 ) > EPS || fabs( mur2[1] / q2 - 1 ) > EPS) {
         //    std::cerr<<"[nnlo::fill_fastnlo]. Error. This is a flexible scale table, but the definition of the scales does not agree with the convention [1]"<< std::endl;
         //    exit(3);
         // }
         // if ( fabs( muf2[2] / p2 - 1 ) > EPS || fabs( mur2[2] / p2 - 1 ) > EPS ) {
         //    std::cerr<<"[nnlo::fill_fastnlo]. Error. This is a flexible scale table, but the definition of the scales does not agree with the convention [2]"<< std::endl;
         //    exit(3);
         // }
         // if ( fabs( muf2[3] / p2 - 1 ) > EPS || fabs( mur2[3] / q2 - 1 ) > EPS) {
         //    std::cerr<<"[nnlo::fill_fastnlo]. Error. This is a flexible scale table, but the definition of the scales does not agree with the convention [3]"<< std::endl;
         //    exit(3);
         // }
         // if ( fabs( muf2[4] / q2 - 1 ) > EPS || fabs( mur2[4] / (p2*p2) - 1 ) > EPS) {
         //    std::cerr<<"[nnlo::fill_fastnlo]. Error. This is a flexible scale table, but the definition of the scales does not agree with the convention [4]"<< std::endl;
         //    exit(3);
         // }
         // if ( fabs( muf2[5] / (q2*q2) - 1 ) > EPS || fabs( mur2[5] / p2 - 1 ) > EPS ) {
         //    std::cerr<<"[nnlo::fill_fastnlo]. Error. This is a flexible scale table, but the definition of the scales does not agree with the convention [5]"<< std::endl;
         //    exit(3);
         // }
         // // TEST
         // if ( fabs( muf2[6] / (q2*q2) - 1 ) > EPS || fabs( mur2[6] / (p2*p2) - 1 ) > EPS ) {
         //    std::cerr<<"[nnlo::fill_fastnlo]. Error. This is a flexible scale table, but the definition of the scales does not agree with the convention [6]"<< std::endl;
         //    exit(3);
         // }
      }

      else if ( !sc.FlexibleScaleTable ) {
         // TODO Improve. Only 0,1,2 of NNLOJET runcard SCALES are accessed.
         // SCALES must be defined to fulfil the following conditions!
         for ( int is = 0 ; is< 3 && is<(int)sc.ScaleVariationFactors.size() ; is++ ){
            // mur(1) is base scale
            double mur0  = mur2_hook(1);
            //double muf0  = muf2_hook(1);
            // This loops over 1,2,3
            double mur2  = mur2_hook(is+1);
            double muf2  = muf2_hook(is+1);
            if ( is==0 && mur2!=muf2 ) { // mur==muf for first "scale variation"
               std::cerr<<"[nnlo::fill_fastnlo]. Error. First scale setting: muf and mur must be identical, but muf2/mur2="<<muf2/mur2<< std::endl;
               exit(2);
            }
            // else if ( is==1 && fabs(muf2/mur2/0.25-1) > 1.e-6 ) {
            //    std::cerr<<"[nnlo::fill_fastnlo]. Error. Second scale variation must be muf*0.5=mur*0.5, but muf2/mur2="<<muf2/mur2<< std::endl;
            //    exit(2);
            // }
            // else if ( is==2 && fabs(muf2/mur2/4-1) > 1.e-6 ) {
            //    std::cerr<<"[nnlo::fill_fastnlo]. Error. Third scale variation must be muf*2=mur*2, but muf2/mur2="<<muf2/mur2<< std::endl;
            //    exit(2);
            // }
            if ( fabs(muf2/mur2-1) > 1.e-6 ) { // mur/muf = 1 for all "scale variations" ==> symmetric variation of muf, mur versus base scale mur(1)
               std::cerr<<"[nnlo::fill_fastnlo]. Error. Scale variation is="<<is<<" must be muf*0.5=mur*0.5, but muf2/mur2="<<muf2/mur2<< std::endl;
               exit(2);
            }
            if ( is==1 && fabs(mur2/mur0/0.25-1) > 1.e-6 ) { // mur/mu0 = 2 && muf/mu0 = 2
               std::cerr<<"[nnlo::fill_fastnlo]. Error. Second scale variation must be muf*0.5=mur*0.5, but mur2/mur2(0)="<<mur2/mur0<< std::endl;
               exit(2);
            }
            if ( is==2 && fabs(mur2/mur0/4.-1) > 1.e-6 ) { // mur/mu0 = 1/2 && muf/mu0 = 1/2
               std::cerr<<"[nnlo::fill_fastnlo]. Error. Second scale variation must be muf*2=mur*2, but mur2/mur2(0)="<<mur2/mur0<< std::endl;
               exit(2);
            }
         }
         //         if ( fabs(wt[3])>TINY ) {
         // KR: Use C++ macro constants like DBL_MIN instead
         if ( fabs(wt[3])>DBL_MIN ) {
            std::cerr<<"[nnlo::fill_fastnlo]. Error. fastNLO only support scale settings of [0] mur=muf; [1] muf*0.5=mur*0.5; [2] muf*2=muf*2"<< std::endl; // Yes.
            exit(2);
         }
      }
//}

   // ---- initialize fastNLOCreate
   // --- Collect info for naming
   std::string wrmfile;
   std::string steerfile;
   std::string tablename;
   std::string runid = runid_.srunid;
   std::string pname = nnlo::GetProcessName();
   std::string jname = jobname_.sjobname;
   jname.erase(remove(jname.begin(), jname.end(), ' '),jname.end());
   std::string fname = outfile_.fname;
   fname.erase(remove(fname.begin(), fname.end(), ' '),fname.end());
   int seed = fnloUtils::GetSeed();
   std::cout << "[fnloUtils::InitFastNLO] Gridname: " << gridname << std::endl;
   std::cout << "[fnloUtils::InitFastNLO] Integration input: " << fname << std::endl;
   std::cout << "[fnloUtils::InitFastNLO] Process name: " << pname << ", job name: " << jname << std::endl;
   // If unitphase, i.e. grid warmup, then name warmfile individually for each subprocess for later check & merge
   if ( getunitphase_() ) {
     wrmfile = pname+"."+jname+"."+gridname;
   }
   // If NOT unitphase, i.e. production, then expect identical warmfile for all subprocesses
   else {
     wrmfile = pname+"."+gridname;
   }
   size_t ipos = gridname.find_first_of("_");
   std::string scenname = gridname.substr(0,ipos);
   steerfile = pname+"."+scenname;
   tablename = pname+"."+jname+"."+runid+"."+gridname+".s"+std::to_string(seed);
   wrmfile   += ".wrm";
   steerfile += ".str";
   std::cout << "[fnloUtils::InitFastNLO] Warmup filename: " << wrmfile << std::endl;
   std::cout << "[fnloUtils::InitFastNLO] Steering filename: " << steerfile << std::endl;
   std::cout << "[fnloUtils::InitFastNLO] Tablename: " << tablename << std::endl;

   sc.OutputFilename = tablename;
   fnloUtils::ftablename[id] = tablename;
   fnloUtils::ftable[id] = new fastNLOCreate(gc,pc,sc,wrmfile,steerfile);
   fnloUtils::fIsInclJet[id] = (gridname.find("ji") != std::string::npos);
   fnloUtils::fIs820GeV[id] = (gridname.find("820") != std::string::npos);
   fnloUtils::fIs137[id]    = (gridname.find("137") != std::string::npos);
   std::cout << "[fnloUtils::InitFastNLO] Power of alpha_s(): " << nnlo::GetNPow() << ", order of calculation: " << nnlo::GetOrder() << std::endl;
   fnloUtils::ftable[id]->SetOrderOfAlphasOfCalculation(nnlo::GetNPow());
   fnloUtils::ftable[id]->SetCacheSize(30);
   //fnloUtils::ftable[id]->fastNLOTable::SetNcontrib(1); // to be removed later!

   if ( fnloUtils::ftable[id]->GetIsWarmup() && !getunitphase_() ) {
      std::cerr<<"[nnlo::init_fastnlo] ERROR! No warmup-file found (hence this cannot be a production run), but the job is NOT running with UNIT_PHASE. Aborted!"<< std::endl;
      exit(5);
   } else if ( !fnloUtils::ftable[id]->GetIsWarmup() && getunitphase_() ) {
      std::cerr<<"[nnlo::init_fastnlo] ERROR! Warmup-file found (hence this is a production run), but the job is running WITH UNIT_PHASE. Aborted!"<< std::endl;
      exit(6);
   }

   // Read fastNLO steering file???
   //    const std::string str = config.subprocConfig.fileName;
   //    const std::string steeringNameSpace = phasespaceFilePath();
   //    READ_NS(str, steeringNameSpace);

   std::cout << "[fnloUtils::InitFastNLO] fastNLO initialized. " << std::endl;
   return ;
}



// _____________________________________________________________________ //
namespace fnloUtils {

   // _____________________________________________________________________ //
   fnloSuprocesses::fnloSuprocesses(const std::string& config){
      //!< constructor reading config file with subprocesses
      ReadConfigFile(config);
   }

   // _____________________________________________________________________ //
  void fnloSuprocesses::ReadConfigFile(const std::string& config) {
      fConfigFile = config;
      this->resize(1000);
      fInitialValues.resize(1000);
      std::ifstream cf(config);
      if ( !cf.is_open() ){
         std::cerr<<"[fnloSubprocesses::ReadConfigFile()] Error. Could not open file: "<<config<< std::endl;
         exit(2);
      }
      int id, np, v1, v2;
      int nl = 0;
      cf>>id; // first line is '0'
      while ( !cf.eof() ) {
         cf>>id>>np;
         if(cf.fail()) break;
         nl++;
         //cout<<"nl: "<<nl<<"\tid="<<id<<"\tnp="<<np<< std::endl;
         for ( int i = 0 ; i<np ; i++ ) {
            cf>>v1>>v2;
            //cout<<"\t"<<v1<<"\t"<<v2<< std::endl;
            fInitialValues[id].push_back(v1);
            fInitialValues[id].push_back(v2);
         }
      }
      std::cout<<"[fnloSuprocesses::ReadConfigFile()] found "<<nl<<" subprocess definitions."<< std::endl;

      SetUniqueValues();
   }


   // _____________________________________________________________________ //
   bool fnloSuprocesses::CompareVectors(const std::vector<int>& v1, const std::vector<int>& v2 ) const {
      //!< compare to vectors if all elements are equal
      if ( v1.size() != v2.size() ) return false;
      for ( unsigned int i=0 ; i<v1.size() ; i++ ) {
         if ( v1[i] != v2[i] ) return false;
      }
      return true;
   }


   // _____________________________________________________________________ //
   void fnloSuprocesses::SetUniqueValues() {
      //!< init member SetUniqueValues and 'mother' vector
      for ( unsigned int i = 0 ; i<fInitialValues.size() ; i++ ) {
         if ( !fInitialValues[i].empty() ){
            // look if values are already present ?
            unsigned int umax = fUniqueValues.size();
            for ( unsigned int u = 0 ; u <= umax ; u++ ) {
               if ( u==fUniqueValues.size() ) fUniqueValues.push_back(fInitialValues[i]);
               else if ( CompareVectors(fInitialValues[i],fUniqueValues[u]) )  break;
            }
         }
      }

      // sort a bit ?!
      // put g-g to 0 if available
      vector<int> gg(2);//{0,0}; gg -> X
      for ( unsigned int u = 0 ; u < fUniqueValues.size() ; u++ ) {
         if ( CompareVectors(gg,fUniqueValues[u])) {
            std::swap(fUniqueValues[u],fUniqueValues[0]);
         }
      }
      gg[1]=99;//{0,99}; // DIS eg -> X
      for ( unsigned int u = 0 ; u < fUniqueValues.size() ; u++ ) {
         if ( CompareVectors(gg,fUniqueValues[u])) {
            std::swap(fUniqueValues[u],fUniqueValues[0]);
         }
      }

      // define our map (which is actually a vector)
      for ( unsigned int i = 0 ; i<fInitialValues.size() ; i++ ) {
         (*this)[i] = -1; // init values
         if ( !fInitialValues[i].empty() ){
            for ( unsigned int u = 0 ; u < fUniqueValues.size() ; u++ ) {
               if ( CompareVectors(fInitialValues[i],fUniqueValues[u]) )
                  (*this)[i] = u; // define our final map
            }
         }
      }
      std::cout<<"[fnloSuprocesses::SetUniqueValues()] found "<<fUniqueValues.size()<<" unique subprocess definitions."<<std::endl;
      std::cout<<"[fnloSuprocesses::SetUniqueValues()] Initial values from config file:"<<std::endl;
      for ( unsigned int i = 0 ; i<fInitialValues.size() ; i++ ) if ( fInitialValues[i].size()) std::cout<<i<<"\t"<<fInitialValues[i] << std::endl;
      std::cout<<"[fnloSuprocesses::SetUniqueValues()] Processes definitions:"<<std::endl;
      for ( unsigned int i = 0 ; i<fUniqueValues.size() ; i++ )  std::cout<<i<<"\t"<<fUniqueValues[i] << std::endl<< std::endl;
      std::cout<<"[fnloSuprocesses::SetUniqueValues()] Mapping: "<<(*this)<< std::endl<< std::endl;
   }

   // _____________________________________________________________________ //
   std::vector<std::vector<int > >  fnloSuprocesses::GetPartonCombinations() const {
      auto out = fUniqueValues;
      int n = 0;
      for (auto& i : out ) i.emplace(i.begin(),n++);
      //for ( auto& i : PartonCombinationsLO ) std::cout<<i<< std::endl;
      return out;
   }

} // end namespace


// _____________________________________________________________________ //
// Temporary pointer to APPLGRID lumi function for fastNLO
//const lumi_pdf *fnlopdf;



//-----------------------------------------------------------------------------
//
// --- --- --- fastNLO_nnlo --- --- --
//
//-----------------------------------------------------------------------------

// _____________________________________________________________________ //
void nnlo::init_fastnlo( const int& id, const std::string& gridname, const int& nbins, const double& lo ) {
   //!<
   //!< Remember init variables for later initialization
   //!<
   if ( id >= (int)fnloUtils::fInits.size() ) {
      fnloUtils::ftable.resize(id+1);
      fnloUtils::fIsInclJet.resize(id+1);
      fnloUtils::fIs820GeV.resize(id+1);
      fnloUtils::fIs137.resize(id+1);
      fnloUtils::ftablename.resize(id+1);
      fnloUtils::fInits.resize(id+1);
   }
   fnloUtils::fInits[id] = fnloUtils::fnloInit{int(id), std::string(gridname), int(nbins), &lo };
   fnloUtils::ftable[id] = NULL;
}



// _____________________________________________________________________ //
void nnlo::fill_fastnlo( const int& id, const double& obs, const double* wt, const double* wtref ) {
   //!<
   //!< Fill fastNLO table
   //!<
   //   say::SetGlobalVerbosity(say::INFO);

   // { // --- debugging
   //    double x1 = parfrac_.x1;
   //    double x2 = parfrac_.x2;
   //    double Q2 = mur2_hook(1);
   //    double xfx[13] = {0};
   //    evolvepdf_(x1,sqrt(Q2),xfx);
   //    double as = alphaspdf_(sqrt(Q2));

   //    std::cout<<"id="<<id<<"\tfnloProc="<<fnloUtils::sp[currentprocess_.iproc]
   //          <<"\tx1="<<x1
   //          <<"\tx2="<<x2
   //          <<"\tQ2="<<Q2
   //          <<"\tobs="<<obs
   //          <<"\tgluon-x1="<<xfx[6]
   //          <<"\tas="<<as
   //          <<"\twt="<<wt[0]<<std::endl;

   // }

   static unsigned int nfills = 0;
   static unsigned int nmess  = 0;
   static unsigned int nfmod  = 1;
   nfills++;
   if ( nmess < 10 && nfills % nfmod == 0 ) {
      nmess++;
      std::cout << "[nnlo::fill_fastnlo]: Call no. " << nfills << " of fastnlo fill routine." << std::endl;
   } else if ( nmess == 10 ) {
      nmess = 0;
      nfmod *= 10;
   }

   static const int obs_id_ptj1 = nnlo::GetNLegs() >= 1 ? obs_id_hook("ptj1") : -1;
   static const int obs_id_ptj2 = nnlo::GetNLegs() >= 2 ? obs_id_hook("ptj2") : -1;
   static const int obs_id_ptj3 = nnlo::GetNLegs() >= 3 ? obs_id_hook("ptj3") : -1;
   static const int obs_id_ptj4 = nnlo::GetNLegs() >= 4 ? obs_id_hook("ptj4") : -1;
   //static const int obs_id_Q2("q2");

   // --- init fastNLO
   if ( fnloUtils::ftable[id]==NULL ) fnloUtils::InitFastNLO(id, wt);
   // What is the following line for?
   if ( fnloUtils::ftable[id]==NULL ) return;
   // Check UNIT_PHASE setting of NNLOJET
   // If yes ==> weights nonsense, only phase space counts for grid warmup
   static const bool unitphase = getunitphase_();
   // Check sum of all absolute(!) weights
   double wtsum = 0.;
   for ( int i = 0; i < 7; i++ ) {
     wtsum += std::abs(wt[i]);
   }
   // Tiny weights and not grid warmup ==> nothing to do
   static unsigned int nzwgt = 0;
   static unsigned int nzmss = 0;
   static unsigned int nzmod = 1;
   // if ( fabs(wt[0]) < DBL_MIN && !unitphase ) { // WRONG!!!
   if ( wtsum < DBL_MIN && !unitphase ) {
      nzwgt++;
      if ( nzmss < 10 && nzwgt % nzmod == 0 ) {
         nzmss++;
         std::cout << "[nnlo::fill_fastnlo]: " << nzwgt << "th occurrence of return because of tiny scale weight sum = " << wtsum << std::endl;
      } else if ( nzmss == 10 ) {
         nzmss = 0;
         nzmod *= 10;
      }
      return;
   }

   // --- Keep track on what we have calculated
   fnloUtils::RememberContribution();

   /*
   //static const int& ieorder = order_.ieorder;
   static const int& nloops = ndimcurrent_.nv;

   int nord   = nnlo::GetOrder();//ndimcurrent_.norder;
   int njet     = inppar_.njets;
   if ( nnlo::IsDIS() ) njet -= 1;
   int ne = nord - nloops;
   //std::cout<<"nl="<<nloops<<"\tne="<<ne<<"\t\tnord="<<nord<<"\tIsDIS="<<nnlo::IsDIS()<<"\tnjet="<<njet<<std::endl;
   fnloUtils::fNlNe.insert({ nloops,ne });
   */

   // --- calculate weights and x-values
   double __w2[maxscl] = {0};
   double  __x1=1, __x2=1;
   for ( int is = 0 ; is<maxscl ;is++ ) {//maxscl is a global var in nnlo_common.h
      if ( wt[is] != 0 ) {
         __w2[is] = wt[is];
         if ( fnloUtils::fIs137[id] ) { // zeus is using 1./137
            __w2[is] *= 1./eweakz_.amzZ/eweakz_.amzZ/137./137.;
         }
         if ( fnloUtils::fIs820GeV[id] ) { //
            __w2[is] *= 820./920.;
         }
         nnlo::GetFillWeights(__w2[is],__x1,__x2,is);
         if ( __x1 < 0 ) std::cerr<<"fastnlo_utils. ERROR. x1 from GetFillWeights smaller than 0!! x="<<__x1<<std::endl;

         if ( unitphase ) __w2[is]=0;
      }
      else {
         __w2[is]=0;
      }
   }


   const bool DoTest = false;
   if ( DoTest ) {
      // --- set reference weight for cross check with NNLOJET; ignore during warmup
      if ( ! fnloUtils::ftable[id]->GetIsWarmup() && !unitphase ) {
         fnloUtils::CalculateReferenceWeight(__w2[0],__x1,__x2);
      }
   }

   // --- sanity check
   if ( fnloUtils::sp[currentprocess_.iproc]==-1 ) std::cout<<"Error! iproc="<<currentprocess_.iproc<<"\tProcId="<<fnloUtils::sp[currentprocess_.iproc]<<std::endl;

   // std::cout<<"Obs: "<<obs<<"\tid="<<id<<"\tx1="<<__x1<<std::endl;
   // for ( int is = 0 ; is<maxscl ; is++ ) std::cout<<"\t"<<__w2[is];
   // std::cout<<std::endl;
   // for ( int is = 0 ; is<maxscl ; is++ ) std::cout<<"\t"<<sqrt(mur2_hook(is+1));
   // std::cout<<std::endl;
   // for ( int is = 0 ; is<maxscl ; is++ ) std::cout<<"\t"<<sqrt(muf2_hook(is+1));
   // std::cout<<std::endl;
   // std::cout<<std::endl;

   // --- pass variables to fastNLO
   //double muf2  = muf2_hook(is+1); // muf may change, but we have to fill the value from the
   std::vector<double> vObs;
   if ( nnlo::IsDIS() && fnloUtils::fIsInclJet[id] ) {
      vObs.resize(nnlo::GetNLegs());
      if ( vObs.size()>=1 ) vObs[0] = obs_val_hook(obs_id_ptj1);
      if ( vObs.size()>=2 ) vObs[1] = obs_val_hook(obs_id_ptj2);
      if ( vObs.size()>=3 ) vObs[2] = obs_val_hook(obs_id_ptj3);
      if ( vObs.size()>=4 ) vObs[3] = obs_val_hook(obs_id_ptj4);
      for ( auto& iv : vObs ) {
         // KR: Use C++ macro constants like DBL_MIN instead
         //         if ( fabs(iv) < TINY ) iv=0;
         if ( fabs(iv) < DBL_MIN ) iv=0;
      }
      //if ( fabs(vObs[3]) > TINY ) std::cout<<"obs4="<<vObs[3]<<std::endl;
   }
   else {
      vObs.push_back(obs);
   }

   for ( unsigned int iObs = 0 ; iObs<vObs.size() ; iObs++ ) {
      double thisobs = vObs[iObs];
      if ( thisobs == 0 ) {
         //std::cout<<"Filling: IsIncl="<<fnloUtils::fIsInclJet[id]<<"\tobs="<<thisobs<<"\tcontinue"<<std::endl;
         continue;
      }
      // KR: Use C++ macro constants like DBL_MIN instead
      //      else if ( fabs(thisobs) < TINY ) continue;
      else if ( fabs(thisobs) < DBL_MIN ) continue;
      if ( fnloUtils::ftable[id]->GetIsFlexibleScale() ) {
         double mf  = sqrt(muf2_hook(1));
         double mr  = sqrt(mur2_hook(1));
         std::vector<double> wscl = fnloUtils::SclIndepWgt(__w2,mf,mr);
         double s1 = 0;
         double s2 = 0;
         if ( fnloUtils::fIsInclJet[id] ) {
            // DIS: s1=Q2, s2=pt, pp: s1=pt, s2=XY
            s1 = nnlo::IsDIS() ? mf : thisobs;
            s2 = nnlo::IsDIS() ? thisobs : mf;
         }
         else {
            // id_ptavg_j12 = 24
            static int id_obs_mf = fabs(mf/(int)(mf+0.5)-1) < 1.e-10 ? (int)(mf+0.5) : -1; // fabs(mf-round(mf)) < 1.e-10
            static int id_obs_mr = fabs(mr/(int)(mr+0.5)-1) < 1.e-10 ? (int)(mr+0.5) : -1;
            if ( id_obs_mf != -1 ) mf = obs_val_hook(id_obs_mf);
            if ( id_obs_mr != -1 ) mr = obs_val_hook(id_obs_mr);
            s1 = nnlo::IsDIS() ? mf : mr; // muf: scale1 (pp: scale2)
            s2 = nnlo::IsDIS() ? mr : mf; // mur: scale2 (pp: scale1)
            if ( s1==0 || s2==0 ) {
               std::cout<<"WARNING! scale is zero (this happens if a 2-jet observable is taken as scale for a 1-jet observable): s1="<<s1<<"\ts2="<<s2<<"\tobs="<<thisobs<<"\tIsIncl="<<fnloUtils::fIsInclJet[id]<<"\tid="<<id<<"\t"<<fnloUtils::ftable[id]->GetFilename()<<std::endl;
            }
         }
         //std::cout<<"Filling: IsIncl="<<fnloUtils::fIsInclJet[id]<<"\tobs="<<thisobs<<"\ts1="<<s1<<"\ts2="<<s2<<std::endl;
         fnloUtils::ftable[id]->fEvent.SetProcessId(fnloUtils::sp[currentprocess_.iproc]);
         fnloUtils::ftable[id]->fEvent.SetX1(__x1);
         if ( __x1 < 0 ) std::cerr<<"fastnlo_utils. ERROR. x1 smaller than 0!! x="<<__x1<<std::endl;
         fnloUtils::ftable[id]->fEvent.SetX2(__x2);
         fnloUtils::ftable[id]->fEvent.SetSigma(wtref[0]);
         if ( fnloUtils::fIs820GeV[id] ) { //
            double x820 = __x1 * 920./820.;
            fnloUtils::ftable[id]->fEvent.SetX1(x820);
            if ( x820 > 1 ) continue;
         }
         fnloUtils::ftable[id]->fScenario.SetObservable0(thisobs);
         fnloUtils::ftable[id]->fScenario.SetObsScale1(s1); // in DIS, this must be Q2. In pp we take 'mur' from the runcard
         fnloUtils::ftable[id]->fScenario.SetObsScale2(s2); // in DIS, this is mur from the runcard. In pp we take 'muf' from the runcard
         // --- fill event within fastNLO
         fnloUtils::ftable[id]->fEvent.SetWeight_MuIndependent(wscl[0]);
         if ( nnlo::GetOrder() >= 1 ) {
            if ( nnlo::IsDIS() ){
               //fnloUtils::ftable[id]->fEvent.SetWeight_MuIndependent(wscl[0] + log(s1)*(wscl[1]+wscl[2]));  // DIS needs log(mu^2/q2) factors for ln(mur) and ln(muf) terms
               fnloUtils::ftable[id]->fEvent.SetWeight_MuIndependent(wscl[0]);  // DIS: new convention: no log(q2) terms any longer
               fnloUtils::ftable[id]->fEvent.SetWeight_log_mur(wscl[1]);
               fnloUtils::ftable[id]->fEvent.SetWeight_log_muf(wscl[2]);
               fnloUtils::ftable[id]->fEvent.SetWeight_log_murr(wscl[3]); // DIS: no log(q2) factor here!
            } else {
               fnloUtils::ftable[id]->fEvent.SetWeight_log_mur(wscl[1]);
               fnloUtils::ftable[id]->fEvent.SetWeight_log_muf(wscl[2]);
               fnloUtils::ftable[id]->fEvent.SetWeight_log_murr(wscl[3]);
            }
            if ( nnlo::GetOrder() >= 2 ) {
               fnloUtils::ftable[id]->fEvent.SetWeight_log_muff(wscl[4]);
               fnloUtils::ftable[id]->fEvent.SetWeight_log_murf(wscl[5]);
            }
         }
         fnloUtils::ftable[id]->Fill();
      }
      else { // tables with pre-calculated muf-values
         int nScl = nnlo::GetOrder()>0 ? 3 : 1; // fastNLO convention
         if ( unitphase ) nScl = 1;
         if ( fnloUtils::ftable[id]->GetIsWarmup() ) nScl = 1;
         for ( int is = 0 ; is<nScl ; is++ ) { //maxscl is a global var in nnlo_common.h
            // --- ignore unfilled weights during warmup
            double fillweight = __w2[is];
            if (fnloUtils::ftable[id]->GetIsWarmup()) fillweight = 1;
            // KR: Use C++ macro constants like DBL_MIN instead
            //            if ( fabs(fillweight) > TINY || unitphase) {
            if ( fabs(fillweight) > DBL_MIN || unitphase) {
               double mur2  = mur2_hook(is+1); // mur is defined to be always constant is=0 !
               if ( fnloUtils::fIsInclJet[id] ) {
                  mur2 = thisobs*thisobs * mur2_hook(is+1)/mur2_hook(1);
               }
               fnloUtils::ftable[id]->fEvent.SetProcessId(fnloUtils::sp[currentprocess_.iproc]);
               fnloUtils::ftable[id]->fEvent.SetX1(__x1);
               fnloUtils::ftable[id]->fEvent.SetX2(__x2);
               fnloUtils::ftable[id]->fScenario.SetObservable0(thisobs);
               fnloUtils::ftable[id]->fScenario.SetObsScale1(sqrt(mur2));
               //fnloUtils::ftable[id]->fScenario.SetObsScale2(sqrt(muf2));
               // --- fill event within fastNLO
               fnloUtils::ftable[id]->fEvent.SetWeight(fillweight);
               //std::cout << "Filled weight = " << fillweight << "\t is="<<is<<std::endl;
               fnloUtils::ftable[id]->Fill(is);
            }
         }
      }
   }
   return ;

}



// _____________________________________________________________________ //
void nnlo::fill_ref_fastnlo( const int& id, const double& obs, const double* wt ) {
   // Do nothing for the moment
   // test weights for real radiation
   // if ( wt[0]!=0 ) std::cout<<" wt/ref="<<fnloUtils::_myfillwgt/wt[0]<<"\t\ttest wt="<<wt[0]<<"\twtfill="<<fnloUtils::_myfillwgt<<std::endl;
   // if ( wt[0]!=0 ) std::cout<<"                  "<<"\t\ttest wt="<<wt[0]<<"\twtfill="<<fnloUtils::_myfillwgt<<std::endl;

   // KR: Use C++ macro constants like DBL_MIN instead
   //   if (  fabs(wt[0])> TINY && fnloUtils::_myfillwgt!=0 && fabs(1-fnloUtils::_myfillwgt/wt[0])>1.e-5 ) {
   if (  fabs(wt[0])> DBL_MIN && fnloUtils::_myfillwgt!=0 && fabs(1-fnloUtils::_myfillwgt/wt[0])>1.e-5 ) {
      std::cout<<"ERROR !! Your fill weight is not correctly calculated!"<<std::endl;
      std::cout<<"  wt/ref="<<fnloUtils::_myfillwgt/wt[0]<<"\t\ttest wt="<<wt[0]<<"\twtfill="<<fnloUtils::_myfillwgt<<std::endl;
      std::cout<<"  currentprocess_.iproc:   "<<currentprocess_.iproc<< std::endl;
      std::cout<<"  xregions_.xreg1      :   "<<xregions_.xreg1<<"\t1./xr1="<<1./xregions_.xreg1<< std::endl;
      std::cout<<"  xregions_.xreg2      :   "<<xregions_.xreg2<<"\t1./xr1="<<1./xregions_.xreg2<< std::endl;
      std::cout<<"  parfrac_.x1          :   "<<parfrac_.x1<<"\t1./x1="<<1./parfrac_.x1<<"\t(1-x1)="<<1-parfrac_.x1<<"\t1./(1-x1)="<<1./(1-parfrac_.x1)<< std::endl;
      std::cout<<"  parfrac_.x2          :   "<<parfrac_.x2<<"\t1./x2="<<1./parfrac_.x2<<"\t(1-x2)="<<1-parfrac_.x2<<"\t1./(1-x2)="<<1./(1-parfrac_.x2)<< std::endl;
      //cout<<"  ndim                 :   "<<ndimcurrent_.ndim<< std::endl;
      std::cout<<"  nv                   :   "<<ndimcurrent_.nv<< std::endl;
      std::cout<<"  iproc                :   "<<currentprocess_.iproc<< std::endl;
      std::cout<<"  igx                  :   "<<gridix_.igx<< std::endl;
      std::cout<<"  colflag              :   "<<colourflag_.colflag<< std::endl;
      std::cout<<"  xcoljac              :   "<<xcoljac_.xcoljac<<"\t\t1./xcoljac="<<1./xcoljac_.xcoljac<< std::endl;
      std::cout<<"  jacobian_.jflag      :   "<<jacobian_.jflag<< std::endl;

      static int ii=0;
      if ( ii++ >50 ) {
         std::cout <<" ende test weight and 'x' value for real radiation. use gluon channels only! " <<std::endl;
         exit(1);
      }
    }


}



// _____________________________________________________________________ //
void nnlo::term_fastnlo( const int& id, int nweights ) {
   //!<
   //!< Terminate fastNLO: Write table to disk
   //!<

   if ( vegasiterweight_.it != runinfo_.itmax2 ) return;
   std::cout << "[term_fastnlo()] Write fastNLO table id="<<id<<" with nweights " <<nweights << std::endl;

   //! --- Set fastNLO verbosity level
   //say::SetGlobalVerbosity(say::DEBUG);

   // --- set number of event
   double totwgt = vegasiterweight_.swgt;
   //double totwgt = runinfo_.nevents;
   double wgtfac = 1./runinfo_.itmax2;
   std::cout<<"[nnlo::term_fastnlo] VegasIterwgt:   "<<vegasiterweight_.swgt<<"\tniave="<<vegasiterweight_.iave<<"\tit="<<vegasiterweight_.it<<std::endl;
   std::cout<<"[nnlo::term_fastnlo] Runinfo itmax1: "<<runinfo_.itmax1<<"\titmax2="<<runinfo_.itmax2<<"\tnevents="<<runinfo_.nevents<<std::endl;
   std::cout<<"[nnlo::term_fastnlo] Total weight of table: "<<totwgt<<std::endl;
   std::cout<<"[nnlo::term_fastnlo] Weight factor: "<<wgtfac<<std::endl;

   if ( fnloUtils::ftable[id]==NULL ) {
      std::cout<<"[nnlo::term_fastnlo] Warning! No event was passed into fastNLO table! Do not write empty file!"<<std::endl;
      return;
   }
   //std::string runid = runid_.srunid;
   fnloUtils::ftable[id]->MultiplyCoefficientsByConstant(totwgt * wgtfac);
   fnloUtils::ftable[id]->SetNumberOfEvents(totwgt); // *runinfo_.nshota

   // std::string filename = fnloUtils::ftablename[id] +"."+fnloUtils::GetOrderName( &fnloUtils::fNlNe )+"."+runid+".tab";
   // std::string filename = fnloUtils::ftablename[id] +"."+fnloUtils::GetOrderName( &fnloUtils::fNlNe )+".tab.gz";
   std::string filename = fnloUtils::ftablename[id]+".tab.gz";
   fnloUtils::ftable[id]->SetFilename(filename);


   // --- write table
   if (fnloUtils::ftable[id]->GetIsWarmup())
      std::cout << "[nnlo::term()] fastNLO warmup run finished ..." << std::endl;
   else
      std::cout << "[nnlo::term()] fastNLO production run finished, writing table '" << fnloUtils::ftablename[id] << "'." << std::endl;

//   fnloUtils::ftable[id]->WriteTable(filename); // buggy. This needs a fix!
   fnloUtils::ftable[id]->WriteTable();

   return;
}


#endif

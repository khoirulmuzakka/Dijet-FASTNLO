//
//   @file    fastnlo_utils.cxx
//
//
//   @author D. Britzger, K. Rabbertz
//
//
//

#include "amconfig.h"

#ifdef HAVE_FASTNLO

#include <iostream>
#include <fstream>
#include <algorithm> //c++98
#include <utility> //c++11
#include <sys/stat.h> //c++98

/// fastnlo headers
#include "fastnlo_utils.h"
//#include "fastnlotk/fastNLOCreate.h"

#include "nnlo_common.h"
#include "nnlo_utils.h"



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
   gc.Name = "NNLOJET"; //!< Name and version of generator
   gc.References.push_back("JHEP 1607 (2016) 133 [arXiv:1605.04295]"); //!< References for generator
   gc.References.push_back("Phys. Rev. Lett. 117, 042001 (2016) [arXiv:1606.03991]"); //!< References for generator
   gc.UnitsOfCoefficients = 15;
   return gc;
}



// _____________________________________________________________________ //
//fastNLO::ProcessConstants  fnloUtils::GetNNLOJET_ProcConstsPP(int LOord, const std::vector<std::vector<int> >& partLiCos) {
fastNLO::ProcessConstants fnloUtils::GetNNLOJET_ProcConstsPP() {
   //!< Get default process constants for pp processes

   // --- external input from NNLOJET
   std::vector<std::vector<int> > partLiCos = fnloUtils::sp.GetPartonCombinations();
   if ( partLiCos.empty() ) {
      std::cerr<<"[fnloUtils::GetNNLOJET_ProcConstsPP] Error. No parton combinations found. "
               << "Please read config file into Subprocess-class: fnloUtils::sp .Exiting."<<std::endl;
      exit(5);
   }


   fastNLO::ProcessConstants pc;
   pc.LeadingOrder     = nnlo::GetPowLO(); // Order in alpha_s of leading order process
   pc.NPDF             = 2; // No. of PDFs involved
   pc.NPDFDim          = 2; // 2: full-matrix (then we don't need 'asymmetric processes

   int nprocs = partLiCos.size();
   pc.IPDFdef1         = 3; // Flag 1 to define PDF linear combinations
   pc.IPDFdef2         = 0; // Flag 2 to define PDF linear combinations
   pc.IPDFdef3LO       = nprocs; // No. of  subprocesses
   pc.IPDFdef3NLO      = nprocs;
   pc.IPDFdef3NNLO     = nprocs;

   pc.NSubProcessesLO  = nprocs; //No. of LO   subprocesses
   pc.NSubProcessesNLO = nprocs;
   pc.NSubProcessesNNLO= nprocs;


   // parton combinations in two different formats
   pc.PDFCoeffLO       = VectorToPairs(partLiCos);//!< PDF Linear combinations
   pc.PDFCoeffNLO      = pc.PDFCoeffLO;
   pc.PDFCoeffNNLO     = pc.PDFCoeffLO;
   pc.PDFLiCoInLO      = partLiCos; //!< PDF Linear combinations
   pc.PDFLiCoInNLO     = pc.PDFLiCoInLO;
   pc.PDFLiCoInNNLO    = pc.PDFLiCoInLO;


   //AsymmetricProcesses = ; // empty
   pc.Name             = "Name of process not set and only specified in NNLOJET runcard.";
   pc.References.push_back("Please see NNLOJET manual for more reference."); // references for processes
   return pc;
}



// _____________________________________________________________________ //
fastNLO::ProcessConstants fnloUtils::GetNNLOJET_ProcConstsDIS() {
   //!< Get default process constants for pp processes

   // --- external input from NNLOJET
   //int LOord = inppar_.njets - 1; // DIS: minus 1
   std::cout<< "\n\n LOord = "<<nnlo::GetPowLO()<<std::endl<<std::endl;
   std::vector<std::vector<int> > partLiCos = fnloUtils::sp.GetPartonCombinations();
   if ( partLiCos.empty() ) {
      std::cerr<<"[fnloUtils::GetNNLOJET_ProcConstsDIS] Error. No parton combinations found. "
               << "Please read config file into Subprocess-class: fnloUtils::sp .Exiting."<<std::endl;
      exit(5);
   }

   fastNLO::ProcessConstants pc;
   pc.LeadingOrder     = nnlo::GetPowLO(); // Order in alpha_s of leading order process
   pc.NPDF             = 1; // No. of PDFs involved
   pc.NPDFDim          = 0; // 2: full-matrix (then we don't need 'asymmetric processes

   int nprocs = partLiCos.size();
   pc.IPDFdef1         = 2; // Flag 1 to define PDF linear combinations
   pc.IPDFdef2         = 0; // Flag 2 to define PDF linear combinations
   pc.IPDFdef3LO       = nprocs; // No. of  subprocesses
   pc.IPDFdef3NLO      = nprocs;
   pc.IPDFdef3NNLO     = nprocs;

   pc.NSubProcessesLO  = nprocs; //No. of LO   subprocesses
   pc.NSubProcessesNLO = nprocs;
   pc.NSubProcessesNNLO= nprocs;


   // parton combinations in two different formats
   pc.PDFCoeffLO       = VectorToPairs(partLiCos);//!< PDF Linear combinations
   pc.PDFCoeffNLO      = pc.PDFCoeffLO;
   pc.PDFCoeffNNLO     = pc.PDFCoeffLO;
   pc.PDFLiCoInLO      = partLiCos; //!< PDF Linear combinations
   pc.PDFLiCoInNLO     = pc.PDFLiCoInLO;
   pc.PDFLiCoInNNLO    = pc.PDFLiCoInLO;

   //AsymmetricProcesses = ; // empty
   pc.Name             = "Name of process not set and only specified in NNLOJET runcard.";
   pc.References.push_back("Please see NNLOJET manual for more reference."); // references for processes
   return pc;
}



// _____________________________________________________________________ //
fastNLO::ScenarioConstants fnloUtils::GetNNLOJET_ScenConsts() {
   //!< Get default process constants for pp processes

   fastNLO::ScenarioConstants sc;

   sc.ScenarioName = "fnlXYZ"; // no white spaces here!
   sc.ScenarioDescription.push_back("test for nnlojet"); //< Description of the scenario
   sc.PublicationUnits = 12; //< we are usually working in units of [pb]

   // Dimensionality of binning
   //   also decides if SingleDifferentialBinning or DoubleDifferentialBinning is used
   //   1: single-differential,
   //   2: double-differential;
   sc.DifferentialDimension = 1;

   // Labels (symbol and unit) for the measurement dimensions (from outer to inner "loop"),
   // e.g. "|y|" and "p_T [GeV]".
   // This may also help to define the observables to be calculated in an automatized way!
   sc.DimensionLabels.push_back("varName");

   // Specify:
   //  0 : the cross section is NOT differential,
   //      i.e. there are two bin borders (but NO division (normalization) by bin width);
   //  1 : the cross section is point-wise differential, i.e. only one point is given;
   //  2 : the cross section is bin-wise differential,   i.e. there are two bin borders and division by bin width
   sc.DimensionIsDifferential.push_back(2);

   sc.CalculateBinSize = true; //< Calculate bin width from lower and upper bin boundaries
   sc.BinSizeFactor = 1; //< Possibility to provide additional normalization factor, e.g. of 2 for bins in |y|
   //sc.BinSize = ; //< If 'CalculateBinSize' is 'false' provide table with bin widths for normalization.

   sc.ScaleDescriptionScale1 = "scale1"; //< "<pT_1,2>_[GeV]" # This defines the scale to be used (Note: The 1st scale should always be in units of [GeV]!)
   sc.ScaleDescriptionScale2 = "scale2"; //< "pT_max_[GeV]"   # Specify 2nd scale name and unit (ONLY for flexible-scale tables)

   //DB: binning remains empty here and has to be set later.
   // Observable binning Use either 'SingleDifferentialBinning' or
   //    'DoubleDifferentialBinning' or 'TripleDifferentialBinning' in
   //    accordance with 'DifferentialDimension' above
   // sc.SingleDifferentialBinning;
   // sc.DoubleDifferentialBinning; //< Observable binning
   // sc.TripleDifferentialBinning; //< Observable binning

   sc.CenterOfMassEnergy = inputpar_.roots; //todo //< Center-of-mass energy in GeV. LHC Next Run II: 13000
   sc.PDF1 = 2212; //< PDF of 1st hadron (following PDG convention: proton 2212).
   sc.PDF2 = 2212; //< PDF of 2nd hadron (following PDG convention: proton 2212).

   sc.OutputFilename = "test.tab";//< Filename of fastNLO output table
   sc.OutputPrecision = 8; //< Number of decimal digits to store in output table (def.=8).

   // Create table fully flexible in mu_f
   //    larger size, and requires scale independent weights during creation
   //    or table with fixed number of mu_f scale factors, def.=false.
   sc.FlexibleScaleTable = false;

   // Factorization scale variations
   //   only needed for fixed-scale tables,
   //   List of scale factors must include factor '1',
   //   Scale factors will be ordered according to fastNLO convention:
   //   (1, min, ... , max).
   sc.ScaleVariationFactors = {1.0, 0.5, 2.0};

   sc.ReadBinningFromSteering = true; // Specify if binning is read from fScenConst or from warmup
   sc.ApplyPDFReweighting  = true; //  Apply reweighting of pdfs for an optimized interpolation, def.=true.

   // For warmup-run! Set limits for scale nodes to bin borders, if possible
   sc.CheckScaleLimitsAgainstBins = true;

   //   Choose fastNLO interpolation kernels and distance measures
   sc.X_Kernel = "Lagrange";
   sc.X_DistanceMeasure = "sqrtlog10"; //"3rdrtlog10"
   sc.X_NNodes = 30;
   //sc.X_NNodeCounting = "NodesMax";
   sc.X_NNodeCounting = "NodesPerBin";

   if ( nnlo::IsDIS() ) {
      sc.X_Kernel = "Catmull";
      sc.X_DistanceMeasure = "log10"; // we like to have many nodes at LOW-x !
      sc.X_NNodes = 18;
   }

   sc.Mu1_Kernel = "Lagrange";
   sc.Mu1_DistanceMeasure = "loglog025";//"loglog025";
   sc.Mu1_NNodes = 6;
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
   sc.Mu2_NNodes = 6;

   return sc;
}



// _____________________________________________________________________ //
void fnloUtils::SetNNLOJETDefaultBinning(fastNLO::ScenarioConstants& sc, const int& nbins, const double& lo, const double& hi) {
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
   std::string fname = nnlo::GetProcessName();
   std::cout<<"[fnloUtils::GetConfigFile()] The process name is reported to be: "<<fname<<std::endl;
   fname += "-" + fnloUtils::GetOrderName() + ".config";
   fname = fnloUtils::GetConfigDir() + "/" + nnlo::GetProcessName() + "/" + fname;
   std::cout<<"[fnloUtils::GetConfigFile()] Using config file: "<<fname<<std::endl;
   return fname;
}



// _____________________________________________________________________ //
std::vector<std::string> fnloUtils::GetConfigFiles() {
    //!< get name of the process
    //  "ZJ-LO.config"
   std::vector<std::string> fnames;// = fnloUtils::GetProcessName();
   std::string pref = fnloUtils::GetConfigDir()+"/"+nnlo::GetProcessName()+"/"+ nnlo::GetProcessName();
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
   std::cout<<"[fnloUtils::GetConfigFile()]  Found order "<<nnlo::GetOrder()<<". Using config files: "<<std::endl;
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



std::vector<double> fnloUtils::SclIndepWgt(const double* wt, double muf, double mur){
   //! calculate scale independent weights for fastNLO


   std::vector<double> ret{0,0,0,0,0,0};;
   if ( nnlo::GetOrder()==0 ) {
      //! LO does not has scale dependent terms
      ret[0] = wt[0];
      return ret;
   }


//                                    14          11       5       14       11       5        1          1      1    1       1      1
// Out[4]= {{13, 4, -9, -9, 1, 1}, {-(--), -2, 2, --, 0, -(-)}, {-(--), -2, --, 2, -(-), 0}, {-, 0, 0, -(-), 0, -}, {-, 0, -(-), 0, -, 0}, {1, 1, -1, -1, 0, 0}}
//                                    3           2        6       3        2        6        3          2      6    3       2      6

   double EPS = 1.e-6;
   //unsigned int nSclV = nnlo::GetOrder()==1  ? 3 : 6;
   const unsigned int nSclV =  6;

   if ( nSclV==3 ) {
      // -- not working, but I don't know why!
      static const double SclInvNLO[3][3] = {{3, -2, 0},{0, 1, -1},{-1, 0, 1}}; // i+1
      for ( unsigned int j = 0 ; j<nSclV ; j++ ){
         for ( unsigned int i = 0 ; i<nSclV ; i++ ){
            ret[j] += wt[i+1]*SclInvNLO[j][i];
         }
         //if ( fabs(ret[j]/wt[1]) < EPS ) ret[j] = 0;
      }
      // -- not working, but I don't know why!
      // const double SclInvNLO[3][3] = {{5.5,-2,-2.5},{-1.5,1,0.5},{-0.5,0,0.5}}; // i+3
      // for ( unsigned int j = 0 ; j<nSclV ; j++ ){
      //         for ( unsigned int i = 0 ; i<nSclV ; i++ ){
      //            ret[j] += wt[i+3]*SclInvNLO[j][i];
      //         }
      // }
   }
   else {
      static const double SclInvNNLO[6][6] = {
         {13,4,-9,-9,1,1},
         {-14./3.,-2,2,5.5,0,-5./6.},
         {-14./3.,-2,5.5,2,-5./6.,0},
         {1./3.,0,0,-0.5,0,1./6.},
         {1./3.,0,-0.5,0,1./6.,0},
         {1,1,-1,-1,0,0} };
      int nRet = ( nnlo::GetOrder()==1 ) ? 4 : nSclV;
      for ( int j = 0 ; j<nRet ; j++ ){
         for ( unsigned int i = 0 ; i<nSclV ; i++ ){
            ret[j] += wt[i+1]*SclInvNNLO[j][i];
         }
         if ( fabs(ret[j]/wt[1]) < EPS ) ret[j] = 0; // remove very small weights
      }
      // if ( nnlo::GetOrder()==1 ) {
      //    // ret[3]=0;
      //    ret[4]=0;// must be exact zero for fastNLO
      //    ret[5]=0;// must be exact zero for fastNLO
      // }
   }

//nnlo
   // std::vector<double> ret2(6);
   // const double SclInvNNLO[6][6] = {
   //    {13,4,-9,-9,1,1},
   //    {-14./3.,-2,2,5.5,0,-5./6.},
   //    {-14./3.,-2,5.5,2,-5./6.,0},
   //    {1./3.,0,0,-0.5,0,1./6.},
   //    {1./3.,0,-0.5,0,1./6.,0},
   //    {1,1,-1,-1,0,0} };
   //    for ( unsigned int j = 0 ; j<nSclV ; j++ ){
   //    for ( unsigned int i = 0 ; i<nSclV ; i++ ){
   //       ret2[j] += wt[i+1]*SclInvNNLO[j][i];
   //    }
   //    //if ( fabs(ret[j]/wt[1]) < EPS ) ret[j] = 0;
   //    }

   // double r = log(mur*mur);
   // double f = log(muf*muf);
   // #include "ScaleIndepWgt.C"

   // double EPS = 1.e-6;
   // for ( unsigned int j = 0 ; j<6 ; j++ ){
   //    ret[j] = 0;
   //    for ( unsigned int i = 0 ; i<6 ; i++ ){
   //    ret[j] += wt[i]*SclInv[j][i];
   //    }
   //    if ( fabs(ret[j]/wt[0]) < EPS ) ret[j] = 0;
   // }

   const bool DoDbgTest = false;
   if ( DoDbgTest ) {
      double r = log(mur*mur);
      double f = log(muf*muf);
      double MyWgt0     = ret[0] + r*ret[1] + f*ret[2] + r*r*ret[3] + f*f*ret[4] + r*f*ret[5];
      //double MyWgtNLO   = ret[0] + r*ret[1] + f*ret[2];
      //double MyRefWgt = ret[0] + 2*r*ret[1] + 2*f*ret[2] + 4*r*r*ret[3] + 4*f*f*ret[4] + 4*r*f*ret[5];
      EPS = 1.e-2;
      if  ( fabs(MyWgt0/wt[0]-1) > EPS || !std::isfinite(MyWgt0)) {
         std::cerr<<"[fnloUtils::SclIndepWgt]. Warning!"<<std::endl;
         std::cout<<"  Cond1: "<<(fabs(MyWgt0/wt[0]-1) > 1.e-8)<<"\tvalCond1="<<fabs(MyWgt0/wt[0]-1)<<std::endl;
         //std::cout<<"  Cond1: "<<(fabs(MyWgtNLO/wt[0]-1) > 1.e-8)<<"\tvalCond1="<<fabs(MyWgtNLO/wt[0]-1)<<std::endl;
         std::cout<<"  muf="<<muf<<"\tmur="<<mur<<std::endl;
         std::cout<<"  wt []";
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
         // std::cout<<"  ret2[]";
         // for ( unsigned int j = 0 ; j<6 ; j++ )std::cout<<"\t"<<ret2[j];
         // std::cout<<std::endl;

         std::cout<<"  [TEST] Wgt0  ="<<wt[0]<<"\tMyWgt="<<MyWgt0<<"\tMy/Ref: "<<MyWgt0/wt[0]<<std::endl; //MyRefWeight!
         printf ("\tMy0/Ref= %e  \n",MyWgt0/wt[0]);
         std::cout<<std::endl;
         //exit(3);
      }
   }
   return ret;
}



// _____________________________________________________________________ //
void fnloUtils::InitFastNLO( const int& id, const double* wt ) {
   using namespace std;
   // static const int& iproc = currentprocess_.iproc;
   // int nloops = channel_.iloops[iproc-1]; /// loop order
   // int ipars  = channel_.ipars[iproc-1];  /// number of additional particles + 1
   // cout<<"        nloops          "<<nloops<<endl;
   // cout<<"        ipars           "<<ipars<<endl;
   // cout<<"        ndim.norder     "<<ndimcurrent_.norder<<endl;
   // cout<<"        inppar.njets    "<<inppar_.njets<<endl;
   // cout<<"        order_.ieorder  "<<order_.ieorder<<endl;
   // cout<<"        GetNPos         "<<nnlo::GetNPow()<<endl;
   // cout<<"        GetOrder()      "<<nnlo::GetOrder()<<endl;
   // cout<<endl;

   // --- welcome
   cout << endl;
   cout << "[fnloUtils::InitFastNLO] fastNLO initialisation"  << endl;
   if ( fnloUtils::ftable[id] != NULL )  return;
   cout << "[fnloUtils::InitFastNLO] NNLOJET provides process: " << nnlo::GetProcessName() << endl;

   // --- get settings of histogram [id] from NNLOJET
   const string& gridname = fnloUtils::fInits[id].gridname;
   const int& nbins       = fnloUtils::fInits[id].nbins;
   const double& lo       = *(fnloUtils::fInits[id].lo);
   const double& hi       = *(fnloUtils::fInits[id].hi);
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
   if ( gridname.find("iff") != string::npos ) {
      sc.X_NNodes = 50;
      sc.X_DistanceMeasure = "sqrtlog10"; // we like to have increasingly many nodes at high-x
   }
   if ( nnlo::IsDIS() && gridname.find("q2") == string::npos ) {
      std::cout<<"[fnloUtils::InitFastNLO] IsDIS, but not differential as function of q2, using more q2 nodes, less x-nodes"<<std::endl;
      sc.Mu1_NNodes = 14; // more q2 nodes
      sc.Mu2_NNodes = 6;
      if ( gridname.find("iff") != string::npos ) sc.X_NNodes = 42;
   }
   if ( muf2_hook(2) == mur2_hook(2) &&  fabs(sqrt(mur2_hook(2))/exp(1)-1)<1.e-6) sc.FlexibleScaleTable = true;
   else if ( muf2_hook(1) == mur2_hook(1) ) sc.FlexibleScaleTable = false;
   else {
      cout<<"[fnloUtils::InitFastNLO] Unrecognized scale settings. Please choose either a fixed-scale or a flexible-scale table."<<endl;
      cout<<"                       muf2(1)="<<muf2_hook(1)<<"\tmur2(1)="<<mur2_hook(1)<<endl;
      cout<<"                       muf2(2)="<<muf2_hook(2)<<"\tmur2(2)="<<mur2_hook(2)<<endl;
      cout<<" Exiting."<<endl;
      exit(4);
   }
   fnloUtils::SetNNLOJETDefaultBinning(sc, nbins, lo, hi);
   // if ( nnlo::IsDIS() ){
   //    sc.CheckScaleLimitsAgainstBins = sc.FlexibleScaleTable;
   // }

   // --- adapt for muf-variations
   static const double sclfac[maxscl] = {1,0.5,2.,0.505, 202,0.005,200}; // ugly convention: muf*0.5, muf*2, mur+muf * 0.5, mur+muf* 2, mur*0.5 mur*2
   cout<<"[fnloUtils::InitFastNLO] wt  []:";
   for ( int is = 0 ; is<maxscl ; is++ ) printf("\t%5.4e",wt[is]);
   cout<<endl;
   cout<<"[fnloUtils::InitFastNLO] muf2[]:";
   for ( int is = 0 ; is<maxscl ; is++ ) printf("\t%5.4e",muf2_hook(is+1));
   cout<<endl;
   cout<<"[fnloUtils::InitFastNLO] mur2[]:";
   for ( int is = 0 ; is<maxscl ; is++ ) printf("\t%5.4e",mur2_hook(is+1));
   cout<<endl;

   //if ( !getunitphase_() ) {
      for ( int is = 1 ; is<maxscl && fabs(wt[is])> TINY; is++ )  {
         cout<<"[fnloUtils::InitFastNLO] wt["<<is<<"]="<<wt[is]<<endl;
         sc.ScaleVariationFactors.push_back(sclfac[is]);
      }
      cout<<"[fnloUtils::InitFastNLO] Number of scale factors found: "<<sc.ScaleVariationFactors.size()<<"\t Please respect the convention how to specify those in the run-card."<<endl;
      cout<<"[fnloUtils::InitFastNLO] Scale variation factors: "<<sc.ScaleVariationFactors<<endl;

      if ( sc.FlexibleScaleTable ) {
         // double q2 = muf2_hook(1);
         // double p2 = mur2_hook(1);
         std::cout<<"mur: ";
         for ( int is = 0 ; is<maxscl ; is++ ) std::cout<<"\t"<<sqrt(mur2_hook(1+is)) ;
         std::cout<<std::endl;
         std::cout<<"muf: ";
         for ( int is = 0 ; is<maxscl ; is++ ) std::cout<<"\t"<<sqrt(muf2_hook(1+is)) ;
         std::cout<<std::endl;

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
         //    cerr<<"[nnlo::fill_fastnlo]. Error. This is a flexible scale table, but the definition of the scales does not agree with the convention [1]"<<endl;
         //    exit(3);
         // }
         // if ( fabs( muf2[2] / p2 - 1 ) > EPS || fabs( mur2[2] / p2 - 1 ) > EPS ) {
         //    cerr<<"[nnlo::fill_fastnlo]. Error. This is a flexible scale table, but the definition of the scales does not agree with the convention [2]"<<endl;
         //    exit(3);
         // }
         // if ( fabs( muf2[3] / p2 - 1 ) > EPS || fabs( mur2[3] / q2 - 1 ) > EPS) {
         //    cerr<<"[nnlo::fill_fastnlo]. Error. This is a flexible scale table, but the definition of the scales does not agree with the convention [3]"<<endl;
         //    exit(3);
         // }
         // if ( fabs( muf2[4] / q2 - 1 ) > EPS || fabs( mur2[4] / (p2*p2) - 1 ) > EPS) {
         //    cerr<<"[nnlo::fill_fastnlo]. Error. This is a flexible scale table, but the definition of the scales does not agree with the convention [4]"<<endl;
         //    exit(3);
         // }
         // if ( fabs( muf2[5] / (q2*q2) - 1 ) > EPS || fabs( mur2[5] / p2 - 1 ) > EPS ) {
         //    cerr<<"[nnlo::fill_fastnlo]. Error. This is a flexible scale table, but the definition of the scales does not agree with the convention [5]"<<endl;
         //    exit(3);
         // }
         // // TEST
         // if ( fabs( muf2[6] / (q2*q2) - 1 ) > EPS || fabs( mur2[6] / (p2*p2) - 1 ) > EPS ) {
         //    cerr<<"[nnlo::fill_fastnlo]. Error. This is a flexible scale table, but the definition of the scales does not agree with the convention [6]"<<endl;
         //    exit(3);
         // }
      }

      else if ( !sc.FlexibleScaleTable ) {
         for ( int is = 0 ; is< 3 && is<(int)sc.ScaleVariationFactors.size() ; is++ ){
            double mur0  = mur2_hook(1);
            //double muf0  = muf2_hook(1);
            double mur2  = mur2_hook(is+1);
            double muf2  = muf2_hook(is+1);
            if ( is==0 && mur2!=muf2 ) {
               cerr<<"[nnlo::fill_fastnlo]. Error. First scale setting: muf and mur must be identical, but muf2/mur2="<<muf2/mur2<<endl;
               exit(2);
            }
            // else if ( is==1 && fabs(muf2/mur2/0.25-1) > 1.e-6 ) {
            //    cerr<<"[nnlo::fill_fastnlo]. Error. Second scale variation must be muf*0.5=mur*0.5, but muf2/mur2="<<muf2/mur2<<endl;
            //    exit(2);
            // }
            // else if ( is==2 && fabs(muf2/mur2/4-1) > 1.e-6 ) {
            //    cerr<<"[nnlo::fill_fastnlo]. Error. Third scale variation must be muf*2=mur*2, but muf2/mur2="<<muf2/mur2<<endl;
            //    exit(2);
            // }
            if ( fabs(muf2/mur2-1) > 1.e-6 ) {
               cerr<<"[nnlo::fill_fastnlo]. Error. Scale variation is="<<is<<" must be muf*0.5=mur*0.5, but muf2/mur2="<<muf2/mur2<<endl;
               exit(2);
            }
            if ( is==1 && fabs(mur2/mur0/0.25-1) > 1.e-6 ) {
               cerr<<"[nnlo::fill_fastnlo]. Error. Second scale variation must be muf*0.5=mur*0.5, but mur2/mur2(0)="<<mur2/mur0<<endl;
               exit(2);
            }
            if ( is==2 && fabs(mur2/mur0/4.-1) > 1.e-6 ) {
               cerr<<"[nnlo::fill_fastnlo]. Error. Second scale variation must be muf*2=mur*2, but mur2/mur2(0)="<<mur2/mur0<<endl;
               exit(2);
            }
         }
         if ( fabs(wt[3])>TINY ) {
            cerr<<"[nnlo::fill_fastnlo]. Error. fastNLO only support scale settings of [0] mur=muf; [1] muf*0.5=mur*0.5; [2] muf*2=muf*2"<<endl;
            exit(2);
         }
      }
//}



   // ---- initialize fastNLOCreate
   // --- Collect info for naming
   string wrmfile   = gridname;
   string tablename = gridname;
   string runid = runid_.srunid;
   string pname = nnlo::GetProcessName();
   string jname = jobname_.sjobname;
   jname.erase(remove(jname.begin(), jname.end(), ' '),jname.end());
   string fname = outfile_.fname;
   fname.erase(remove(fname.begin(), fname.end(), ' '),fname.end());
   int seed = fnloUtils::GetSeed();
   cout << "[fnloUtils::InitFastNLO] Input gridname: " << gridname << endl;
   cout << "[fnloUtils::InitFastNLO] Initial output filename: " << fname << endl;
   cout << "[fnloUtils::InitFastNLO] Process name (NNLOJET): " << pname << ", job name: " << jname << endl;
   size_t iposu = gridname.find("_");
   size_t ipos1 = gridname.find("-");
   size_t iposn = gridname.find_last_of("-");
   // Daniel's naming (in case gridname contains an underscore '_')
   if ( iposu != string::npos ) {
      wrmfile   = pname+"."+jname+"."+tablename;
      tablename = pname+"."+jname+"."+runid+"."+tablename+".s"+to_string(seed);
   }
   // Klaus' naming (in case gridname contains a minus '-')
   else if ( ipos1 != string::npos ) {
      string procname = gridname.substr(0,ipos1);
      string histname = gridname.substr(iposn+1,string::npos);
      cout << "[fnloUtils::InitFastNLO] Process name (gridname): " << procname << ", histogram name: " << histname << endl;
      // If unitphase, i.e. grid warmup, then name wrmfile individually for later check & merge
      // If NOT unitphase, i.e. production, then set wrmfile name from gridname identical for all subprocs
      wrmfile = gridname;
      if ( getunitphase_() ) wrmfile = pname+"."+jname+"."+histname;
      tablename = pname+"."+jname+"."+histname;
   }
   else {
      wrmfile   = tablename;
      tablename = tablename;
   }
   wrmfile += ".wrm";

   sc.OutputFilename = tablename;
   fnloUtils::ftablename[id] = tablename;
   fnloUtils::ftable[id] = new fastNLOCreate(wrmfile,gc,pc,sc);
   fnloUtils::fIsInclJet[id] = (gridname.find("ji") != string::npos);
   fnloUtils::fIs820GeV[id] = (gridname.find("820") != string::npos);
   fnloUtils::fIs137[id]    = (gridname.find("137") != string::npos);
   //cout<<"\n\n\n ORDER OF alpha_s(): "<<pc.LeadingOrder+nnlo::GetOrder()<<"\t\t ORDER OF CALCULATION: "<<nnlo::GetOrder()<<"\n\n"<<endl;
   cout<<"\n\n\n POWER OF alpha_s(): "<<nnlo::GetNPow()<<"\t\t ORDER OF CALCULATION: "<<nnlo::GetOrder()<<"\n\n"<<endl;
   fnloUtils::ftable[id]->SetOrderOfAlphasOfCalculation(nnlo::GetNPow());
   fnloUtils::ftable[id]->SetCacheSize(30);
   //fnloUtils::ftable[id]->fastNLOTable::SetNcontrib(1); // to be removed later!

   if ( fnloUtils::ftable[id]->GetIsWarmup() && !getunitphase_() ) {
      cerr<<"[nnlo::init_fastnlo] ERROR! No warmup-file found (hence this cannot be a production run), but the job is NOT running with UNIT_PHASE. Aborted!"<<endl;
      exit(5);
   } else if ( !fnloUtils::ftable[id]->GetIsWarmup() && getunitphase_() ) {
      cerr<<"[nnlo::init_fastnlo] ERROR! Warmup-file found (hence this is a production run), but the job is running WITH UNIT_PHASE. Aborted!"<<endl;
      exit(6);
   }

   cout << "[fnloUtils::InitFastNLO] fastNLO initialized. " << endl;
   return ;
}



// _____________________________________________________________________ //
namespace fnloUtils {
   using namespace std;

   // _____________________________________________________________________ //
   fnloSuprocesses::fnloSuprocesses(const string& config){
      //!< constructor reading config file with subprocesses
      ReadConfigFile(config);
   }

   // _____________________________________________________________________ //
   void fnloSuprocesses::ReadConfigFile(const string& config) {
      fConfigFile = config;
      this->resize(1000);
      fInitialValues.resize(1000);
      ifstream cf(config);
      if ( !cf.is_open() ){
         cerr<<"[fnloSubprocesses::ReadConfigFile()] Error. Could not open file: "<<config<<endl;
         exit(2);
      }
      int id, np, v1, v2;
      int nl = 0;
      cf>>id; // first line is '0'
      while ( !cf.eof() ) {
         cf>>id>>np;
         if(cf.fail()) break;
         nl++;
         //cout<<"nl: "<<nl<<"\tid="<<id<<"\tnp="<<np<<endl;
         for ( int i = 0 ; i<np ; i++ ) {
            cf>>v1>>v2;
            //cout<<"\t"<<v1<<"\t"<<v2<<endl;
            fInitialValues[id].push_back(v1);
            fInitialValues[id].push_back(v2);
         }
      }
      cout<<"[fnloSuprocesses::ReadConfigFile()] found "<<nl<<" subprocess definitions."<<endl;

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
      for ( unsigned int i = 0 ; i<fInitialValues.size() ; i++ ) if ( fInitialValues[i].size()) cout<<i<<"\t"<<fInitialValues[i] <<endl;
      std::cout<<"[fnloSuprocesses::SetUniqueValues()] Processes definitions:"<<std::endl;
      for ( unsigned int i = 0 ; i<fUniqueValues.size() ; i++ )  cout<<i<<"\t"<<fUniqueValues[i] <<endl<<endl;
      cout<<"[fnloSuprocesses::SetUniqueValues()] Mapping: "<<(*this)<<endl<<endl;
   }

   // _____________________________________________________________________ //
   std::vector<std::vector<int > >  fnloSuprocesses::GetPartonCombinations() const {
      auto out = fUniqueValues;
      int n = 0;
      for (auto& i : out ) i.emplace(i.begin(),n++);
      //for ( auto& i : PartonCombinationsLO ) cout<<i<<endl;
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
void nnlo::init_fastnlo( const int& id, const std::string& gridname, const int& nbins, const double& lo, const double& hi ) {
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
   fnloUtils::fInits[id] = fnloUtils::fnloInit{int(id), std::string(gridname), int(nbins), &lo, &hi };
   fnloUtils::ftable[id] = NULL;
}



// _____________________________________________________________________ //
void nnlo::fill_fastnlo( const int& id, const double& obs, const double* wt, const double* wtref ) {
   //!<
   //!< Fill fastNLO table
   //!<
   //say::SetGlobalVerbosity(say::INFO);

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

   //  }

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
   // Tiny weight and not grid warmup ==> nothing todo
   if ( fabs(wt[0]) < TINY && !unitphase ) return;

   // Check weights of other scales, not yet implemented/tested
   // if ( fabs(wt[0])<1.e-20 && fabs(wt[1])>1.e-10 ) {
   //    std::cout<<"Found wt[0]=0, but wt[1]!=0: ";
   //    for ( int is = 0 ; is<maxscl ;is++ ) {//maxscl is a global var in nnlo_common.h
   //    std::cout<<"\t"<<wt[is];
   //    }
   //    std::cout<<std::endl;
   // }

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
            __w2[is] *= 920./820.;
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
         if ( fabs(iv) < TINY ) iv=0;
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
      else if ( fabs(thisobs) < TINY ) continue;
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
            s1 = nnlo::IsDIS() ? mf : mr; // muf: scale1
            s2 = nnlo::IsDIS() ? mr : mf; // mur: scale2
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
            if ( fabs(fillweight) > TINY || unitphase) {
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
   using namespace std;
   // test weights for real radiation
   // if ( wt[0]!=0 ) std::cout<<" wt/ref="<<fnloUtils::_myfillwgt/wt[0]<<"\t\ttest wt="<<wt[0]<<"\twtfill="<<fnloUtils::_myfillwgt<<std::endl;
   // if ( wt[0]!=0 ) std::cout<<"                  "<<"\t\ttest wt="<<wt[0]<<"\twtfill="<<fnloUtils::_myfillwgt<<std::endl;

   if (  fabs(wt[0])> TINY && fnloUtils::_myfillwgt!=0 && fabs(1-fnloUtils::_myfillwgt/wt[0])>1.e-5 ) {
      std::cout<<"ERROR !! Your fill weight is not correctly calculated!"<<std::endl;
      cout<<"  wt/ref="<<fnloUtils::_myfillwgt/wt[0]<<"\t\ttest wt="<<wt[0]<<"\twtfill="<<fnloUtils::_myfillwgt<<std::endl;
      cout<<"  currentprocess_.iproc:   "<<currentprocess_.iproc<<endl;
      cout<<"  xregions_.xreg1      :   "<<xregions_.xreg1<<"\t1./xr1="<<1./xregions_.xreg1<<endl;
      cout<<"  xregions_.xreg2      :   "<<xregions_.xreg2<<"\t1./xr1="<<1./xregions_.xreg2<<endl;
      cout<<"  parfrac_.x1          :   "<<parfrac_.x1<<"\t1./x1="<<1./parfrac_.x1<<"\t(1-x1)="<<1-parfrac_.x1<<"\t1./(1-x1)="<<1./(1-parfrac_.x1)<<endl;
      cout<<"  parfrac_.x2          :   "<<parfrac_.x2<<"\t1./x2="<<1./parfrac_.x2<<"\t(1-x2)="<<1-parfrac_.x2<<"\t1./(1-x2)="<<1./(1-parfrac_.x2)<<endl;
      //cout<<"  ndim                 :   "<<ndimcurrent_.ndim<<endl;
      cout<<"  nv                   :   "<<ndimcurrent_.nv<<endl;
      cout<<"  iproc                :   "<<currentprocess_.iproc<<endl;
      cout<<"  igx                  :   "<<gridix_.igx<<endl;
      cout<<"  colflag              :   "<<colourflag_.colflag<<endl;
      cout<<"  xcoljac              :   "<<xcoljac_.xcoljac<<"\t\t1./xcoljac="<<1./xcoljac_.xcoljac<<endl;
      cout<<"  jacobian_.jflag      :   "<<jacobian_.jflag<<endl;

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
      std::cout << "[nnlo::term()] fastNLO run finished, writing table '" << fnloUtils::ftablename[id] << "'." << std::endl;

//   fnloUtils::ftable[id]->WriteTable(filename); // buggy. This needs a fix!
   fnloUtils::ftable[id]->WriteTable();

   return;
}


#endif

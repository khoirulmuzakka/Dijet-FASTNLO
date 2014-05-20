///@file********************************************************************
///
///     fastNLO_reader: FNLOCPPREAD
///     Program to read fastNLO v2 tables and derive
///     QCD cross sections using PDFs e.g. from LHAPDF
///
///     D. Britzger, K. Rabbertz
///
///********************************************************************
#include <cfloat>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <utility>
#include "fastnlotk/fastNLOInterpolCatmullRom.h"
#include "fastnlotk/fastNLOTable.h"
#include "fastnlotk/fastNLOCreate.h"
#include "fastnlotk/fastNLOReader.h"

#include "fastnlotk/fastNLOLHAPDF.h"

#include "fastnlotk/Alphas.h"
#include "fastnlotk/fastNLOAlphas.h"
#include "fastnlotk/fastNLOCRunDec.h"
#include "fastnlotk/fastNLOLHAPDF.h"
//#include "fastnlotk/fastNLOUser.h"
//#include "fastnlotk/fastNLODiffUser.h"
// Optional:
// If compiled --with-qcdnum support the following line(s) can be commented in
//#include "fastnlotk/fastNLOQCDNUMAS.h"
// End --with-qcdnum
// If compiled --with-hoppet support the following line(s) can be commented in
//#include "fastnlotk/fastNLOHoppet.h"
// End --with-hoppet

/// Function prototype for flexible-scale function
double Function_Mu(double s1, double s2);

//__________________________________________________________________________________________________________________________________
int main(int argc, char** argv) {

   // namespaces
   using namespace std;
   using namespace say;          // namespace for 'speaker.h'-verbosity levels
   using namespace fastNLO;      // namespace for fastNLO constants

   //---  Initialization for nice printing
   const string CSEPS = "##################################################################################\n";
   const string LSEPS = "#---------------------------------------------------------------------------------\n";
   const string CSEPL = "####################################################################################################################################################################\n";

   //---  Parse commmand line
   // Don't miss help output, default is INFO
   SetGlobalVerbosity(MANUAL);
   shout>>"\n";
   shout>>" "<<CSEPS;
   shout<<" fnlo-read: Program Steering"<<endl;
   shout>>" "<<LSEPS;
   string tablename = "table.tab";
   if (argc <= 1) {
      shout<<" fnlo-read: WARNING! No table name given,"<<endl;
      shout<<" taking the default table.tab instead!"<<endl;
      shout<<"   For an explanation of command line arguments type:"<<endl;
      shout<<"   ./fnlo-tk-cppread -h"<<endl;
   } else {
      tablename = (const char*) argv[1];
      if (tablename == "-h") {
         man<<""<<endl;
         man<<"Usage: ./fnlo-tk-cppread [arguments]"<<endl;
         man<<"Table input file, def. = table.tab"<<endl;
         man<<"PDF set, def. = cteq6m.LHpdf"<<endl;
         man<<"   For LHAPDF5: Give full path(s), if the PDF file is not in the cwd."<<endl;
         man<<"   For LHAPDF6: Drop filename extensions and give PDF directory instead."<<endl;
         man<<"Number of mu_r, mu_f scale settings to investigate, if possible, def. = 1, max. = 7"<<endl;
         man<<"Name of desired alpha_s evolution code, def. = GRV."<<endl;
         man<<"   Alternatives are: LHAPDF, RUNDEC, and"<<endl;
         man<<"                     QCDNUM, or HOPPET, IF compiled with these options!"<<endl;
         man<<"                     In the latter case please comment in the optional lines in main.cc!"<<endl;
         man<<""<<endl;
         man<<"Use \"_\" to skip changing a default argument."<<endl;
         man<<""<<endl;
         return 0;
      } else if (tablename == "_") {
         tablename = "table.tab";
         shout<<" fnlo-read: WARNING! No table name given,"<<endl;
         shout<<" taking the default table.tab instead!"<<endl;
      } else {
         shout<<" fnlo-read: Evaluating table: " << tablename << endl;
      }
   }

   //---  PDF set
   string PDFFile = "X";
   if (argc > 2) {
      PDFFile = (const char*) argv[2];
   }
   if (argc <= 2 || PDFFile == "_") {
      PDFFile = "cteq6m.LHpdf";
      shout<<" fnlo-read: WARNING! No PDF set given,"<<endl;
      shout<<" taking cteq6m.LHpdf instead!"<<endl;
   } else {
      shout<<" fnlo-read: Using PDF set   : " << PDFFile << endl;
   }

   //--- Number of scale settings
   unsigned int nscls = 1;
   const unsigned int nsclmax = 7;
   const double xmur[] = { 1.0, 0.5, 2.0, 0.5, 1.0, 1.0, 2.0 };
   const double xmuf[] = { 1.0, 0.5, 2.0, 1.0, 0.5, 2.0, 1.0 };
   string ch2tmp = "X";
   if (argc > 3) {
      ch2tmp = (const char*) argv[3];
   }
   if (argc <= 3 || ch2tmp == "_") {
      shout<<" fnlo-read: No request given for number of scale settings,"<<endl;
      shout<<"            investigating primary scale only."<<endl;
   } else {
      nscls = atoi(argv[3]);
      if (nscls < 1) {
         printf(" # fnlo-read: ERROR! No scale setting or even less??? Aborting! nscls = %i\n",nscls);
         exit(1);
      } else if (nscls > nsclmax) {
         printf(" # fnlo-read: ERROR! Too many scale settings requested, aborting! nscls = %i\n",nscls);
         exit(1);
      } else {
         shout<<" fnlo-read: If possible, will try to do "<<nscls<<" scale setting(s)."<<endl;
      }
   }

   //--- alpha_s evolution code
   string AsEvolCode = "GRV";
   if (argc > 4) {
      AsEvolCode = (const char*) argv[4];
   }
   if (argc <= 4 || AsEvolCode == "_") {
      shout<<" fnlo-read: No request given for alpha_s evolution code,"<<endl;
      shout<<"            using GRV default."<<endl;
   } else {
      shout<<" fnlo-read: Using alpha_s evolution code: " << AsEvolCode << endl;
   }

   //---  Too many arguments
   if (argc > 5) {
      printf("fnlo-read: ERROR! Too many arguments, aborting!\n");
      return 1;
   }
   shout>>" "<<CSEPS;
   //---  End of parsing arguments


   // ************************** fastNLO and example documentation starts here ****************************
   // --- fastNLO user: Hello!
   //     If you use fastNLO for the first time, please read through the
   //     documentation and comments carefully in order to calculate
   //     a reasonable cross section.
   //     All comments that start with '--- fastNLO user:' are intended as a
   //     short documentation for various options, that can be changed by you.
   //
   //     In fastNLO version 2, there are two different types of tables.
   //     Although internally they are implemented slightly differently, both are called
   //     v2 for their larger flexiblity compared to version 1.4.
   //     The simpler ones, v2.0, are extended versions of this previous format
   //     v1.4 from which a conversion into v2.0 is possible, but without profiting
   //     of the improvements, of course.
   //     The second type of tables, v2.1, are called 'flexible-scale' tables
   //     which store the matrix elements in a scale independent way and
   //     also the scale variables are stored more generally.
   //     These tables give you the possibility to change in addition to the renormalization
   //     also the factorization scale by arbitrary factors and have the possiblity to
   //     change the formula according to which the scale is derived.
   //
   //     Please check (see point 2 below), which type of table you are using and
   //     then refer to the comments and functions suitable for this fastNLO table.
   //
   //
   //      0.  This Introduction
   //      1.  Instantiation of fastNLO classes
   //           a. basic FastNLOUser class
   //           b. FastNLOLHAPDF
   //           c. fastNLOAlphas
   //           d. FastNLOCRunDec
   //      2.  Print table information
   //      3.  Calculate cross sections
   //      4.  Modify PDF settings for LHAPDF-interfaces
   //      5.  Change alphas values
   //      6.  Units of the calculation
   //      7.  Contributions and order of calculation
   //      8.  Scale settings
   //      9.  Flexible-scale concept
   //     10.  Access cross section and k-factor
   //     11.  Print-out cross section
   //     12.  Verbosity level
   //     13.  FastNLO for jet-production diffractive DIS
   //           a. Introduction
   //           b. The FastNLOUser.h class
   //           c. Calculate diffractive cross sections.
   //           d. Diffractive DIS example code
   //     14.  Example code
   //     15.  Example analysis


   // 1a.
   // ------- Initialize table for fastNLOReader ------- //
   // --- fastNLO user:
   //     In addition to a fastNLO table two additional ingredients are required:
   //     - the PDF set and
   //     - the alpha_s evolution to be used
   //     These can be freely defined by the user by making an instance of your class
   //     that derives from the FastNLOReader class and passing the name of the
   //     fastNLO table as an argument, e.g.:
   //        FastNLOUser* fnlo = new FastNLOUser( tablename );
   //
   //     To facilitate using fastNLOReader a number of predefined user classes
   //     of FastNLOUser exist, interfacing to
   //     LHAPDF (PDF and alpha_s, see M. Whalley, D. Bourilkov, R. Group, hep-ph/0508110),
   //     GRV Alphas (default alpha_s evolution used in fastNLO for crosschecks,
   //                 based on M. Glueck, E. Reya, A. Vogt, Eur.Phys.J.C5:461-470,1998, hep-ph/9806404;
   //                 PDF from LHAPDF),
   //     CRunDec (alpha_s evolution up to 4 loops, see B. Schmidt, M. Steinhauser,
   //              Comput.Phys.Commun. 183 (2012) 1845-1848, arXiv:1201.6149;
   //              PDF from LHAPDF).
   //
   //     Their use is explained in the following.
   //
   //
   // 1b.
   //     Initialize with PDF from LHAPDF and corresponding alphas value and
   //     evolution for this PDF set. A change of the alpha_s value is only
   //     possible through the choice of the PDF file/set and member, e.g. CT10as.LHgrid
   //         FastNLOLHAPDF fnlolhapdf( tablename , PDFFile , PDFMember );
   //
   //     Print information from LHAPDF
   //         fnlolhapdf.PrintPDFInformation();
   //         int npdf = fnlolhapdf.GetNPDFMembers();
   //         int imaxpdf = fnlolhapdf.GetNPDFMaxMember(); // imaxpdf = npdf - 1
   //
   //     ( Please note that because of a feature in gfortran the output via your LHAPDF
   //       installation may be asynchronous to the C++ output. Usually, the gfortran
   //       output comes at the end after all C++ output, but this depends on your actual system.
   //       You can try to set the environment variable GFORTRAN_UNBUFFERED_ALL to yes
   //       in your shell to get it synchronized. Keep your fingers crossed. )
   //
   //
   // 1c.
   //     Initialize with PDF from LHAPDF and GRV alphas evolution (default)
   //         fastNLOAlphas fnlo( tablename , PDFFile , PDFMember );
   //
   //     Change the alpha_s value through
   //         fnlo.SetAlphasMz(0.1179);
   //     Change values of the alpha_s evolution code through:
   //         Alphas::SetNf(5);
   //         Alphas::SetMz(91.70);
   //     For all options see Alphas.h
   //
   //
   // 1d.
   //     Initialize with PDF from LHAPDF and RunDec alpha_s evolution.
   //         FastNLOCRunDec fnlo( tablename , PDFFile , PDFMember );
   //     Change the alpha_s value for all instances, by:
   //         fnlo.SetAlphasMz(0.1179);
   //     Change values of the alpha_s evolution code through:
   //         fnlo.SetNf(5);
   //         fnlo.SetNloop(4);
   //         fnlo.SetMz(91.70);
   //
   //     (Note: CTEQ6M:   M_Z = 91.70,   alpha_s(M_Z) = 0.1179;
   //            PDG 2012: M_Z = 91.1876, alpha_s(M_Z) = 0.1184)
   //


   // 2.
   // ---- Table information ---- //
   // --- fastNLO user: For a comprehensive insight into the fastNLO variables
   //     you can use:
   //             fnlo.PrintFastNLOTableConstants(0);
   //


   // 3.
   // ---- (Re-)calculate cross sections ---- //
   // --- fastNLO user: Before you can access the fastNLO computed
   //     cross sections, you always have to call CalcCrossSection()!
   //     So, before accessing the cross sections, please call:
   //             fnlo.CalcCrossSection();


   // 4.
   // ------- Select another PDF set and member ------- //
   // --- fastNLO user: You can select another PDF set and member here.
   //     With LHAPDF, you can set the PDF set and member using e.g.:
   //           fnlo.SetLHAPDFFilename( string PDFFile );
   //           fnlo.SetLHAPDFMember( int PDFMember );
   //


   // 5.
   // ------- Changing the alpha_s(M_Z) value and/or evolution ------- //
   // --- fastNLO user:
   //     The alpha_s evolution is provided by the code of the chosen
   //     interface, e.g. GRV alpha_s for the fnlo instance here.
   //     The value of alpha_s(M_Z) can be changed from its default PDG 2012 values
   //     like this:
   //
   //            fnlo.SetAlphasMz(0.1179);
   //
   //     (Note: CTEQ6M:   M_Z = 91.70,   alpha_s(M_Z) = 0.1179;
   //            PDG 2012: M_Z = 91.1876, alpha_s(M_Z) = 0.1184)
   //
   //     To use a different alpha_s evolution code one has to interface it.
   //     Here, for example, we use the above-mentioned CRunDec code:
   //
   //  FastNLOCRunDec fnlocrundec( tablename , PDFFile , 0 );
   //  fnlocrundec.SetMz(91.1876);
   //  fnlocrundec.SetAlphasMz(0.1184);
   //  fnlocrundec.CalcCrossSection();


   // 6.
   // ------- Set the units of your calculation (kPublicationUnits or kAbsoluteUnits) ------- //
   // --- fastNLO user: You can choose the units in which you want
   //     to access (or print) your cross-section results.
   //     There are two possibilites:
   //       - The default option is 'publication units', i.e. divided by
   //         bin widths if done so in the relevant publication
   //            fnlo.SetUnits(fastNLO::kPublicationUnits);
   //       - The other option is 'absolute' units in barn, but still in
   //         the same magnitude as in the publication (e.g. pb, fb, nb, etc.)
   //
   //       fnlo.SetUnits(kAbsoluteUnits); // in namespace fastNLO
   //     or
   //       fnlo.SetUnits(kPublicationUnits); // in namespace fastNLO


   // 7.
   // ------- Set the calculation order (if available) ------- //
   // --- fastNLO user: Each fastNLO table comes typically with
   //     various contributions.
   //     Currently, five different types of contributions have been tested.
   //     Three can be combined to give a scale, PDF and alpha_s dependent
   //     cross-section, one is a fixed multiplicative correction and, at last,
   //     also data points with uncertainties might be included in a table.
   //     For calculating a cross section, by default only the LO & NLO contributions
   //     are used. However, each contribution can be swiched on or off separately.
   //     Please make sure to avoid combinations that do not make sense,
   //     e.g. 2-loop threshold corrections with LO pQCD.
   //
   //     For switching a contribution on/off, its type must be known:
   //       - kFixedOrder                  -> Fixed order calculation (in alpha_s)
   //       - kThresholdCorrection         -> Threshold corrections
   //       - kElectroWeakCorrection       -> Electroweak corrections (not derived yet)
   //       - kNonPerturbativeCorrections  -> Non-perturbative corrections|Hadronisation corrections
   //     plus one must know the 'Id' of this contribution, which can be printed e.g.
   //     by calling
   //        fnlo.PrintTableInfo();
   //
   //     To switch a contribution on/off please use:
   //            bool SetOn = fnlo.SetContributionON( contrib, Id, on/off )
   //     and in particular for switching on check on the return value SetOn that it actually worked.
   //     Here, 'contrib' is not the contribution number, but the type
   //     as given above: kFixedOrder, ...
   //     Within each type the contributions are counted separately starting with Id=0.
   //     The total number of contributions then counts all contributions of all types.


   // 8.
   // ------- Selecting the scale treatment ------- //
   // --- fastNLO user: The simplest way to modify the predefined renormalization and
   //     factorization scales is to provide a scale factor by which the default scale
   //     is multiplied. These factors must be positive and not too small (> 1.e-6).
   //     Otherwise they can in principal (within reason) be set arbitrarily for
   //     flexible-scale tables. For the normal v2 tables the choice of factors for the
   //     factorization scale is limited to some fixed values, usually 0.5, 1.0, and 2.0
   //     plus sometimes also 0.25, see the respective table information.
   //     Note: If threshold corrections are available and switched on for evaluation,
   //     the scale factors for the renormalization and factorization scale must be identical.
   //
   //     The function call to set the scale factors is:
   //         bool SetScales = fnlo.SetScaleFactorsMuRMuF(xmur, xmuf);
   //     where xmur, xmuf are the scale factors. Check the return value in order to verify
   //     that the selected scale factors could actually be activated.
   //
   //     The return value of this function call is boolean and returns false, if the
   //     the requested scale factors can not be chosen. In this case, the last legal
   //     values remain unchanged.


   // 9.
   // ----- Additional possibilities for scales in 'flexible-scale' tables (v2.1) ----- //
   //     First check, if your table is a flexible-scale table or not
   //          bool IsFlex = fnlo.GetIsFlexibleScaleTable()
   //     You can choose a function to define how
   //     to compute the renormalization and factorization scale.
   //     Each 'flexible-scale' table comes with two variables that can be used
   //     for calculating the scales. They are called scale1 and scale2 and
   //     at least one needs to have a dimension in "GeV".
   //     DIS tables have typically stored scale1 = Q and scale2 = pt, while
   //     hadron-hadron tables might have for example scale1 = pt and scale2 = y.
   //     Other settings are imaginable. Please check, which obervables exactly
   //     are stored as scale variables!
   //
   //     There are two possibilities, how you can define your scale now:
   //
   //       - use predefined functions using e.g.
   //            fnlo.SetMuRFunctionalForm(fastNLO::EScaleFunctionalForm);
   //         for changing the calculation of the renormalizatoin scale.
   //         Please refer to FastNLOReader.h for all options of EScaleFunctionalForm.
   //
   //       - or you can pass a function pointer to FastNLOReader using
   //            fnlo.SetExternalFuncForMuR( double (*Func)(double,double) );
   //         to pass any function using scale1 and scale2 to fastNLO.
   //
   //     WARNING: Some choice had to be made for the default settings. Please think
   //     carefully about the choice of the scales ...
   //     Default setting for DIS tables:
   //       - mu_r:  kQuadraticMean      -> mu_r = sqrt( (Q^2 + scale2^2)/2. ) // because scale1=Q!
   //       - mu_f:  kScale1             -> mu_f = Q
   //     Default setting for pp and ppbar tables:
   //       - mu_r:  kScale1             -> mu_r = scale1
   //       - mu_f:  kScale1             -> mu_f = scale1
   //
   //     Valid calls are e.g.:
   //     fnlo.SetMuRFunctionalForm(fastNLO::kScale1);        // set function how to calculate mu_r from scale1 and scale2
   //     fnlo.SetMuFFunctionalForm(fastNLO::kScale1);        // set function how to calculate mu_f from scale1 and scale2
   //     fnlo.SetMuRFunctionalForm(fastNLO::kQuadraticMean); // set function how to calculate mu_r from scale1 and scale2
   //     fnlo.SetMuFFunctionalForm(fastNLO::kScale1);        // set function how to calculate mu_f from scale1 and scale2
   //     fnlo.SetExternalFuncForMuR( &Function_Mu );         // set external function to calculate mu_r from scale1 and scale2
   //     fnlo.SetMuRFunctionalForm(fastNLO::kExpProd2);      // set function how to calculate mu_f from scale1 and scale2
   //     fnlo.SetMuFFunctionalForm(fastNLO::kExpProd2);      // set function how to calculate mu_f from scale1 and scale2
   //
   // INFO: All above-mentioned scale changing functions automatically perform a refilling of the
   //       fastNLO internal PDF cache. To switch it off you can use a boolean, like:
   //       fnlo.SetMuFFunctionalForm(fastNLO::kScale1 , false );


   // 10.
   // ---- Access cross sections ---- //
   // --- fastNLO user: To access the cross section from fastNLO
   //     you should use:
   //           vector < double > xs = fnlo.GetCrossSection();
   //     If you want to have a pointer to an array of numbers you might use
   //           vector < double > xs = fnlo.GetCrossSection();
   //           double* cs = &xs[0];
   //
   //     Further you can access the "k-factor", which is calculated with all
   //     'contributions' that are switched on (e.g. non-perturbative corrections)
   //     against the LO fixed-order contribution.
   //     Remark:
   //          - the proverbial k-factor is NLO vs. LO
   //          - 1-loop threshold corrections are vs. LO
   //          - 2-loop threshold corrections are vs. NLO
   //          - non-perturbative corrections usually are vs. NLO
   //
   //           vector < double > kFactors = fnlo.GetKFactors();


   // 11.
   // ---- Printing ---- //
   // --- fastNLO user: For an easy overview of your cross section calculation
   //     you might use the following print methods:
   //             fnlo.PrintCrossSections();
   //
   //     Or print it like the Fortran reader code:
   //             fnlo.PrintCrossSectionsDefault();


   // 12.
   // ------- Set fastNLOReader verbosity ------- //
   // --- fastNLO user:
   //     The following line sets the verbosity level of fastNLOReader
   //     Six different levels are implemented, the default is INFO:
   //     DEBUG, MANUAL, INFO, WARNING, ERROR, SILENT
   //         SetGlobalVerbosity(WARNING);
   //     Alternatively, a specific verbosity level can be set
   //     to any instance:
   //         fnlo.SetVerbosity(level);


   // 13.
   // ------- FastNLO for jets in diffractive DIS ------- //
   // 13a.
   //  FastNLO is also applicable to jets in diffractive DIS.
   //  The calculation of jet cross sections in diffractive
   //  DIS is performed by adapting the slicing method,
   //  where the xpom integration is performed during the evaluation
   //  of the fastNLO table. The differential cross section
   //  in xpom is calcualted by a rescaling of the center-of-mass
   //  energy of the incident hadron.
   //  The boundaries of the integration interval are automatically
   //  smoothed out.
   //  More details on the applied method can be found on the
   //  website, i.e.
   //  http://fastnlo.hepforge.org/docs/talks/20120912_fastNLOv2_DBritzger_DiffractiveFastNLO.pdf
   //
   // 13b.
   // --- fastNLO user:
   //  In order to calculate diffractive DIS processes, the user
   //  has to provide a diffractive PDF, as well as an alpha_s
   //  evolution code. Both pieces have to be implemented in the
   //  FastNLODiffUser.h file, where the functions
   //     double FastNLODiffUser::EvolveAlphas(double Q)
   //     bool FastNLODiffUser::InitPDF()
   //     vector<double> FastNLODiffUser::GetDiffXFX(double xpom, double zpom, double muf)
   //  have to be implemented in a reasonable way.
   //  Some examples and more help on this, can provide the authors.
   //  The implementation of the alpha_s evolution code can also be
   //  adapted e.g. from fastNLOAlphas.h or FastNLOCRunDec.h.
   //
   // 13c.
   //  The calculation of diffractive cross sections performs
   //  an integration of xpom. This is done by a simple Riemann integration.
   //  Four possibilities to define the slicing are implemented.
   //  1. Use a logarithmic xpom slicing
   //     Set the number of slices, the xpom_min and xpom_max range, e.g.:
   //       fnlodiff->SetXPomLogSlicing( 12, pow(10.,-2.3), pow(10.,-1) );
   //  2. Use a linear xpom slicing
   //       fnlodiff->SetXPomLinSlicing( 12, 0.0, 0.1 );
   //  3. Use an exponential xpom slicing
   //  4. Set your individual xpom slicing. This basically also allows
   //     to implement a MC integration.
   //         nStep:     number of slices
   //         xpom[]:    central value of each slice
   //         dxpom[]:   width of each slice
   //     fnlodiff->SetXPomSlicing(int nStep, double* xpom, double* dxpom);
   //
   //  To calculate and access the cross sections use:
   //        vector<double> xs = fnlodiff->GetDiffCrossSection();
   //
   //  If you want to calculate cross sections as fucntion of xpom,
   //  you have to calculate each xpom bin by setting the 'xpomslicing', and
   //  summing all bins by yourself.
   //  WARNING:
   //  In this case, one always have to call SetUnits(fastNLO::kAbsoluteUnits) !
   //
   //  Tipp 1: Some brief studies showed, that already with ca. 10 slices, the
   //  cross section converges sufficiently fast. The linear slicing is
   //  preferred over the logarithmic slicing.
   //  Tipp 2:
   //  Choosing Q2 (or pT) as factorization scale increases the speed significantly.
   //
   // 13d.
   //  In the following example code the class description FastNLODiffUser may be
   //  replaced by a specific interface class to a diffractive PDF (see 13b). This class
   //  has to be added to the include statements above, e.g.:
   //     #include "fastnlo/FastNLODiffUser.h"
   //  Some example code could look like (uncomment the following lines,
   //  comment out the other examples under 14. and 15., and recompile):
   //
   //    // ---- Example code for jet-cross sections in diffractive DIS ---- //
   //    //  we setup an instance of the FastNLODiffUser class
   //    FastNLODiffUser fnlodiff( tablename );
   //
   //    //  If you want to receive your cross section in
   //    //   pb/GeV or in pb. Here we choose pb/GeV
   //    fnlodiff.SetUnits(fastNLO::kPublicationUnits);
   //
   //    // Set the xpom integration interval and method
   //    fnlodiff.SetXPomLinSlicing( 12, 0.0, 0.1 );
   //
   //    // Optional:
   //    // make your scale definition (see above)
   //    fnlodiff.SetMuFFunctionalForm(kQuadraticSum);
   //    fnlodiff.SetMuRFunctionalForm(kQuadraticSum);
   //    fnlodiff.SetScaleFactorsMuRMuF(1.0,1.0);
   //
   //    // calculate and access the cross section
   //    vector<double>  xs = fnlodiff.GetDiffCrossSection();
   //    // Print it
   //    fnlodiff.PrintCrossSections();
   //    // ------------------------------------------------------------------ //


   // 14.
   // ---- Example of a cross section calculation with some nice standardized output

   // Select verbosity level
   SetGlobalVerbosity(WARNING);

   // Instead of instantiating a class via e.g.
   //   FastNLOAlphas fnlo(tablename, PDFFile, PDFMember);
   // and accessing member functions via fnlo.memberfunction()
   // we create a pointer to the yet unfilled fnlo table using the FastNLOLHAPDF class.
   // (The other classes inherit the interface structure from FastNLOLHAPDF!)
   // As a difference to previous examples we now have to dereference the pointer via
   //   fnlo->memberfunction()!
   //
   fastNLOLHAPDF* fnlo = NULL;
   if (AsEvolCode == "GRV") {
      fnlo = new fastNLOAlphas(tablename);
   } else if (AsEvolCode == "LHAPDF") {
      fnlo = new fastNLOLHAPDF(tablename);
   } else if (AsEvolCode == "RUNDEC") {
      fnlo = new fastNLOCRunDec(tablename);
      // Optional:
      // If compiled --with-qcdnum support the following line(s) can be commented in
      //   } else if (AsEvolCode == "QCDNUM") {
      //      fnlo = new fastNLOQCDNUMAS(tablename);
      // End --with-qcdnum
      // If compiled --with-hoppet support the following line(s) can be commented in
      //   } else if (AsEvolCode == "HOPPET") {
      //      fnlo = new fastNLOHoppetAs(tablename);
      // End --with-hoppet
   } else {
      printf("fnlo-read: ERROR! Unknown alpha_s evolution code %s!\n",AsEvolCode.c_str());
      printf("           If you compiled with optional QCDNUM or HOPPET support, please\n");
      printf("           do not forget to comment in the marked lines in main.cc!\n");
      exit(1);
   }

   // Print some fastNLO table info
   fnlo->PrintTableInfo();
   fnlo->PrintFastNLOTableConstants(0);

   // Define the PDF set and member
   fnlo->SetLHAPDFFilename(PDFFile);
   fnlo->SetLHAPDFMember(0);

   // The table and PDF initialization could also be done in one step, e.g.:
   //   FastNLOAlphas fnlo(tablename, PDFFile, 0);
   //
   // From here on parameter settings can be read from LHAPDF, e.g.
   // to check the upper limit of the PDF member numbering do
   //   int imaxpdf = fnlo->GetNPDFMaxMember();
   // Note: Usually there is a member no. 0 corresponding to the central result,
   //       so the total number of members (fnlo->GetNPDFMembers()) is imaxpdf + 1
   //
   // Some remarks:
   // 1) Values returned via LHAPDF are not always coherent or consistent
   //    between the different PDF sets. Take care if you want to use them
   //    for more than just information.
   //
   // 2) Astonishingly, there are no functions to access the value used for M_Z or
   //    directly the value for alpha_s(M_Z).
   //
   // 3) The PDF sets within LHAPDF are usually available in the form of
   //    precalculated interpolation grids. Changing of parameters therefore
   //    is not possible.
   //    However, with fastNLO that separates matrix-element coefficients,
   //    PDF weights, and factors of alpha_s^n it is possible to replace
   //    the alpha_s evolution and the parameters therein. By setting one of
   //    the parameters below this effects ONLY the alpha_s evolution, if
   //    allowed by the corresponding code. For LHAPDF this has no effect.

   // The number of loops used for the alpha_s evolution; the LO is nloop = 1, i.e
   // this number corresponds to LHAPDF: getOrderAlphaS + 1
   // Order n PDFs usually should be accompanied by an alpha_s evolution of the same order.
   // int nloop = fnlo->GetNLoop();
   // cout << "Read from LHAPDF: Number of loops = " << nloop << endl;
   fnlo->SetNLoop(2);// NLO

   // int nflavor = fnlo->GetNFlavor();
   // cout << "Read from LHAPDF: Number of flavors = " << nflavor << endl;
   fnlo->SetNFlavor(5);// CTEQ
   //   fnlo->SetNFlavor(0);// NNPDF

   // Unfortunately, LHAPDF5 has no function to access M_Z!
   //   double Mz = 91.174;  // ABM11
   //   double Mz = 91.188;  // CTEQ
   //   double Mz = 91.187;  // HERAPDF
   double Mz = 91.1876; // MSTW, PDG 2013
   //   double Mz = 91.2;    // NNPDF
   fnlo->SetMz(Mz);

   // Unfortunately, LHAPDF5 neither has a function to access alphas(M_Z) directly!
   // double asmz = fnlo->GetAlphasMz(Mz);
   // cout << "Read from LHAPDF: alpha_s at M_Z = " << asmz << endl;
   fnlo->SetAlphasMz(0.1184);// PDG 2013
   //   fnlo->SetAlphasMz(0.1180);// CT10-NLO
   //   fnlo->SetAlphasMz(0.1190);// NNPDF21-NLO

   // Read quark masses
   // for (int iq=1;iq<7;iq++) {
   //    double mq = fnlo->GetQMass(iq);
   //    cout << "Read from LHAPDF: For quark PDG code " << iq << " the quark mass is = " << mq << endl;
   // }
   //   double mt = fnlo->GetQMass(6);
   //   fnlo->SetQMass(6,mt);

   // Calculate cross sections
   fnlo->InitEvolveAlphas();
   fnlo->CalcCrossSection();
   // Uncomment this to actually print out the result
   //   fnlo->PrintCrossSectionsDefault();
   //
   // Example code to print out data points (if available)
   //   fnlo->PrintCrossSectionsData();
   //
   // Example code to print out alpha_s(Q) values
   // for (int iq = 10; iq < 2010; iq = iq + 10) {
   //    double mu = iq;
   //    double as = fnlo->CalcAlphas(mu);
   //    printf("%#18.11E %#18.11E\n",mu,as);
   // }
   // ************************************************************************************************


   // 15.
   // ---- Example to do some cross section analysis ---- //
   // Some initialization
   string CSEP41("#########################################");
   string DSEP41("=========================================");
   string SSEP41("-----------------------------------------");
   string CSEP = CSEP41 + CSEP41 + CSEP41 + CSEP41;
   string DSEP = DSEP41 + DSEP41 + DSEP41 + DSEP41;
   string SSEP = SSEP41 + SSEP41 + SSEP41 + SSEP41;
   printf("\n");
   printf("%s",CSEPL.c_str());
   printf("fnlo-read: Calculate my cross sections\n");
   printf("%s",CSEPL.c_str());

   // Instance fastNLO (For this example we assume fnlo was instantiated already above ...)
   // fastNLOAlphas fnlo( tablename , PDFFile , 0 );

   // Check on existence of LO (Id = -1 if not existing)
   int ilo   = fnlo->ContrId(kFixedOrder, kLeading);
   if (ilo < 0) {
      error["fnlo-read"] << "LO not found, nothing to be done!" << endl;
      exit(1);
   } else {
      info["fnlo-read"] << "The LO contribution has Id: " << ilo << endl;
   }
   // Check on existence of NLO (Id = -1 if not existing)
   int inlo  = fnlo->ContrId(kFixedOrder, kNextToLeading);
   if (inlo < 0) {
      info["fnlo-read"] << "No NLO contribution found!" << endl;
   } else {
      info["fnlo-read"] << "The NLO contribution has Id: " << inlo << endl;
   }
   // Check on existence of threshold corrections
   int ithc1 = fnlo->ContrId(kThresholdCorrection, kLeading);
   int ithc2 = fnlo->ContrId(kThresholdCorrection, kNextToLeading);
   if (ithc1 < 0) {
      info["fnlo-read"] << "1-loop threshold corrections not found!" << endl;
   } else {
      info["fnlo-read"] << "1-loop threshold corrections have Id: " << ithc1 << endl;
   }
   if (ithc2 < 0) {
      info["fnlo-read"] << "2-loop threshold corrections not found!" << endl;
   } else {
      info["fnlo-read"] << "2-loop threshold corrections have Id: " << ithc2 << endl;
   }
   // Check on existence of non-perturbative corrections from LO MC
   int inpc1 = fnlo->ContrId(kNonPerturbativeCorrection, kLeading);
   if (inpc1 < 0) {
      info["fnlo-read"] << "Non-perturbative corrections not found!" << endl;
   } else {
      info["fnlo-read"] << "Non-perturbative corrections have Id: " << inpc1 << endl;
   }

   // Run over all pre-defined scale settings xmur, xmuf
   for (unsigned int iscls=0; iscls<nscls; iscls++) {

      // Switch on LO & NLO, switch off anything else
      if (!(ilo   < 0)) {
         bool SetOn = fnlo->SetContributionON(kFixedOrder, ilo, true);
         if (!SetOn) {
            error["fnlo-read"] << "LO not found, nothing to be done!" << endl;
            error["fnlo-read"] << "This should have been caught before!" << endl;
            exit(1);
         }
      }
      if (!(inlo  < 0)) {
         bool SetOn = fnlo->SetContributionON(kFixedOrder, inlo, true);
         if (!SetOn) {
            error["fnlo-read"] << "NLO not found, nothing to be done!" << endl;
            error["fnlo-read"] << "This should have been caught before!" << endl;
            exit(1);
         }
      }
      if (!(ithc1 < 0)) {
         fnlo->SetContributionON(kThresholdCorrection, ithc1, false);
      }
      if (!(ithc2 < 0)) {
         fnlo->SetContributionON(kThresholdCorrection, ithc2, false);
      }
      if (!(inpc1 < 0)) {
         fnlo->SetContributionON(kNonPerturbativeCorrection, inpc1, false);
      }

      // Define result vectors
      bool lscvar = false;
      bool lthcvar = false;
      vector < double > qscl;
      vector < double > xslo;
      vector < double > xsnlo;
      vector < double > kfac;
      vector < double > xsthc1;
      vector < double > kthc1;
      vector < double > xsthc2;
      vector < double > kthc2;
      vector < double > xsnpc1;
      vector < double > knpc1;
      vector < double > xsnpc2;
      vector < double > knpc2;
      vector < double > xsewk1;
      vector < double > kewk1;

      // Set MuR and MuF scale factors for pQCD cross sections and test availability
      lscvar = fnlo->SetScaleFactorsMuRMuF(xmur[iscls], xmuf[iscls]);
      if (!lscvar) {
         warn["fnlo-read"] << "The selected scale variation (xmur, xmuf) = ("
                           << fnlo->GetScaleFactorMuR() << ","
                           << fnlo->GetScaleFactorMuF() << ") is not possible with this table, skipped completely!" << endl;
         //         printf("fnlo-read: WARNING! The selected scale variation (xmur, xmuf) = (% #10.3f, % #10.3f) is not possible with this table, skipped completely!\n",fnlo->GetScaleFactorMuR(),fnlo->GetScaleFactorMuF());
         // skip completely
         continue;
      }
      if (fnlo->GetIsFlexibleScaleTable()) {
         fnlo->SetMuFFunctionalForm(kScale1);
         fnlo->SetMuRFunctionalForm(kScale1);
         //      fnlo->SetMuFFunctionalForm(kScale2);
         //      fnlo->SetMuRFunctionalForm(kScale2);
      }

      // Calculate cross section
      fnlo->CalcCrossSection();

      // Get LO & NLO results
      xsnlo = fnlo->GetCrossSection();
      kfac  = fnlo->GetKFactors();
      // Set order for Q scale determination, rel. to LO: 0 --> LO, 1 --> NLO
      int irelord = 1;
      qscl  = fnlo->GetQScales(irelord);
      xslo  = xsnlo;
      for (unsigned int i=0; i<xslo.size(); i++) {
         if (abs(kfac[i]) > DBL_MIN) {
            xslo[i] = xslo[i]/kfac[i];
         } else {
            xslo[i] = -1.;
         }
      }

      // Get threshold corrections
      if ( !(inlo < 0 || ithc2 < 0) ) {
         bool SetOn = fnlo->SetContributionON(kThresholdCorrection, ithc2, true);
         if (!SetOn) {
            warn["fnlo-read"] << "2-loop threshold corrections could not be switched on, skip threshold correction factors!" << endl;
            //            printf("fnlo-read: WARNING! 2-loop threshold corrections could not be switched on, skip threshold correction factors!\n");
         }

         // Set MuR and MuF scale factors for pQCD + THC cross sections and test availability
         lthcvar = SetOn ? fnlo->SetScaleFactorsMuRMuF(xmur[iscls], xmuf[iscls]) : SetOn;
         if (!lthcvar) {
            warn["fnlo-read"] << "The selected scale variation (xmur, xmuf) = ("
                              << fnlo->GetScaleFactorMuR() << ","
                              << fnlo->GetScaleFactorMuF() << ") is not possible with this table, skip threshold correction factors!" << endl;
            //            printf("fnlo-read: WARNING! The selected scale variation (xmur, xmuf) = (% #10.3f, % #10.3f) is not possible with this table, skip threshold correction factors!\n",fnlo->GetScaleFactorMuR(),fnlo->GetScaleFactorMuF());
            // skip this part, check on lthcvar later for proper, i.e. no printout
         } else {
            fnlo->CalcCrossSection();
            xsthc2 = fnlo->GetCrossSection();
            kthc2  = fnlo->GetKFactors();
            for (unsigned int i=0; i<xsnlo.size(); i++) {
               if (abs(xsnlo[i]) > DBL_MIN) {
                  kthc2[i] = xsthc2[i]/xsnlo[i];
               } else {
                  kthc2[i] = -1.;
               }
            }
         }
      } else if ( !(ilo < 0 || ithc1 < 0)) {
         if ( !(inlo < 0) ) fnlo->SetContributionON(kFixedOrder, inlo, false);
         bool SetOn = fnlo->SetContributionON(kThresholdCorrection, ithc1, true);
         if (!SetOn) {
            warn["fnlo-read"] << "1-loop threshold corrections could not be switched on, skip threshold correction factors!" << endl;
            //            printf("fnlo-read: WARNING! 1-loop threshold corrections could not be switched on, skip threshold correction factors!\n");
         }

         // Set MuR and MuF scale factors for pQCD + THC cross sections and test availability
         lthcvar = SetOn ? fnlo->SetScaleFactorsMuRMuF(xmur[iscls], xmuf[iscls]) : SetOn;
         if (!lthcvar) {
            warn["fnlo-read"] << "The selected scale variation (xmur, xmuf) = ("
                              << fnlo->GetScaleFactorMuR() << ","
                              << fnlo->GetScaleFactorMuF() << ") is not possible with this table, skip threshold correction factors!" << endl;
            //            printf("fnlo-read: WARNING! The selected scale variation (xmur, xmuf) = (% #10.3f, % #10.3f) is not possible with this table, skip threshold correction factors!\n",fnlo->GetScaleFactorMuR(),fnlo->GetScaleFactorMuF());
            // skip this part, check on lthcvar later for proper, i.e. no printout
         } else {
            fnlo->CalcCrossSection();
            xsthc1 = fnlo->GetCrossSection();
            kthc1  = fnlo->GetKFactors();
            for (unsigned int i=0; i<xslo.size(); i++) {
               if (abs(xslo[i]) > DBL_MIN) {
                  kthc1[i] = xsthc1[i]/xslo[i];
               } else {
                  kthc1[i] = -1.;
               }
            }
         }
      }
      // Get non-perturbative corrections
      if ( !(inpc1 < 0) ) {
         if ( !(inlo < 0) ) {
            bool SetOn = fnlo->SetContributionON(kFixedOrder, inlo, true);
            if (!SetOn) {
               error["fnlo-read"] << "NLO not found, nothing to be done!" << endl;
               error["fnlo-read"] << "This should have been caught before!" << endl;
               exit(1);
            }
         }
         if ( !(ithc1 < 0) ) fnlo->SetContributionON(kThresholdCorrection, ithc1, false);
         if ( !(ithc2 < 0) ) fnlo->SetContributionON(kThresholdCorrection, ithc2, false);
         bool SetOn = fnlo->SetContributionON(kNonPerturbativeCorrection, inpc1, true);
         if (!SetOn) {
            error["fnlo-read"] << "NPC1 not found, nothing to be done!" << endl;
            error["fnlo-read"] << "This should have been caught before!" << endl;
            exit(1);
         }
         fnlo->CalcCrossSection();
         xsnpc1 = fnlo->GetCrossSection();
         knpc1  = fnlo->GetKFactors();
         for (unsigned int i=0; i<kfac.size(); i++) {
            if (abs(kfac[i]) > DBL_MIN) {
               knpc1[i] = knpc1[i]/kfac[i];
            } else {
               knpc1[i] = -1.;
            }
         }
      }

      // Start print out
      cout << DSEP << endl;
      printf(" My Cross Sections\n");
      printf(" The scale factors xmur, xmuf chosen here are: % #10.3f, % #10.3f\n",fnlo->GetScaleFactorMuR(),fnlo->GetScaleFactorMuF());
      cout << SSEP << endl;

      // Get table constants relevant for print out
      const int NDim = fnlo->GetNumDiffBin();
      unsigned int NDimBins[NDim];
      vector < string > DimLabel = fnlo->GetDimLabels();
      const int NObsBin = fnlo->GetNObsBin();
      vector < vector < double > > LoBin(NObsBin);
      vector < vector < double > > UpBin(NObsBin);
      for (int i=0; i<NObsBin; i++) {
         LoBin[i].resize(2);
         UpBin[i].resize(2);
         for (int j=0; j<2; j++) {
            LoBin[i][j]= fnlo->GetLoBin(i,j);
            UpBin[i][j]= fnlo->GetUpBin(i,j);
         }
      }
      vector < double > BinSize = fnlo->GetBinSize();

      // Print
      string header0 = "  IObs  Bin Size IODimO ";
      string header1 = " IODimI ";
      string header2 = " LO cross section";
      if (inlo>-1) {
         header2 += "   NLO cross section   KNLO";
      }
      if (ithc2>-1 && lthcvar) {
         header2 += "      KTHC2";
      } else if (ithc1>-1 && lthcvar) {
         if (inlo>-1) {
            header2 += "      KTHC1";
         } else {
            header2 += "    KTHC1";
         }
      }
      if (inpc1>-1) {
         header2 += "     KNPC1";
      }
      if (NDim == 1) {
         printf("%s [ %-17s ]  <%-12.12s> %s\n",
                // TODO: ScaleDescription Murks
                header0.c_str(),DimLabel[0].c_str(),"pT_jet_[GeV]",header2.c_str());
         cout << SSEP << endl;
         NDimBins[0] = 0;
         for (unsigned int i=0; i<xslo.size(); i++) {
            NDimBins[0]++;
            if (ilo > -1 && inlo > -1 && ithc2 > -1 && lthcvar && inpc1 > -1 ) {
               printf(" %5.i % -#10.4g %5.i  % -#10.4g  % -#10.4g % -#10.4g     %#18.11E %#18.11E %#9.5F %#9.5F %#9.5F",
                      i+1,BinSize[i],NDimBins[0],LoBin[i][0],UpBin[i][0],
                      qscl[i],xslo[i],xsnlo[i],kfac[i],kthc2[i],knpc1[i]);
            } else if (ilo > -1 && inlo > -1 && ithc2 > -1 && lthcvar) {
               printf(" %5.i % -#10.4g %5.i  % -#10.4g  % -#10.4g % -#10.4g     %#18.11E %#18.11E %#9.5F %#9.5F",
                      i+1,BinSize[i],NDimBins[0],LoBin[i][0],UpBin[i][0],
                      qscl[i],xslo[i],xsnlo[i],kfac[i],kthc2[i]);
            } else if (ilo > -1 && inlo > -1 && inpc1 > -1) {
               printf(" %5.i % -#10.4g %5.i  % -#10.4g  % -#10.4g % -#10.4g     %#18.11E %#18.11E %#9.5F %#9.5F",
                      i+1,BinSize[i],NDimBins[0],LoBin[i][0],UpBin[i][0],
                      qscl[i],xslo[i],xsnlo[i],kfac[i],knpc1[i]);
            } else if (ilo > -1 && inlo > -1 && ithc1 > -1 && lthcvar) {
               printf(" %5.i % -#10.4g %5.i  % -#10.4g  % -#10.4g % -#10.4g     %#18.11E %#18.11E %#9.5F %#9.5F",
                      i+1,BinSize[i],NDimBins[0],LoBin[i][0],UpBin[i][0],
                      qscl[i],xslo[i],xsnlo[i],kfac[i],kthc1[i]);
            } else if (ilo > -1 && inlo > -1) {
               printf(" %5.i % -#10.4g %5.i  % -#10.4g  % -#10.4g % -#10.4g     %#18.11E %#18.11E %#9.5F",
                      i+1,BinSize[i],NDimBins[0],LoBin[i][0],UpBin[i][0],
                      qscl[i],xslo[i],xsnlo[i],kfac[i]);
            } else if (ilo > -1 && ithc1 > -1 && lthcvar) {
               printf(" %5.i % -#10.4g %5.i  % -#10.4g  % -#10.4g % -#10.4g     %#18.11E %#9.5F",
                      i+1,BinSize[i],NDimBins[0],LoBin[i][0],UpBin[i][0],
                      qscl[i],xslo[i],kthc1[i]);
            } else if (ilo > -1) {
               printf(" %5.i % -#10.4g %5.i  % -#10.4g  % -#10.4g % -#10.4g     %#18.11E",
                      i+1,BinSize[i],NDimBins[0],LoBin[i][0],UpBin[i][0],
                      qscl[i],xslo[i]);
            } else {
               printf("fnlo-read: Nothing to report!\n");
               continue;
            }
            printf("\n");
         }
      } else if (NDim == 2) {
         printf("%s [ %-17s ] %s [ %-17s ]  <%-12.12s> %s\n",
                // header0.c_str(),DimLabel[0].c_str(),header1.c_str(),DimLabel[1].c_str(),fnlo->GetScaleDescription(0).c_str(),header2.c_str());
                // TODO: ScaleDescription Murks
                // Invert dimension numbering to from outer to inner
                header0.c_str(),DimLabel[0].c_str(),header1.c_str(),DimLabel[1].c_str(),"pT_jet_[GeV]",header2.c_str());
         cout << SSEP << endl;
         for (unsigned int i=0; i<xslo.size(); i++) {
            for (int j=0; j<NDim; j++) {
               if (i==0) {
                  NDimBins[j] = 1;
               } else if (LoBin[i-1][j] < LoBin[i][j]) {
                  NDimBins[j]++;
               } else if (LoBin[i][j] < LoBin[i-1][j]) {
                  NDimBins[j] = 1;
               }
            }
            if (ilo > -1 && inlo > -1 && ithc2 > -1 && lthcvar && inpc1 > -1 ) {
               printf(" %5.i % -#10.4g %5.i  % -#10.4g  % -#10.4g  %5.i  % -#10.4g  % -#10.4g % -#10.4g     %#18.11E %#18.11E %#9.5F %#9.5F %#9.5F",
                      i+1,BinSize[i],NDimBins[0],LoBin[i][0],UpBin[i][0],
                      NDimBins[1],LoBin[i][1],UpBin[i][1],qscl[i],xslo[i],xsnlo[i],kfac[i],kthc2[i],knpc1[i]);
            } else if (ilo > -1 && inlo > -1 && ithc2 > -1 && lthcvar) {
               printf(" %5.i % -#10.4g %5.i  % -#10.4g  % -#10.4g  %5.i  % -#10.4g  % -#10.4g % -#10.4g     %#18.11E %#18.11E %#9.5F %#9.5F",
                      i+1,BinSize[i],NDimBins[0],LoBin[i][0],UpBin[i][0],
                      NDimBins[1],LoBin[i][1],UpBin[i][1],qscl[i],xslo[i],xsnlo[i],kfac[i],kthc2[i]);
            } else if (ilo > -1 && inlo > -1 && inpc1 > -1) {
               printf(" %5.i % -#10.4g %5.i  % -#10.4g  % -#10.4g  %5.i  % -#10.4g  % -#10.4g % -#10.4g     %#18.11E %#18.11E %#9.5F %#9.5F",
                      i+1,BinSize[i],NDimBins[0],LoBin[i][0],UpBin[i][0],
                      NDimBins[1],LoBin[i][1],UpBin[i][1],qscl[i],xslo[i],xsnlo[i],kfac[i],knpc1[i]);
            } else if (ilo > -1 && inlo > -1 && ithc1 > -1 && lthcvar) {
               printf(" %5.i % -#10.4g %5.i  % -#10.4g  % -#10.4g  %5.i  % -#10.4g  % -#10.4g % -#10.4g     %#18.11E %#18.11E %#9.5F %#9.5F",
                      i+1,BinSize[i],NDimBins[0],LoBin[i][0],UpBin[i][0],
                      NDimBins[1],LoBin[i][1],UpBin[i][1],qscl[i],xslo[i],xsnlo[i],kfac[i],kthc1[i]);
            } else if (ilo > -1 && inlo > -1) {
               printf(" %5.i % -#10.4g %5.i  % -#10.4g  % -#10.4g  %5.i  % -#10.4g  % -#10.4g % -#10.4g     %#18.11E %#18.11E %#9.5F",
                      i+1,BinSize[i],NDimBins[0],LoBin[i][0],UpBin[i][0],
                      NDimBins[1],LoBin[i][1],UpBin[i][1],qscl[i],xslo[i],xsnlo[i],kfac[i]);
            } else if (ilo > -1 && ithc1 > -1 && lthcvar) {
               printf(" %5.i % -#10.4g %5.i  % -#10.4g  % -#10.4g  %5.i  % -#10.4g  % -#10.4g % -#10.4g     %#18.11E %#9.5F",
                      i+1,BinSize[i],NDimBins[0],LoBin[i][0],UpBin[i][0],
                      NDimBins[1],LoBin[i][1],UpBin[i][1],qscl[i],xslo[i],kthc1[i]);
            } else if (ilo > -1) {
               printf(" %5.i % -#10.4g %5.i  % -#10.4g  % -#10.4g  %5.i  % -#10.4g  % -#10.4g % -#10.4g     %#18.11E",
                      i+1,BinSize[i],NDimBins[0],LoBin[i][0],UpBin[i][0],
                      NDimBins[1],LoBin[i][1],UpBin[i][1],qscl[i],xslo[i]);
            } else {
               printf("fnlo-read: Nothing to report!\n");
               continue;
            }
            printf("\n");
         }
      } else {
         printf("fnlo-read: WARNING! Print out optimized for up to two dimensions. No output for %1.i dimensions.\n",NDim);
      }
   }

   // Print data if available, checks on availability internally
   //   fnlo->PrintCrossSectionsData();

   return 0;
}



//__________________________________________________________________________________________________________________________________


double Function_Mu(double s1, double s2) {
   // --- fastNLO user: This is an example function
   //     to demonstrate how you might perform the
   //     definition of the scales using a
   //     'flexible-scale'-table
   double mu = s1*exp(0.3*s2);
   return mu;
}

//__________________________________________________________________________________________________________________________________

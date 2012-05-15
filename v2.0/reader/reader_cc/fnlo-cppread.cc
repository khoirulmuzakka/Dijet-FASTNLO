//********************************************************************
//     
//     fastNLO_reader: FNLOCPPREAD
//     Program to read fastNLO v2 tables and derive
//     QCD cross sections using PDFs from LHAPDF
//     
//     D. Britzger, K. Rabbertz
//
//     Contains:
//
//********************************************************************
#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include <cstdlib>
#include <cfloat>
#include "FastNLOUser.h"
#include "FastNLODiffUser.h"
#include "Alphas.h"

// Function prototype for flexible-scale function 
double Function_Mu(double s1, double s2 );

//__________________________________________________________________________________________________________________________________

int fnlocppread(int argc, char** argv){

  using namespace std;

  //---  Initialization for nice printing
  const string CSEPS = "##################################################################################\n";
  const string LSEPS = "#---------------------------------------------------------------------------------\n";
  const string CSEPL = "####################################################################################################################################################################\n";

  //---  Parse commmand line
  printf("\n");
  printf(" %s",CSEPS.c_str());
  printf(" # fnlo-read: Program Steering\n");
  printf(" %s",LSEPS.c_str());
  string tablename = "table.tab";
  if ( argc <= 1 ){
    printf(" # fnlo-read: WARNING! No table name given,\n");
    printf(" # taking the default table.tab instead!\n");
    printf(" #   For an explanation of command line arguments type:\n");
    printf(" #   ./fnlo-cppread -h\n");
  } else {
    tablename = (const char*) argv[1];
    if (tablename == "-h"){ 
      printf(" #\n");
      printf(" # Usage: ./fnlo-cppread [arguments]\n");
      printf(" # Table input file, def. = table.tab\n");
      printf(" # PDF set, def. = cteq6mE.LHgrid\n");
      printf(" # Number of mu_r, mu_f scale settings to ");
      printf("investigate, if possible, def. = 1, max. = 7\n");
      printf(" #\n");
      printf(" # Give full path(s) if these are not in the cwd.\n");
      printf(" # Use \"_\" to skip changing a default argument.\n");
      printf(" #\n");
      return 0;
    } else if (tablename == "_") {
      tablename = "table.tab";
      printf("\n # fnlo-read: WARNING! No table name given,\n");
      printf(" # taking the default table.tab instead!\n");
    } else {
      cout << " # fnlo-read: Evaluating table: " << tablename << endl;
    }
  }

  //---  PDF set
  string PDFFile = "X";
  if ( argc > 2 ){
    PDFFile = (char*) argv[2];
  }
  if ( argc <= 2 || PDFFile == "_"){
    PDFFile = "cteq6m.LHpdf";
    printf(" # fnlo-read: WARNING! No PDF set given,\n");
    printf(" # taking cteq6mE.LHgrid instead!\n");
  } else {
    cout << " # fnlo-read: Using PDF set   : " << PDFFile << endl;
  }

  //--- Number of scale settings
  unsigned int nscls = 1;
  const unsigned int nsclmax = 7;
  const double xmur[] = { 1.0, 0.5, 2.0, 0.5, 1.0, 1.0, 2.0 }; 
  const double xmuf[] = { 1.0, 0.5, 2.0, 1.0, 0.5, 2.0, 1.0 }; 
  string ch2tmp = "X";
  if ( argc > 3 ){
    ch2tmp = (char*) argv[3];
  }
  if ( argc <= 2 || ch2tmp == "_" ){
    printf(" # fnlo-read: No request given for number of scale settings,\n");
    printf(" #            investigating primary scale only.\n");
  } else {
    nscls = atoi(argv[3]);
    if ( nscls < 1 ) {
      printf(" # fnlo-read: ERROR! No scale setting or even less??? Aborting! nscls = %i\n",nscls);
      exit(1);
    } else if ( nscls > nsclmax ){
      printf(" # fnlo-read: ERROR! Too many scale settings requested, aborting! nscls = %i\n",nscls);
      exit(1);
    } else {
      printf(" # fnlo-read: If possible, will try to do %i scale setting(s).\n",nscls);
    }
  }

  //---  Too many arguments
  if ( argc > 4 ){
    printf("fnlo-read: ERROR! Too many arguments, aborting!\n");
    return 1;
  }
  printf(" %s",CSEPS.c_str());
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
  //     which have encoded an advanced storage of matrix elements and scale variables. 
  //     These tables give you the possibility to change in addition to the renormalization
  //     also the factorization scale by arbitrary factors and have the possiblity to
  //     change the formula according to which the scale is derived.
  //
  //     Please check, which type of table you are using and then refer to the comments and
  //     functions suitable for this fastNLO table.
  


  // ------- Initialize fastNLOReader ------- //
  // --- fastNLO user: Make an instance of your class that derives
  //     from the FastNLOReader class and
  //     pass the name of the fastNLO table as an argument.
  //
  //        FastNLOUser* fnloreader = new FastNLOUser( tablename );
  //
  //     The example class for LHAPDF has overwriten the constructor
  //     and thus takes also the PDF filename (and PDF set).
  //     TBD: ??? C++ chinese ???
  FastNLOUser* fnloreader = new FastNLOUser( tablename , PDFFile , 0 );
  


  // ------- Select a PDF set and member ------- //
  // --- fastNLO user: You can select the PDF here.
  //     With LHAPDF, you can set the PDF set and member using e.g.:
  //           fnloreader->SetLHAPDFfilename( PDFFile );
  //           fnloreader->SetLHAPDFset( int PDFSet );
  //     Or you talk directly to LHAPDF using LHAPDF::LHAPDF().
  //
  //     >>>>   WARNING  <<<< 
  //     fastNLO communictes to LHAPDF using the LHAPDFwrap class.
  //     This class is a global singleton class, which means, that if
  //     it is e.g. changed to another PDF-name (e.g. by a second
  //     instance of FastNLOReader), then all FastNLOReader 
  //     instances access this 'new' LHAPDF file or member set (of course
  //     only after calling FillPDFCache() and CalcCrossSection()).
  //     
  //     You can print some information about the currently initialized
  //     LHAPDF using
  //		fnloreader->PrintCurrentLHAPDFInformation();
  //	
  //	 If you are not really sure, which pdf set is currently used, then
  //     please reset the LHAPDF filename and memberset id.
  //
  //     Since the PDF evolution is typically provided in some external
  //     code, you can choose the interface to the evolution routine.
  //     By default LHAPDF is used:
  //           fnloreader->SetPDFInterface(FastNLOReader::kLHAPDF);
  //	
  //   TBD: What about the PDFFile above? What to use?
  //   TBD: Which alternatives for SetPDFInterface(FastNLOReader::kLHAPDF) ?
  //
  //   fnloreader->SetLHAPDFfilename( PDFFile );
  //   fnloreader->SetLHAPDFset( 0 );
  //  TBD: The following doesn't work ...???
  //  fnloreader->PrintCurrentLHAPDFInformation();
  
  
  
  // --- fastNLO user: After defining the PDF or changing the factorization
  //     scale you always have to refill the fastNLO internal PDF cache
  //     calling:
  //            fnloreader->FillPDFCache();
  //
  fnloreader->FillPDFCache();



  // ------- Setting Alpha_s value ------- //
  // --- fastNLO user: With fastNLO, the user can choose the value of alpha_s(Mz).
  //     Furthermore, there are multiple options for the evolution
  //     code to calculate alpha_s(mu_r) according to the RGE.
  //     To set the value of alpha_s(Mz) at the Z0 mass, please use:
  //            fnloreader->SetAlphasMz(0.1185);
  //     You always have to specify an alpha_s(Mz) value.
  //    
  //     Further you can specify the desired evolution code using
  //            fnloreader->SetAlphasEvolution(FastNLOReader::EAlphasEvolution);
  //
  //     The following options can be chosen for the evolution:
  //      - kGRV		The default option in fastNLO.
  //                            Here you can select 2-,3- or 4-loop iterative solutions of the RGE and
  //                            set any value for M_Z, the number of flavors, the flavor matching, etc...
  //                            Please see the documentation in Alphas.h and Alphas.cc.
  //                            For usage you can access the static class directly like e.g.:
  //                                Alphas::SetNLoop(2);
  //                                Alphas::SetMz(91.1786);
  //				INFO: When setting kGRV, you initialize Alphas::Alphas() with
  //				some reasonable values and overwrite previous settings in this class!
  //      - kNLOJET		This is the alpha_s evolution, which is used by NLOJet++ as default
  //      - kCTEQpdf		This alpha_s evolution is used in the CTEQ6 PDFs. [Sure?]
  //      - kLHAPDFInternal	With this option, you access the alpha_s evolution, which is defined
  //                            within the LHAPDF-file. You cannot change alphas(Mz) here!
  //      - kQCDNUMInternal	Using kQCDNUM as PDF evolution code, you can make use of the
  //				QCDNUM alpha_s evolution. Please see the QCDNUM manual for options.
  //      - kFixed		Take always the fixed value of SetAlphasMz() without any evolution code.
  //
  //     INFO: The choice of the 'number of flavors' in your alpha_s evolution might be 
  //     inconsistent with the number of flavors that is used for calculating the matrix elements.
  //     NLOJet++ typically uses 5 massless flavors.
  //
  fnloreader->SetAlphasEvolution(FastNLOReader::kGRV);
  //fnloreader->SetAlphasEvolution(FastNLOReader::kNLOJET); 
  //fnloreader->SetAlphasEvolution(FastNLOReader::kLHAPDFInternal); 
  //fnloreader->SetAlphasEvolution(FastNLOReader::kCTEQpdf);
  //fnloreader->SetAlphasMz(0.1168);
  //fnloreader->SetAlphasEvolution(FastNLOReader::kNLOJET);
  //  fnloreader->SetAlphasMz(0.1179);
  // TBD: Use again value in released reader code for comparison
  fnloreader->SetAlphasMz(0.1185);
  
  // TBD: What about alpha_s cache? Fill here also?



  // ------- Set the units of your calculation (kPublicationUnits or kAbsoluteUnits) ------- //
  // --- fastNLO user: You can choose the units in which you want
  //     to access (or print) your cross-section results.
  //     There are two possibilites:
  //       - The default option is 'publication units', i.e. divided by 
  //         bin widths if done so in the relevant publication
  //            fnloreader->SetUnits(FastNLOReader::kPublicationUnits);
  //       - The other option is 'absolute' units in barn, but still in
  //         the same magnitude as in the publication (e.g. pb, fb, nb, etc.)
  //            fnloreader->SetUnits(FastNLOReader::kAbsoluteUnits);
  fnloreader->SetUnits(FastNLOReader::kPublicationUnits);
  //fnloreader->SetUnits(FastNLOReader::kAbsoluteUnits);



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
  //       - kFixedOrder		  -> Fixed order calculation (in alpha_s)
  //       - kThresholdCorrection	  -> Threshold corrections
  //       - kElectroWeakCorrection	  -> Electroweak corrections (not derived yet)
  //       - kNonPerturbativeCorrections  -> Non-perturbative corrections|Hadronisation corrections
  //     plus one must know the 'Id' of this contribution, which can be printed e.g.
  //     by calling      fnloreader->PrintTableInfo();
  //
  //     To switch a contribution on/off please use:
  //            fnloreader->SetContributionON( contrib, Id, on/off ) 
  //     Here, 'contrib' is not the contribution number, but the type which
  //     is encoded in an enum list with names as given above: kFixedOrder, ...
  //     Within each type the contributions are counted separately starting with Id=0.
  //     The total number of contributions then counts all contributions of all types.
  fnloreader->PrintTableInfo();
  


  if ( !fnloreader->GetIsFlexibleScaleTable() ) {
     // ---- options for scales in v2.0 tables ---- //
     // --- fastNLO user: Here you can specify the options if you use a v2.0 fastNLO table
     //     (normal case in pp/ppbar) but NOT a 'flexible-scale' table (v2.1):
     //     A fastNLO table comes usually with various precalculated scale variations.
     //     There, the renormalization and the factorization scale are varied simultaneously
     //     by a series of predefined factors, usually 1, 2, 1/2, 1/4.
     //     The welcome message should show you the information about the available
     //     scale variations. For accessing a table with a certain factor you need to
     //     know the scale Id (int) and then use:
     //            fnloreader->SetScaleVariation(0);
     //
     //     Furthermore you have the possiblity to vary the renormalization scale by any
     //     factor. For a scale variation factor of 0.5, please use e.g.:
     //            fnloreader->SetScaleFactorMuR(0.5);
     //     WARNING: This option is inconsistent with threshold corrections. These require
     //     in the present implementation to always have mu_r = mu_f.
  }  

  
  if ( fnloreader->GetIsFlexibleScaleTable() ) {
     // ---- options for scales in 'flexible-scale' tables (v2.1) ---- //
     // --- fastNLO user: You can choose a function to define how
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
     //            fnloreader->SetMuRFunctionalForm(FastNLOReader::EScaleFunctionalForm);
     //         for changing the calculation of the renormalizatoin scale.
     //         Please refer to FastNLOReader.h for all options of EScaleFunctionalForm.
     //
     //       - or you can pass a function pointer to FastNLOReader using
     //            fnloreader->SetExternalFuncForMuR( double (*Func)(double,double) );
     //         to pass any function using scale1 and scale2 to fastNLO.
     //    
     //     Further you can define a scale factor (one for the renormalization
     //     and one for the factorization scale) which is multiplied to the scale
     //     independently from the predefined function, like e.g. for mu_r:
     //         mu_r = scalefac * f(scale1,scale2)
     //     Please use e.g. a factor of 1.5 for mu_r:
     //             fnloreader->SetScaleFactorMuR(1.5);
     //
     //     For changing the factorization scale, replace all 'MuR' by 'MuF' in the function calls.
     //  
     //     INFO: All above-mentioned functions automatically perform a refilling of the
     //     fastNLO internal PDF-cache. To switch it off you can use a boolean, like:
     //             fnloreader->SetMuFFunctionalForm(FastNLOReader::kScale1 , false );
     //
     //     WARNING: Some choice had to be made for the default settings. Please think
     //     carefully about the choice of the scales ...
     //     Default setting for DIS tables:
     //       - mu_r:  kQuadraticMean	-> mu_r = sqrt( (Q^2 + scale2^2)/2. ) // because scale1=Q!
     //       - mu_f:  kScale1		-> mu_f = Q
     //     Default setting for pp and ppbar tables:
     //       - mu_r:  kScale1		-> mu_r = scale1
     //       - mu_f:  kScale1		-> mu_f = scale1
     //
     //  Valid calls are e.g.:
     //     fnloreader->SetMuFFunctionalForm(FastNLOReader::kScale1);	 // set function how to calculate mu_f from scale1 and scale2
     //     fnloreader->SetMuRFunctionalForm(FastNLOReader::kScale1);	 // set function how to calculate mu_f from scale1 and scale2
     //     fnloreader->SetMuRFunctionalForm(FastNLOReader::kQuadraticMean); // set function how to calculate mu_r from scale1 and scale2
     //     fnloreader->SetMuFFunctionalForm(FastNLOReader::kScale1);	 // set function how to calculate mu_f from scale1 and scale2
     //     fnloreader->SetExternalFuncForMuR( &Function_Mu );		 // set external function to calculate mu_r from scale1 and scale2
     //     fnloreader->SetMuFFunctionalForm(FastNLOReader::kExpProd2);	 // set function how to calculate mu_f from scale1 and scale2
     //     fnloreader->SetMuRFunctionalForm(FastNLOReader::kExpProd2);	 // set function how to calculate mu_f from scale1 and scale2
     fnloreader->SetMuFFunctionalForm(FastNLOReader::kScale1);	 // set function how to calculate mu_f from scale1 and scale2
     fnloreader->SetMuRFunctionalForm(FastNLOReader::kScale1);	 // set function how to calculate mu_f from scale1 and scale2

 

     // ---- (Re-)calculate cross sections ---- //
     // --- fastNLO user: Before you can access the fastNLO computed
     //     cross sections, you always have to call CalcCrossSection()!
     //     If you are not sure, whether you have recalculated the internal
     //     PDF cache with your current scale choice you further can call 
     //     FillPDFCache() before calling CalcCrossSection().
     //     So, before accessing the cross sections, please call:
     //             fnloreader->CalcCrossSection();
  
  }


  
  // ---- fastNLO user: choice of scale factor
  //    fnloreader->SetScaleFactorMuR(1.5);				 // set scale factor for mu_r
  //    fnloreader->SetScaleFactorMuF(0.66);				 // set scale factor for mu_f



  // ---- Get cross sections ---- //
  // --- fastNLO user: To access the cross section from fastNLO
  //     you should use:
  //           vector < double > xs = fnloreader->GetCrossSection();
  //     If you want to have a pointer to an array of numbers you might use
  //           vector < double > xs = fnloreader->GetCrossSection();
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
  //           vector < double > kFactors = fnloreader->GetKFactors();



  // ---- Printing ---- //
  // --- fastNLO user: For an easy overview of your cross section calculation
  //     you might use the following print methods:
  //             fnloreader->PrintCrossSections();
  //
  //     Or print it (almost exactly) like the Fortran reader code:
  //             fnloreader->PrintCrossSectionsDefault();
  //     WARNING: The latter call makes implicit changes to the settings to
  //              reproduce the 'default' results. Call it BEFORE you make
  //              your settings in case.



  // ---- Information ---- //
  // --- fastNLO user: For a comprehensive insight into the fastNLO variables
  //     you can use:
  //             fnloreader->PrintFastNLOTableConstants(0);
  //     
  //     For a comparision with a Reference cross section calculated with
  //     NLOJet++ you might use:
  //             fnloreader->PrintCrossSectionsWithReference();
  //     WARNING: 
  //     In v2.0 tables there are not necessarily any reference cross sections stored.
  //     Usually, these would come as extra tables.
  //     In v2.1 there are always three reference tables stored. Please check
  //     which one you access with this method. These reference cross sections also might
  //     be empty. The three references are calculated usually with different scales,
  //     which might be very different from your choice!
  //     Therefore this comparison DOES not provide any information about
  //     the precision of this fastNLO table.
  //
  //     INFO: NLOJet++ typically uses the NLOJet++-like alpha_s evolution with a value of 0.1179
  //     and a PDF similar to the cteq6m.LHgrid pdf-file.
  // ************************************************************************************************
  fnloreader->PrintFastNLOTableConstants(0);
  
  // The presented example is done automatically for print out here  
  //  fnloreader->PrintCrossSectionsDefault(0);
  
  // Example code to print out data points (if available)
  //  fnloreader->PrintCrossSectionsData();
  
  // Example code to access cross sections and K factors:
  //  fnloreader->PrintFastNLODemo();
  


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

  // Check on existence of LO and NLO (Id = -1 if not existing)
  int ilo   = fnloreader->ContrId(FastNLOReader::kFixedOrder, FastNLOReader::kLeading); 
  int inlo  = fnloreader->ContrId(FastNLOReader::kFixedOrder, FastNLOReader::kNextToLeading);
  if ( ilo < 0 || inlo < 0 ){
    printf("fnlo-read: ERROR! LO and/or NLO not found, nothing to be done!\n");
    exit(1);
    //  } else {
    //    printf("fnlo-read: LO and NLO contributions have Id's: %i and %i\n",ilo,inlo);
  }
  // Check on existence of 2-loop threshold corrections
  int ithc2 = fnloreader->ContrId(FastNLOReader::kThresholdCorrection, FastNLOReader::kNextToLeading);
  // if ( ithc2 < 0 ){
  //   printf("fnlo-read: 2-loop threshold corrections not found!\n");
  // } else {
  //   printf("fnlo-read: 2-loop threshold corrections have Id: %i\n",ithc2);
  // }
  // Check on existence of non-perturbative corrections from LO MC
  int inpc1 = fnloreader->ContrId(FastNLOReader::kNonPerturbativeCorrection, FastNLOReader::kLeading);
  // if ( inpc1 < 0 ){
  //   printf("fnlo-read: Non-perturbative corrections not found!\n");
  // } else {
  //   printf("fnlo-read: Non-perturbative corrections have Id: %i\n",inpc1);
  // }

  // Switch on LO & NLO, switch off anything else (to be verbose about it set last entry to "true")
  if (! (ilo   < 0)) {fnloreader->SetContributionON( FastNLOReader::kFixedOrder, 0, true, false );} 
  if (! (inlo  < 0)) {fnloreader->SetContributionON( FastNLOReader::kFixedOrder, 1, true, false );}
  if (! (ithc2 < 0)) {fnloreader->SetContributionON( FastNLOReader::kThresholdCorrection, ithc2, false, false );}
  if (! (inpc1 < 0)) {fnloreader->SetContributionON( FastNLOReader::kNonPerturbativeCorrection, inpc1, false, false );}
  // Temporary: Also don't print the cross sections out even when existing for this example
  ithc2 = -1;
  inpc1 = -1;

  // Find scale variation for intended MuF scale factor value
  const int nsclsfmax = 10;
  double fxmuf[nsclsfmax];
  if ( !fnloreader->GetIsFlexibleScaleTable() ) {
    // Get number of available scale variations for selected set of contributions and
    // check on available scale factors, in particular for MuF
    int nsclsf = fnloreader->GetNScaleVariations();
    if ( nsclsf > nsclsfmax ) {
      printf("fnlo-read: WARNING! Found more scale variations than I can deal with: %i\n",nsclsf); 
      printf("           Using only the first 10 of them.\n");
      nsclsf = nsclsfmax;
    }
    // With active threshold corrections, only default scale (0) usable for the moment!
    for (int iscls=0; iscls<nsclsf; iscls++){
      fnloreader->SetScaleVariation(iscls);
      fxmuf[iscls] = fnloreader->GetScaleFactorMuF();
      //      printf("fnlo-read: MuF scale factor for scale variation no. %i is: %7.3f\n",iscls,fxmuf[iscls]);
    }
    
    // Run over all requested scale settings xmur, xmuf if possible
    for (int iscls=0; iscls<nscls; iscls++){
      int isclf = -1;
      for (int jscls=0; jscls<nsclsf; jscls++){
	if ( abs(xmuf[iscls]-fxmuf[jscls]) < DBL_MIN ){
	  isclf = jscls;
	  break;
	}
      }
      if ( isclf > -1 ) {
	// Select MuF scale variation
	//   Do not refill PDF cache (2nd arg. = false) --> will be done explicitly, be verbose (3rd arg. = true)
	fnloreader->SetScaleVariation(isclf, false, false);
	
	// Set MuR scale factor
	//   Do not refill PDF & alpha_s caches (2nd arg. = false) --> will be done explicitly, be verbose (3rd arg. = true)
	fnloreader->SetScaleFactorMuR(xmur[iscls], false, false);

	// This refers to flex scales ...
	// } else {
	//   // Set MuF scale factor, only usable with flexible-scale tables
	//   //   Do not refill PDF cache (2nd arg. = false) --> will be done explicitly, be verbose (3rd arg. = true)
	//   //  double fxmuf = 1.0; 
	//   //  fnloreader->SetScaleFactorMuF(fxmuf, false, true);
	// }

	// If the MuR scale was changed the alpha_s cache MUST be refilled
	fnloreader->FillAlphasCache();
	// If the PDF or a scale was changed the PDF cache MUST be refilled
	fnloreader->FillPDFCache();
  
	// Calculate cross section
	fnloreader->CalcCrossSection();

	// Get results
	vector < double > xsnlo = fnloreader->GetCrossSection();
	vector < double > kfac  = fnloreader->GetKFactors();
	vector < double > xslo  = xsnlo;
	for (unsigned int i=0;i<xslo.size();i++){
	  if ( abs(kfac[i]) > DBL_MIN ){
	    xslo[i] = xslo[i]/kfac[i];
	  } else {
	    xslo[i] = -1.;
	  }
	}
	vector < double > xsthc2;
	vector < double > kthc;
	vector < double > xsnpc;
	vector < double > knpc;
	
	// Start print out
	cout << DSEP << endl;
	printf(" My Cross Sections\n");
	printf(" The scale factors chosen here are: % #10.3f, % #10.3f\n",fnloreader->GetScaleFactorMuR(),fnloreader->GetScaleFactorMuF());
	cout << SSEP << endl;
    
	// Get table constants relevant for print out
	// TBD: This Getter should be renamed!!!
	int NDim = fnloreader->GetNDiffBin();
	unsigned int NDimBins[NDim];
	vector < string > DimLabel = fnloreader->GetDimensionLabel();
	vector < vector < double > > LoBin = fnloreader->GetLowBinEdge();
	vector < vector < double > > UpBin = fnloreader->GetUpBinEdge();
	vector < double > BinSize = fnloreader->GetBinSize();
	
	// Print
	if ( NDim == 2 ){
	  string header0 = "  IObs  Bin Size IODim1 "; 
	  string header1 = "   IODim2 ";
	  string header2 = " LO cross section   NLO cross section   K NLO";
	  if ( ithc2>-1 ){
	    header2 += "     K THC";
	  }
	  if ( inpc1>-1 ){
	    header2 += "     K NPC";
	  }
	  printf("%s [ %-12s ] %s [  %-12s  ] %s\n",
		 header0.c_str(),DimLabel[0].c_str(),header1.c_str(),DimLabel[1].c_str(),header2.c_str());
	  cout << SSEP << endl;
	  for ( unsigned int i=0; i<xslo.size(); i++ ){ 
	    for ( int j=0; j<NDim; j++ ){ 
	      if ( i==0 ){
		NDimBins[j] = 1;
	      } else if ( LoBin[i-1][j] < LoBin[i][j]){
		NDimBins[j]++;
	      } else if ( LoBin[i][j] < LoBin[i-1][j]){
		NDimBins[j] = 1;
	      }
	    }
	    if ( ithc2<0 && inpc1<0 ) {
	      printf(" %5.i % -#10.4g %5.i % -#10.4g % -#10.4g %5.i  %-#8.2E  %-#8.2E %#18.11E %#18.11E %#9.5F",
		     i+1,BinSize[i],NDimBins[0],LoBin[i][0],UpBin[i][0],
		     NDimBins[1],LoBin[i][1],UpBin[i][1],xslo[i],xsnlo[i],kfac[i]);
	    } else if ( inpc1<0 ) {
	      printf(" %5.i % -#10.4g %5.i % -#10.4g % -#10.4g %5.i  %-#8.2E  %-#8.2E %#18.11E %#18.11E %#9.5F %#9.5F",
		     i+1,BinSize[i],NDimBins[0],LoBin[i][0],UpBin[i][0],
		     NDimBins[1],LoBin[i][1],UpBin[i][1],xslo[i],xsnlo[i],kfac[i],kthc[i]);
	    } else if ( ithc2<0 ) {
	      printf(" %5.i % -#10.4g %5.i % -#10.4g % -#10.4g %5.i  %-#8.2E  %-#8.2E %#18.11E %#18.11E %#9.5F %#9.5F",
		     i+1,BinSize[i],NDimBins[0],LoBin[i][0],UpBin[i][0],
		     NDimBins[1],LoBin[i][1],UpBin[i][1],xslo[i],xsnlo[i],kfac[i],knpc[i]);
	    } else {
	      printf(" %5.i % -#10.4g %5.i % -#10.4g % -#10.4g %5.i  %-#8.2E  %-#8.2E %#18.11E %#18.11E %#9.5F %#9.5F %#9.5F",
		     i+1,BinSize[i],NDimBins[0],LoBin[i][0],UpBin[i][0],
		     NDimBins[1],LoBin[i][1],UpBin[i][1],xslo[i],xsnlo[i],kfac[i],kthc[i],knpc[i]);
	    }
	    printf("\n");
	  }
	} else {
	  printf("fnlo-read: WARNING! Print out optimized for two dimensions. No output for %1.i dimensions.\n",NDim);
	}
      }
    }
  }
  return 0;
}



//__________________________________________________________________________________________________________________________________


double Function_Mu(double s1, double s2 ){
  // --- fastNLO user: This is an example function
  //     to demonstrate how you might perform the
  //     definition of the scales using a 
  //     'flexible-scale'-table
  double mu = s1*exp(0.3*s2);
  return mu;
}

//__________________________________________________________________________________________________________________________________

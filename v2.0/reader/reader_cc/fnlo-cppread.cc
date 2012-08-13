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
#include "FastNLOLHAPDF.h"
#include "FastNLOAlphas.h"
#include "FastNLONLOJETLIKE.h"
#include "FastNLODiffUser.h"
#include "FastNLOCRunDec.h"
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
    PDFFile = (const char*) argv[2];
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
    ch2tmp = (const char*) argv[3];
  }
  if ( argc <= 3 || ch2tmp == "_" ){
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


  say::SetGlobalVerbosity(say::INFO);

  // ------- Initialize fastNLOReader ------- //
  // --- fastNLO user: Make an instance of your class that derives
  //     from the FastNLOReader class and
  //     pass the name of the fastNLO table as an argument.
  //
  //        FastNLOUser* fnloreader = new FastNLOUser( tablename );
  //
  //     The example class for LHAPDF has overwritten the constructor
  //     and thus takes also the PDF filename (and PDF set).
  //     TBD: ??? C++ chinese ???

  
  
  cout<<" ================ CRunDec ================== " << endl;
  FastNLOCRunDec* fnlocrundec = new FastNLOCRunDec( tablename , PDFFile , 0 );
  fnlocrundec->SetAlphasMz(0.1182);
  fnlocrundec->FillPDFCache();
  fnlocrundec->CalcCrossSection();
  fnlocrundec->PrintCrossSectionsDefault();

  cout<<" ================ LHAPDF ================== " << endl;
  FastNLOLHAPDF* fnlo = new FastNLOLHAPDF( tablename , PDFFile , 0 );
  fnlo->PrintCrossSectionsDefault();
  

  cout<<" ================ Alphas ================== " << endl;
  FastNLOAlphas* fnloreader = new FastNLOAlphas( tablename , PDFFile , 0 );
  fnloreader->SetAlphasMz(0.1179);
  fnloreader->PrintCrossSectionsDefault();

  cout<<" ================ NLOJETLIKE ================== " << endl;
  FastNLONLOJETLIKE* nloj = new FastNLONLOJETLIKE(tablename);
  nloj->SetMuFFunctionalForm(fastNLO::kQuadraticMean); // set function how to calculate mu_r from scale1 and scale2
  nloj->SetMuRFunctionalForm(fastNLO::kQuadraticMean); // set function how to calculate mu_r from scale1 and scale2
  nloj->FillPDFCache();
  nloj->CalcCrossSection();
  nloj->PrintCrossSectionsWithReference();


  // ------- Set the units of your calculation (kPublicationUnits or kAbsoluteUnits) ------- //
  // --- fastNLO user: You can choose the units in which you want
  //     to access (or print) your cross-section results.
  //     There are two possibilites:
  //       - The default option is 'publication units', i.e. divided by 
  //         bin widths if done so in the relevant publication
  //            fnloreader->SetUnits(fastNLO::kPublicationUnits);
  //       - The other option is 'absolute' units in barn, but still in
  //         the same magnitude as in the publication (e.g. pb, fb, nb, etc.)
  //fnloreader->SetUnits(fastNLO::kAbsoluteUnits);
  fnloreader->SetUnits(fastNLO::kPublicationUnits);
  
  
  
  // ------- Select a PDF set and member ------- //
  // --- fastNLO user: You can select the PDF here.
  //     With LHAPDF, you can set the PDF set and member using e.g.:
  //           fnloreader->SetLHAPDFfilename( PDFFile );
  //           fnloreader->SetLHAPDFset( int PDFSet );
  //     Or you talk directly to LHAPDF using LHAPDF::LHAPDF().
  //
  //     >>>>   WARNING  <<<< 
  //     fastNLO communicates to LHAPDF using the LHAPDFwrap class.
  //     This class is a global singleton class, which means, that if
  //     it is e.g. changed to another PDF name (e.g. by a second
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
  //           fnloreader->SetPDFInterface(fastNLO::kLHAPDF);
  //	
  //   TBD: What about the PDFFile above? What to use?
  //   TBD: Which alternatives for SetPDFInterface(fastNLO::kLHAPDF) ?
  //
  //   fnloreader->SetLHAPDFfilename( PDFFile );
  //   fnloreader->SetLHAPDFset( 0 );
  //  TBD: The following doesn't work ...???
  //  fnloreader->PrintCurrentLHAPDFInformation();


  
  // ------- Setting Alpha_s value ------- //
  // --- fastNLO user: 
  //     The alpha_s evolution is provided in the user code
  // TBD: Use again value in released reader code for comparison
  //fnloreader->SetAlphasMz(0.1185);



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
  
  
  
  // ------- Selecting the scale treatment ------- //
  // --- fastNLO user: The simplest way to modify the predefined renormalization and
  //     factorization scales is to provide a scale factor by which the default scale
  //     is multiplied. These factors must be positive and not too small (> 1.E-6).
  //     Otherwise they can in principal (within reason) be set arbitrarily for
  //     flexible-scale tables. For the normal v2 tables the choice of factors for the
  //     factorization scale is limited to some fixed values, usually 0.5, 1.0, and 2.0
  //     plus sometimes also 0.25, see the respective table information.
  //     Note: If threshold corrections are available and switched on for evaluation,
  //     the scale factors for the renormalization and factorization scale must be identical. 
  //
  //     The function call to set the scale factors is:
  //         fnloreader->SetScaleFactorsMuRMuF(xmur, xmuf, ReFillCache);
  //     where xmur, xmuf are the scale factors, and
  //     ReFillCache, by default true, should not be changed.
  //
  //     The return value of this function call is boolean and returns false, if the
  //     the requested scale factors can not be chosen. In this case, the last legal
  //     values remain unchanged.
  
  // ----- Additional possibilities for scales in 'flexible-scale' tables (v2.1) ----- //
  if ( fnloreader->GetIsFlexibleScaleTable() ) {
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
    //     WARNING: Some choice had to be made for the default settings. Please think
    //     carefully about the choice of the scales ...
    //     Default setting for DIS tables:
    //       - mu_r:  kQuadraticMean	-> mu_r = sqrt( (Q^2 + scale2^2)/2. ) // because scale1=Q!
    //       - mu_f:  kScale1		-> mu_f = Q
    //     Default setting for pp and ppbar tables:
    //       - mu_r:  kScale1		-> mu_r = scale1
    //       - mu_f:  kScale1		-> mu_f = scale1
    //
    //     Valid calls are e.g.:
    //     fnloreader->SetMuRFunctionalForm(fastNLO::kScale1);	 // set function how to calculate mu_r from scale1 and scale2
    //     fnloreader->SetMuFFunctionalForm(fastNLO::kScale1);	 // set function how to calculate mu_f from scale1 and scale2
    //     fnloreader->SetMuRFunctionalForm(fastNLO::kQuadraticMean); // set function how to calculate mu_r from scale1 and scale2
    //     fnloreader->SetMuFFunctionalForm(fastNLO::kScale1);	 // set function how to calculate mu_f from scale1 and scale2
    //     fnloreader->SetExternalFuncForMuR( &Function_Mu );		 // set external function to calculate mu_r from scale1 and scale2
    //     fnloreader->SetMuRFunctionalForm(fastNLO::kExpProd2);	 // set function how to calculate mu_f from scale1 and scale2
    //     fnloreader->SetMuFFunctionalForm(fastNLO::kExpProd2);	 // set function how to calculate mu_f from scale1 and scale2
  }
  // INFO: All above-mentioned scale changing functions automatically perform a refilling of the
  //       fastNLO internal PDF cache. To switch it off you can use a boolean, like:
  //       fnloreader->SetMuFFunctionalForm(fastNLO::kScale1 , false );
  
  
  
  // ---- (Re-)calculate cross sections ---- //
  // --- fastNLO user: Before you can access the fastNLO computed
  //     cross sections, you always have to call CalcCrossSection()!
  //     If you are not sure, whether you have recalculated the internal
  //     PDF cache with your current scale choice you further can call 
  //     FillPDFCache() before calling CalcCrossSection().
  //     So, before accessing the cross sections, please call:
  //             fnloreader->CalcCrossSection();



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
  int ilo   = fnloreader->ContrId(fastNLO::kFixedOrder, fastNLO::kLeading); 
  int inlo  = fnloreader->ContrId(fastNLO::kFixedOrder, fastNLO::kNextToLeading);
  if ( ilo < 0 || inlo < 0 ){
    printf("fnlo-read: ERROR! LO and/or NLO not found, nothing to be done!\n");
    exit(1);
    //  } else {
    //    printf("fnlo-read: LO and NLO contributions have Id's: %i and %i\n",ilo,inlo);
  }
  // Check on existence of 2-loop threshold corrections
  int ithc2 = fnloreader->ContrId(fastNLO::kThresholdCorrection, fastNLO::kNextToLeading);
  // if ( ithc2 < 0 ){
  //   printf("fnlo-read: 2-loop threshold corrections not found!\n");
  // } else {
  //   printf("fnlo-read: 2-loop threshold corrections have Id: %i\n",ithc2);
  // }
  // Check on existence of non-perturbative corrections from LO MC
  int inpc1 = fnloreader->ContrId(fastNLO::kNonPerturbativeCorrection, fastNLO::kLeading);
  // if ( inpc1 < 0 ){
  //   printf("fnlo-read: Non-perturbative corrections not found!\n");
  // } else {
  //   printf("fnlo-read: Non-perturbative corrections have Id: %i\n",inpc1);
  // }

  // Switch on LO & NLO, switch off anything else
  if (! (ilo   < 0)) {fnloreader->SetContributionON( fastNLO::kFixedOrder, 0, true );} 
  if (! (inlo  < 0)) {fnloreader->SetContributionON( fastNLO::kFixedOrder, 1, true );}
  if (! (ithc2 < 0)) {fnloreader->SetContributionON( fastNLO::kThresholdCorrection, ithc2, false );}
  if (! (inpc1 < 0)) {fnloreader->SetContributionON( fastNLO::kNonPerturbativeCorrection, inpc1, false );}
  // Temporary: Also don't print the cross sections out even when existing for this example
  ithc2 = -1;
  inpc1 = -1;

  // Run over all pre-defined scale settings xmur, xmuf
  for (int iscls=0; iscls<nscls; iscls++){
    // Set MuR and MuF scale factors
    bool lscvar = fnloreader->SetScaleFactorsMuRMuF(xmur[iscls], xmuf[iscls], true);
    if ( !lscvar ) {
      printf("fnlo-cppread: WARNING! The selected scale variation (xmur, xmuf) = (% #10.3f, % #10.3f) is not possible, skipped!\n");
      continue;
    }
    if ( fnloreader->GetIsFlexibleScaleTable() ) {
      fnloreader->SetMuFFunctionalForm(fastNLO::kScale1);
      fnloreader->SetMuRFunctionalForm(fastNLO::kScale1);
      //      fnloreader->SetMuFFunctionalForm(fastNLO::kScale2);
      //      fnloreader->SetMuRFunctionalForm(fastNLO::kScale2);
    }
    
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
    printf(" The scale factors xmur, xmuf chosen here are: % #10.3f, % #10.3f\n",fnloreader->GetScaleFactorMuR(),fnloreader->GetScaleFactorMuF());
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

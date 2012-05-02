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

  //---  Too many arguments
  if ( argc > 3 ){
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
  //     The second type of tables, v2.1, are called 'flexible-scale'-tables
  //     which have encoded an advanced storage of matrix elements and scale variables. 
  //     These tables give you the possibility to change the renormalization and
  //     the factorization scales independently and also have the possiblity to
  //     change the calculation of the scale.
  //
  //     Please check, which type of table you are using and then
  //     refer only to the comments and functions that are suitable for this
  //     fastNLO table.
  



  // -------- Initialize fastNLOReader --------- //
  // --- fastNLO user: Make an instance of your class that derived 
  //     from the FastNLOReader class and
  //     pass the name of the fastNLO table as an argument.
  //
  //        FastNLOUser* fnloreader = new FastNLOUser( tablename );
  //
  //     The example class for LHAPDF has overwriten the constructor
  //     and thus takes also the PDF-filename (and PDFset).
  FastNLOUser* fnloreader = new FastNLOUser( tablename , PDFFile , 0 );
  


  // ---- 'Setting'/init pdf ---- //
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
  //		fnloreader->PrintCurrentLHAPDFInformation().
  //	
  //	 If you are not really sure, which pdf set is currently used, then
  //     please reset the LHAPDF filename and memberset id.
  //
  //     Since the PDF evolution is typically provided in some external
  //     code, you can choose the interface to the evolution routine.
  //     By default LHAPDF is used:
  //           fnloreader->SetPDFInterface(FastNLOReader::kLHAPDF);
  //	
  //   fnloreader->SetLHAPDFfilename( PDFFile );
  //   fnloreader->SetLHAPDFset( 0 );



  // --- fastNLO user: After defining the PDF or changing the factorization
  //     scale you always have to refill the fastNLO internal PDF cache
  //     calling:
  //            fnloreader->FillPDFCache();
  //
  fnloreader->FillPDFCache();



  // ---- Setting Alpha_s value ---- //
  // --- fastNLO user: With fastNLO, the user can choose the value of
  //     alpha_s(Mz). Further, there are multiple options for the evolution
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
  fnloreader->SetAlphasMz(0.1179);

  //fnloreader->SetAlphasEvolution(FastNLOReader::kLHAPDFInternal);

  
  // ---- Set the units of your calculation (kPublicationUnits or kAbsoluteUnits) ---- //
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



  // ---- Set the calculation order (if available)---- //
  // --- fastNLO user: Each fastNLO table comes typically with
  //     various contributions.
  //     Currently, five different types of contributions have been tested.
  //     Three can be combined to give a scale, PDF and alpha_s dependent
  //     cross-section, one is a fixed multiplicative correction and, at last,
  //     also data points with uncertainties might be included in a table.
  //     For calculating a cross-section, by default only the LO & NLO contributions
  //     are used. However, each contribution can be swiched on or off separately.
  //     Please make sure to avoid combinations that do not make sense,
  //     e.g. 2-loop threshold corrections with LO pQCD.
  //     
  //     For switching a contribution on/off, its type must be known:
  //       - kFixedOrder		  -> Fixed order calculation (in alpha_s)
  //       - kThresholdCorrection	  -> Threshold corrections
  //       - kElectroWeakCorrection	  -> Electroweak corrections (not derived yet)
  //       - kNonPerturbativeCorrections  -> Non-perturbative corrections|Hadronisation corrections
  //     plus one must know the 'Id' of this contribution, which is typically
  //     printed when reading a table. To switch contribution on/off please use:
  //            fnloreader->SetContributionON( contrib, Id, on/off ) 
  //     To show the Id's of each contribution please call:
  fnloreader->PrintTableInfo();

  fnloreader->PrintFastNLOTableConstants(0);


  //****************************************************
  // 
  //   Some example options/uses are shown below
  // 
  //****************************************************



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

  
  if ( !fnloreader->GetIsFlexibleScaleTable() ) {
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
     //    fnloreader->SetMuRFunctionalForm(FastNLOReader::kQuadraticMean); // set function how to calculate mu_r from scale1 and scale2
     //    fnloreader->SetMuFFunctionalForm(FastNLOReader::kScale1);	 // set function how to calculate mu_f from scale1 and scale2
     //    fnloreader->SetExternalFuncForMuR( &Function_Mu );		 // set external function to calculate mu_r from scale1 and scale2
 

     // ---- (Re-)calculate cross sections ---- //
     // --- fastNLO user: Before you can access the fastNLO computed
     //     cross sections, you always have to call CalcCrossSection()!
     //     If you are not sure, whether you have recalculated the internal
     //     PDF-cache with your current scale choice you further can call 
     //     FillPDFCache() before calling CalcCrossSection().
     //     So, before accessing the cross sections, please call:
     //             fnloreader->CalcCrossSection();
  
     fnloreader->SetMuFFunctionalForm(FastNLOReader::kScale1);	 // set function how to calculate mu_f from scale1 and scale2
     fnloreader->SetMuRFunctionalForm(FastNLOReader::kScale1);	 // set function how to calculate mu_f from scale1 and scale2
     //   fnloreader->SetMuFFunctionalForm(FastNLOReader::kExpProd2);	 // set function how to calculate mu_f from scale1 and scale2
     //   fnloreader->SetMuRFunctionalForm(FastNLOReader::kExpProd2);	 // set function how to calculate mu_f from scale1 and scale2
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
  //     you might use following print methods:
  //             fnloreader->PrintCrossSections();
  //     Or print it (almost exaclty) like the Fortran reader code:
  //             fnloreader->PrintCrossSectionsDefault();
  //
  //  fnloreader->PrintCrossSectionsDefault();



  // ---- Information ---- //
  // --- fastNLO user: For a comprehensive insight into the fastNLO variables
  //     you can use:
  fnloreader->PrintFastNLOTableConstants(0);
  //     
  //     For a comparision with a Reference cross section calculated with
  //     NLOJet++ you might use:
  //             fnloreader->PrintCrossSectionsWithReference();
  //     WARNING: in v2.1 there are always three reference tables stored. Please check
  //     which one you access with this method. These reference cross sections also might
  //     be empty. In v2.0-tables, there are not necessarily any reference cross sections stored.
  //     The three references are calculated usually with different scales, that might be very 
  //     different from your choice. Therefore this comparison DOES not provide any information
  //     about the precision of this fastNLO table.
  //     INFO: NLOJet++ typically uses the NLOJet++-like alpha_s evolution with a value of 0.1179
  //     and a PDF similar to the cteq6m.LHgrid pdf-file.
  // ************************************************************************************************



  // ---- do sth. useful ---- //
  printf("\n");
  printf("%s",CSEPL.c_str());
  printf("fnlo-read: Calculate cross sections\n");
  printf("%s",CSEPL.c_str());
  
  // Give some info on contribution Ids
  // fnloreader->PrintTableInfo();
  
  // Example code to access cross sections and K factors:
  fnloreader->PrintFastNLODemo();
  
  // The presented example is done automatically for print out here  
  fnloreader->PrintCrossSectionsDefault();

  // Example code to print out data points (if available)
  fnloreader->PrintCrossSectionsData();

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

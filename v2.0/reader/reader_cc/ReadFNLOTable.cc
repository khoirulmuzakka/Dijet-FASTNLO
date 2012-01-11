#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include "FastNLOReader.h"
#include "Alphas.h"

//__________________________________________________________________________________________________________________________________
//
//
//__________________________________________________________________________________________________________________________________

double Function_Mu(double s1, double s2 );
void PrintWelcomeMessage( );

int main(int argc, char** argv){

  using namespace std;

  //---  Initialization for nice printing
  string cseps = " ##################################################################################\n";
  string lseps = " # --------------------------------------------------------------------------------\n";
  string csepl = "##################################################################";

  //---  Display welcome message
  PrintWelcomeMessage();
  
  // ---------------------------- parse commmand line ---------------------------- //
  printf(" # ReadFNLOTable: Program Steering\n");
  cout << lseps;
  string tablename = "table.tab";
  if ( argc <= 1 ){
    printf(" # ReadFNLOTable: WARNING! No table name given,\n");
    printf(" # taking the default table.tab instead!\n");
    printf(" #   For an explanation of command line arguments type:\n");
    printf(" #   ./ReadFNLOTable -h\n");
    // printf("Please specify filename.\n");
    // printf("Usage:\n  %s [tablename] <[pdfname]>\n",argv[0]);
    // return 0;
  } else {
    tablename = (const char*) argv[1];
    if (tablename == "-h"){ 
      printf(" #\n");
      printf(" # Usage: ./ReadFNLOTable [arguments]\n");
      printf(" # Table input file, def. = table.tab\n");
      printf(" # PDF set, def. = cteq6mE.LHgrid\n");
      printf(" #\n");
      printf(" # Give full path(s) if these are not in the cwd.\n");
      printf(" # Use \"_\" to skip changing a default argument.\n");
      printf(" #\n");
      return 0;
    } else if (tablename == "_") {
      tablename = "table.tab";
      printf("\n # ReadFNLOTable: WARNING! No table name given,\n");
      printf(" # taking the default table.tab instead!\n");
    } else {
      cout << " # ReadFNLOTable: Evaluating table: " << tablename << endl;
    }
  }

  //---  PDF set
  string PDFFile = "X";
  if ( argc > 2 ){
    PDFFile = (char*) argv[2];
  }
  if ( argc <= 2 || PDFFile == "_"){
    PDFFile = "cteq6m.LHpdf";
    printf(" # ReadFNLOTable: WARNING! No PDF set given,\n");
    printf(" # taking cteq6mE.LHgrid instead!\n");
  } else {
    cout << " # ReadFNLOTable: Using PDF set   : " << PDFFile << endl;
  }

  //---  Too many arguments
  if ( argc > 3 ){
    printf("ReadFNLOTable: ERROR! Too many arguments, aborting!\n");
    return 1;
  }

  // ---------------------------- end of parsing arguments ---------------------------- //



  // ************************** FastNLO and example documentation starts here ****************************

  // --- fastNLO user: Hello!
  //     If you use FastNLO for the first time, please read through the
  //     documentation and comments carefully in order to calculate
  //     a reasonable cross section.
  //
  //     In FastNLO version 2 (v2.0), there are two different types of tables, while
  //     those also have a little bit a different concept. We call those tables
  //     v2.0 as for tables similar (but with much larger flexiblity) to version 1.4.
  //     Further there are tables available that are called 'flexible-scale'-tables
  //     or sometimes v2.1, which have encoded a slightly advanced storage of matrix 
  //     elements and also of the scale variables. Those tables give you the possibility
  //     to change the renormalization and the factorization scale independently and
  //     also make a choice of how you calculate those scale variables a-posteriori.
  //
  //     You must check, which type of table you are currently using !!!
  //     Then please refer only to the comments and functions that are suitable for your
  //     FastNLO table.
  

  // -------- initialize FastNLOReader --------- //
  // --- fastNLO user: Make an instance of the FastNLO reader and
  //     pass the name of the FastNLO table as an argument.
  FastNLOReader* fnloreader = new FastNLOReader( tablename );
  


  // ---- 'Setting'/init pdf ---- //
  // --- fastNLO user: You can set the PDF information here.
  //     Since the PDF-evolution is typically provided in some external
  //     code, the user can choose the interface to the PDF routine. By
  //     default LHAPDF is set:
  //           fnloreader->SetPDFInterface(FastNLOReader::kLHAPDF);
  //     
  //     If you are using LHAPDF, you can set the PDF-name and the
  //     PDF-set using e.g.:
  //           fnloreader->SetLHAPDFfilename(PDFFile);
  //           fnloreader->SetLHAPDFset( int PDFSet );
  //     Or you directly talk to LHAPDF using LHAPDF::LHAPDF().
  //	
  fnloreader->SetLHAPDFfilename(PDFFile);
  fnloreader->SetLHAPDFset( 0 );



  // --- fastNLO user: After defining the pdf (or you change the
  //     factorization scale), you always have to refill the FastNLO-internal
  //     PDF cache calling:
  //            fnloreader->FillPDFCache();
  //
  fnloreader->FillPDFCache();



  // ---- Setting Alpha_s value ---- //
  // --- fastNLO user: With FastNLO, the user can choose the value of
  //     alpha_s(Mz). Further, there are multiple options for the evolution
  //     code to calculate alpha_s(mu_r) according to the RGE.
  //     To set the value of alpha_s(Mz) at the Z0-mass, please use:
  //            fnloreader->SetAlphasMz(0.1185);
  //     You always have to specify an alpha_s(Mz) value.
  //    
  //     Further you can specify the desired evolution code using
  //            fnloreader->SetAlphasEvolution(FastNLOReader::EAlphasEvolution);
  //     Following options can be choosen:
  //      - kGRV		The standard option in FastNLO. Here you can set
  //                            any value for Mz, and 2-/3-/4- loop iterative solution
  //                            of the RGE, number of flavors, flovor matching, etc...
  //                            Please see documentation in Alphas.h and Alphas.cc. For usage
  //                            you can access the static-class directly like e.g.:
  //                                Alphas::SetNLoop(2);
  //                                Alphas::SetMz(91.1786);
  //				INFO: When setting kGRV, you initialize Alphas::Alphas() with
  //				some reasonable values and overwrite previous settings to this class!!!
  //      - kNLOJET		This is the alpha_s evolution, which is used by nlojet++ as default
  //      - kCTEQpdf		This alpha_s evolution is used for the cteq6 PDFs. [Sure?]
  //      - kLHAPDFInternal	With this option, you access the alhpa_s evolution, which is defined
  //                            within the LHAPDF-file. You cannot change alphas(Mz) here !!!
  //      - kQCDNUMInternal	Using kQCDNUM as PDF-evolution code, you can further make use of the
  //				QCDNUM alpha_s evolution. Please see the QCDNUM manual for options then.
  //      - kFixed		Take always the fixed value of SetAlphasMz() without any evolution code.
  //
  //     INFO: The choice of the 'number of flavors' in your alhpa_s evolution might be 
  //     inconsistent with the number of flavors that are used for calculating the matrix elements.
  //     nlojet++ typically uses 5 massless flavors.
  //
  fnloreader->SetAlphasEvolution(FastNLOReader::kGRV);
  fnloreader->SetAlphasMz(0.1185);


  
  // ---- Set the units of your calculation (kPublicationUnits or kAbsoluteUnits) ---- //
  // --- fastNLO user: You can choose the dimension in which
  //     you want to access (or print) your cross section results.
  //     There are two possibilites:
  //       - The default option are 'publication units', i.e. divided by 
  //         bin widths if done so in publication
  //            fnloreader->SetUnits(FastNLOReader::kPublicationUnits);
  //       - The other option is 'absolute' units in barn, but still in
  //         the same magnitude as in the publication (e.g. pb, fb, nb, etc...)
  //            fnloreader->SetUnits(FastNLOReader::kAbsoluteUnits);
  fnloreader->SetUnits(FastNLOReader::kPublicationUnits);



  // ---- Set the calculation order (if available)---- //
  // --- fastNLO user: Each FastNLO table comes typically with
  //     various contributions.
  //     There are five different types of information,
  //     while four of them can contribute to the cross section and
  //     one can be the data card.
  //     For calculating the cross section, by default all contributions
  //     are used. However, you can switch on or off each single
  //     contribution seperately.
  //     
  //     For switching a contribution oon/off, you must know the type of it, e.g.
  //       - kFixedOrder		-> Fixed order calculation 
  //       - kThresholdCorrection	-> Threshold corrections
  //       - kElectroWeakCorrection	-> Electro-weak corrections
  //     And you must know the 'Id' of this contribution, which is typically
  //     printed when reading a table. Please use:
  //            fnloreader->SetContributionON( contrib, Id, on/off) 
  //


  // ****************************************************
  // 
  //   Several example options are shown here
  // 
  // ****************************************************


  // ---- options for scales in v2.0 tables ---- //
  // --- fastNLO user: Here you can specify the options if you use a v2.0 fastNLO-table, 
  //     which is NOT a 'flexible-scale'-table (v2.1):
  //     A fastNLO table comes usually with various precalculated scale-variations.
  //     There, the renormalization and the factorization scale are varied simultaneously
  //     by a certain factor. The welcome screen should show you the information about
  //     the available scale variations. For accessing a table with a certain factor you
  //     need to know the id (int) and then use:
  //            fnloreader->SetScaleVariation(0);
  //
  //     Further you have the possiblity to vary the renormalization scale by any
  //     factor. Please use e.g.:
  //            fnloreader->SetScaleFactorMuR(0.5);
  //     WARNING: This option is inconsistent with threshold corrections. Those are
  //     automatically switched OFF to preserve consistency.
  

  
  // ---- options for scales in 'flexible-sclale' tables (v2.1) ---- //
  // --- fastNLO user: You can choose a function, how you would like
  //     to compute the renormalization and factorization scale. 
  //     Each 'flexible scale'-table comes with two variables that can be used 
  //     for calculating the scales. Those are called scale1 and scale2, while
  //     at least one is in dimension GeV.
  //     DIS tables have typically stored scale1 = Q and scale2 = pt, while
  //     HHC tables have typcially scale1 = pt and scale2 = y, however this
  //     can be any other observable. 
  //     Please check, which obervables are stored as scale variables !!!
  //
  //     There are two possibilities, how you can define your scale now:
  //       - use predefined functions using e.g.
  //            fnloreader->SetMuRFunctionalForm(FastNLOReader::EScaleFunctionalForm);
  //         for changing the calculation of the renormalizatoin scale.
  //         Please refer to FastNLOReader.h for all options of EScaleFunctionalForm.
  //       - or you can pass a function pointer to FastNLOReader using
  //            fnloreader->SetExternalFuncForMuR( double (*Func)(double,double) );
  //         to pass any function using scale1 and scale2 to FastNLO.
  //    
  //     Further you can define a scale-factor (one for the renormalization
  //     and one for the factorization scale) which is independently from the
  //     predefined function multiplied to the scale, like e.g. for mu_r
  //         mu_r = scalefac * f(scale1,scale2)
  //     Please use e.g. for mu_r:
  //             fnloreader->SetScaleFactorMuR(1.5);
  //
  //     For changing the factorization scale, replace all 'MuR' by 'MuF' in the function calls.
  //  
  //     INFO: All above mentioned functions perform automatically a refilling of the
  //     FastNLO-internal PDF-cache. To switch it off you can use a boolean, like:
  //             fnloreader->SetMuFFunctionalForm(FastNLOReader::kScale1 , false );
  //
  //     WARNING: There is some choice made for a default calculation. Please always think
  //     about the choice of the scales and make a definition.
  //     Default setting for DIS tables:
  //       - mu_r:  kQuadraticMean	-> mu_r = sqrt( (Q^2 + scale2^2)/2.) // because scale1=Q!
  //       - mu_f:  kScale1		-> mu_f = Q
  //     Default setting for pp and ppbar tables:
  //       - mu_r:  kScale1		-> mu_r = scale1
  //       - mu_f:  kScale1		-> mu_f = scale1
  //
  //  Possible calls are e.g.:
  //    fnloreader->SetMuRFunctionalForm(FastNLOReader::kQuadraticMean);	// set function how to calculate mu_r from scale1 and scale2
  //    fnloreader->SetMuFFunctionalForm(FastNLOReader::kScale1);		// set function how to calculate mu_f from scale1 and scale2
  //    fnloreader->SetExternalFuncForMuR( &Function_Mu );			// set external function to calculate mu_r from scale1 and scale2
  //    fnloreader->SetScaleFactorMuR(1.5);					// set scale factor for mu_r
  //    fnloreader->SetScaleFactorMuF(0.66);					// set scale factor for mu_f
  fnloreader->SetExternalFuncForMuR( &Function_Mu );
  fnloreader->SetExternalFuncForMuF( &Function_Mu );
    


  // ---- (Re-)calculate cross sections ---- //
  // --- fastNLO user: Before you can access the FastNLO computed
  //     cross sections, you always have to call CalcCrossSection() !!!
  //     If you are not sure, if you have recalculated the internal
  //     alphas-cache and PDF-cache with your current scale choice,
  //     you further can call FillAlphasCache() and FillPDFCache() before.
  //     So please call:
  //             fnloreader->CalcCrossSection();
  fnloreader->CalcCrossSection();


  
  // ---- Get cross sections ---- //
  // --- fastNLO user: To access the cross section from FastNLO
  //     you should use:
  //           vector < double > xs = fnloreader->GetXSection();
  //     If you want to have a pointer to an array of numbers you might use
  //           double* cs = &fnloreader->GetXSection()[0];
  //     
  //     Further you can access the k-factor, which is calculated with all
  //     'contributions' (e.g. no Non-perturabtive corrections) against the
  //     LO fixed-order contribution.
  //           vector < double > kFactors = fnloreader->GetKFactors();
  //



  // ---- Printing ---- //
  // --- fastNLO user: For an easy overview of your cross section calculation
  //     you might use following print methods:
  //             fnloreader->PrintCrossSections();
  //     Or print it (almost exaclty) like the fortran reading code:
  //             fnloreader->PrintCrossSectionsLikeFreader();



  // ---- Information ---- //
  // --- fastNLO user: For a comprehensive insight into the FastNLO variables
  //     you can use:
  //             fnloreader->Print();
  //     
  //     For a comparision with a Reference cross section calculated with
  //     nlojet++ you might use:
  //             fnloreader->PrintCrossSectionsWithReference();
  //     WARNING: in v2.1 there are always three reference tables stored. Please check
  //     which one you access with this method. Those reference cross sections also might
  //     be empty. In v2.0-tables, there are not necessarily any reference cross sections stored.
  //     Those three references are calculated usually with different scales, that might be very 
  //     different from your choice. Therefore this comparision DOES not provide any information
  //     about the precision of this FastNLO table.
  //     INFO: NLOJET++ typically uses the nlojet-like alpha_s evolution with a value of 0.1179
  //     and a PDF similar to the cteq6m.LHgrid pdf-file.
  //


  // ************************************************************************************************


  // ---- do sth. useful ---- //
  cout << csepl;
  printf("ReadFNLOTable: Calculate cross sections\n");
  cout << csepl;
  int nscale = fnloreader->GetNScaleVariations();
  for (int iscale = 0;iscale < nscale ;iscale++) {
    fnloreader->SetScaleVariation(iscale);
    fnloreader->CalcCrossSection();
    fnloreader->PrintCrossSectionsLikeFreader();
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


void PrintWelcomeMessage( ){
   // Say hello to the FastNLO user.
   string cseps = " ##################################################################################\n";
   string lseps = " # --------------------------------------------------------------------------------\n";
   string csepl = "##################################################################";
   csepl = csepl + csepl + "\n";
   printf("\n");
   cout << cseps;
   printf(" # ReadFNLOTable\n");
   cout << cseps;
   printf(" # Program to read fastNLO v2 tables and derive\n");
   printf(" # QCD cross sections using PDFs from LHAPDF\n");
   cout << cseps;
}

//__________________________________________________________________________________________________________________________________

#include <iostream>
#include <iomanip>
#include "FastNLOReader.h"
#include <string>

//__________________________________________________________________________________________________________________________________
//
//
//__________________________________________________________________________________________________________________________________


int main(int argc, char** argv){

  using namespace std;

  //---  Initialization
  string cseps = " ##################################################################\n";
  string lseps = " # ----------------------------------------------------------------\n";
  string csepl = "##################################################################";
  csepl = csepl + csepl + "\n";

  //---  Initial output
  printf("\n");
  cout << cseps;
  printf(" # ReadFNLOTable\n");
  cout << cseps;
  printf(" # Program to read fastNLO v2 tables and derive\n");
  printf(" # QCD cross sections using PDFs from LHAPDF\n");
  cout << cseps;
  
  //---  Parse command line
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
    cout << " # ReadFNLOTable: Using PDF set: " << PDFFile << endl;
  }

  //---  Too many arguments
  if ( argc > 3 ){
    printf("ReadFNLOTable: ERROR! Too many arguments, aborting!\n");
    return 1;
  }

  // if ( argc >= 3 ){
  //   string ThePDF = (const char*) argv[2];
  //   if ( ThePDF == "ct10" || ThePDF == "CT10" ){
  //     PDFFile	= "CT10.LHgrid";		
  //   } else if ( ThePDF == "cteq66" || ThePDF=="CTEQ66" ){
  //     PDFFile	= "cteq66.LHgrid";		
  //   } else if ( ThePDF == "cteq6m" || ThePDF=="CTEQ6M" ){
  //     cout << "ToDo: Cteq6m not yet proper implemented. (PDFEigenvektorer)"<<endl;
  //     PDFFile	= "cteq6m.LHpdf";		
  //   } else if ( ThePDF == "cteq65" || ThePDF=="CTEQ65" ){
  //     PDFFile	= "cteq65.LHgrid";		
  //   } else if ( ThePDF == "herapdf10" || ThePDF=="HERAPDF10" ){
  //     PDFFile	= "HERAPDF10_EIG.LHgrid";	
  //   } else if ( ThePDF == "MSTW2008as" || ThePDF == "MSTW2008AS" || ThePDF=="mstw2008as" ){
  //     PDFFile	= "MSTW2008nlo90cl_asmz+90cl.LHgrid";	
  //   } else if ( ThePDF == "MSTW2008" || ThePDF=="mstw2008" ){
  //     PDFFile	= "MSTW2008nlo90cl.LHgrid";	
  //   } else if ( ThePDF == "nnpdf20" || ThePDF=="NNPDF20" ){
  //     PDFFile	= "NNPDF20_100.LHgrid";		
  //   } else if ( ThePDF == "nnpdf21" || ThePDF=="NNPDF21" ){
  //     PDFFile	= "NNPDF21_100.LHgrid";		
  //   } else if ( ThePDF == "herapdf15" || ThePDF=="HERAPDF15" ){
  //     PDFFile	= "HERAPDF1.5_EIG.LHgrid";	
  //   } else{
  //     PDFFile	= ThePDF;
  //     printf("Warning. I try to set the pdf-file directly to %s.\n",PDFFile.c_str());
  //   }
  // }
  // else {
  //   printf("Using default pdf file %s.\n",PDFFile.c_str());
  // }  

  // -------- initialize FastNLOReader --------- //
  FastNLOReader* fnloreader = new FastNLOReader( tablename );
  
  // ---- 'Setting'/init pdf ---- //
  //  fnloreader->SetPDFInterface(FastNLOReader::kLHAPDF);			// Interfaces to different pdf routines
  fnloreader->SetLHAPDFfilename(PDFFile);
  fnloreader->SetLHAPDFset(0);
  fnloreader->FillPDFCache();	// pdf is 'external'! you always have to call FillPDFCache();

  // ---- Setting Alpha_s value ---- //
  //   fnloreader->SetAlphasEvolution(FastNLOReader::kNLOJET);			// set the precoded alpha_s evolution codes
  fnloreader->SetAlphasMz(0.1179);						// you MUST specify some alpha_s value
  
  // ---- Set the units of your calculation (kPublicationUnits or kAbsoluteUnits) ---- //
  //  fnloreader->SetUnits(FastNLOReader::kAbsoluteUnits);
  // ---- The default should be publication units, i.e. divided by bin widths if done so in publication
  fnloreader->SetUnits(FastNLOReader::kPublicationUnits);

  // ---- Set the calculation order (if available)---- //
  fnloreader->SetCalculationOrder(FastNLOReader::kAllAvailableOrders);		// Set the order of your calculation (Mind: k-factor is always calculated as ratio to LO calcuation)

  // ---- Set the calculation order (if available)---- //
  fnloreader->SetCalculationOrder(FastNLOReader::kAllAvailableOrders);		// Set the order of your calculation (Mind: k-factor is always calculated as ratio to LO calcuation)



  // ****************************************************
  // 
  //   Several example options are shown here
  // 
  // ****************************************************

  // ---- options ---- //
  //fnloreader->Print();		// print FastNLO internal variables


  // ---- options for v2.0 tables ---- //
  //fnloreader->SetScaleVariation(0);	// set scale variation table
  

  
  // ---- options for v2.1 tables ---- //
  //   fnloreader->SetMuRFunctionalForm(FastNLOReader::kQuadraticMean);		// set function how to calculate mu_f from scale1 and scale2
  //   fnloreader->SetMuFFunctionalForm(FastNLOReader::kScale1);		// set function how to calculate mu_f from scale1 and scale2
  //   fnloreader->SetScaleFactorMuR(1.5);					// set scale factor for mu_r
  //   fnloreader->SetScaleFactorMuF(0.66);					// set scale factor for mu_f
    

  // ---- (Re-)calculate cross sections ---- //
  //  fnloreader->CalcCrossSection();
  

  // ****************************************************

  // ---- do sth. useful ---- //
  cout << csepl;
  printf("ReadFNLOTable: Calculate cross sections\n");
  cout << csepl;
  int iscalevar = 3;
  fnloreader->SetScaleVariation(iscalevar);
  fnloreader->CalcCrossSection();
  fnloreader->PrintCrossSectionsLikeFreader(iscalevar);
  iscalevar = 2;
  fnloreader->SetScaleVariation(iscalevar);
  fnloreader->CalcCrossSection();
  fnloreader->PrintCrossSectionsLikeFreader(iscalevar);
  iscalevar = 0;
  fnloreader->SetScaleVariation(iscalevar);
  fnloreader->CalcCrossSection();
  fnloreader->PrintCrossSectionsLikeFreader(iscalevar);
  iscalevar = 1;
  fnloreader->SetScaleVariation(iscalevar);
  fnloreader->CalcCrossSection();
  fnloreader->PrintCrossSectionsLikeFreader(iscalevar);




  // ---- get cross sections ---- //
  //   vector < double > xs = fnloreader->GetXSection();
  //   vector < double > xsref = fnloreader->GetReferenceXSection();
  //   vector < double > kFactors = fnloreader->GetKFactors();
  

  return 0;

}

//__________________________________________________________________________________________________________________________________

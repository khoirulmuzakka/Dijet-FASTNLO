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

  // -------- parse input --------- //
  if ( argc <= 1 ){
    printf("Please specify filename.\n");
    printf("Usage:\n  %s [tablename] <[pdfname]>\n",argv[0]);
    return 0;
  }
  
  string tablename = (const char*) argv[1];

  string PDFFile	= "cteq6m.LHpdf";
  
  if ( argc >= 3 ){
    string ThePDF = (const char*) argv[2];
    if ( ThePDF == "ct10" || ThePDF == "CT10" ){
      PDFFile	= "CT10.LHgrid";		
    } else if ( ThePDF == "cteq66" || ThePDF=="CTEQ66" ){
      PDFFile	= "cteq66.LHgrid";		
    } else if ( ThePDF == "cteq6m" || ThePDF=="CTEQ6M" ){
      cout << "ToDo: Cteq6m not yet proper implemented. (PDFEigenvektorer)"<<endl;
      PDFFile	= "cteq6m.LHpdf";		
    } else if ( ThePDF == "cteq65" || ThePDF=="CTEQ65" ){
      PDFFile	= "cteq65.LHgrid";		
    } else if ( ThePDF == "herapdf10" || ThePDF=="HERAPDF10" ){
      PDFFile	= "HERAPDF10_EIG.LHgrid";	
    } else if ( ThePDF == "MSTW2008as" || ThePDF == "MSTW2008AS" || ThePDF=="mstw2008as" ){
      PDFFile	= "MSTW2008nlo90cl_asmz+90cl.LHgrid";	
    } else if ( ThePDF == "MSTW2008" || ThePDF=="mstw2008" ){
      PDFFile	= "MSTW2008nlo90cl.LHgrid";	
    } else if ( ThePDF == "nnpdf20" || ThePDF=="NNPDF20" ){
      PDFFile	= "NNPDF20_100.LHgrid";		
    } else if ( ThePDF == "nnpdf21" || ThePDF=="NNPDF21" ){
      PDFFile	= "NNPDF21_100.LHgrid";		
    } else if ( ThePDF == "herapdf15" || ThePDF=="HERAPDF15" ){
      PDFFile	= "HERAPDF1.5_EIG.LHgrid";	
    } else{
      PDFFile	= ThePDF;
      printf("Warning. I try to set the pdf-file directly to %s.\n",PDFFile.c_str());
    }
  }
  else {
    printf("Using default pdf file %s.\n",PDFFile.c_str());
  }  

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

  
  // ---- Set the units of your calculation ---- //
  fnloreader->SetUnits(FastNLOReader::kAbsoluteUnits);				// Set the unit of your calculation (kPublicationUnits or kAbsoluteUnits)

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
  // fnloreader->SetScaleVariation(0);
  // fnloreader->CalcCrossSection();
  // fnloreader->PrintCrossSectionsLikeFreader();
  fnloreader->SetScaleVariation(1);
  fnloreader->CalcCrossSection();
  fnloreader->PrintCrossSectionsLikeFreader();
  // fnloreader->SetScaleVariation(2);
  // fnloreader->CalcCrossSection();
  // fnloreader->PrintCrossSectionsLikeFreader();
  // fnloreader->SetScaleVariation(3);
  // fnloreader->CalcCrossSection();
  // fnloreader->PrintCrossSectionsLikeFreader();




  // ---- get cross sections ---- //
  //   vector < double > xs = fnloreader->GetXSection();
  //   vector < double > xsref = fnloreader->GetReferenceXSection();
  //   vector < double > kFactors = fnloreader->GetKFactors();
  

  return 0;

}

//__________________________________________________________________________________________________________________________________

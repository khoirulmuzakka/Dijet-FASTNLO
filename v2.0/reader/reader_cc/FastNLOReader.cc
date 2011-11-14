// Author: Daniel Britzger
// DESY, 23/07/2011

//////////////////////////////////////////////////////////////////////////
//                                                                      //
//  FastNLOReader                                                       //
//                                                                      //
//  FastNLOReader is a standalone code for reading                      //
//  FastNLO tables of version 2.0 and v2.1 for DIS processes            //
//  It is also optimized for an integration into                        //
//  the H1Fitter project.                                               //
//                                                                      //
//  FastNLO is developed by                                             //
//    D. Britzger, T. Kluge, K. Rabbertz, F. Stober, M. Wobisch         //
//    (publication in preparation)                                      //
//    http://projects.hepforge.org/fastnlo                              //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "FastNLOReader.h"

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <LHAPDF/LHAPDF.h>
#include "Alphas.h"

using namespace std;


//______________________________________________________________________________


extern "C"{
  void fpdfxq_(int *iset, const double *x, const double *q2, double *pdfs,int *ichk);
  void evolution_();
  //int getalf_( double* alfs, double* r2 );
}


//______________________________________________________________________________


FastNLOReader::FastNLOReader(void)
{
   //
   // do not call the standard constructor
   // 
   BlockB_LO		= NULL;
   BlockB_NLO		= NULL;
   BlockB_LO_Ref	= NULL;
   BlockB_NLO_Ref	= NULL;
   fUnits		= kPublicationUnits;
   fOrder		= kAllAvailableOrders;
   cout << "FastNLOReader::FastNLOReader. Please set a filename using SetFilename(<name>)! "<<endl;
}


FastNLOReader::FastNLOReader(string filename)
{
   BlockB_LO		= NULL;
   BlockB_NLO		= NULL;
   BlockB_LO_Ref	= NULL;
   BlockB_NLO_Ref	= NULL;
   fUnits		= kPublicationUnits;
   fOrder		= kAllAvailableOrders;
   
   SetFilename(filename);
}


//______________________________________________________________________________


FastNLOReader::~FastNLOReader(void)
{
   if ( BlockB_LO )		delete BlockB_LO;
   if ( BlockB_NLO )		delete BlockB_NLO;
   if ( BlockB_LO_Ref )		delete BlockB_LO_Ref;
   if ( BlockB_NLO_Ref )	delete BlockB_NLO_Ref;
}


//______________________________________________________________________________



void FastNLOReader::SetAlphasEvolution(EAlphasEvolution AlphasEvolution){
   fAlphasEvolution = AlphasEvolution; 
   if (AlphasEvolution==kLHAPDFInternal || AlphasEvolution==kQCDNUMInternal ) {
      cout << "Warning. You cannot change the Alpha_s(Mz) value."<<endl; 
   }

   if ( AlphasEvolution == kGRV ){
      Alphas::SetMz(91.1876);
      //Alphas::SetAlphasMz(ALPSMZ);
      Alphas::SetNf(5);
      Alphas::SetNLoop(2);
      Alphas::SetFlavorMatchingOn(true);
   }
}


//______________________________________________________________________________



void FastNLOReader::SetFilename(string filename){
   ffilename	= filename;
   Init();
}


//______________________________________________________________________________



void FastNLOReader::Init(){
  printf(" ***************************************************************** \n");
  printf(" *  \n");
  printf(" *  FastNLO Reader - version 0.6\n");
  printf(" *  \n");
  printf(" *  This code is a reader for FastNLO tables, which were\n");
  printf(" *  calculated using nlojet++ v4.1.3 and FastNLO v2.0 (or higher).\n"); 
  printf(" *  \n"); 
  printf(" *  The FastNLO package by\n"); 
  printf(" *    D. Britzger, T. Kluge, K. Rabbertz, F. Stober, M. Wobisch\n");
  printf(" *    (publication in preparation)\n");
  printf(" *  \n");
  printf(" *  for NLOJET++ please cite \n");
  printf(" *    Z. Nagy, Phys. Rev. Lett. 88, 122003 (2002),\n");
  printf(" *    Z. Nagy, Phys. Rev. D68, 094002 (2003)\n");
  printf(" *  \n");
  printf(" ***************************************************************** \n");

  ReadTable();
  PrintFastNLOTableConstants();

  SetPDFInterface(FastNLOReader::kLHAPDF);
  SetAlphasEvolution(FastNLOReader::kGRV);
  
  InitScalevariation();

}


//______________________________________________________________________________



void FastNLOReader::InitScalevariation(){
  
  fScaleFacMuR	= 1.;
  fScaleFacMuF	= 1.;
  fScalevar	= 0;

  if ( BlockB_NLO->NScaleDep != 3 ){
    // this is an 'original' v2.0 table
    printf (" *  This table has following %d scale variations for 'theory-error' determination.\n",BlockB_NLO->Nscalevar[0]);
    printf (" *    scalevar #n -> scalefactor\n");
    for ( int i = 0 ; i<BlockB_NLO->Nscalevar[0]; i++ ){
      printf (" *         '%d'    ->    %4.2f .\n", i, BlockB_NLO->ScaleFac[0][i]);
    }
    printf (" *    Setting scale factor to %4.2f and varying mu_f and mu_r simultaneously.\n",BlockB_NLO->ScaleFac[0][0]);
    fScalevar	= 0;
  }

  else if ( BlockB_NLO->NScaleDep == 3 ){
    // this is a MuVar table. You can vary mu_f and mu_r independently by any factor
    // and you can choose the functional form of mu_f and mu_r as functions of
    // scale1 and scale1 (called partly scaleQ2 and scalePt).
    
    if ( BlockB_NLO->ScaleDescript[0].size() <0 ) {
      printf("Error. No scaledescription available.\n"); // the code will crash soon.
      fMuFFunc	= kScale1;
      fMuRFunc	= kScale1;
      return;
    }

    // ---- DIS ---- //
    if ( BlockB_LO->NPDFDim == 0 ) {
      fMuRFunc	= kQuadraticMean;
      fMuFFunc	= kScale1;
      printf (" *    Setting factorization scale to mu_f^2 = %s^2 .\n", BlockB_NLO->ScaleDescript[0][0].c_str() );
      if ( BlockB_NLO->ScaleDescript[0].size() == 2 ){
	printf (" *    Setting renormalization scale to mu_r^2 = (%s^2 + %s^2)/2 .\n", BlockB_NLO->ScaleDescript[0][0].c_str() , BlockB_NLO->ScaleDescript[0][1].c_str() );
      }
      else if ( BlockB_NLO->ScaleDescript[0].size() == 1 &&  BlockB_LO->ScaleNodeScale2[0].size() > 3 ){
	printf("FastNLOReader::InitScalevariation. Warning. Could not find description for scale variables.\n");
	printf (" *    Setting renormalization scale to mu_r^2 = (scale1^2 + scale2^2)/2 .\n" );
      }
      else if ( BlockB_NLO->ScaleDescript[0].size() == 1 &&  BlockB_LO->ScaleNodeScale2[0].size() <= 3  ){
	printf("FastNLOReader::InitScalevariation. Warning. This table has only one scale variable %s stored.\n", BlockB_NLO->ScaleDescript[0][0].c_str() );
	printf (" *    Setting renormalization scale to mu_r^2 = %s^2 .\n", BlockB_NLO->ScaleDescript[0][0].c_str() );
	fMuRFunc	= kScale1;
      }
      else {
	printf("Error. I don't know what to do.\n");
      }
    }
        
    else printf("Error. Unknown process.\n");

  }
  
  else {
    printf("FastNLOReader::InitScalevariation(). ERROR. Could not identify table..\n");
  }
  
  
}


//______________________________________________________________________________



double FastNLOReader::CalcMu( FastNLOReader::EMuX kMuX , double scale1, double scale2, double scalefac ){
   //
   //  Calculate the scales with the defined function and the 
   //  corresponding prefactor.
   //
  
   if ( kMuX == kMuR && fScaleFacMuR != scalefac ) printf("Error. Sth. went wrong with the scales.\n");
   if ( kMuX == kMuF && fScaleFacMuF != scalefac ) printf("Error. Sth. went wrong with the scales.\n");
  
   EScaleFunctionalForm Func;
   if ( kMuX  == FastNLOReader::kMuR )		Func	= fMuRFunc;    // return renormalization scale
   else						Func	= fMuFFunc;    // return factorization scale
   //   else if ( kMuX  == FastNLOReader::kMuF )	Func	= fMuFFunc;    // return factorization scale
   //   else printf( "I dont know what to do.\n");
  
   double mu = 0;

   if		( Func == kScale1 )		mu	= scale1 ;
   else if	( Func == kScale2 )		mu	= scale2 ;
   else if	( Func == kQuadraticSum )	mu	= FuncMixedOver1(scale1,scale2) ;
   else if	( Func == kQuadraticMean )	mu	= FuncMixedOver2(scale1,scale2) ;
   else if	( Func == kQuadraticSumOver4 )	mu	= FuncMixedOver4(scale1,scale2) ;
   else if	( Func == kScaleMax )		mu	= FuncMax(scale1,scale2);
   else if	( Func == kScaleMin )		mu	= FuncMin(scale1,scale2);
   else printf( "Error. could not identify functional form for scales calculation.\n");
  
   return scalefac * mu;

}


//______________________________________________________________________________
double FastNLOReader::FuncMixedOver1(double scale1 , double scale2 ){
  return ( sqrt( (pow(scale1,2) + pow(scale2,2))  / 1. ) ) ;
}

//______________________________________________________________________________
double FastNLOReader::FuncMixedOver2(double scale1 , double scale2 ){
  return ( sqrt( (pow(scale1,2) + pow(scale2,2))  / 2. ) ) ;
}

//______________________________________________________________________________
double FastNLOReader::FuncMixedOver4(double scale1 , double scale2 ){
  return ( sqrt( (pow(scale1,2) + pow(scale2,2))  / 4. ) ) ;
}

//______________________________________________________________________________
double FastNLOReader::FuncLinearMean(double scale1 , double scale2 ){
  return ( scale1 + scale2 ) / 2. ;
}

//______________________________________________________________________________
double FastNLOReader::FuncLinearSum(double scale1 , double scale2 ){
  return scale1 + scale2;
}

//______________________________________________________________________________
double FastNLOReader::FuncMax(double scale1 , double scale2 ){
  if ( scale1 > scale2 ) return scale1;
  else return scale2;
}

//______________________________________________________________________________
double FastNLOReader::FuncMin(double scale1 , double scale2 ){
  if ( scale1 < scale2 ) return scale1;
  else return scale2;
}




//______________________________________________________________________________



double FastNLOReader::SetScaleVariation(int scalevar , bool ReFillCache ){ 
  
  // ------------------------------------------------
  //   Set the scalevariation factor for detemining the
  //   'theory'-error. Usually, you have tables stored with
  //   factors of 0.5, 1 and 2 times the nominal scale.
  //     corresponding to:
  //     scalevar -> scalefactor
  //        '0'   ->   1
  //        '1'   ->   0.5
  //        '2'   ->   2
  //   This method returns the scalefactor correspoding to
  //   the choosen 'scalevar'.
  // ------------------------------------------------

  if ( BlockB_NLO->NScaleDep == 3 ){
    printf("FastNLOReader::SetScaleVariation(). Info: This is a v2.1 table, therefore, you can choose all possible scale variations. Your Scalevar has to be '0'.\n");
    printf("    Please use SetScaleFacMuR(double) and SetScaleFacMuF(double).\n");
    return 0;
  }
  

  if (  scalevar >= BlockB_NLO->Nscalevar[0]  ){
     printf("Warning in FastNLOReader::SetScaleVariation. This table has only %d scalevariations stored. You wanted to acces number %d. Using '0' instead.\n", BlockB_NLO->Nscalevar[0] ,scalevar );
     fScalevar	= 0;
     return BlockB_NLO->ScaleFac[0][0] ;
  }
  
  fScalevar	= scalevar;
  printf(" * FastNLOReader::SetScaleVariation. Scalefactor of %4.2f for the nominal scale is chosen.\n",BlockB_NLO->ScaleFac[0][fScalevar]);

  if ( ReFillCache ){
    FillAlphasCache();
    FillPDFCache();
  }

  return BlockB_NLO->ScaleFac[0][fScalevar];
  
}




//______________________________________________________________________________



void FastNLOReader::SetFunctionalForm( EScaleFunctionalForm func , FastNLOReader::EMuX kMuX  ){
  //
  //  For MuVar tables this method sets the functional form of
  //  the renormalization or the factorization scale.
  //     func:  Choose a pre-defined function
  //     kMuX:  is it for mu_r or for mu_f
  //

  if ( BlockB_NLO->NScaleDep != 3 ) {
    printf("FastNLOReader::SetFunctionalForm. Warning. This is not a MuVar table.\n");
    printf("      SetFunctionalForm has no impact.\n");
    printf("      Please use another file, if you want to change your scale-definition.\n");
    return;
  }

  if ( kMuX == kMuR ) fMuRFunc = func;
  else fMuFFunc = func;

  if	( func == kScale2 || func == kQuadraticSum||  func == kQuadraticMean ||  func == kQuadraticSumOver4 ||  func == kScaleMax|| func == kScaleMin ) {
     if ( BlockB_NLO->ScaleNodeScale2[0].size() <= 3){
	printf("FastNLOReader::SetFunctionalForm. Error. There is no second scale variable available in this table.\n");
	printf("      Please use FastNLOReader::kScale1 only.\n");
	if ( kMuX == kMuR ) fMuRFunc = kScale1;
	else fMuFFunc = kScale1;
     }
     for(int i=0;i<NObsBin;i++){
	if ( BlockB_NLO->ScaleNodeScale2[i].size() < 6 ){
	   printf("FastNLOReader::SetFunctionalForm. Warning. Scale2 has only very little nodes (n=%d) in bin %d.\n",BlockB_LO->ScaleNodeScale2[i].size(),i);
	}
     }
  }
}


//______________________________________________________________________________


void FastNLOReader::SetMuRFunctionalForm( EScaleFunctionalForm func , bool ReFillCache  ){
   SetFunctionalForm(func,kMuR);

  if ( ReFillCache ){
    FillAlphasCache();
    FillPDFCache();
  }

}


//______________________________________________________________________________


void FastNLOReader::SetMuFFunctionalForm( EScaleFunctionalForm func , bool ReFillCache  ){
  SetFunctionalForm(func,kMuF);

  if ( ReFillCache ){
    FillPDFCache();
  }

}

//______________________________________________________________________________



void FastNLOReader::SetScaleFactorMuR(double fac , bool ReFillCache ){
  // 
  // Set scale factor for scale variations in MuVar tables
  // You have to ReFill your cache!
  // This is done automatically, but if you want to do it by yourself
  // set ReFillCache=false
  //

  if ( BlockB_NLO->NScaleDep != 3 ) {
    printf("FastNLOReader::SetScaleFactorMuR. Warning. This is not a MuVar table.\n");
    printf("      SetScaleFactorMuR has no impact.\n");
    printf("      Please use SetScaleVariation(int) instead.\n");
  }
  fScaleFacMuR = fac;

  if ( ReFillCache ){
    FillAlphasCache();
    FillPDFCache();
  }
}


//______________________________________________________________________________


void FastNLOReader::SetScaleFactorMuF(double fac , bool ReFillCache ){
  // 
  // Set scale factor for scale variations in MuVar tables
  // You have to ReFill your cache.
  // This is done automatically, but if you want to do it by yourself
  // set ReFillCache=false
  //
  
  if ( BlockB_NLO->NScaleDep != 3 ) {
    printf("FastNLOReader::SetScaleFactorMuF. Warning. This is not a MuVar table.\n");
    printf("      SetScaleFactorMuF has no impact.\n");
    printf("      Please use SetScaleVariation(int) instead.\n");
  }
  fScaleFacMuF = fac;

  if ( ReFillCache ){
    FillPDFCache();
  }
}


//______________________________________________________________________________



void FastNLOReader::ReadTable(void)
{
  //
  // Read in the FastNLO Table
  //
  
  // ---- check if file exists ----- //
  FILE* fp = fopen(ffilename.c_str(), "r");
  if (fp) {
    fclose(fp);
  } else {
    printf("Error. FastNLO table file does not exists. Was looking for: %s. Exiting.\n",ffilename.c_str());
    exit(1);
  } 
  
  // open stream
  ifstream* instream = new ifstream(ffilename.c_str(),ios::in);

  ReadBlockA1(instream);
  ReadBlockA2(instream);

  PrintBlockA1();
  PrintBlockA2();
  
  int nblocks	= Ncontrib + Ndata;
  printf(" * FastNLOReader::ReadTable(). Reading %d B-Blocks.\n",nblocks);
  for(int i=0;i<nblocks;i++){
    FastNLOBlockB* blockb	= new FastNLOBlockB( "ReadingBlockB", NObsBin , instream );
    if ( blockb->IContrFlag1==1 && blockb->IContrFlag2==1 && blockb->IRef==0 ){
      if ( blockb->NScaleDep != 3 )   blockb->SetName("BlockB. LO. v2.0.");
      if ( blockb->NScaleDep == 3 )   blockb->SetName("BlockB. LO. v2.1.");
      BlockB_LO		= blockb;
    }
    else if ( blockb->IContrFlag1==1 && blockb->IContrFlag2==2 && blockb->IRef==0 ){
      if ( blockb->NScaleDep != 3 )   blockb->SetName("BlockB. NLO. v2.0.");
      if ( blockb->NScaleDep == 3 )   blockb->SetName("BlockB. NLO. v2.1.");
      BlockB_NLO	= blockb;
    }
    else if ( blockb->IContrFlag1==1 && blockb->IContrFlag2==1 && blockb->IRef==1 ){
      if ( blockb->NScaleDep != 3 )   blockb->SetName("BlockB. LO Reference. v2.0.");
      if ( blockb->NScaleDep == 3 )   blockb->SetName("BlockB. LO Reference. v2.1.");
      BlockB_LO_Ref		= blockb;
    }
    else if ( blockb->IContrFlag1==1 && blockb->IContrFlag2==2 && blockb->IRef==1 ){
      if ( blockb->NScaleDep != 3 )   blockb->SetName("BlockB. NLO Reference. v2.0.");
      if ( blockb->NScaleDep == 3 )   blockb->SetName("BlockB. NLO Reference. v2.1.");
      BlockB_NLO_Ref	= blockb;
    }
    else {
      printf("FastNLOReader::ReadTable(). Error in initializing the 'B'-Blocks. Your FastNLO Table file might be incompatible with this reader. Exiting.\n");
      printf("IContrFlag1 = %d , IContrFlag2 = %d , IRef = %d\n", blockb->IContrFlag1==1, blockb->IContrFlag2==2, blockb->IRef==1);
    exit(1);      
    }
  }

  if ( (BlockB_LO_Ref == NULL || BlockB_NLO_Ref == NULL) && BlockB_LO->NScaleDep!=3 ){
    printf("No Reference Tables were found.\n");
  }
  
  if ( BlockB_LO == NULL ){
    printf("ERROR. Could not find any LO Calculation (BlockB_LO).\n");exit(1);
  }


  if ( BlockB_LO )  BBlocks.push_back(BlockB_LO);
  if ( BlockB_NLO ) BBlocks.push_back(BlockB_NLO);

  // todo add NNLO block

  
  //NPDFDim	= BlockB_LO->NPDFDim;

}

//______________________________________________________________________________


void FastNLOReader::ReadBlockA1(istream *table){
  //
  //  Read in information which is called Block A1 from file 
  //

  table->peek();
  if (table->eof()){
    printf("FastNLOReader::Read: Cannot read from file.\n");
    return;
  }
   
  int key = 0;
  *table >> key;
  if(key != tablemagicno){
    printf("FastNLOReader::Read: At beginning of block found %d instead of %d.\n",key,tablemagicno);
    return;
  };
  *table >> Itabversion;
  if ( Itabversion < 20000 ){
    printf("fnloBlockA1::Read. ERROR. This reader is only compatible with FastNLO v2.0 tables and higher.\n");  
    printf("       This FastNLO-table (file) is of version %6.4f.\n",Itabversion/10000.);
    printf("       Please download a compatible reader from the website or use the APPL_grid interface.\n");
    printf("       Exiting.\n");
    exit(1);
  }
  *table >> ScenName;
  *table >> Ncontrib;
  *table >> Nmult;
  *table >> Ndata;
  *table >> NuserString;
  *table >> NuserInt;
  *table >> NuserFloat;
  *table >> Imachine;
  key=0;
  *table >> key;
  if(key != tablemagicno){
    printf("FastNLOReader::Read: At end of block found %d instead of %d.\n",key,tablemagicno);
    return;
  };
  // Put magic number back
  for(int i=0;i<(int)(log10((double)key)+1);i++){
    table->unget();
  }
}


//______________________________________________________________________________


void FastNLOReader::ReadBlockA2(istream *table){
  //
  //  Read in information which is called Block A2 from file 
  //
  
  table->peek();
   if (table->eof()){
      printf("FastNLOReader::Read: Cannot read from file.\n");
      return;
   }

   int key = 0;
   *table >> key;
   if(key != tablemagicno){
      printf("FastNLOReader::Read: At beginning of block found %d instead of %d.\n",key,tablemagicno);
      return;
   };

   *table >> Ipublunits;
   int NScDescript;
   *table >> NScDescript;
   ScDescript.resize(NScDescript);
   char buffer[257];
   table->getline(buffer,256);
   for(int i=0;i<NScDescript;i++){
      table->getline(buffer,256);
      ScDescript[i] = buffer;
      //      StripWhitespace(ScDescript[i]);
   }

   *table >> Ecms;
   *table >> ILOord;
   *table >> NObsBin;
   *table >> NDim;
   DimLabel.resize(NDim);
   table->getline(buffer,256);
   for(int i=0;i<NDim;i++){
      table->getline(buffer,256);
      DimLabel[i] = buffer;
      //      StripWhitespace(DimLabel[i]);
   }

   IDiffBin.resize(NDim);
   for(int i=0;i<NDim;i++){
      *table >>  IDiffBin[i];
   }
   LoBin.resize(NObsBin);
   UpBin.resize(NObsBin);
   //KR: Set rapidity index also when reading a table
   RapIndex.push_back(0);
   //   int irap = 0;
   for(int i=0;i<NObsBin;i++){
      LoBin[i].resize(NDim);
      UpBin[i].resize(NDim);
      for(int j=0;j<NDim;j++){
         *table >>  LoBin[i][j];
         if(IDiffBin[j]==2) *table >>  UpBin[i][j];
      }
      //      cout << "iobs1: " << i << ", LoBin i: " << LoBin[i][1] << endl;
      if ( i > 0 ) {
        if ( LoBin[i][1] != LoBin[i-1][1] ) {
          //      cout << "iobs2: " << i << ", LoBin i-1: " << LoBin[i-1][1] << ", LoBin i: " << LoBin[i][1] << endl;
          RapIndex.push_back(i);
	  //irap++;
	  //cout << "irap: " << irap << ", RapIndex: " << RapIndex[irap] << endl;
        }
      }
   }

   BinSize.resize(NObsBin);
   for(int i=0;i<NObsBin;i++){
     *table >> BinSize[i];
   }

   *table >> INormFlag;
   if(INormFlag>1){
      *table >> DenomTable;
   }
   if(INormFlag>0){
      IDivLoPointer.resize(NObsBin);
      IDivUpPointer.resize(NObsBin);
      for(int i=0;i<NObsBin;i++){
         *table >> IDivLoPointer[i];
         *table >> IDivUpPointer[i];
      }
   }

   key=0;
   *table >> key;
   if(key != tablemagicno){
      printf("FastNLOReader::Read: At end of block found %d instead of %d.\n",key,tablemagicno);
      return;
   };
   // Put magic number back
   for(int i=0;i<(int)(log10((double)key)+1);i++){
      table->unget();
   }
}



//______________________________________________________________________________


void FastNLOReader::PrintFastNLOTableConstants(){
  //
  // Print FastNLO(Reader) internal variables
  //

  PrintBlockA1();
  PrintBlockA2();
  if ( BlockB_LO      ) BlockB_LO	->Print();
  if ( BlockB_NLO     ) BlockB_NLO	->Print();
  if ( BlockB_LO_Ref  ) BlockB_LO_Ref	->Print();
  if ( BlockB_NLO_Ref ) BlockB_NLO_Ref	->Print();
}



//______________________________________________________________________________


void FastNLOReader::PrintBlockA1(){
  printf("\n **************** FastNLO Table: BlockA1 ****************\n\n");
  printf(" A1  tablemagicno                  %d\n",tablemagicno);
  printf(" A1  Itabversion                   %d\n",Itabversion);
  printf(" A1  ScenName                       %s\n",ScenName.data());
  printf(" A1  Ncontrib                      %d\n",Ncontrib);
  printf(" A1  Nmult                         %d\n",Nmult);
  printf(" A1  Ndata                         %d\n",Ndata);
  printf(" A1  NuserString                   %d\n",NuserString);
  printf(" A1  NuserInt                      %d\n",NuserInt);
  printf(" A1  NuserFloat                    %d\n",NuserFloat);
  printf(" A1  Imachine                      %d\n",Imachine);
  printf("\n ********************************************************\n\n");
}



//______________________________________________________________________________


void FastNLOReader::PrintBlockA2(){
  printf("\n **************** FastNLO Table: BlockA2 ****************\n\n");
  printf(" A2  Ipublunits                    %d\n",Ipublunits);
  for(int i=0;i<ScDescript.size();i++){
    printf(" A2  ScDescript[%d]                 %s\n",i,ScDescript[i].data());
  }
  printf(" A2  Ecms                          %7.4f\n",Ecms);
  printf(" A2  ILOord                        %d\n",ILOord);
  printf(" A2  NDim                          %d\n",NDim);
  for(int i=0;i<NDim;i++){
    printf(" A2   - DimLabel[%d]                %s\n",i,DimLabel[i].data());
  }
  for(int i=0;i<NDim;i++){
    printf(" A2   - IDiffBin[%d]               %d\n",i,IDiffBin[i]);
  }
  printf(" A2  NObsBin                       %d\n",NObsBin);
  for(int i=0;i<NObsBin;i++){
    for(int j=0;j<NDim;j++){
      printf(" A2   -  - LoBin[%d][%d]             %7.4f\n", i,j,LoBin[i][j]);
      if(IDiffBin[j]==2)
        printf(" A2   -  - UpBin[%d][%d]             %7.4f\n", i,j,UpBin[i][j]);
    }
   }
   for(int i=0;i<NObsBin;i++){
     printf(" A2   - BinSize[%d]                %7.4f\n", i,BinSize[i]);
   }
   printf(" A2  INormFlag                     %d\n",INormFlag);

   if(INormFlag>1){
     printf(" A2  DenomTable                    %s\n",DenomTable.data());
   }
   if(INormFlag>0){
      for(int i=0;i<NObsBin;i++){
        printf(" A2   - IDivLoPointer[%d]               %d\n",i,IDivLoPointer[i]);
        printf(" A2   - IDivUpPointer[%d]               %d\n",i,IDivUpPointer[i]);
      }
   }
   printf("\n ********************************************************\n\n");
}


//______________________________________________________________________________


void FastNLOReader::PrintCrossSections( ){
  //
  // Print Cross sections in NLO, k-factors and Reference table cross sections
  //
  
  if ( XSection.empty() ){
     CalcCrossSection();
  }
  if ( XSectionRef.empty() && XSectionRef_s1.empty() ){
    CalcReferenceCrossSection();
  }

  vector < double > xs = XSection;

  printf(" *  \n");
  printf(" *  FastNLO Cross sections for\n");
  for ( int i = 0 ; i < ScDescript.size() ; i++ ){
    printf(" *     %s\n",ScDescript[i].c_str());
  }
  printf(" *  at sqrt(s) = %8.2f GeV\n", Ecms);
  printf(" *  \n");
  printf(" *  This is a %s-differential table in %s", ( (NDim==1)?"single":"double"),DimLabel[0].c_str());
  if ( NDim==2 ) printf(" and in %s",DimLabel[1].c_str());
  printf(".\n");
  printf(" *\n");

  string Aunits[16]  = { "[b] --   ","","","[mb] --  ","","","[mu b] --","","","[nb] --  ","","","[pb] --  ","","","[fb] --  "};
  string Nounits[16] = { " --      ","",""," --      ","",""," --      ","",""," --      ","",""," --      ","",""," --      "};
  string* unit = fUnits==kAbsoluteUnits ? Aunits : Nounits;
  

  if ( NDim == 2 ){
    double lobindim2 = -42;
    printf(" *  - Bin - |   ---  %5s  ---        -- XS-FNLO %s -- k-factor -- |\n",DimLabel[0].c_str(),unit[Ipublunits].c_str());
    printf(" *  --------------------------------------------------------------------\n");
    for ( unsigned int i=0;i<xs.size();i++){
      if ( LoBin[i][1] != lobindim2 ){
	printf(" *                  ---->  from %9.3f to %9.3f in %s  <----\n",LoBin[i][1],UpBin[i][1],DimLabel[1].c_str());
	lobindim2 = LoBin[i][1];
      }
      printf(" *   %4.0f   | %9.3f - %9.3f       % 9.4e           % 5.2f      |\n",i*1.,LoBin[i][0],UpBin[i][0],xs[i],kFactor[i]);
    }
  }

  else {
    printf("   ---  %5s  ---        - Bin -       -- XS-FNLO --  \n",DimLabel[0].c_str());
    for ( unsigned int i=0;i<xs.size();i++){
      printf("  %9.3f - %9.3f   %3.0f         % 9.4e\n",LoBin[i][0],UpBin[i][0],i*1.,xs[i]);
    }
  }
  printf(" *  --------------------------------------------------------------------\n");
}


//______________________________________________________________________________


void FastNLOReader::PrintCrossSectionsLikeFreader( ){
  //
  // Print observable binnning and cross sections at
  // LO, NLO and k-factors like in Fortran Reader for comparison
  //
  
  if ( XSection.empty() ){
    CalcCrossSection();
  }
  
  vector < double > xs = XSection;
  
  string CSEP50("##################################################");
  string DSEP50("==================================================");
  string SSEP50("--------------------------------------------------");
  string CSEP = CSEP50 + CSEP50 + CSEP50;
  string DSEP = DSEP50 + DSEP50 + DSEP50;
  string SSEP = SSEP50 + SSEP50 + SSEP50;
  
  cout << DSEP << endl;
  printf(" Cross Sections\n");
  printf(" The scale factor is: \n");
  cout << SSEP << endl;
  
  if ( NDim == 2 ){
    string header[3] = { "  IObs  Bin Size IODim1 ", 
			 "      IODim2",
			 "              LO cross section    NLO cross section   K factor"};
    string label[2] = { "[ " + DimLabel[0] + "     ]", "[ " + DimLabel[1] + "          ]"};
    unsigned int NDimBins[NDim];
    printf("%s %s %s %s %s\n",header[0].c_str(),label[0].c_str(),header[1].c_str(),label[1].c_str(),header[2].c_str());
    cout << SSEP << endl;
    for ( unsigned int i=0; i<xs.size(); i++ ){ 
      for ( unsigned int j=0; j<NDim; j++ ){ 
	if ( i==0 ){
	  NDimBins[j] = 1;
	} else if ( LoBin[i-1][j] < LoBin[i][j]){
	  NDimBins[j]++;
	} else if ( LoBin[i][j] < LoBin[i-1][j]){
	  NDimBins[j] = 1;
	}
      }
      printf(" %5.i  %-#10.4g%5.i  %-#10.4g %-#10.4g %5.i %-#10.4g %-#10.4g %-#18.11E %-#18.11E %-#18.11E\n",
	     i+1,BinSize[i],NDimBins[0],LoBin[i][0],UpBin[i][0],
	     NDimBins[1],LoBin[i][1],UpBin[i][1],xs[i]/kFactor[i],xs[i],kFactor[i]);
    }
  } else {
  }
}


//______________________________________________________________________________


void FastNLOReader::PrintCrossSectionsWithReference( ){
   //
   //  Print Cross sections in NLO, k-factors and Reference table cross sections
   //
   //  Please mention, that the reference cross section can be easily deviating
   //  more than 20% (scales, pdfs, alpha_s, etc...). This does not mean that
   //  the table is wrong!
   //

  
  if ( XSection.empty() ){
     CalcCrossSection();
  }
  if ( XSectionRef.empty() && XSectionRef_s1.empty() ){
    CalcReferenceCrossSection();
  }

  vector < double > xs = XSection;
  vector < double > xsref;
  
  if ( BlockB_LO->NScaleDep == 3 ){
    if ( fMuFFunc == kScale1 && fMuRFunc == kScale1 )			xsref = XSectionRef_s1;
    else if ( fMuFFunc == kScale2 && fMuRFunc == kScale2 )		xsref = XSectionRef_s2;
    else if ( fMuFFunc == kQuadraticMean && fMuRFunc == kQuadraticMean )xsref = XSectionRefMixed;
    else xsref = XSectionRefMixed;
  }
  else xsref = XSectionRef;


  printf(" *  \n");
  printf(" *  FastNLO Cross sections for\n");
  for ( int i = 0 ; i < ScDescript.size() ; i++ ){
    printf(" *     %s\n",ScDescript[i].c_str());
  }
  printf(" *  at sqrt(s) = %8.2f GeV\n", Ecms);
  printf(" *  \n");
  printf(" *  This is a %s-differential table in %s", ( (NDim==1)?"single":"double"),DimLabel[0].c_str());
  if ( NDim==2 ) printf(" and %s",DimLabel[1].c_str());
  printf(" *  \n");
  printf(" *  Please mention, that the reference cross section can easily deviating up to more\n *  than 20%% due to different scale choices, alhpa_s value/evolution, PDFs, etc.");
  printf(" *  This does not mean, that this FastNLO table is wrong!\n\n");
  printf(" *  There are three reference cross sections stored for different scale choices.\n");
  printf(" *  If you have choosen mu_r=mu_f=%s, or mu_r=mu_f=%s or mu_r=mu_f=sqrt((%s^2+%s^2)/2), then you access automatically the corresponding reference cross section.\n",BlockB_NLO->ScaleDescript[0][0].c_str(),BlockB_NLO->ScaleDescript[0][1].c_str(),BlockB_NLO->ScaleDescript[0][0].c_str(),BlockB_NLO->ScaleDescript[0][1].c_str());
  printf(" *  In any other case your reference cross section is calculated using mu_r=mu_f=sqrt((%s^2+%s^2)/2).\n",BlockB_NLO->ScaleDescript[0][0].c_str(),BlockB_NLO->ScaleDescript[0][1].c_str());
  printf(" *  To be fully consistent with the nlojet++ reference cross section, you also have to adjust alpha_s and the alpha_s evolution accordingly.\n\n");

  printf(".\n");
  printf(" *\n");

  string Aunits[16]  = { "[b] --   ","","","[mb] --  ","","","[mu b] --","","","[nb] --  ","","","[pb] --  ","","","[fb] --  "};
  string Nounits[16] = { " --      ","",""," --      ","",""," --      ","",""," --      ","",""," --      ","",""," --      "};
  string* unit = fUnits==kAbsoluteUnits ? Aunits : Nounits;

  if ( NDim == 2 ){
    double lobindim2 = -321312;
    printf(" *  - Bin - |   ---  %5s  ---        -- XS-FNLO %s -- k-factor -- |  -- XS-ref (NLOJET++) --    Diff [%%]\n",DimLabel[0].c_str(),unit[Ipublunits].c_str());
    printf(" *  -----------------------------------------------------------------------------------------------------------\n");
    for ( unsigned int i=0;i<xs.size();i++){
      if ( LoBin[i][1] != lobindim2 ){
	printf(" *                    ---->  from %9.3f to %9.3f in %s  <----\n",LoBin[i][1],UpBin[i][1],DimLabel[1].c_str());
	lobindim2 = LoBin[i][1];
      }
      printf(" *   %4.0f   | %9.3f - %9.3f      % 9.4e           % 5.3f      |     % 9.4e            % 5.4f\n",i*1.,LoBin[i][0],UpBin[i][0],xs[i],kFactor[i],xsref[i],(xs[i]-xsref[i])/xsref[i]*100.);
    }
  }

  else {
    printf("FastNLOReader::PrintCrossSections( ). Info. Single differential printing of cross sections not yet nicely implemented.\n");
    printf("   ---  %s  ---        - Bin -    -- XS-FNLO  --       -- XS-ref (NLOJET++) --    Diff [%%]\n",DimLabel[0].c_str());
    for ( unsigned int i=0;i<xs.size();i++){
      printf("  %9.3f - %9.3f   %3.0f         % 9.4e           % 9.4e          % 5.4f\n",LoBin[i][0],UpBin[i][0],i*1.,xs[i],xsref[i],(xs[i]-xsref[i])/xsref[i]*100.);
    }
  }
  printf(" *  ------------------------------------------------------------------------------------------------------------\n");
}


//______________________________________________________________________________


vector < double > FastNLOReader::GetCrossSection( ){
  // Get fast calculated NLO cross section
  
  if ( XSection.empty() ){
    CalcCrossSection();
  }
  
  return XSection;
}


//______________________________________________________________________________


vector < double > FastNLOReader::GetReferenceCrossSection( ){
  // Get reference cross section from direct nlojet++ calculation
  
  if ( XSectionRef.empty() && XSectionRef_s1.empty() ){
    CalcReferenceCrossSection();
  }
  
  if ( BlockB_LO->NScaleDep == 3 ){
    if ( fMuFFunc == kScale1 && fMuRFunc == kScale1 )			return XSectionRef_s1;
    else if ( fMuFFunc == kScale2 && fMuRFunc == kScale2 )		return XSectionRef_s2;
    else if ( fMuFFunc == kQuadraticMean && fMuRFunc == kQuadraticMean )return XSectionRefMixed;
    else return XSectionRefMixed;
  }
  else return XSectionRef; // XSectionRef from BlockB-Ref
    
}


//______________________________________________________________________________


void FastNLOReader::CalcReferenceCrossSection( ){
  //
  //  Initialize the internal arrays for the reference cross
  //  sections with the information from the FastNLO file
  //
  
  XSectionRef.clear();
  XSectionRef.resize(NObsBin);

  XSectionRefMixed.clear();		  
  XSectionRef_s1.clear();
  XSectionRef_s2.clear();
  XSectionRefMixed.resize(NObsBin);		  
  XSectionRef_s1.resize(NObsBin);
  XSectionRef_s2.resize(NObsBin);

  if ( BlockB_LO_Ref && BlockB_NLO_Ref ){

    for(int i=0;i<NObsBin;i++){
       double unit = fUnits==kAbsoluteUnits ? BinSize[i] : 1.;
       for(int l=0;l<BlockB_LO_Ref->NSubproc;l++){ 
	  XSectionRef[i] +=  BlockB_LO_Ref->SigmaTilde[i][0][0][0][l] * unit; // no scalevariations in LO tables
       }
       for(int l=0;l<BlockB_NLO_Ref->NSubproc;l++){ 
	  XSectionRef[i] +=  BlockB_NLO_Ref->SigmaTilde[i][fScalevar][0][0][l] * unit;
       }
    }
    
  }
  
  if ( BlockB_LO->NScaleDep == 3 ){
    for(int i=0;i<NObsBin;i++){
       double unit = fUnits==kAbsoluteUnits ? BinSize[i] : 1.;
       for(int n=0;n<BlockB_NLO->NSubproc;n++) {
	  XSectionRefMixed[i]		+= BlockB_LO ->SigmaRefMixed[i][n] * unit;
	  XSectionRef_s1[i]		+= BlockB_LO ->SigmaRef_s1[i][n] * unit;
	  XSectionRef_s2[i]		+= BlockB_LO ->SigmaRef_s2[i][n] * unit;
       }
       for(int n=0;n<BlockB_NLO->NSubproc;n++) {
	  XSectionRefMixed[i]		+= BlockB_NLO->SigmaRefMixed[i][n] * unit;
	  XSectionRef_s1[i]		+= BlockB_NLO->SigmaRef_s1[i][n] * unit;
	  XSectionRef_s2[i]		+= BlockB_NLO->SigmaRef_s2[i][n] * unit;
       }
    }
  }
  
  if ( BlockB_LO->NScaleDep != 3 && ( BlockB_NLO_Ref==NULL ) )
    printf("FastNLOReader::CalcReferenceXSection( ). Warning. No reference cross sections available.\n");
       
}


//______________________________________________________________________________


void FastNLOReader::CalcCrossSection( ){
   //
   //  Initialize the internal arrays with the NLO cross sections
   //  with the information from the FastNLO file, the pdf and
   //  the defined alpha_s
   //
   
   XSection_LO.clear();
   XSection.clear();
   XSection_LO.resize(NObsBin);
   XSection.resize(NObsBin);
   kFactor.clear();
   kFactor.resize(NObsBin);
   
   int iLOBs = 0;
   for ( unsigned int i = 0 ; i<BBlocks.size() ; i++ ){
      int tableorder = BBlocks[i]->Npow - ILOord; // 0:LO, 1:NLO, 2:NNLO
      // 			    kAllAvailableOrders       = 0,    // calculate all available orders in this table
      // 			    kLO                       = 1,    // return only LO calculation
      // 			    kFullNLO                  = 2,    // calculate full NLO cross sectoin
      // 			    kNLOOnly                  = 3,    // calculate NLO correction
      // 			    kApproxNNLO               = 4,    // calculate LO+NLO+NNLO(threshold)
      // 			    kNNLOOnly                 = 5,    // calculate only corrections to NLO
      // 			    kHigherOrderCorr          = 6     // calculate higher order corrections to LO
      if ( ( tableorder == 0 && ( fOrder==kLO || fOrder==kFullNLO || fOrder==kApproxNNLO ) ) ||
	   ( tableorder == 1 && ( fOrder==kFullNLO || fOrder==kNLOOnly || fOrder==kApproxNNLO || fOrder==kHigherOrderCorr ) ) || 
	   ( tableorder == 2 && ( fOrder==kNNLOOnly || fOrder==kApproxNNLO ) )||
	   fOrder==kAllAvailableOrders ) {

		// ---- DIS ---- //
	        if ( BBlocks[i]->IPDFdef1 == 2 ){
		   if ( BBlocks[i]->NPDFDim == 0 ) {
		      if ( BlockB_LO->NScaleDep != 3 ){ // v2.0
			 CalcCrossSectionDISv20(BBlocks[i]);
		      }
		      else if ( BlockB_LO->NScaleDep == 3 ){ // v2.1
			 CalcCrossSectionDISv21(BBlocks[i]);
		      }
		   }
		}
		// ---- pp ---- //
		else if (  BBlocks[i]->IPDFdef1 == 3 ){
		   if ( BBlocks[i]->NPDFDim == 1 ) {
		      CalcCrossSectionHHCv20(BBlocks[i]);
		   }
		   else {
		      printf("CalcCrossSection(). only half matrices for hh is implemented.\n"); exit(1);
		   }
		}
		else {
		   printf("CalcCrossSection(). tables not yet implemented.\n");
		}
      }
      
      // calculate LO cross sections
      if ( tableorder == 0 ){
		// ---- DIS ---- //
		if ( BBlocks[i]->IPDFdef1 == 2 ){
		   if ( BBlocks[i]->NPDFDim == 0 ) {
		      iLOBs++;
		      if ( BlockB_LO->NScaleDep != 3 ){
			 CalcCrossSectionDISv20(BBlocks[i],true);
		      }
		      else if ( BlockB_LO->NScaleDep == 3 ){
			 CalcCrossSectionDISv21(BBlocks[i],true);
		      }
		   }
		}
		// ---- pp ---- //
		else if (  BBlocks[i]->IPDFdef1 == 3 ){
		   if ( BBlocks[i]->NPDFDim == 1 ) {
		      iLOBs++;
		      CalcCrossSectionHHCv20(BBlocks[i],true);
		   }
		   else {
		      printf("CalcCrossSection(). only half matrices for hh is implemented.\n"); exit(1);
		   }
		}
      }
   }
   if ( iLOBs != 1 ) printf("CalcCrossSection(). Warning. There are %d LO-tables instead of only one.\n",iLOBs);
   
   // ---- k-factor calculation ---- //
   for(int i=0;i<NObsBin;i++){
      kFactor[i]	= XSection[i] / XSection_LO[i];
   }

}

//______________________________________________________________________________


void FastNLOReader::CalcCrossSectionDISv20( FastNLOBlockB* B , bool IsLO ){
   //
   //  Cross section calculation for DIS tables in v2.0 format
   //
   
   vector<double>* XS = IsLO ? &XSection_LO : &XSection;
   for(int i=0;i<NObsBin;i++){
      int nxmax = B->GetNxmax(i);
      double unit = fUnits==kAbsoluteUnits ? BinSize[i] : 1.;
      for(int j=0;j<B->GetTotalScalenodes();j++){
	 int scalenode1 = j;
	 int scalenode2 = j;
	 if (B->NScaleDim>1){          
	    scalenode1 = j / B->Nscalenode[1];
	    scalenode2 = j % B->Nscalenode[1];
	 }
	 for(int k=0;k<nxmax;k++){ 
	    // LO Block
	    for(int l=0;l<B->NSubproc;l++){ 
	       XS->at(i)		+=  B->SigmaTilde[i][fScalevar][j][k][l] *  B->AlphasTwoPi_v20[i][scalenode2]  *  B->PdfLc[i][scalenode2][k][l] * unit;
	    }
	 }
      }
   }
}

//______________________________________________________________________________


void FastNLOReader::CalcCrossSectionDISv21( FastNLOBlockB* B , bool IsLO){
   //
   //  Cross section calculation for DIS tables in v2.1 format
   //

   vector<double>* XS = IsLO ? &XSection_LO : &XSection;
   for(int i=0;i<NObsBin;i++){
      int nxmax = B->GetNxmax(i);
      double unit = fUnits==kAbsoluteUnits ? BinSize[i] : 1.;
      for(int jS1=0;jS1<B->ScaleNodeScale1[i].size();jS1++){
	 for(int kS2=0;kS2<B->ScaleNodeScale2[i].size();kS2++){
	    double Q2   = B->ScaleNodeScale1[i][jS1]*B->ScaleNodeScale1[i][jS1];
	    
	    double mur	= CalcMu( kMuR , B->ScaleNodeScale1[i][jS1] ,  B->ScaleNodeScale2[i][kS2] , fScaleFacMuR );
	    double muf	= CalcMu( kMuF , B->ScaleNodeScale1[i][jS1] ,  B->ScaleNodeScale2[i][kS2] , fScaleFacMuF );
	    double mur2 = pow(mur,2);
	    double muf2 = pow(muf,2);
	    
	    for(int x=0;x<nxmax;x++){ 
	       for(int n=0;n<B->NSubproc;n++){ 
		  double as	= B->AlphasTwoPi[i][jS1][kS2];
		  double pdflc	= B->PdfLcMuVar[i][x][jS1][kS2][n];
		  XS->at(i)	+=  B->SigmaTildeMuIndep[i][x][jS1][kS2][n] *                     as * pdflc * unit;
		  XS->at(i)	+=  B->SigmaTildeMuFDep [i][x][jS1][kS2][n] * std::log(muf2/Q2) * as * pdflc * unit;
		  XS->at(i)	+=  B->SigmaTildeMuRDep [i][x][jS1][kS2][n] * std::log(mur2/Q2) * as * pdflc * unit;
	       }
	    }
	 }
      }
   }
}

//______________________________________________________________________________


void FastNLOReader::CalcCrossSectionHHCv20( FastNLOBlockB* B , bool IsLO ){
   //
   //  Cross section calculation for HH tables in v2.0 format
   //
   
   int scaleVar = B->Npow == ILOord ? 0 : fScalevar;
   
   vector<double>* XS = IsLO ? &XSection_LO : &XSection;
   for(int i=0;i<NObsBin;i++){
      int nxmax = B->GetNxmax(i);
      double unit = fUnits==kAbsoluteUnits ? BinSize[i] : 1.;
      for(int j=0;j<B->GetTotalScalenodes();j++){
	 int scalenode1 = j;
	 int scalenode2 = j;
	 if (B->NScaleDim>1){          
	    scalenode1 = j / B->Nscalenode[1];
	    scalenode2 = j % B->Nscalenode[1];
	 }
	 for(int k=0;k<nxmax;k++){ 
	    for(int l=0;l<B->NSubproc;l++){ 
	       XS->at(i)		+=  B->SigmaTilde[i][scaleVar][j][k][l] *  B->AlphasTwoPi_v20[i][scalenode2]  *  B->PdfLc[i][scalenode2][k][l] * unit;
// 	       static int count = 0;
// 	       count++;
// 	       printf("%3d  %2d  %3d  %2d  %16.12e    %16.12e    %16.12e\n",i,j,k,l,B->SigmaTilde[i][fScalevar][j][k][l],B->AlphasTwoPi_v20[i][scalenode2],B->PdfLc[i][scalenode2][k][l]);
// 	       if ( count >20 ) exit(1);
	    }
	 }
      }
   }
}

//______________________________________________________________________________


void FastNLOReader::SetUnits( EUnits Unit ){
   if ( fUnits != Unit ){
      fUnits  = Unit;
      CalcCrossSection();
   }
   else {
      // nothing todo
   }
}



//______________________________________________________________________________


void FastNLOReader::SetAlphasMz( double AlphasMz , bool ReCalcCrossSection ){
  //
  //  Set the alpha_s value at M_Z
  //
  
  if ( AlphasMz != fAlphasMz ){
    fAlphasMz	= AlphasMz;		// new alpha_s value
    FillAlphasCache();
    if ( ReCalcCrossSection ) CalcCrossSection(); 
  }
  else {
    // nothing to do!
  }
  
}

//______________________________________________________________________________


void FastNLOReader::FillAlphasCache(){
  //
  //  Fill the internal alpha_s cache.
  //  This is usally called automatically. Only if you
  //  make use of ReFillCache==false options, you have
  //  to take care of this filling by yourself.
  //

   for ( unsigned int i = 0 ; i<BBlocks.size() ; i++ ){
      if ( BBlocks[i]->NScaleDep != 3 ){
	 FillAlphasCacheInBlockBv20( BBlocks[i]  );
      }
      else if ( BBlocks[i]->NScaleDep == 3 ){
	 FillAlphasCacheInBlockBv21( BBlocks[i]  );
      }
   }

}


//______________________________________________________________________________


void FastNLOReader::FillAlphasCacheInBlockBv20( FastNLOBlockB* B ){
   // 
   //  Internal method for filling alpha_s cache
   //
   
   int scaleVar = B->Npow == ILOord ? 0 : fScalevar;

   for(int i=0;i<NObsBin;i++){
      for(int j=0;j<B->GetTotalScalenodes();j++){
	 int scalenode1 = j;
	 int scalenode2 = j;
	 if (B->NScaleDim>1){          
	    scalenode1 = j / B->Nscalenode[1];
	    scalenode2 = j % B->Nscalenode[1];
	 }
	 double mur	= B->ScaleNode[i][0][scaleVar][scalenode1];
	 double as		= GetAlphas(mur);
	 double alphastwopi = pow( as/TWOPI , B->Npow );
	 B->AlphasTwoPi_v20[i][j] = alphastwopi;
      }
   }
}


//______________________________________________________________________________


void FastNLOReader::FillAlphasCacheInBlockBv21( FastNLOBlockB* B ){
   // 
   //  Internal method for filling alpha_s cache
   //

   for(int i=0;i<NObsBin;i++){
      for(int jS1=0;jS1<B->ScaleNodeScale1[i].size();jS1++){
	 for(int kS2=0;kS2<B->ScaleNodeScale2[i].size();kS2++){
	    // 	    double Q2   = B->ScaleNodeScale1[i][jS1]*B->ScaleNodeScale1[i][jS1];
	    // 	    double Pt   = B->ScaleNodeScale2[i][kS2];
	    double mur		= CalcMu( kMuR , BlockB_LO->ScaleNodeScale1[i][jS1] ,  BlockB_LO->ScaleNodeScale2[i][kS2] , fScaleFacMuR );
	    double as		= GetAlphas(mur);
	    double alphastwopi	= pow( as/TWOPI, B->Npow );
	    B->AlphasTwoPi[i][jS1][kS2] = alphastwopi;
	 }
      }
   }
}


//______________________________________________________________________________


double FastNLOReader::GetAlphas( double Q ){
  // 
  //  Internal method for caluclating the alpha_s(mu)
  //
  
  //switch (AlphasEvolution )
  if ( fAlphasEvolution == kGRV )			return GetAlphasGRV	( Q , fAlphasMz );
  else if ( fAlphasEvolution == kNLOJET )		return GetAlphasNLOJET	( Q , fAlphasMz );
  else if ( fAlphasEvolution == kCTEQpdf )		return GetAlphasCTEQpdf	( Q , fAlphasMz );
  else if ( fAlphasEvolution == kFastNLO )		return GetAlphasFastNLO	( Q , fAlphasMz );
  else if ( fAlphasEvolution == kLHAPDFInternal )	return GetAlphasLHAPDF	( Q );
  else if ( fAlphasEvolution == kQCDNUMInternal )	return GetAlphasQCDNUM	( Q );
  else return 0;
}



//______________________________________________________________________________


double FastNLOReader::GetAlphasLHAPDF(double Q){
  //
  // Implementation of Alpha_s evolution as function of Mu_r only.
  //
  // the alpha_s evolution is done within LHAPDF.
  // 
  // WARNING: You cannot change alpha_s(Mz), but is is
  // defined with the pdf.
  //
  
  return LHAPDF::alphasPDF(Q); 
}


//______________________________________________________________________________


double FastNLOReader::GetAlphasQCDNUM(double Q){
   //
   // Sorry! This function is still under developement.
   printf("FastNLOReader::GetAlphasQCDNUM. ERROR. Not yet implemented.\n");
   return 0.;
}


//______________________________________________________________________________


double FastNLOReader::GetAlphasNLOJET(double Q, double alphasMZ){
  //
  // Implementation of Alpha_s evolution as function of Mu_r.
  //
  // this is the evolution, which is used by nlojet++ and cteq6m.
  // Be aware of the Mz-value from 2001 of 91.70.
  //
  // Values used in nlojet++ and cteq6m:
  //   alphasmz   = 0.1179  
  //   double b0  = 1.2202
  //   double b1  = 0.4897
  //   double Mz  = 91.70
  //
  // as evolution by lhpdf.c 
  //

  // the original parameters from the cteq6 pdf
  //   // Alpha QCD //
  //   1, 1, 0, -1, 0.1179, 91.70, 1.3, 4.5, 180.0,

  // #define b0 1.2202
  // #define b1 0.4897
  // #define b2 0.1913
  //   double BETA0 =  (11. - 2./3.*NF); // The beta coefficients of the QCD beta function
  //   double BETA1 =  (51. - 19./3.*NF);

  double b0  = 1.2202;
  double b1  = 0.4897;

  //double Mz	= 91.187;
  double Mz	= 91.70;
  double L = log(Q/Mz);
  L = (b0 + alphasMZ*b1)*L;

  return alphasMZ/(1.0 + alphasMZ*L);

}


//______________________________________________________________________________


double FastNLOReader::GetAlphasGRV(double MU, double ALPSMZ){

   return Alphas::GetAlphasMu(MU,ALPSMZ);

}


//______________________________________________________________________________


double FastNLOReader::GetAlphasCTEQpdf(double Q, double alphasMZ){
  //
  // Implementation of Alpha_s evolution as function of Mu_r.
  //
  // alpha_s evolution as it is within by cteq-pdf-1.0.4 and used in nlojet 4.1.3 without crack

  double as_twopi = alphasMZ/TWOPI;
  double Mz	= 91.187;
  //double Mz	= 91.70;
  //   double Mz	= MZ;

  //int ord=2;
  int nf=5;

  const int NF	= 5;

  double b0 =  (11. - 2./3.*NF); // The beta coefficients of the QCD beta function
  double b1 =  (51. - 19./3.*NF);

  //double astolmd(unsigned int ord, unsigned int nf, double as, double q)
  // ----------------------------------
  double t8 = 1.0/(b0*as_twopi);
  //       /*  it is a leading order evolution  */
  //       if(ord <= 1) return q*exp(-t8);
  /*  at NLO or higer order level returns with the NLO alpha_s  */
  double as0, as1, ot, lt, br = (51.0 - 19.0/3.0*nf)/(b0*b0);
  do {
    lt = log(2.0*t8)/t8;
    ot = t8;

    as0 = (1.0 - br*lt)/(b0*t8);
    as1 = (-1.0 - br*(1.0/t8-2.0*lt))/(b0*t8*t8);
    t8 += (as_twopi - as0)/as1;
  } while(fabs(ot-t8)/ot > 1e-5);
  double lmd = Mz*exp(-t8);
  // ----------------------------------

  double t = log(Q/lmd);
  double asMz = 1.0/(b0*t);

  return asMz*(1.0-b1/b0*asMz*log(2.0*t)) *TWOPI;

}



//______________________________________________________________________________


double FastNLOReader::GetAlphasFastNLO(double Q, double alphasMZ){
  //
  // Implementation of Alpha_s evolution as function of Mu_r.
  //
  // This is the 'original' FNLOv2.0 implementation.
  //

  const int NF	= 5;
  double BETA0 =  (11. - 2./3.*NF); // The beta coefficients of the QCD beta function
  double BETA1 =  (51. - 19./3.*NF);

  //    // This is from NLOJET++, alpha.cc
  double Mz	= 91.187;
  double res	= alphasMZ;
  double b0	= BETA0/TWOPI;
  double w = 1.0 + b0*alphasMZ*log(Q/Mz);
  res /= w;
  double b1 = BETA1/TWOPISQR;
  res *= 1.0 - alphasMZ*b1/b0*log(w)/w;
  return res;

}




//______________________________________________________________________________


void FastNLOReader::FillPDFCache( bool ReCalcCrossSection ){
   //
   //  Fill the internal pdf cache.
   //  This function has to be called by the user, since the 
   //  pdf parameters and evolutions are calculated externally.
   //
   
   if ( fPDFInterface == kLHAPDF ){
      if ( fLHAPDFfilename == ""){
	 printf("FastNLOReader::FillPDFCache(). ERROR. You must specify a LHAPDF filename first or you have to specify kH1FITTER..\n"); exit(1);
      }
      InitLHAPDF();
   }
   else if ( fPDFInterface == kH1FITTER ){
     evolution_();
   }
   
   for ( unsigned int i = 0 ; i<BBlocks.size() ; i++ ){
      if (BBlocks[i]->NScaleDim>1){
	 printf("FastNLOReader::FillBlockBPDFLCsWithLHAPDF. WOW! NScaleDim>1! This is usually not the case!\n");
	 //scaleindex2 = 1; // If we use multiple scales, then mu_f is by convention the second scale -> index=1 
	 //fScalevar2 = fScalevar % NfScalevar[1]; 
      }
      
      // linear: DIS-case
      // ---- DIS ---- //
      if ( BBlocks[i]->IPDFdef1 == 2 ){
	 if ( BBlocks[i]->NPDFDim == 0 ) {
	    if	 ( BBlocks[i]->NScaleDep != 3 )		FillBlockBPDFLCsDISv20(BBlocks[i]);
	    else if ( BBlocks[i]->NScaleDep == 3 )	FillBlockBPDFLCsDISv21(BBlocks[i]);
	 }
      }
      // ---- pp ---- //
      else if (  BBlocks[i]->IPDFdef1 == 3 ){
	 if ( BBlocks[i]->NPDFDim == 1 ) {
	    FillBlockBPDFLCsHHCv20(BBlocks[i]);
	 }
	 else {
	    printf("FastNLOReader::FillBlockBPDFLCs(). only half matrices for hh is implemented.\n"); exit(1);
	 }
      }
      else {
	 printf("FastNLOReader::FillBlockBPDFLCs(). tables not yet implemented.\n");
      }
   }   
   if ( ReCalcCrossSection ) CalcCrossSection();

}



//______________________________________________________________________________


void FastNLOReader::InitLHAPDF(){
  //
  //  Initalize some necessary LHAPDF parameters
  //
    
  if ( fLHAPDFfilename == ""){
    printf("FastNLOReader::FillPDFCacheLHAPDF(). ERROR. You must specify a LHAPDF filename first.\n"); exit(1);
  }

  //string LHAPDFfile = fLHAPDFpath+"/"+fLHAPDFfilename;

  // ---- check if file exists ----- //
  //   FILE* fp = fopen(LHAPDFfile.c_str(), "r");
  //   if (fp) {
  //     fclose(fp);
  //   } else {
  //     printf("Error. LHAPDF file does not exists. Was looking in/for: %s.\n Exiting.\n",LHAPDFfile.c_str());
  //     exit(1);
  //   } 

  //LHAPDF::setVerbosity(LHAPDF::SILENT);
  LHAPDF::setVerbosity(LHAPDF::LOWKEY);
  //cout << " LHAPDF version: " << LHAPDF::getVersion() <<endl;
  LHAPDF::initPDFSetByName(fLHAPDFfilename);
  fnPDFs = LHAPDF::numberPDF();
  if ( fnPDFs < fiPDFSet ){
    cout << "Error. There are only " << fnPDFs << " pdf sets within this LHAPDF file. You were looking for set number " << fiPDFSet << endl;
  }

  LHAPDF::initPDF(fiPDFSet);

}


//______________________________________________________________________________


void FastNLOReader::FillBlockBPDFLCsDISv20( FastNLOBlockB* B ){
   int scaleVar = B->Npow == ILOord ? 0 : fScalevar;
   vector<double> xfx(13); // PDFs of all partons
   if ( B->NScaleDep != 3 ){
      for(int i=0;i<NObsBin;i++){
	 int nxmax = B->GetNxmax(i);
	 for(int j=0;j<B->Nscalenode[0];j++){
	    for(int k=0;k<nxmax;k++){ 
	       
	       double xp	= B->XNode1[i][k];
	       double muf	= B->ScaleNode[i][0][scaleVar][j];
	       xfx = GetXFX(xp,muf);
	       //xfx = LHAPDF::xfx(xp,muf); // LHAPDF::xfx_p_(x,muf,0,0)
	       vector < double > buffer = CalcPDFLinearCombDIS( xfx , B->NSubproc );
	       for(int l=0;l<B->NSubproc;l++){ 
		  B->PdfLc[i][j][k][l] = buffer[l];
	       }
	    }
	 }
      }
   }
}

//______________________________________________________________________________


void FastNLOReader::FillBlockBPDFLCsDISv21( FastNLOBlockB* B ){
   
   if ( B->PdfLcMuVar.empty() ) { cout<< "empty."<<endl; exit(1);}// [i][x][jS1][kS2][l]

   vector<double> xfx(13); // PDFs of all partons

   for(int i=0;i<NObsBin;i++){
      int nxmax = B->GetNxmax(i);
      
      // speed up! if mu_f is only dependent on one variable, we can safe the loop over the other one

      for(int x=0;x<nxmax;x++){ 
	 double xp	= B->XNode1[i][x];
	
	 if ( fMuFFunc != kScale1 &&  fMuFFunc != kScale2 ) { // that't the standard case!
	    for(int jS1=0;jS1<B->ScaleNodeScale1[i].size();jS1++){
	       for(int kS2=0;kS2<B->ScaleNodeScale2[i].size();kS2++){
		  double muf = CalcMu( kMuF , BlockB_LO->ScaleNodeScale1[i][jS1] ,  BlockB_LO->ScaleNodeScale2[i][kS2] , fScaleFacMuF );
		  
		  xfx = GetXFX(xp,muf);
		  //xfx = LHAPDF::xfx(xp, muf); // LHAPDF::xfx_p_(x,muf,0,0)
		  vector < double > buffer = CalcPDFLinearCombDIS( xfx , B->NSubproc );
		  for(int l=0;l<B->NSubproc;l++){ 
		     B->PdfLcMuVar[i][x][jS1][kS2][l] = buffer[l];
		  }
	       }
	    }
	 }
	 else if ( fMuFFunc == kScale2 ){	// speed up
	    for(int kS2=0;kS2<B->ScaleNodeScale2[i].size();kS2++){
	       double muf = CalcMu( kMuF , 0 ,  BlockB_LO->ScaleNodeScale2[i][kS2] , fScaleFacMuF );
	       xfx = GetXFX(xp,muf);
	       //xfx = LHAPDF::xfx(xp, muf); // LHAPDF::xfx_p_(x,muf,0,0)
	       vector < double > buffer = CalcPDFLinearCombDIS( xfx , B->NSubproc );
	       for(int jS1=0;jS1<B->ScaleNodeScale1[i].size();jS1++){
		  for(int l=0;l<B->NSubproc;l++){ 
		     B->PdfLcMuVar[i][x][jS1][kS2][l] = buffer[l];
		  }
	       }
	    }
	 }
	 else if ( fMuFFunc == kScale1 ){	// speed up
	    for(int jS1=0;jS1<B->ScaleNodeScale1[i].size();jS1++){
	       double muf = CalcMu( kMuF , BlockB_LO->ScaleNodeScale1[i][jS1] , 0 , fScaleFacMuF );
	       xfx = GetXFX(xp,muf);
	       //xfx = LHAPDF::xfx(xp, muf); // LHAPDF::xfx_p_(x,muf,0,0)
	       vector < double > buffer = CalcPDFLinearCombDIS( xfx , B->NSubproc );
	       for(int kS2=0;kS2<B->ScaleNodeScale2[i].size();kS2++){
		  for(int l=0;l<B->NSubproc;l++){ 
		     B->PdfLcMuVar[i][x][jS1][kS2][l] = buffer[l];
		  }
	       }
	    }
	 }
      }
   }
}


//______________________________________________________________________________


void FastNLOReader::FillBlockBPDFLCsHHCv20( FastNLOBlockB* B ){
   int scaleVar = B->Npow == ILOord ? 0 : fScalevar;
   vector < vector < double > > xfx; // PDFs of all partons
   if ( B->NScaleDep != 3 ){
      for(int i=0;i<NObsBin;i++){
	 int nxmax = B->GetNxmax(i);
         int nxbins1 = B->Nxtot1[i]; // number of columns in half matrix
         xfx.resize(nxbins1);
	 for(int j=0;j<B->Nscalenode[0];j++){
	    //for(int k=0;k<nxmax;k++){ 
	       // determine all pdfs of hadron1
	       for(int k=0;k<nxbins1;k++){ 
		  double xp	= B->XNode1[i][k];
		  double muf	= B->ScaleNode[i][0][scaleVar][j];
		  xfx[k]	= GetXFX(xp,muf);
	       }
	       int x1bin = 0;
	       int x2bin = 0;
	       for(int k=0;k<nxmax;k++){ 
		  // 		  vector < double > buffer  = CalcPDFLinearComb(xfx[x1bin],xfx[x2bin]); //calculate linear combinations
		  // 		  for(int l=0;l<B->NSubproc;l++){ 
		  // 		     PdfLc[i][j][k][l] = buffer[l];
		  // 		  }
		  // ----- if pp ---- //
		  if ( B->NPDFPDG[0] == B->NPDFPDG[1] ){
		     B->PdfLc[i][j][k] = CalcPDFLinearCombHHC( xfx[x1bin], xfx[x2bin], B->NSubproc ) ;//CalcPDFLinearComb(xfx[x1bin],xfx[x2bin],B->IPDFdef1, B->IPDFdef2, B->NSubproc); //calculate linear combinations
		  }
		  // ----- if ppbar ---- //
		  else if ( B->NPDFPDG[0] == -B->NPDFPDG[1] ){
		     vector < double > xfxbar(13);
		     for ( unsigned int p = 0 ; p<13 ; p++ ){
			xfxbar[p] = xfx[x2bin][12-p];
		     }
		     B->PdfLc[i][j][k] = CalcPDFLinearCombHHC( xfx[x1bin], xfxbar, B->NSubproc ) ;
		  }
		  else {
		     printf("FastNLOReader::FillBlockBPDFLCsHHCv20(). This is not pp, nor ppbar, nor pbarpbar!\n"); exit(1);
		  }
		  x1bin++;
		  if(x1bin>x2bin){
		     x1bin = 0;
		     x2bin++;
		  }
	       }
// 	       //xfx = LHAPDF::xfx(xp,muf); // LHAPDF::xfx_p_(x,muf,0,0)
// 	       vector < double > buffer = CalcPDFLinearComb(xfx,xfx,B->IPDFdef1, B->IPDFdef2, B->NSubproc ); //calculate linear combinations
// 	       for(int l=0;l<B->NSubproc;l++){ 
// 		  B->PdfLc[i][j][k][l] = buffer[l];
// 	       }
//	    }
	 }
      }
   }
}


//______________________________________________________________________________



vector<double> FastNLOReader::GetXFX(double xp, double muf){
   //
   //  Internal method.
   //  GetXFX is used to get the parton array from the
   //  pre-defined pdf-interface.
   // 
  
   if ( fPDFInterface == kLHAPDF ){
      return LHAPDF::xfx(xp,muf);
   }
   else if ( fPDFInterface == kH1FITTER ){
      int iqnset = 1;
      int iqnchk = 0;
      double muf2	= muf*muf;
      vector < double > a;
      a.resize(13);
      fpdfxq_(&iqnset, &xp, &muf2, &a[0], &iqnchk); 
      return a;
   }
   else {
      vector < double > a;
      return a;
   }
}


//______________________________________________________________________________


vector<double> FastNLOReader::CalcPDFLinearCombDIS( vector<double> pdfx1 , int NSubproc){
   //
   //  Internal method.
   //  CalcPDFLinearComb is used to calculate the 
   //  linear combinations of different parton flavors
   //  according to the used subprocesses of your
   //  calculation/table.
   //
  
   vector < double > pdflc;
   pdflc.resize(3);

   pdflc[1] = pdfx1[6]; //gluon
  
   for(int l=0;l<13;l++){
      double temp = (l==6 ? 0.0 : pdfx1[l]);
      if (!(l&1)) temp *= 4.;
      pdflc[0] += temp; // delta
   }
   pdflc[0] /= 9.;
  
   if(NSubproc>2){ // only from NLO
      for(int l=0;l<6;l++){
	 pdflc[2] += pdfx1[5-l] + pdfx1[l+7]; // sigma
      }
   }

   return pdflc;
    
}



//______________________________________________________________________________


vector<double> FastNLOReader::CalcPDFLinearCombHHC( vector<double> pdfx1 , vector<double> pdfx2 , int NSubproc){
   //
   //  Internal method.
   //  CalcPDFLinearComb is used to calculate the 
   //  linear combinations of different parton flavors
   //  according to the used subprocesses of your
   //  calculation/table.
   //
  

   vector < double > njpdflc;
   njpdflc.resize(7);
  
   // ---------------------------------------------------
   // remember from nlojet++ proc-hhc/weight.cc
   //   const char *weight_label_hhc[7] = {"gg", "qg", "gq", "qr", "qq", "qqb", "qrb"};
   // ---------------------------------------------------
   // remember from nlojet++ proc-hhc/process.cc
   //    njpdflc[0] = A0*B0;
   //    njpdflc[1] = (A + Ab)*B0;
   //    njpdflc[2] = A0*(B + Bb);
   //    njpdflc[3] = A*B + Ab*Bb - D;
   //    njpdflc[4] = D;
   //    njpdflc[5] = Db;
   //    njpdflc[6] = A*Bb +Ab*B - Db;
   // ---------------------------------------------------
   unsigned int nu	= 2;
   unsigned int nd	= 3;
   
   int ia, iq;
   //weight_hhc njpdflc;
   //   static double __f1[13], __f2[13];
   static const int iu[3] = {2,4,6}, id[3] = {1,3,5};
   
   //   //----- calculat the pdfs -----
   //   double *f1 = __f1+6, *f2 = __f2+6;
   
   //   this -> hadronA(x1, mf2, nu, nd, f1);
   //   this -> hadronB(x2, mf2, nu, nd, f2);
   
   //----- gluon pdfs -----
   double A0 = pdfx1[0+6];
   double B0 = pdfx2[0+6];
   
   
   //---- up type quarks -----
   double q1, q2, a1, a2;
   double A = 0.0, B = 0.0, Ab = 0.0, Bb = 0.0, D = 0.0, Db = 0.0; 
   
   for(unsigned int u = 0; u < nu && u < 3; u++) {
      ia = -(iq = iu[u]);
      q1 = pdfx1[iq+6]; q2 = pdfx2[iq+6];
      a1 = pdfx1[ia+6]; a2 = pdfx2[ia+6];
      
      A += q1; Ab += a1; B += q2; Bb += a2;
      D += q1*q2 + a1*a2; Db += q1*a2 + a1*q2;
   }
    
   //----- down type quarks -----
   for(unsigned int d = 0; d < nd && d < 3; d++) {
      ia = -(iq = id[d]);
      //cout << "ia="<<ia<<"\tiq="<<iq<<"\td="<<d<< endl;
      q1 = pdfx1[iq+6]; q2 = pdfx2[iq+6];
      a1 = pdfx1[ia+6]; a2 = pdfx2[ia+6];
      
      A += q1; Ab += a1; B += q2; Bb += a2;
      D += q1*q2 + a1*a2; Db += q1*a2 + a1*q2;
   }
    
   njpdflc[0] = A0*B0;
   njpdflc[1] = (A + Ab)*B0;
   njpdflc[2] = A0*(B + Bb);
   njpdflc[3] = A*B + Ab*Bb - D;
   njpdflc[4] = D;
   njpdflc[5] = Db;
   njpdflc[6] = A*Bb +Ab*B - Db;
   
   
   vector < double > retval(7);
   
   retval[0] = njpdflc[0];
   retval[1] = njpdflc[3];
   retval[2] = njpdflc[4];
   retval[3] = njpdflc[5];
   retval[4] = njpdflc[6];
   retval[5] = njpdflc[1];
   retval[6] = njpdflc[2];
   
   if (NSubproc == 6) {
      retval[5] += retval[6];    // -- sum is correct in ref-mode
      retval.resize(6);
   }

   return retval;
    
}



//______________________________________________________________________________



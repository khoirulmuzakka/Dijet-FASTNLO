// Author: Daniel Britzger
// DESY, 23/07/2011

//////////////////////////////////////////////////////////////////////////
//                                                                      //
//  fastNLO_reader_2.1.0                                                //
//  D. Britzger, T. Kluge, K. Rabbertz, F. Stober, M. Wobisch           //
//                                                                      //
//  The projects web page can be found at:                              // 
//    http://projects.hepforge.org/fastnlo                              //
//                                                                      //
//  If you use this code, please cite:                                  //
//    T. Kluge, K. Rabbertz and M. Wobisch, hep-ph/0609285              //
//    D. Britzger, T. Kluge, K. Rabbertz, F. Stober, M. Wobisch,        //
//       arXiv:1109.1310                                                //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "FastNLOReader.h"

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <cfloat>
#include "Alphas.h"

using namespace std;


//______________________________________________________________________________


// some names for nice output
const string FastNLOReader::fContrName[20] = {
  "Fixed order calculation", "Threshold corrections", "Electroweak corrections", "Non-perturbative corrections",
  "Undefined", "Undefined", "Undefined", "Undefined", "Undefined", "Undefined", "Undefined",
  "Quark compositeness", "ADD-LED", "TeV 1-ED", "Unknown", "Unknown", "Unknown", "Unknown", "Unknown", "Unknown" };
const string FastNLOReader::fOrdName[4][4] = { { "LO",     "NLO",    "NNLO"   , "N3LO"    },
					       { "1-loop", "2-loop", "3-loop" , "4-loop"  },
					       { "Undef" , "Undef" , "Undef"  , "Undef"   },
					       { "LO MC" , "NLO MC", "NNLO MC", "N3LO MC" } };
const string FastNLOReader::fNSDep[4] = {"v2.0","v2.0","v2.0","v2.1"};
int FastNLOReader::WelcomeOnce = 0;

//______________________________________________________________________________


FastNLOReader::FastNLOReader(string filename)
{
  InitMembers();
  SetFilename(filename);
}


//______________________________________________________________________________


FastNLOReader::~FastNLOReader(void)
{
  for ( unsigned int j = 0 ; j<BBlocksSMCalc.size() ; j++ ){
    if ( !BBlocksSMCalc[j].empty() ){
      for ( unsigned int i = 0 ; i<BBlocksSMCalc[j].size() ; i++ ){
	delete BBlocksSMCalc[j][i];
      }
      BBlocksSMCalc.clear();
    }
  }
}


//______________________________________________________________________________



void FastNLOReader::InitMembers(){
  BlockB_Data		= NULL;
  BlockB_LO_Ref		= NULL;
  BlockB_NLO_Ref	= NULL;
  fUnits		= kPublicationUnits;
  fAlphasMz		= 0.118500001;
  fMuRFunc		= FastNLOReader::kScale1;
  fMuFFunc		= FastNLOReader::kScale1;
}


//______________________________________________________________________________



void FastNLOReader::SetAlphasEvolution(EAlphasEvolution AlphasEvolution){
  if (AlphasEvolution==kExternAs ) {
    cout << "FastNLOReader::SetAlphasEvolution. Info. Alphas(Mz) is received from some derive class."<<endl; 
    cout << "    SetAlphasMz(double) and thus the internal Alphas(Mz) value might have no more influence."<<endl;
  }

  fAlphasEvolution = AlphasEvolution; 
  FillAlphasCache();
}


//______________________________________________________________________________



void FastNLOReader::SetGRVtoPDG2011_2loop(bool Print){
   printf("FastNLOReader::SetAlphasEvolution. Info. Resetting GRV Alphas::Alphas evolution.\n");
   Alphas::SetMz(91.1876); // PDG 2011
   Alphas::SetNf(5);
   Alphas::SetNLoop(2);
   Alphas::SetFlavorMatchingOn(true);
   if ( Print )  Alphas::PrintInfo();
}


//______________________________________________________________________________



void FastNLOReader::SetFilename(string filename){
  ffilename	= filename;
  Init();
}


//______________________________________________________________________________



void FastNLOReader::Init(){
  ReadTable();
  //int iprint = 2;
  //PrintFastNLOTableConstants(iprint);
  InitScalevariation();
  //SetAlphasEvolution(FastNLOReader::kGRV);
}


//______________________________________________________________________________



void FastNLOReader::InitScalevariation(){
  
  fScaleFacMuR	= 1.;
  fScaleFacMuF	= 1.;
  fScalevar	= -1;

  if ( !GetIsFlexibleScaleTable() ){
     for (int iscls=0; iscls<GetNScaleVariations(); iscls++){
	const double muFac = BBlocksSMCalc[0][1]->ScaleFac[0][iscls];
 	if (abs(muFac-1.0) < 1.e-7){
	   SetScaleVariation(iscls,false);
	   break;
	}
     }
     if ( fScalevar == -1 ){
	printf(" FastNLOReader::InitScalevariation(). Error. Could not found scale variation with scale factor 1.0. Exiting.\n");
	exit(1);
     }
  }
  else {
    // this is a MuVar table. You can vary mu_f and mu_r independently by any factor
    // and you can choose the functional form of mu_f and mu_r as functions of
    // scale1 and scale1 (called partly scaleQ2 and scalePt).
    
    if ( BBlocksSMCalc[0][0]->ScaleDescript[0].size() <0 ) {
      printf("Error. No scaledescription available.\n"); // the code will crash soon.
      SetFunctionalForm( kScale1 , kMuR );
      SetFunctionalForm( kScale1 , kMuF );
      return;
    }

    // ---- DIS ---- //
    if ( BBlocksSMCalc[0][0]->NPDFDim == 0 ) {
      SetFunctionalForm( kQuadraticMean , kMuR );
      SetFunctionalForm( kScale1 , kMuF );
    }
    // ---- HHC --- //
    else if (  BBlocksSMCalc[0][0]->NPDFDim == 1 ) {
      SetFunctionalForm( kScale1 , kMuR );
      SetFunctionalForm( kScale1 , kMuF );
    }
    else {
      printf("Error. Unknown process.\n");
      exit(1);
    }
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
  
  EScaleFunctionalForm Func = (kMuX == FastNLOReader::kMuR) ? fMuRFunc : fMuFFunc;
  
  double mu = 0;

  if		( Func == kScale1 )		mu	= scale1;
  else if	( Func == kScale2 )		mu	= scale2;
  else if	( Func == kQuadraticSum )	mu	= FuncMixedOver1(scale1,scale2);
  else if	( Func == kQuadraticMean )	mu	= FuncMixedOver2(scale1,scale2);
  else if	( Func == kQuadraticSumOver4 )	mu	= FuncMixedOver4(scale1,scale2);
  else if	( Func == kLinearMean )		mu	= FuncLinearMean(scale1,scale2);
  else if	( Func == kLinearSum )		mu	= FuncLinearSum(scale1,scale2);
  else if	( Func == kScaleMax )		mu	= FuncMax(scale1,scale2);
  else if	( Func == kScaleMin )		mu	= FuncMin(scale1,scale2);
  else if	( Func == kExpProd2 )		mu	= FuncExpProd2(scale1,scale2);
  else if	( Func == kExtern  )		mu	= (kMuX==FastNLOReader::kMuR) ? (*Fct_MuR)(scale1,scale2) : (*Fct_MuF)(scale1,scale2);
  else printf( "Error. could not identify functional form for scales calculation.\n");
  
  return scalefac * mu;

}


//______________________________________________________________________________
double FastNLOReader::FuncMixedOver1(double scale1 , double scale2 ){
  return ( sqrt( (pow(scale1,2) + pow(scale2,2))  / 1. ) );
}

//______________________________________________________________________________
double FastNLOReader::FuncMixedOver2(double scale1 , double scale2 ){
  return ( sqrt( (pow(scale1,2) + pow(scale2,2))  / 2. ) );
}

//______________________________________________________________________________
double FastNLOReader::FuncMixedOver4(double scale1 , double scale2 ){
  return ( sqrt( (pow(scale1,2) + pow(scale2,2))  / 4. ) );
}

//______________________________________________________________________________
double FastNLOReader::FuncLinearMean(double scale1 , double scale2 ){
  return ( scale1 + scale2 ) / 2.;
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
double FastNLOReader::FuncExpProd2(double scale1 , double scale2 ){
  return ( scale1 * exp(0.3*scale2) );
}




//______________________________________________________________________________



double FastNLOReader::SetScaleVariation(int scalevar , bool ReFillCache , bool Verbose ){ 
  
  // ------------------------------------------------
  //   Set the scalevariation factor for determining the
  //   'theory'-error. Usually, you have tables stored with
  //   factors of 0.5, 1 and 2 times the nominal scale.
  //     corresponding to:
  //     scalevar -> scalefactor
  //        '0'   ->   1.00
  //        '1'   ->   0.50
  //        '2'   ->   2.00
  //   This method returns the scalefactor correspoding to
  //   the chosen 'scalevar'.
  // ------------------------------------------------

   if ( GetIsFlexibleScaleTable() ){
      printf("FastNLOReader::SetScaleVariation(). Info: This is a flexible-scale table!\n");
      printf("  You can choose freely (within reason) a factorization scale factor. Your Scalevar has to be '0'.\n");
      printf("  Please use SetScaleFacMuR(double) and SetScaleFacMuF(double) to set scale factors.\n");
      return 0;
   }

  // Check for maximal scale variation of all rel. and active SM calcs
  int scalevarmax = GetNScaleVariations(); 
  if ( scalevar >= scalevarmax ){
    printf("FastNLOReader::SetScaleVariation. WARNING! This table has only %d scale variation(s) stored!\n",scalevarmax);
    printf("  for the currently active contributions. You wanted to access the non-existing number %d.\n",scalevar);
    printf("  Using '0' instead.\n");
    fScalevar = 0;
    return B_NLO()->ScaleFac[0][0];
  }
  
  fScalevar	= scalevar;
  fScaleFacMuF	= B_NLO()->ScaleFac[0][fScalevar];
  if (Verbose) {printf("FastNLOReader::SetScaleVariation: Selecting MuF table according to a multiplicative scale factor of the factorization scale of %4.2f times the nominal scale.\n",fScaleFacMuF);}
  
  if ( !BBlocksSMCalc[kThresholdCorrection].empty() ){
    bool lkth = false;
    for ( unsigned int i = 0 ; i <BBlocksSMCalc[kThresholdCorrection].size() ; i++ ){
      if ( bUseSMCalc[kThresholdCorrection][i] ) {lkth = true;}
    }
    if ( lkth && abs(fScaleFacMuR-fScaleFacMuF) > DBL_MIN ) {
      fScaleFacMuR = fScaleFacMuF;
      printf("FastNLOReader::SetScaleVariation. WARNING! Threshold corrections do not allow\n");
      printf("  a posteriori variations of the renormalization scale!\n");
      printf("  The scale factor for MuR has been set equal to the one for MuF = %7.3f\n",fScaleFacMuF);
      printf("  Either select a different simultaneous scale variation i, if possible, via\n");
      printf("  FastNLOReader::SetScaleVariation(i);\n");
      printf("  or deactivate first all threshold corrections using\n");
      printf("  FastNLOReader::SetContributionON(kTresholdCorrections,Id,false).\n");
    }
  }

  FillAlphasCache();
  if ( ReFillCache ){
    FillPDFCache();
  }

  return B_NLO()->ScaleFac[0][fScalevar];
  
}



//______________________________________________________________________________



void FastNLOReader::SetFunctionalForm( EScaleFunctionalForm func , FastNLOReader::EMuX kMuX , bool Verbose ){
  //
  //  For MuVar tables this method sets the functional form of
  //  the renormalization or the factorization scale.
  //     func:  Choose a pre-defined function
  //     kMuX:  is it for mu_r or for mu_f ?
  //

  if ( BBlocksSMCalc[0][0]->NScaleDep != 3 ) {
    printf("FastNLOReader::SetFunctionalForm. WARNING! This is not a flexible-scale table.\n");
    printf("  SetFunctionalForm has no impact.\n");
    printf("  Please use a flexible-scale table, if you want to change your scale-definition.\n");
    return;
  }


  // ---- prepare printout ---- //
  const string sname[2] = {"renormalization","factorization"};
  const string smu[2] = {"mu_r","mu_f"};
  const int isc = kMuX==kMuR ? 0 : 1;
  char fname[100];

  switch (func){
  case kScale1: sprintf(fname,"%s^2",BBlocksSMCalc[0][0]->ScaleDescript[0][0].c_str());
     break;
  case kScale2: sprintf(fname,"%s^2",BBlocksSMCalc[0][0]->ScaleDescript[0][1].c_str());
     break;
  case kQuadraticSum: sprintf(fname,"%s^2 + %s^2",BBlocksSMCalc[0][0]->ScaleDescript[0][0].c_str(),BBlocksSMCalc[0][0]->ScaleDescript[0][1].c_str());
     break;
  case kQuadraticMean: sprintf(fname,"(%s^2 + %s^2)/2",BBlocksSMCalc[0][0]->ScaleDescript[0][0].c_str(),BBlocksSMCalc[0][0]->ScaleDescript[0][1].c_str());
     break;
  case kQuadraticSumOver4: sprintf(fname,"(%s^2 + %s^2)/4",BBlocksSMCalc[0][0]->ScaleDescript[0][0].c_str(),BBlocksSMCalc[0][0]->ScaleDescript[0][1].c_str());
     break;
  case kLinearMean: sprintf(fname,"((%s+%s)/2)^2",BBlocksSMCalc[0][0]->ScaleDescript[0][0].c_str(),BBlocksSMCalc[0][0]->ScaleDescript[0][1].c_str());
     break;
  case kLinearSum: sprintf(fname,"(%s+%s)^2",BBlocksSMCalc[0][0]->ScaleDescript[0][0].c_str(),BBlocksSMCalc[0][0]->ScaleDescript[0][1].c_str());
     break;
  case kScaleMax: sprintf(fname,"max(%s^2,%s^2)",BBlocksSMCalc[0][0]->ScaleDescript[0][0].c_str(),BBlocksSMCalc[0][0]->ScaleDescript[0][1].c_str());
      break;
  case kScaleMin: sprintf(fname,"min(%s^2,%s^2)",BBlocksSMCalc[0][0]->ScaleDescript[0][0].c_str(),BBlocksSMCalc[0][0]->ScaleDescript[0][1].c_str());
     break;
  case kExpProd2: sprintf(fname,"(%s*exp(0.3*%s)^2)",BBlocksSMCalc[0][0]->ScaleDescript[0][0].c_str(),BBlocksSMCalc[0][0]->ScaleDescript[0][1].c_str());
     break;
  case kExtern: sprintf(fname,"f_ext(%s,%s)",BBlocksSMCalc[0][0]->ScaleDescript[0][0].c_str(),BBlocksSMCalc[0][0]->ScaleDescript[0][1].c_str());
     break;
  default: printf("unknown scale choice.\n");
  }
  
   
  // ---- setting scale ---- //
  if ( Verbose ) {printf (" *    Setting %s scale to %s^2 = %1.2f^2 * %s.\n",
			  sname[isc].c_str(),smu[isc].c_str(), 
			  (kMuX==kMuR?fScaleFacMuR:fScaleFacMuF), fname);}
  if ( kMuX == kMuR ) fMuRFunc = func;
  else fMuFFunc = func;


  // ---- cross check ---- //
  if	( func == kScale2 || func == kQuadraticSum ||  func == kQuadraticMean || func == kQuadraticSumOver4 
	  || func == kLinearMean || func == kLinearSum  ||  func == kScaleMax|| func == kScaleMin ) {
    if ( BBlocksSMCalc[0][1]->ScaleNodeScale2[0].size() <= 3){
      printf("FastNLOReader::SetFunctionalForm. ERROR! There is no second scale variable available in this table.\n");
      printf("      Using FastNLOReader::kScale1 only.\n");
      SetFunctionalForm(kScale1,kMuX);
    }
    for(int i=0;i<NObsBin;i++){
      if ( BBlocksSMCalc[0][1]->ScaleNodeScale2[i].size() < 4 ){
	printf("FastNLOReader::SetFunctionalForm. WARNING! Scale2 has only very little nodes (n=%zd) in bin %d.\n",BBlocksSMCalc[0][0]->ScaleNodeScale2[i].size(),i);
      }
    }
  }
}


//______________________________________________________________________________


void FastNLOReader::SetMuRFunctionalForm( EScaleFunctionalForm func , bool Verbose ){
  SetFunctionalForm(func,kMuR,Verbose);
  FillAlphasCache();
}


//______________________________________________________________________________


void FastNLOReader::SetMuFFunctionalForm( EScaleFunctionalForm func , bool ReFillCache , bool Verbose ){
  SetFunctionalForm(func,kMuF,Verbose);
  if ( ReFillCache ){
    FillPDFCache();
  }
}



//______________________________________________________________________________



bool FastNLOReader::SetScaleFactorsMuRMuF( double xmur, double xmuf, bool ReFillCache, bool Verbose ){
  // 
  // Set renormalization and factorization scale factors simultaneously for scale variations in all v2 tables.
  // You have to ReFill your cache!
  // This is done automatically, but if you want to do it by yourself set ReFillCache = false.
  // 
  // The function aborts the whole program if non-sensical scale factors < 1.E-6 are requested. 
  // The function returns true if the requested scale factors can be used with the available table.
  // If threshold corrections are selected, only xmur / xmuf = 1. is allowed.
  //    If this is not the case, xmur and xmuf are unchanged,
  //    a warning is printed and the function returns false!
  // If it is NOT a flexibleScaleTable and there is no scalevar-table for xmuf,
  //    xmur and xmuf are unchanged, a warning is printed and the function returns false!

  // Check whether xmur and xmuf are positive and at least larger than 1.E-6
  if ( xmur < 1.E-6 || xmuf < 1.E-6 ) {
    printf("FastNLOReader::SetScaleFactorsMuRMuF: ERROR! Selected scale factors too small ( < 1.E-6 )!\n");
    printf("  Program aborted!\n");
    exit(1);
  }
  
  // Check whether threshold corrections exist and are activated
  if ( !BBlocksSMCalc[kThresholdCorrection].empty() ){
    bool lkth = false;
    for ( unsigned int i = 0 ; i <BBlocksSMCalc[kThresholdCorrection].size() ; i++ ){
      if ( bUseSMCalc[kThresholdCorrection][i] ) {
	lkth = true;
	break;
      }
    }
    if ( lkth && abs(xmur-xmuf) > DBL_MIN ) {
      printf("FastNLOReader::SetScaleFactorsMuRMuF: WARNING! Threshold corrections do not allow\n");
      printf("  different scale factors for MuR and MuF, nothing changed!\n");
      printf("  Please do only symmetric scale variation, i.e. xmur = xmuf,\n");
      printf("  or deactivate first all threshold corrections using\n");
      printf("  FastNLOReader::SetContributionON(kTresholdCorrections,Id,false).\n");
      return false;
    }
  }

  // Deal with factorization scale first
  // Check whether corresponding xmuf variation exists in case of v2.0 table 
  if ( !GetIsFlexibleScaleTable() ) {
    const double xmuf0 = B_NLO()->ScaleFac[0][fScalevar];
    const int ns = GetNScaleVariations();
    if ( Verbose ) {
      cout << "FastNLOReader::SetScaleFactorsMuRMuF: NScaleVarMax = " << ns << ", B->ScaleFac[0].size() = " << B_NLO()->ScaleFac[0].size() << endl;
    }
    int sf = -1;
    for ( int is = 0 ; is<ns ; is++ ) {
      if ( abs(B_NLO()->ScaleFac[0][is]-xmuf) < DBL_MIN ){
	sf = is;
	break;
      }
    }
    if ( sf == -1 ) {
      printf("FastNLOReader::SetScaleFactorsMuRMuF: WARNING! Could not find table with given mu_f scale factor of %6.3f.\n",xmuf);
      printf("  Nothing changed!\n");
      return false;
    } else {
      SetScaleVariation( sf , ReFillCache , Verbose );
    }
  } else {
    fScaleFacMuF = xmuf;
    if ( GetIsFlexibleScaleTable() && Verbose ) {SetFunctionalForm( fMuFFunc , kMuF );} // just for printout
    if ( ReFillCache ){
      FillPDFCache();
    }
  }
  if ( Verbose ) {printf("FastNLOReader::SetScaleFactorsMuRMuF: Setting multiplicative scale factor for factorization scale to %4.2f times the nominal scale.\n",xmuf);}
  
  // Now the renormalization scale
  fScaleFacMuR = xmur;
  if ( ReFillCache ){
    FillAlphasCache();
  }
  if ( GetIsFlexibleScaleTable() && Verbose ) {SetFunctionalForm( fMuRFunc , kMuR );} // just for printout
  if ( Verbose ) {printf("FastNLOReader::SetScaleFactorsMuRMuF: Setting multiplicative scale factor for renormalization scale to %4.2f times the nominal scale.\n",xmur);}
  
  return true;
}



//______________________________________________________________________________



void FastNLOReader::ReadTable(void)
{
  //
  // Read in the FastNLO Table
  //
  
  // Check whether file exists
  FILE* fp = fopen(ffilename.c_str(), "r");
  if (fp) {
    fclose(fp);
  } else {
    cout << "FastNLOReader::ReadTable: ERROR! Table file " << ffilename.c_str() << endl;
    cout << "                          not found, exiting!" << endl;
    exit(1);
  } 
  
  // Open stream
  ifstream* instream = new ifstream(ffilename.c_str(),ios::in);
  
  // Read block A1 and A2
  ReadBlockA1(instream);
  ReadBlockA2(instream);

  // Initialize lists for BlockB's
  BBlocksSMCalc.resize(10);
  BBlocksNewPhys.resize(10);
  
  bUseSMCalc.resize(BBlocksSMCalc.size());
  bUseNewPhys.resize(BBlocksNewPhys.size());

  // Initialize BlockB's
  FastNLOBlockB* BlockB_LO   = NULL;
  FastNLOBlockB* BlockB_NLO  = NULL;
  FastNLOBlockB* BlockB_THC1 = NULL;
  FastNLOBlockB* BlockB_THC2 = NULL;
  FastNLOBlockB* BlockB_NPC1 = NULL;
  int nblocks = Ncontrib;
  
  for(int i=0;i<nblocks;i++){
    // Read block
    FastNLOBlockB* blockb = new FastNLOBlockB( "ReadingBlockB", NObsBin , instream );
    blockb->SetIc(i+1);
    char nbuf[400];
    if ( blockb->IDataFlag && !BlockB_Data ) { // Data
      blockb->SetName("Data");
      BlockB_Data = blockb;
    }
    else if ( blockb->IDataFlag && BlockB_Data ) { // Data, but data already initalized
      printf("FastNLOReader::ReadTable(): WARNING! Only one data block is allowed, skipped!\n");
    }
    else if ( blockb->IRef == 1 ) { // Reference table, implemented only for LO or NLO
      if ( blockb->IContrFlag1==1 && blockb->IContrFlag2==1 ){
	if ( blockb->NScaleDep != 3 )   blockb->SetName("BlockB. LO Reference. v2.0.");
	if ( blockb->NScaleDep == 3 )   blockb->SetName("BlockB. LO Reference. v2.1.");
	BlockB_LO_Ref		= blockb;
      }
      else if ( blockb->IContrFlag1==1 && blockb->IContrFlag2==2 ){
	if ( blockb->NScaleDep != 3 )   blockb->SetName("BlockB. NLO Reference. v2.0.");
	if ( blockb->NScaleDep == 3 )   blockb->SetName("BlockB. NLO Reference. v2.1.");
	BlockB_NLO_Ref	= blockb;
      }
      else {
	printf("FastNLOReader::ReadTable(): ERROR! Reference tables are only implemented for fixed order, stopped!\n");
	exit(1);
      }
    }
    else if ( blockb->IRef == 0 && !blockb->IAddMultFlag ) { // Additive corrections
      if ( blockb->IContrFlag1==1 ) { // Fixed order
	if ( blockb->IContrFlag2==1 ){ // LO
	  sprintf(nbuf,"BlockB. %s %s %s",fContrName[blockb->IContrFlag1].c_str(),fOrdName[blockb->IContrFlag1-1][blockb->Npow-ILOord].c_str(),fNSDep[blockb->NScaleDep].c_str());
	  blockb->SetName(nbuf);
	  BlockB_LO  = blockb;
	}
	else if ( blockb->IContrFlag2==2 ){ //NLO
	  sprintf(nbuf,"BlockB. %s %s %s",fContrName[blockb->IContrFlag1].c_str(),fOrdName[blockb->IContrFlag1-1][blockb->Npow-ILOord].c_str(),fNSDep[blockb->NScaleDep].c_str());
	  blockb->SetName(nbuf);
	  BlockB_NLO = blockb;
	}
      }
      else if ( blockb->IContrFlag1==2 ){ // Threshold corrections
	if ( blockb->IContrFlag2==1 ){ // 1-loop
	  sprintf(nbuf,"BlockB. %s %s %s",fContrName[blockb->IContrFlag1].c_str(),fOrdName[blockb->IContrFlag1-1][blockb->Npow-ILOord].c_str(),fNSDep[blockb->NScaleDep].c_str());
	  blockb->SetName(nbuf);
	  BlockB_THC1 = blockb;
	}
	else if ( blockb->IContrFlag2==2 ){ // 2-loop
	  sprintf(nbuf,"BlockB. %s %s %s",fContrName[blockb->IContrFlag1].c_str(),fOrdName[blockb->IContrFlag1-1][blockb->Npow-ILOord].c_str(),fNSDep[blockb->NScaleDep].c_str());
	  blockb->SetName(nbuf);
	  BlockB_THC2 = blockb;
	}
	else {
	  printf("FastNLOReader::ReadTable(): ERROR! Threshold correction implemented only up to 2-loops, stopped!\n");
	  exit(1);
	}
      }
      //	else if ( blockb->IContrFlag1>=3 ){ // 
      // sprintf(nbuf,"BlockB. %s. %s. %s",fNPName[blockb->IContrFlag2].c_str(),fOrdName[blockb->Npow-ILOord].c_str(),fNSDep[blockb->NScaleDep].c_str());
      // blockb->SetName(nbuf);
      // BBlocksNewPhys[blockb->IContrFlag2].push_back(blockb);
      // bUseNewPhys[blockb->IContrFlag2].push_back(true);
      // }
      else {
	printf("FastNLOReader::ReadTable(): ERROR! Further additive corrections not yet implemented, stopped!\n");
	exit(1);
      }
    }
    else if ( blockb->IRef == 0 && blockb->IAddMultFlag ) { // Multiplicative corrections
      if ( blockb->IContrFlag1==4 ) { // Non-perturbative corrections
	sprintf(nbuf,"BlockB. %s %s %s",fContrName[blockb->IContrFlag1].c_str(),fOrdName[blockb->IContrFlag1-1][blockb->IContrFlag2-1].c_str(),fNSDep[blockb->NScaleDep].c_str());
	blockb->SetName(nbuf);
	BlockB_NPC1 = blockb;
      }
      
    }
    else {
      printf("FastNLOReader::ReadTable(): ERROR! Further multiplicative corrections not yet implemented, stopped!\n");
      exit(1);
    }
  }
  
  // Assign NPC, switch off by default
  if ( BlockB_NPC1 ) {
    BBlocksSMCalc[BlockB_NPC1->IContrFlag1-1].push_back(BlockB_NPC1);
    bUseSMCalc[BlockB_NPC1->IContrFlag1-1].push_back(false);
  }

  // Assign THC, switch off by default
  if ( BlockB_THC2 ) {
    BBlocksSMCalc[BlockB_THC2->IContrFlag1-1].push_back(BlockB_THC2);
    bUseSMCalc[BlockB_THC2->IContrFlag1-1].push_back(false);
  }
  if ( BlockB_THC1 ) {
    BBlocksSMCalc[BlockB_THC1->IContrFlag1-1].push_back(BlockB_THC1);
    bUseSMCalc[BlockB_THC1->IContrFlag1-1].push_back(false);
  }
  
  // Assign fixed order calculations (LO must be [0]), switch on by default
  if ( BlockB_LO )  {
    BBlocksSMCalc[0].push_back(BlockB_LO);
    bUseSMCalc[0].push_back(true);
  } else {
    printf("FastNLOReader::ReadTable(): ERROR! Could not find any LO Calculation, stopped!\n");
    exit(1);
  }
  if ( BlockB_NLO ) {
    BBlocksSMCalc[0].push_back(BlockB_NLO);
    bUseSMCalc[0].push_back(true);
  }
  
  // Some printout?
  // PrintTableInfo();
  
}

//______________________________________________________________________________



void FastNLOReader::PrintTableInfo(const int iprint) const {

  //---  Initialization for nice printing
  string CSEPS = "##################################################################################\n";
  string LSEPS = "#---------------------------------------------------------------------------------\n";
  printf("\n");
  printf(" %s",CSEPS.c_str());
  printf(" # Overview on contribution types and numbers contained in table:\n");
  printf(" %s",LSEPS.c_str());
  printf(" # Number of contributions: %2i\n",Ncontrib);
  
  unsigned int icnt = 0;
  for ( unsigned int j = 0 ; j<BBlocksSMCalc.size() ; j++ ){
    if ( !BBlocksSMCalc[j].empty() ){
      for ( unsigned int i = 0 ; i<BBlocksSMCalc[j].size() ; i++ ){
	icnt++;
	cout << " # "<< "  No.: " << icnt << ", type: " << fContrName[j] <<", Id: " << i <<
	  ", order: " << BBlocksSMCalc[j][i]->CtrbDescript[0] <<
	  ", by: " << BBlocksSMCalc[j][i]->CodeDescript[0] << endl;
      }
      if ( iprint > 0 ){
	for ( unsigned int k = 0 ; k<BBlocksSMCalc[j][0]->CodeDescript.size();k++ ) {
	  printf( " # \t\t%s\n",BBlocksSMCalc[j][0]->CodeDescript[k].c_str());
	}
	//BBlocksSMCalc[j][0]->Print(0,0);
      }
    }
  }
  
  for ( unsigned int j = 0 ; j<BBlocksNewPhys.size() ; j++ ){
    if ( !BBlocksNewPhys[j].empty() ){
      for ( unsigned int i = 0 ; i<BBlocksNewPhys[j].size() ; i++ ){
	icnt++;
	cout << " # "<< "  No.: " << icnt << ", type: " << fContrName[j] <<", Id: " << i <<
	  ", order: " << BBlocksNewPhys[j][i]->CtrbDescript[0] <<
	  ", by: " << BBlocksNewPhys[j][i]->CodeDescript[0] << endl;
      }
      printf(" #   -> SM extensions can not be evaluated by this reader! Just skipping those ...\n");
    }
  }

  if ( BlockB_Data ) {
    icnt++;
    cout << " # "<< "  No.: " << icnt << ", type: Data, Id: 0, by: "
	 << BlockB_Data->CodeDescript[0] << endl;
    if ( iprint > 0 ){
      for ( unsigned int k = 0 ; k<BlockB_Data->CodeDescript.size();k++ ) {
	printf( " * \t\t%s\n",BlockB_Data->CodeDescript[k].c_str());
      }
    }
  }
  printf(" %s",CSEPS.c_str());

}

//______________________________________________________________________________


void FastNLOReader::SetContributionON( ESMCalculation eCalc , unsigned int Id , bool SetOn, bool Verbose ){
  if ( bUseSMCalc[eCalc].empty() || BBlocksSMCalc.empty() ){
    printf("FastNLOReader::SetContributionON. Warning. This contribution (%s) does not exist in this table. Cannot switch it On/Off. Ignoring call.\n",fContrName[eCalc].c_str());
    return;
  }
   
  if ( bUseSMCalc[eCalc].size() < Id || BBlocksSMCalc[eCalc].size() < Id || !BBlocksSMCalc[eCalc][Id] ){
    printf("FastNLOReader::SetContributionON. Warning. This Id = %d does not exist for this contribtion. Cannot switch it On/Off. Ignoring call.\n",Id);
    return;
  }
   
  if (Verbose) {printf("fastNLOReader::SetContributionON. %s contribution '%s' with Id = %d.\n",
		       (SetOn?"Activating":"Deactivating"),fContrName[eCalc].c_str(),Id);}
  
  bUseSMCalc[eCalc][Id] = SetOn;   
  
  // if ( eCalc==kThresholdCorrection && SetOn && fScaleFacMuR!=BBlocksSMCalc[0][1]->ScaleFac[0][fScalevar]){
  //   printf("FastNLOReader::SetContributionON. Info. Resetting a-posteriori scale variation factor, since threshold corrections can not be used with an a-posteriori scale variation..\n");
  //   SetScaleFactorMuR(1.,true);
  // }
  
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
  //--- Ncontrib: All contributions including additive, multiplicative or data
  //---           (Then: Nadd = Ncontrib - Nmult - Ndata)
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
    StripWhitespace(&ScDescript[i]);
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
    StripWhitespace(&DimLabel[i]);
  }

  IDiffBin.resize(NDim);
  for(int i=0;i<NDim;i++){
    *table >>  IDiffBin[i];
  }
  LoBin.resize(NObsBin);
  UpBin.resize(NObsBin);
  // Set rapidity index also when reading a table
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



void FastNLOReader::PrintFastNLOTableConstants(const int iprint) const {

  //
  // Define different levels of detail for printing out table content
  // The minimum (iprint = 0) just gives basic scenario information  
  // including the employed code with references. The additional levels
  // are mostly for debugging purposes.
  // (Partially to be implemented!)
  //
  // iprint = 0: No additional printout
  //          1: Print Block A1 & A2 (A1, A2)
  //          2: Also print basic values of Block B (B0)
  //          3: Also print x nodes of Block B for each contribution (BX)
  //          4: Also print scale nodes of Block B for each contribution (BS)
  //          5: Also print sigma tilde of Block B (not implemented yet)
  
  //---  Initialization for nice printing
  string CSEPS = "##################################################################################\n";
  string LSEPS = "#---------------------------------------------------------------------------------\n";
  
  //
  // Print basic scenario information (always)
  //
  printf("\n");
  printf(" %s",CSEPS.c_str());
  printf(" # Information on fastNLO scenario: %s\n",ScenName.data());
  printf(" %s",LSEPS.c_str());
  printf(" # Description:\n");
  for( unsigned int i=0; i<ScDescript.size(); i++ ){
    printf(" #   %s\n",ScDescript[i].data());
  }
  printf(" #\n");
  printf(" # Centre-of-mass energy Ecms: % -#10.4g GeV\n",Ecms);
  printf(" #\n");
  printf(" # Tot. no. of observable bins: %3i in %1i dimensions:\n",NObsBin,NDim);
  printf(" #\n");
  printf(" # No. of contributions: %1i\n",Ncontrib);
  for ( unsigned int j = 0 ; j<BBlocksSMCalc.size() ; j++ ){
    if ( !BBlocksSMCalc.empty() ) {
      for ( unsigned int i = 0 ; i<BBlocksSMCalc[j].size() ; i++ ){
	if ( BBlocksSMCalc[j][i] ){
	  int iContr = BBlocksSMCalc[j][i]->GetIc();
	  BBlocksSMCalc[j][i]->Print(iContr,iprint);
	}
      }
    }
  }
  if ( BlockB_Data  ){
    int iContr = BlockB_Data->GetIc();
    BlockB_Data->Print(iContr,iprint);
  }

  //
  // Print additional info for debugging including internal variables as selected by iprint
  //
  if ( iprint > 0 ) {
    PrintBlockA1();
    PrintBlockA2();
  }
  if ( iprint > 1 ) {
    for ( unsigned int j = 0; j<BBlocksSMCalc.size(); j++ ){
      if ( !BBlocksSMCalc.empty() ) {
	for ( unsigned int i = 0 ; i<BBlocksSMCalc[j].size() ; i++ ){
	  if ( BBlocksSMCalc[j][i] ){
	    int iContr = BBlocksSMCalc[j][i]->GetIc();
	    BBlocksSMCalc[j][i]->Print(iContr,iprint);
	  }
	}
      }
    }
    if ( BlockB_LO_Ref  ){
      int iContr = BlockB_LO_Ref->GetIc();
      BlockB_LO_Ref->Print(iContr,iprint);
    }
    if ( BlockB_NLO_Ref  ){
      int iContr = BlockB_NLO_Ref->GetIc();
      BlockB_NLO_Ref->Print(iContr,iprint);
    }
  }
  
  printf(" #\n");
  printf(" %s",CSEPS.c_str());
}



//______________________________________________________________________________



void FastNLOReader::PrintBlockA1() const{
  printf("\n *****************************************\n");
  printf(" * fastNLO Table: Block A1\n");
  printf(" *****************************************\n");
  printf("  A1  ISep                              %10i\n",tablemagicno);
  printf("  A1  Itabversion                       %10i\n",Itabversion);
  printf("  A1  ScenName                          %s\n",ScenName.data());
  printf("  A1  Ncontrib                          %10i\n",Ncontrib);
  printf("  A1  Nmult                             %10i\n",Nmult);
  printf("  A1  Ndata                             %10i\n",Ndata);
  printf("  A1  NuserString                       %10i\n",NuserString);
  printf("  A1  NuserInt                          %10i\n",NuserInt);
  printf("  A1  NuserFloat                        %10i\n",NuserFloat);
  printf("  A1  Imachine                          %10i\n",Imachine);
  printf(" #########################################\n");
}



//______________________________________________________________________________



void FastNLOReader::PrintBlockA2() const {
  printf("\n *****************************************\n");
  printf(" * fastNLO Table: Block A2\n");
  printf(" *****************************************\n");
  printf("  A2  ISep                              %10i\n",tablemagicno);
  printf("  A2  IpublUnits                        %10i\n",Ipublunits);
  printf("  A2  NscDescript                       %10zi\n",ScDescript.size());
  for(unsigned int i=0;i<ScDescript.size();i++){
    printf("  A2    ScDescript(%1i)                   %s\n",i+1,ScDescript[i].data());
  }
  printf("  A2  Ecms                              % -#10.4g\n",Ecms);
  printf("  A2  ILOord                            %10i\n",ILOord);
  printf("  A2  NobsBin                           %10i\n",NObsBin);
  printf("  A2  NDim                              %10i\n",NDim);
  for(int i=0;i<NDim;i++){
    printf("  A2    DimLabel(%1i)                     %s\n",i+1,DimLabel[i].data());
  }
  for(int i=0;i<NDim;i++){
    printf("  A2    IDiffBin(%1i)                     %10i\n",i+1,IDiffBin[i]);
  }
  for(int i=0;i<NObsBin;i++){
    for(int j=0;j<NDim;j++){
      printf("  A2      LoBin(%3i,%1i)              % #10.4g\n", i+1,j+1,LoBin[i][j]);
      if(IDiffBin[j]==2)
        printf("  A2      UpBin(%3i,%1i)              % #10.4g\n", i+1,j+1,UpBin[i][j]);
    }
  }
  for(int i=0;i<NObsBin;i++){
    printf("  A2    BinSize(%3i)                    % -#10.4g\n", i+1,BinSize[i]);
  }
  printf("  A2  INormFlag                         %10i\n",INormFlag);
  
  if(INormFlag>1){
    printf("  A2  DenomTable                        %s\n",DenomTable.data());
  }
  if(INormFlag>0){
    for(int i=0;i<NObsBin;i++){
      printf("  A2    IDivLoPointer(%3i)              %10i\n",i+1,IDivLoPointer[i]);
      printf("  A2    IDivUpPointer(%3i)              %10i\n",i+1,IDivUpPointer[i]);
    }
  }
  printf(" #########################################\n");
}


//______________________________________________________________________________


void FastNLOReader::PrintCrossSections( ) const {
  //
  // Print Cross sections in NLO, k-factors and Reference table cross sections
  //
  
   //   if ( XSection.empty() )    CalcCrossSection();
   //   if ( XSectionRef.empty() && XSectionRef_s1.empty() )    CalcReferenceCrossSection();}

  vector < double > xs = XSection;

  printf(" *  \n");
  printf(" *  FastNLO Cross sections for\n");
  for ( unsigned int i = 0 ; i < ScDescript.size() ; i++ ){
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


void FastNLOReader::PrintFastNLODemo(){
   //
   // This method prints out cross sections for different scale
   // variation tables. Though it also changes the currently stored
   // settings of this instance!
   //

  // If flexible-scale table, set MuR and MuF functional forms
  if ( BBlocksSMCalc[0][0]->NScaleDep == 3 ){
    SetMuRFunctionalForm(FastNLOReader::kScale1);
    SetMuFFunctionalForm(FastNLOReader::kScale1);
    //SetMuRFunctionalForm(FastNLOReader::kExpProd2);
    //SetMuRFunctionalForm(FastNLOReader::kExpProd2);
  }

  // Check on existence of LO and NLO (Id = -1 if not existing)
  int ilo   = ContrId(FastNLOReader::kFixedOrder, FastNLOReader::kLeading); 
  int inlo  = ContrId(FastNLOReader::kFixedOrder, FastNLOReader::kNextToLeading);
  if ( ilo < 0 || inlo < 0 ){
     printf("FastNLOReader: ERROR! LO and/or NLO not found, nothing to be done!\n");
     exit(1);
  }
  // Check on existence of 2-loop threshold corrections
  int ithc2 = ContrId(FastNLOReader::kThresholdCorrection, FastNLOReader::kNextToLeading);
  // Switched off by default. Don't do scale variations. Not available for the moment.
  //  if ( ithc2 > -1 ) {
  //    SetContributionON( FastNLOReader::kThresholdCorrection, ithc2, false, false );
  //  }
  
  // Pre-define desired order of scale variations
  const int nxmu = 4;
  double xmu[nxmu] = {1.0, 0.25, 0.5, 2.0};
  int   ixmu[nxmu] = { -1,   -1,  -1,  -1};
  // Get number of available scale variations and check on available scale factors,
  // in particular for MuF; set pointers
  int nscls = GetNScaleVariations();
  // With threshold corrections, allow only default scale (0)
  if ( ithc2 > -1 ) {
    nscls = 1;
  }
  for (int iscls=0; iscls<nscls; iscls++){
    SetScaleVariation(iscls);
    double fxmu = fScaleFacMuF;
    for (int i=0; i<nxmu; i++){
      if (abs(xmu[i]-fxmu) < 0.000001){
	ixmu[i] = iscls;
      }
    }
  }

  // Loop over scales
  for (int iscls=0; iscls<nxmu; iscls++){
    // First result is with NLO, LO result via division by K factor
    if (ixmu[iscls] > -1){
      SetScaleVariation(ixmu[iscls]);
      FillAlphasCache();
      FillPDFCache();
      CalcCrossSection();

      // Second result: Include threshold corrections for NLO if available
      vector < double > kthc;
      if ( ithc2 > -1 ) {
	vector < double > stdk = kFactor;
	SetContributionON( FastNLOReader::kThresholdCorrection, ithc2, true, false );
	CalcCrossSection();
	kthc = kFactor;
	// Threshold K factor is NLO including 2-loop vs. NLO
	for (unsigned int i=0;i<kthc.size();i++){
	   if ( abs(kFactor[i]) > DBL_MIN ){
	     kthc[i] = kFactor[i]/stdk[i];
	  } else {
	    kthc[i] = -1.;
	  }
	}
	SetContributionON( FastNLOReader::kThresholdCorrection, ithc2, false, false );
      }
      
      PrintCrossSectionsDefault( kthc );
    }
  }
}


//______________________________________________________________________________


void FastNLOReader::PrintCrossSectionsDefault( const vector <double> kthc ) const {
  //
  // Print observable binnning and cross sections at
  // LO, NLO and K factors like in Fortran Reader for comparison
  //

  // Some initialization
  const string CSEP41("#########################################");
  const string DSEP41("=========================================");
  const string SSEP41("-----------------------------------------");
  const string CSEP = CSEP41 + CSEP41 + CSEP41 + CSEP41;
  const string DSEP = DSEP41 + DSEP41 + DSEP41 + DSEP41;
  const string SSEP = SSEP41 + SSEP41 + SSEP41 + SSEP41;
  

  // Check on existence of 2-loop threshold corrections
  //const int ithc2 = kthc.empty() ? -1 : ContrId( FastNLOReader::kThresholdCorrection, FastNLOReader::kNextToLeading);
  const int ithc2 = kthc.empty() ? -1 : ContrId( kThresholdCorrection,kNextToLeading);


  cout << DSEP << endl;
  printf(" Cross Sections\n");
  printf(" The scale chosen here are: mu_f = % #6.3f * %s, and mu_r = % #6.3f * %s \n",fScaleFacMuF,GetScaleDescription().c_str(),fScaleFacMuR,GetScaleDescription().c_str());
  cout << SSEP << endl;
    
  if ( NDim == 2 ){

     // non-perturbative corrections (just first np correction)
     const int inpc1 = ContrId(FastNLOReader::kNonPerturbativeCorrection, FastNLOReader::kLeading);
     const vector < double > knpc = inpc1>-1 ? BBlocksSMCalc[3][0]->fact : vector<double>(NObsBin);
     

     string header0 = "  IObs  Bin Size IODim1 "; 
     string header1 = "   IODim2 ";
     string header2 = " LO cross section   NLO cross section   K NLO";
     if ( ithc2>-1 )header2 += "     K THC";
     if ( inpc1>-1 )header2 += "     K NPC";
     unsigned int NDimBins[NDim];
     printf("%s [ %-12s ] %s [  %-12s  ] %s\n",
	    header0.c_str(),DimLabel[0].c_str(),header1.c_str(),DimLabel[1].c_str(),header2.c_str());
     cout << SSEP << endl;
     for ( int i=0; i<NObsBin; i++ ){ 
	for ( int j=0; j<NDim; j++ ){ 
	   if ( i==0 )					NDimBins[j] = 1;
	   else if ( LoBin[i-1][j] < LoBin[i][j])	NDimBins[j]++;
	   else if ( LoBin[i][j] < LoBin[i-1][j])	NDimBins[j] = 1;
	}
	if ( ithc2<0 && inpc1<0 ) {
	   printf(" %5.i % -#10.4g %5.i % -#10.4g % -#10.4g %5.i  %-#8.2E  %-#8.2E %#18.11E %#18.11E %#9.5F",
		  i+1,BinSize[i],NDimBins[0],LoBin[i][0],UpBin[i][0],
		  NDimBins[1],LoBin[i][1],UpBin[i][1],XSection_LO[i],XSection[i],kFactor[i]);
	} else if ( inpc1<0 ) {
	   printf(" %5.i % -#10.4g %5.i % -#10.4g % -#10.4g %5.i  %-#8.2E  %-#8.2E %#18.11E %#18.11E %#9.5F %#9.5F",
		  i+1,BinSize[i],NDimBins[0],LoBin[i][0],UpBin[i][0],
		  NDimBins[1],LoBin[i][1],UpBin[i][1],XSection_LO[i],XSection[i],kFactor[i],kthc[i]);
	} else if ( ithc2<0 ) {
	   printf(" %5.i % -#10.4g %5.i % -#10.4g % -#10.4g %5.i  %-#8.2E  %-#8.2E %#18.11E %#18.11E %#9.5F %#9.5F",
		  i+1,BinSize[i],NDimBins[0],LoBin[i][0],UpBin[i][0],
		  NDimBins[1],LoBin[i][1],UpBin[i][1],XSection_LO[i],XSection[i],kFactor[i],knpc[i]);
	} else {
	   printf(" %5.i % -#10.4g %5.i % -#10.4g % -#10.4g %5.i  %-#8.2E  %-#8.2E %#18.11E %#18.11E %#9.5F %#9.5F %#9.5F",
		  i+1,BinSize[i],NDimBins[0],LoBin[i][0],UpBin[i][0],
		  NDimBins[1],LoBin[i][1],UpBin[i][1],XSection_LO[i],XSection[i],kFactor[i],kthc[i],knpc[i]);
	}
	printf("\n");
     }
  } else {
     printf("FastNLOReader: WARNING! Print out optimized for two dimensions. No output for %1.i dimensions.\n",NDim);
  }

}


//______________________________________________________________________________


void FastNLOReader::PrintCrossSectionsData() const{
  //
  // Print summary information on data table if available
  //
  
  if ( !BlockB_Data ){
    return;
  }
  
  vector < double > xs = BlockB_Data->Value;
  
  string CSEP41("#########################################");
  string DSEP41("=========================================");
  string SSEP41("-----------------------------------------");
  string CSEP = CSEP41 + CSEP41 + CSEP41 + CSEP41;
  string DSEP = DSEP41 + DSEP41 + DSEP41 + DSEP41;
  string SSEP = SSEP41 + SSEP41 + SSEP41 + SSEP41;
  
  cout << DSEP << endl;
  printf(" Measurement\n");
  //--- Additionally print out data description
  // for ( unsigned int k = 0 ; k<BlockB_Data->CodeDescript.size();k++ ) {
  //   printf( "\t\t%s\n",BlockB_Data->CodeDescript[k].c_str());
  // }
  cout << SSEP << endl;
  
  if ( NDim == 2 ){
    string header[3] = { "  IObs  Bin Size IODim1 ", 
			 "   IODim2 ",
			 "   X Section QSumUnc.+ QSumUnc.- QSumCor.+ QSumCor.-\n"};
    string label[2] = { "[ " + DimLabel[0] + "     ]", "[  " + DimLabel[1] + "           ]"};
    unsigned int NDimBins[NDim];
    printf("%s %s %s %s %s",header[0].c_str(),label[0].c_str(),header[1].c_str(),label[1].c_str(),header[2].c_str());
    //--- Additionally print out descriptions on all uncertainties
    // for ( int iUnco = 0 ; iUnco<BlockB_Data->Nuncorrel ; iUnco++ ){
    //   printf("%30s   ",BlockB_Data->UncDescr[iUnco].c_str());
    // }
    // for ( int iCorr = 0 ; iCorr<BlockB_Data->Ncorrel ; iCorr++ ){
    //   printf("%30s   ",BlockB_Data->CorDescr[iCorr].c_str());
    // }
    // cout<<endl;

    cout << SSEP << endl;
    for ( unsigned int i=0; i<xs.size(); i++ ){ 
      for ( int j=0; j<NDim; j++ ){ 
	if ( i==0 ){
	  NDimBins[j] = 1;
	} else if ( LoBin[i-1][j] < LoBin[i][j]){
	  NDimBins[j]++;
	} else if ( LoBin[i][j] < LoBin[i-1][j]){
	  NDimBins[j] = 1;
	}
      }
      printf(" %5.i % -#10.4g %5.i % -#10.4g % -#10.4g %5.i  %-#8.2E  %-#8.2E   % -#10.3E",
 	     i+1,BinSize[i],NDimBins[0],LoBin[i][0],UpBin[i][0],
 	     NDimBins[1],LoBin[i][1],UpBin[i][1],xs[i]);
      
      //--- For now add sources quadratically
      double DXsUncorLo = 0.;
      double DXsUncorHi = 0.;
      for ( int iUnco = 0 ; iUnco<BlockB_Data->Nuncorrel ; iUnco++ ){
	DXsUncorLo += BlockB_Data->UncorLo[i][iUnco] * BlockB_Data->UncorLo[i][iUnco];
        DXsUncorHi += BlockB_Data->UncorHi[i][iUnco] * BlockB_Data->UncorHi[i][iUnco];
      }
      DXsUncorLo = -sqrt(DXsUncorLo);
      DXsUncorHi = +sqrt(DXsUncorHi);
      printf(" %+8.2E %+8.2E",DXsUncorHi,DXsUncorLo);
      
      double DXsCorrLo = 0.;
      double DXsCorrHi = 0.;
      for ( int iCorr = 0 ; iCorr<BlockB_Data->Ncorrel ; iCorr++ ){
	DXsCorrLo += BlockB_Data->CorrLo[i][iCorr] * BlockB_Data->CorrLo[i][iCorr];
        DXsCorrHi += BlockB_Data->CorrHi[i][iCorr] * BlockB_Data->CorrHi[i][iCorr];
      }
      DXsCorrLo = -sqrt(DXsCorrLo);
      DXsCorrHi = +sqrt(DXsCorrHi);
      printf(" %+8.2E %+8.2E",DXsCorrHi,DXsCorrLo);
      cout << endl;
    }
  } else {
    printf("WARNING! Sorry, currently no output for more than 2 dimensions!\n");
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
  
  if ( BBlocksSMCalc[0][0]->NScaleDep == 3 ){
    if ( fMuFFunc == kScale1 && fMuRFunc == kScale1 )	{
      printf(" *  FastNLOReader::PrintCrossSectionsWithReference. Info. Taking reference cross sections 's1'\n");
      xsref = XSectionRef_s1;
    }
    else if ( fMuFFunc == kScale2 && fMuRFunc == kScale2 ){
      printf(" *  FastNLOReader::PrintCrossSectionsWithReference. Info. Taking reference cross sections 's2'\n");
      xsref = XSectionRef_s2;
    }
    else if ( fMuFFunc == kQuadraticMean && fMuRFunc == kQuadraticMean ) {
      printf(" *  FastNLOReader::PrintCrossSectionsWithReference. Info. Taking reference cross sections 'mixed'\n");
      xsref = XSectionRefMixed;
    }
    else {
      xsref = XSectionRefMixed;
      printf(" *  FastNLOReader::PrintCrossSectionsWithReference. Info. Taking reference cross sections 'mixed'\n");
    }
  }
  else xsref = XSectionRef;


  printf(" *  \n");
  printf(" *  FastNLO Cross sections for\n");
  for ( unsigned int i = 0 ; i < ScDescript.size() ; i++ ){
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
  printf(" *  If you have choosen mu_r=mu_f=%s, or mu_r=mu_f=%s or mu_r=mu_f=sqrt((%s^2+%s^2)/2), then you access automatically the corresponding reference cross section.\n",
	 BBlocksSMCalc[0][1]->ScaleDescript[0][0].c_str(),BBlocksSMCalc[0][1]->ScaleDescript[0][1].c_str(),BBlocksSMCalc[0][1]->ScaleDescript[0][0].c_str(),BBlocksSMCalc[0][1]->ScaleDescript[0][1].c_str());
  printf(" *  In any other case your reference cross section is calculated using mu_r=mu_f=sqrt((%s^2+%s^2)/2).\n",BBlocksSMCalc[0][1]->ScaleDescript[0][0].c_str(),BBlocksSMCalc[0][1]->ScaleDescript[0][1].c_str());
  printf(" *  To be fully consistent with the nlojet++ reference cross section, you also have to adjust alpha_s and the alpha_s evolution accordingly.\n\n");

  printf("\n");
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


int FastNLOReader::GetNScaleVariations() const {
   if ( GetIsFlexibleScaleTable() ){
      printf("FastNLOReader::GetNScaleVariations(). This is a 'flexible-scale' table, therefore you can choose all desired scale variations.\n");
      return 0;
   }
  
  // Check for maximal scale variation of all rel. and active SM calcs
  // Assume a maximum of 10!
  unsigned int scalevarmax = 10;
  for ( unsigned int j = 0 ; j<BBlocksSMCalc.size() ; j++ ){
    if ( !BBlocksSMCalc.empty() ){
      for ( unsigned int i = 0 ; i<BBlocksSMCalc[j].size() ; i++ ){
	// Do not check pQCD LO or mult. corrections
	if ( bUseSMCalc[j][i] && !BBlocksSMCalc[j][i]->IAddMultFlag &&
	     !( j==kFixedOrder && i==kLeading) ){
	  if ( BBlocksSMCalc[j][i]->Nscalevar[0] < scalevarmax ) {
	    scalevarmax = BBlocksSMCalc[j][i]->Nscalevar[0];
	  }
	}
      }
    }
  }
  return scalevarmax;
}


//______________________________________________________________________________


vector < double > FastNLOReader::GetScaleFactors() const {
  if ( BBlocksSMCalc[0][1]->NScaleDep ==3 ){
    printf("FastNLOReader::GetScaleFactors(). This is a 'flexible scale table', therefore you can choose all desired scale variations.\n");
    return vector<double>();
  }
  return BBlocksSMCalc[0][1]->ScaleFac[0];
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


vector < double > FastNLOReader::GetKFactors( ){
  // Get ratio of fast calculated NLO to LO cross section
  if ( XSection.empty() ){
    CalcCrossSection();
  }
  return kFactor;
}


//______________________________________________________________________________


vector < double > FastNLOReader::GetReferenceCrossSection( ){
  // Get reference cross section from direct nlojet++ calculation
  
  if ( XSectionRef.empty() && XSectionRef_s1.empty() ){
    CalcReferenceCrossSection();
  }
  
  if ( BBlocksSMCalc[0][0]->NScaleDep == 3 ){
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
  
  if ( BBlocksSMCalc[0][0]->NScaleDep == 3 ){
    for(int i=0;i<NObsBin;i++){
      double unit = fUnits==kAbsoluteUnits ? BinSize[i] : 1.;
      for(int n=0;n<BBlocksSMCalc[0][1]->NSubproc;n++) {
	XSectionRefMixed[i]		+= BBlocksSMCalc[0][0] ->SigmaRefMixed[i][n] * unit;
	XSectionRef_s1[i]		+= BBlocksSMCalc[0][0] ->SigmaRef_s1[i][n] * unit;
	XSectionRef_s2[i]		+= BBlocksSMCalc[0][0] ->SigmaRef_s2[i][n] * unit;
      }
      for(int n=0;n<BBlocksSMCalc[0][1]->NSubproc;n++) {
	XSectionRefMixed[i]		+= BBlocksSMCalc[0][1]->SigmaRefMixed[i][n] * unit;
	XSectionRef_s1[i]		+= BBlocksSMCalc[0][1]->SigmaRef_s1[i][n] * unit;
	XSectionRef_s2[i]		+= BBlocksSMCalc[0][1]->SigmaRef_s2[i][n] * unit;
      }
    }
  }
  
  if ( BBlocksSMCalc[0][0]->NScaleDep != 3 && ( BlockB_NLO_Ref==NULL ) )
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

  // perturbative (additive) contributions
  for ( unsigned int j = 0 ; j<BBlocksSMCalc.size() ; j++ ){
    if ( !BBlocksSMCalc[j].empty() ) {
      for ( unsigned int i = 0 ; i<BBlocksSMCalc[j].size() ; i++ ){
	if ( bUseSMCalc[j][i] && !BBlocksSMCalc[j][i]->IAddMultFlag) {
	  if ( BBlocksSMCalc[0][0]->NScaleDep != 3 ){ // v2.0
	    CalcCrossSectionv20(BBlocksSMCalc[j][i]);
	  } else if ( BBlocksSMCalc[0][0]->NScaleDep == 3 ){
	    // ---- v2.1---- //
	    CalcCrossSectionv21(BBlocksSMCalc[j][i]);
	  }
	}
      }
    }
  }


  // contributions from the a-posteriori scale variation
  if ( BBlocksSMCalc[0][0]->NScaleDep!=3 ){
    if ( abs(fScaleFacMuR-BBlocksSMCalc[0][1]->ScaleFac[0][fScalevar]) > DBL_MIN ) {
      CalcAposterioriScaleVariation();
    }
  }
  
  // calculate LO cross sections
  if ( !BBlocksSMCalc[0][0] ) printf("CalcCrossSection(). Warning. There is no LO fixed order calculation.\n");
  else if (BBlocksSMCalc[0][0]->Npow != ILOord)   printf("CalcCrossSection(). Warning. The table, which is supposed to be LO fixed order is of order %d compared to LO order.\n",BBlocksSMCalc[0][0]->Npow);
  else {
    if ( BBlocksSMCalc[0][0]->NScaleDep != 3 )	{ // v2.0
      CalcCrossSectionv20(BBlocksSMCalc[0][0],true);
    }
    else {
      CalcCrossSectionv21(BBlocksSMCalc[0][0],true);
    }
  }
  
  // non-perturbative corrections (multiplicative corrections)
  for ( unsigned int j = 0 ; j<BBlocksSMCalc.size() ; j++ ){
    if ( !BBlocksSMCalc[j].empty() ) {
      for ( unsigned int i = 0 ; i<BBlocksSMCalc[j].size() ; i++ ){
	if ( bUseSMCalc[j][i] && BBlocksSMCalc[j][i]->IAddMultFlag && BBlocksSMCalc[j][i]->IContrFlag1 == 4 ) {
	  for(int iB=0;iB<NObsBin;iB++){
	    if (    BBlocksSMCalc[j][i]->IContrFlag2 == 1 ) {
	      XSection[iB] *= BBlocksSMCalc[j][i]->fact[iB];
	    }
	    //	    XSection_LO[iB]	*= BBlocksSMCalc[j][i]->fact[iB];
	  }
	}
      }
    }
  }
  
  // ---- k-factor calculation ---- //
  for(int i=0;i<NObsBin;i++){
    kFactor[i]	= XSection[i] / XSection_LO[i];
  }
  
}

//______________________________________________________________________________


void FastNLOReader::CalcAposterioriScaleVariation( ){
  double scalefac	= fScaleFacMuR/fScaleFacMuF;
  vector<double>* XS	= &XSection;
  const double n     = B_LO()->Npow;
  const double L     = std::log(scalefac); 
  const double beta0 = (11.*3.-2.*5)/3.;
  for(int i=0;i<NObsBin;i++){
     int nxmax = B_LO()->GetNxmax(i);
     double unit = fUnits==kAbsoluteUnits ? BinSize[i] : 1.;
     for(int j=0;j<B_LO()->GetTotalScalenodes();j++){
	int scalenode1 = j;
	int scalenode2 = j;
	if (BBlocksSMCalc[0][1]->NScaleDim>1){
	   scalenode1 = j / B_LO()->Nscalenode[1];
	   scalenode2 = j % B_LO()->Nscalenode[1];
	}
	double asnp1 = pow(B_LO()->AlphasTwoPi_v20[i][scalenode2],(n+1)/n);//as^n+1
	for(int k=0;k<nxmax;k++){ 
	   for(int l=0;l<B_LO()->NSubproc;l++){ 
	      double clo = B_LO()->SigmaTilde[i][0][scalenode2][k][l] *  B_LO()->PdfLc[i][scalenode2][k][l] * unit;
	      XS->at(i)	+=  asnp1 * clo * n * L * beta0;
	   }
	}
     }
  }
}


//______________________________________________________________________________


void FastNLOReader::CalcCrossSectionv21( FastNLOBlockB* B , bool IsLO){
  //
  //  Cross section calculation for DIS and HHC tables in v2.1 format
  //

  vector<double>* XS = IsLO ? &XSection_LO : &XSection;
  B->fact.resize(NObsBin);
  for(int i=0;i<NObsBin;i++){
    B->fact[i]=0;
    int nxmax = B->GetNxmax(i);
    double unit = fUnits==kAbsoluteUnits ? BinSize[i] : 1.;
    for(unsigned int jS1=0;jS1<B->ScaleNodeScale1[i].size();jS1++){
      for(unsigned int kS2=0;kS2<B->ScaleNodeScale2[i].size();kS2++){
	double Q2   = B->ScaleNodeScale1[i][jS1]*B->ScaleNodeScale1[i][jS1];
	double mur	= CalcMu( kMuR , B->ScaleNodeScale1[i][jS1] ,  B->ScaleNodeScale2[i][kS2] , fScaleFacMuR );
	double muf	= CalcMu( kMuF , B->ScaleNodeScale1[i][jS1] ,  B->ScaleNodeScale2[i][kS2] , fScaleFacMuF );
	double mur2 = pow(mur,2);
	double muf2 = pow(muf,2);
	for(int x=0;x<nxmax;x++){ 
	  for(int n=0;n<B->NSubproc;n++){ 
	    double as	= B->AlphasTwoPi[i][jS1][kS2];
	    double pdflc	= B->PdfLcMuVar[i][x][jS1][kS2][n];
	    if ( pdflc == 0. ) continue;
	    double fac	= as * pdflc * unit;
	    double xsci	=  B->SigmaTildeMuIndep[i][x][jS1][kS2][n] *                  fac;
	    xsci		+= B->SigmaTildeMuFDep [i][x][jS1][kS2][n] * std::log(muf2) * fac;
	    xsci		+= B->SigmaTildeMuRDep [i][x][jS1][kS2][n] * std::log(mur2) * fac;
	    if ( BBlocksSMCalc[0][0]->IPDFdef1 == 2 ) { // DIS tables use log(mu/Q2) instead of log(mu)
	      xsci -= B->SigmaTildeMuFDep [i][x][jS1][kS2][n] * std::log(Q2) * fac;
	      xsci -= B->SigmaTildeMuRDep [i][x][jS1][kS2][n] * std::log(Q2) * fac;
	    }
	    XS->at(i)	+= xsci;
	    B->fact[i]	+= xsci;
	  }
	}
      }
    }
  }
}

//______________________________________________________________________________


void FastNLOReader::CalcCrossSectionv20( FastNLOBlockB* B , bool IsLO ){
  //
  //  Cross section calculation in v2.0 format
  //
   
  int scaleVar		= B->Npow == ILOord ? 0 : fScalevar;
  vector<double>* XS	= IsLO ? &XSection_LO : &XSection;
  B->fact.resize(NObsBin);
  for(int i=0;i<NObsBin;i++){
    B->fact[i] = 0;
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
	  double xsci	= B->SigmaTilde[i][scaleVar][j][k][l] *  B->AlphasTwoPi_v20[i][scalenode2]  *  B->PdfLc[i][scalenode2][k][l] * unit;
	  XS->at(i)	+=  xsci;
	  B->fact[i]	+=  xsci;
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
   fAlphasMz	= AlphasMz;		// new alpha_s value
   FillAlphasCache();
   if ( ReCalcCrossSection ) CalcCrossSection(); 
}

//______________________________________________________________________________


void FastNLOReader::FillAlphasCache(){
  //
  //  Fill the internal alpha_s cache.
  //  This is usally called automatically. Only if you
  //  make use of ReFillCache==false options, you have
  //  to take care of this filling by yourself.
  //
   
  // check if the alpha_s value is somehow reasonable
  double asMz = CalcAlphas(91.18);
  if ( fAlphasEvolution != kFixed && ( asMz > 0.5 || asMz < 0.01 ) ) {
    printf("FastNLOReader. Warning. Your alphas value seems to be unreasonably small/large.\n");
    printf("   Your evolution code calculated alphas(Mz=91.18GeV) = %9.7f.\n",asMz);
  }

  for ( unsigned int j = 0 ; j<BBlocksSMCalc.size() ; j++ ){
    if ( !BBlocksSMCalc.empty() ){
      for ( unsigned int i = 0 ; i<BBlocksSMCalc[j].size() ; i++ ){
	// Check that this contribution type j and no. i should actually be used
	// Otherwise deactivation of e.g. threshold corr. is not respected here
	if ( bUseSMCalc[j][i] && !BBlocksSMCalc[j][i]->IAddMultFlag ){
	  if ( BBlocksSMCalc[j][i]->NScaleDep != 3 ){
	    FillAlphasCacheInBlockBv20( BBlocksSMCalc[j][i]  );
	  }
	  else if ( BBlocksSMCalc[j][i]->NScaleDep == 3 ){
	    FillAlphasCacheInBlockBv21( BBlocksSMCalc[j][i]  );
	  }
	}
      }
    }
  }
}


//______________________________________________________________________________


void FastNLOReader::FillAlphasCacheInBlockBv20( FastNLOBlockB* B ){
  // 
  //  Internal method for filling alpha_s cache
  //
   
  int scaleVar		= B->Npow == ILOord ? 0 : fScalevar;
  // Sanity check that scaleVar is in allowed range
  // For thresh. corr. can otherwise lead to inf and then segfault!
  if ( scaleVar >= GetNScaleVariations() ){
    printf("FastNLOReader::FillAlphasCacheInBlockBv20. ERROR! Trying to refresh\n");
    printf("  cache for non-existing scale variation no. %i while only %i exist in total. Aborted.\n",
	   scaleVar,GetNScaleVariations());
    exit(1);
  }
  double scalefac	= fScaleFacMuR/B->ScaleFac[0][scaleVar];

  for(int i=0;i<NObsBin;i++){
    for(int j=0;j<B->GetTotalScalenodes();j++){
      int scalenode1 = j;
      int scalenode2 = j;
      if (B->NScaleDim>1){          
	scalenode1 = j / B->Nscalenode[1];
	scalenode2 = j % B->Nscalenode[1];
      }
      double mur	= scalefac * B->ScaleNode[i][0][scaleVar][scalenode1];
      double as		= CalcAlphas(mur);
      B->AlphasTwoPi_v20[i][j] = pow( as/TWOPI , B->Npow );
    }
  }
}


//______________________________________________________________________________


void FastNLOReader::FillAlphasCacheInBlockBv21( FastNLOBlockB* B ){
  // 
  //  Internal method for filling alpha_s cache
  //

  for(int i=0;i<NObsBin;i++){
    for(unsigned int jS1=0;jS1<B->ScaleNodeScale1[i].size();jS1++){
      for(unsigned int kS2=0;kS2<B->ScaleNodeScale2[i].size();kS2++){
	double mur		= CalcMu( kMuR , BBlocksSMCalc[0][0]->ScaleNodeScale1[i][jS1] ,  BBlocksSMCalc[0][0]->ScaleNodeScale2[i][kS2] , fScaleFacMuR );
	double as		= CalcAlphas(mur);
	double alphastwopi	= pow( as/TWOPI, B->Npow);
	B->AlphasTwoPi[i][jS1][kS2] = alphastwopi;
      }
    }
  }
}


//______________________________________________________________________________


double FastNLOReader::CalcAlphas( double Q ){
  // 
  //  Internal method for calculating the alpha_s(mu)
  //
  
   switch ( fAlphasEvolution ) {
   case kGRV: return Alphas::CalcAlphasMu( Q , fAlphasMz );
   case kNLOJET: return CalcAlphasNLOJET ( Q , fAlphasMz );
   case kCTEQpdf: return CalcAlphasCTEQpdf        ( Q , fAlphasMz );
   case kExternAs: return EvolveAlphas     ( Q , fAlphasMz );
   case kFixed: return fAlphasMz;
   default: {
      cout << "\nFastNLOReader: ERROR! No alpha_s evolution selected, aborting!\n";
      exit (1);
 
   }
   }
}


//______________________________________________________________________________


double FastNLOReader::CalcAlphasNLOJET(double Q, double alphasMZ){
  //
  // Implementation of Alpha_s evolution as function of Mu_r.
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
  // please cite the nlojet++ references.
  //
  // the original parameters from the cteq6 pdf
  //   // Alpha QCD //
  //   1, 1, 0, -1, 0.1179, 91.70, 1.3, 4.5, 180.0,
  //   double BETA0 =  (11. - 2./3.*NF); // The beta coefficients of the QCD beta function
  //   double BETA1 =  (51. - 19./3.*NF);

  double b0  = 1.2202;
  double b1  = 0.4897;
  // double b2 = 0.1913;

  //double Mz	= 91.187;
  double Mz	= 91.70;
  double L = log(Q/Mz);
  L = (b0 + alphasMZ*b1)*L;

  return alphasMZ/(1.0 + alphasMZ*L);

}


//______________________________________________________________________________


double FastNLOReader::CalcAlphasCTEQpdf(double Q, double alphasMZ){
   //
   // Implementation of Alpha_s evolution as function of Mu_r.
   //
   // alpha_s evolution as it is within by cteq-pdf-1.0.4 and used in nlojet 4.1.3
   // please notice the nlojet++ reference

   double as_twopi = alphasMZ/TWOPI;
   double Mz	= 91.187;
   //double Mz	= 91.70;
   //int ord=2;
   int nf=5;

   const int NF	= 5;
   double b0 =  (11. - 2./3.*NF); // The beta coefficients of the QCD beta function
   double b1 =  (51. - 19./3.*NF);
   double t8 = 1.0/(b0*as_twopi);
   double as0, as1, ot, lt, br = (51.0 - 19.0/3.0*nf)/(b0*b0);
   do {
      lt = log(2.0*t8)/t8;
      ot = t8;

      as0 = (1.0 - br*lt)/(b0*t8);
      as1 = (-1.0 - br*(1.0/t8-2.0*lt))/(b0*t8*t8);
      t8 += (as_twopi - as0)/as1;
   } while(fabs(ot-t8)/ot > 1e-5);
   double lmd = Mz*exp(-t8);
   double t = log(Q/lmd);
   double asMz = 1.0/(b0*t);

   return asMz*(1.0-b1/b0*asMz*log(2.0*t)) *TWOPI;

}



//______________________________________________________________________________


void FastNLOReader::FillPDFCache( bool ReCalcCrossSection ){
  //
  //  Fill the internal pdf cache.
  //  This function has to be called by the user, since the 
  //  pdf parameters and evolutions are calculated externally.
  //
   
  InitPDF();

  // check if the pdf is somehow reasonable
  // vector<double> pdftest = GetXFX(1.e-2,10);
  // if ( pdftest.size() != 13) {
  //   printf("FastNLOReader. Error. The pdf array must have the size of 13 flavors. Exiting.\n");
  //   exit(1);
  // }
  
  // if ( pdftest[6] == 0. ) {
  //   printf("FastNLOReader. Warning. There seems to be no gluon in the pdf.\n");
  // }
  // double sum = 0;
  // for ( int i = 0 ; i<13 ; i++ ) sum+=fabs(pdftest[i]);
  // if ( sum== 0. ) {
  //   printf("FastNLOReader. Error. All 13 pdf probabilities are 0. There might be sth. wrong in your pdf interface. Please check FastNLOUser::GetXFX().\n");
  //   exit(1);
  // }
  // for ( int i = 0 ; i<13 ; i++ ){
  //   if ( pdftest[i] > 1.e30 || ( pdftest[i] < 1.e-10 && pdftest[i] != 0. )) {
  //     printf("FastNLOReader. Warning. The pdf probability of your %d's flavor seeems to be unreasonably large/small (pdf=%8.2e).\n",i,pdftest[i]);
  //   }
  // }
  
  for ( unsigned int j = 0 ; j<BBlocksSMCalc.size() ; j++ ){
    if (  !BBlocksSMCalc.empty() ){
      for ( unsigned int i = 0 ; i<BBlocksSMCalc[j].size() ; i++ ){
	// Check that this contribution type j and no. i should actually be used
	// Otherwise deactivation of e.g. threshold corr. is not respected here
	if ( bUseSMCalc[j][i] && !BBlocksSMCalc[j][i]->IAddMultFlag ){
	  if (BBlocksSMCalc[j][i]->NScaleDim>1){
	    printf("FastNLOReader::FillPDFCache. WOW! NScaleDim>1! This is usually not the case!\n");
	    //scaleindex2 = 1; // If we use multiple scales, then mu_f is by convention the second scale -> index=1 
	    //fScalevar2 = fScalevar % NfScalevar[1]; 
	  }
	       
	  // linear: DIS-case
	  // ---- DIS ---- //
	  if ( BBlocksSMCalc[j][i]->IPDFdef1 == 2 ){
	    if ( BBlocksSMCalc[j][i]->NPDFDim == 0 ) {
	      if	 ( BBlocksSMCalc[j][i]->NScaleDep != 3 )	FillBlockBPDFLCsDISv20(BBlocksSMCalc[j][i]);
	      else if ( BBlocksSMCalc[j][i]->NScaleDep == 3 )	FillBlockBPDFLCsDISv21(BBlocksSMCalc[j][i]);
	    }
	  }
	  // ---- pp ---- //
	  else if (  BBlocksSMCalc[j][i]->IPDFdef1 == 3 ){
	    if ( BBlocksSMCalc[j][i]->NPDFDim == 1 ) {
	      if	 ( BBlocksSMCalc[j][i]->NScaleDep != 3 )	FillBlockBPDFLCsHHCv20(BBlocksSMCalc[j][i]);
	      else						FillBlockBPDFLCsHHCv21(BBlocksSMCalc[j][i]);
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
    }
  }
}     


//______________________________________________________________________________


void FastNLOReader::FillBlockBPDFLCsDISv20( FastNLOBlockB* B ){
  int scaleVar		= B->Npow == ILOord ? 0 : fScalevar;
  double scalefac	= B->ScaleFac[0][scaleVar] == fScaleFacMuF ? 1. : fScaleFacMuF;
  vector<double> xfx(13); // PDFs of all partons
  if ( B->NScaleDep != 3 ){
    for(int i=0;i<NObsBin;i++){
      int nxmax = B->GetNxmax(i);
      for(int j=0;j<B->Nscalenode[0];j++){
	for(int k=0;k<nxmax;k++){ 
	  double xp	= B->XNode1[i][k];
	  double muf	= scalefac * B->ScaleNode[i][0][scaleVar][j];
	  xfx = GetXFX(xp,muf);
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
   
  if ( B->PdfLcMuVar.empty() ) { cout<< "empty."<<endl; exit(1);}

  for(int i=0;i<NObsBin;i++){
    // speed up! if mu_f is only dependent on one variable, we can safe the loop over the other one
    for(int x=0;x<B->GetNxmax(i);x++){ 
      double xp	= B->XNode1[i][x];
      if ( fMuFFunc != kScale1 &&  fMuFFunc != kScale2 ) { // that't the standard case!
	for(unsigned int jS1=0;jS1<B->ScaleNodeScale1[i].size();jS1++){
	  for(unsigned int kS2=0;kS2<B->ScaleNodeScale2[i].size();kS2++){
	    double muf = CalcMu( kMuF , BBlocksSMCalc[0][0]->ScaleNodeScale1[i][jS1] ,  BBlocksSMCalc[0][0]->ScaleNodeScale2[i][kS2] , fScaleFacMuF );
	    B->PdfLcMuVar[i][x][jS1][kS2] = CalcPDFLinearCombDIS( GetXFX(xp,muf) , B->NSubproc );
	  }
	}
      }
      else if ( fMuFFunc == kScale2 ){	// speed up
	for(unsigned int kS2=0;kS2<B->ScaleNodeScale2[i].size();kS2++){
	  double muf = CalcMu( kMuF , 0 ,  BBlocksSMCalc[0][0]->ScaleNodeScale2[i][kS2] , fScaleFacMuF );
	  vector < double > buffer = CalcPDFLinearCombDIS( GetXFX(xp,muf) , B->NSubproc );
	  for(unsigned int jS1=0;jS1<B->ScaleNodeScale1[i].size();jS1++){
	     B->PdfLcMuVar[i][x][jS1][kS2] = buffer;
	  }
	}
      }
      else if ( fMuFFunc == kScale1 ){	// speed up
	for(unsigned int jS1=0;jS1<B->ScaleNodeScale1[i].size();jS1++){
	  double muf = CalcMu( kMuF , BBlocksSMCalc[0][0]->ScaleNodeScale1[i][jS1] , 0 , fScaleFacMuF );
	  vector < double > buffer = CalcPDFLinearCombDIS( GetXFX(xp,muf) , B->NSubproc );
	  for(unsigned int kS2=0;kS2<B->ScaleNodeScale2[i].size();kS2++){
	     B->PdfLcMuVar[i][x][jS1][kS2] = buffer;
	  }
	}
      }
    }
  }
}


//______________________________________________________________________________


void FastNLOReader::FillBlockBPDFLCsHHCv20( FastNLOBlockB* B ){
  int scaleVar		= B->Npow == ILOord ? 0 : fScalevar;
  double scalefac	= fScaleFacMuF/B->ScaleFac[0][scaleVar];

  vector < vector < double > > xfx; // PDFs of all partons
  if ( B->NScaleDep != 3 ){
    for(int i=0;i<NObsBin;i++){
      int nxmax = B->GetNxmax(i);
      int nxbins1 = B->Nxtot1[i]; // number of columns in half matrix
      xfx.resize(nxbins1);
      for(int j=0;j<B->Nscalenode[0];j++){
	// determine all pdfs of hadron1
	for(int k=0;k<nxbins1;k++){ 
	  double xp	= B->XNode1[i][k];
	  double muf	= scalefac * B->ScaleNode[i][0][scaleVar][j];
	  xfx[k]	= GetXFX(xp,muf);
	}
	int x1bin = 0;
	int x2bin = 0;
	for(int k=0;k<nxmax;k++){ 
	  // ----- if pp ---- //
	  if ( B->NPDFPDG[0] == B->NPDFPDG[1] ){
	    B->PdfLc[i][j][k] = CalcPDFLinearCombHHC( xfx[x2bin], xfx[x1bin], B->NSubproc );
	  }
	  // ----- if ppbar ---- //
	  else if ( B->NPDFPDG[0] == -B->NPDFPDG[1] ){
	    vector < double > xfxbar(13);
	    for ( unsigned int p = 0 ; p<13 ; p++ ){
	      xfxbar[p] = xfx[x1bin][12-p];
	    }
	    B->PdfLc[i][j][k] = CalcPDFLinearCombHHC( xfx[x2bin], xfxbar, B->NSubproc );
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
      }
    }
  }
}


//______________________________________________________________________________



void FastNLOReader::FillBlockBPDFLCsHHCv21( FastNLOBlockB* B ){
  if ( B->PdfLcMuVar.empty() ) { cout<< "empty."<<endl; exit(1);}
  vector < vector < double > > xfx; // PDFs of all partons
  for(int i=0;i<NObsBin;i++){
    int nxmax = B->GetNxmax(i);
    int nxbins1 = B->Nxtot1[i]; // number of columns in half matrix
    xfx.resize(nxbins1);

    if ( fMuFFunc != kScale1 &&  fMuFFunc != kScale2 )  { // that't the standard case!
      for(unsigned int jS1=0;jS1<B->ScaleNodeScale1[i].size();jS1++){
	for(unsigned int kS2=0;kS2<B->ScaleNodeScale2[i].size();kS2++){
	  // determine all pdfs of hadron1
	  for(int k=0;k<nxbins1;k++){ 
	    double muf = CalcMu( kMuF , BBlocksSMCalc[0][0]->ScaleNodeScale1[i][jS1] ,  BBlocksSMCalc[0][0]->ScaleNodeScale2[i][kS2] , fScaleFacMuF );
	    double xp	= B->XNode1[i][k];
	    xfx[k] = GetXFX(xp,muf);
	  }
	  int x1bin = 0;
	  int x2bin = 0;
	   
	  for(int x=0;x<nxmax;x++){ 
	    // ----- if pp ---- //
	    if ( B->NPDFPDG[0] == B->NPDFPDG[1] ){
	      B->PdfLcMuVar[i][x][jS1][kS2] = CalcPDFLinearCombHHC( xfx[x2bin], xfx[x1bin], B->NSubproc );
	    }
	    // ----- if ppbar ---- //
	    else if ( B->NPDFPDG[0] == -B->NPDFPDG[1] ){
	      vector < double > xfxbar(13);
	      for ( unsigned int p = 0 ; p<13 ; p++ ){
		xfxbar[p] = xfx[x1bin][12-p];
	      }
	      B->PdfLcMuVar[i][x][jS1][kS2] = CalcPDFLinearCombHHC( xfx[x2bin], xfxbar, B->NSubproc );
	    }
	    else {
	      printf("FastNLOReader::FillBlockBPDFLCsHHCv21(). This is not pp, nor ppbar, nor pbarpbar!\n"); exit(1);
	    }
	    x1bin++;
	    if(x1bin>x2bin){
	      x1bin = 0;
	      x2bin++;
	    }
	  }
	}
      }
    }
    else if ( fMuFFunc == kScale2 ){	// speed up
      for(unsigned int kS2=0;kS2<B->ScaleNodeScale2[i].size();kS2++){
	// determine all pdfs of hadron1
	for(int k=0;k<nxbins1;k++){ 
	  double muf = CalcMu( kMuF , 0 ,  BBlocksSMCalc[0][0]->ScaleNodeScale2[i][kS2] , fScaleFacMuF );
	  double xp	= B->XNode1[i][k];
	  xfx[k] = GetXFX(xp,muf);
	}
	for(unsigned int jS1=0;jS1<B->ScaleNodeScale1[i].size();jS1++){
	  int x1bin = 0;
	  int x2bin = 0;
	  for(int x=0;x<nxmax;x++){ 
	    // ----- if pp ---- //
	    if ( B->NPDFPDG[0] == B->NPDFPDG[1] ){
	      B->PdfLcMuVar[i][x][jS1][kS2] = CalcPDFLinearCombHHC( xfx[x2bin], xfx[x1bin], B->NSubproc );
	    }
	    // ----- if ppbar ---- //
	    else if ( B->NPDFPDG[0] == -B->NPDFPDG[1] ){
	      vector < double > xfxbar(13);
	      for ( unsigned int p = 0 ; p<13 ; p++ ){
		xfxbar[p] = xfx[x1bin][12-p];
	      }
	      B->PdfLcMuVar[i][x][jS1][kS2] = CalcPDFLinearCombHHC( xfx[x2bin], xfxbar, B->NSubproc );
	    }
	    else {
	      printf("FastNLOReader::FillBlockBPDFLCsHHCv21(). This is not pp, nor ppbar, nor pbarpbar!\n"); exit(1);
	    }
	    x1bin++;
	    if(x1bin>x2bin){
	      x1bin = 0;
	      x2bin++;
	    }
	  }
	}
      }
    }
    else if ( fMuFFunc == kScale1 ){	// speed up
      for(unsigned int jS1=0;jS1<B->ScaleNodeScale1[i].size();jS1++){
	// determine all pdfs of hadron1
	for(int k=0;k<nxbins1;k++){ 
	  double muf = CalcMu( kMuF , BBlocksSMCalc[0][0]->ScaleNodeScale1[i][jS1] , 0 , fScaleFacMuF );
	  double xp	= B->XNode1[i][k];
	  xfx[k] = GetXFX(xp,muf);
	}
	for(unsigned int kS2=0;kS2<B->ScaleNodeScale2[i].size();kS2++){
	  int x1bin = 0;
	  int x2bin = 0;
	  for(int x=0;x<nxmax;x++){ 
	    // ----- if pp ---- //
	    if ( B->NPDFPDG[0] == B->NPDFPDG[1] ){
	      B->PdfLcMuVar[i][x][jS1][kS2] = CalcPDFLinearCombHHC( xfx[x2bin], xfx[x1bin], B->NSubproc );
	    }
	    // ----- if ppbar ---- //
	    else if ( B->NPDFPDG[0] == -B->NPDFPDG[1] ){
	      vector < double > xfxbar(13);
	      for ( unsigned int p = 0 ; p<13 ; p++ ){
		xfxbar[p] = xfx[x1bin][12-p];
	      }
	      B->PdfLcMuVar[i][x][jS1][kS2] = CalcPDFLinearCombHHC( xfx[x2bin], xfxbar, B->NSubproc );
	    }
	    else {
	      printf("FastNLOReader::FillBlockBPDFLCsHHCv21(). This is not pp, nor ppbar, nor pbarpbar!\n"); exit(1);
	    }
	    x1bin++;
	    if(x1bin>x2bin){
	      x1bin = 0;
	      x2bin++;
	    }
	  }
	}
      }
    }
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

  double SumQ1  = 0;
  double SumQB1 = 0;
  double SumQ2  = 0;
  double SumQB2 = 0;
  vector <double> Q1(6);
  vector <double> QB1(6);
  vector <double> Q2(6);
  vector <double> QB2(6);
  for ( int k = 0 ; k<6 ; k++ ){
    Q1[k]  = pdfx1[k+7];  //! read 1st PDF at x1
    QB1[k] = pdfx1[5-k];
    SumQ1  += Q1[k];
    SumQB1 += QB1[k];
    Q2[k]  = pdfx2[k+7];//  ! read 2nd PDF at x2
    QB2[k] = pdfx2[5-k];
    SumQ2  += Q2[k];
    SumQB2 += QB2[k];
  }
  double G1     = pdfx1[6];
  double G2     = pdfx2[6]; 

  //   - compute S,A
  double S = 0;
  double A = 0;
  for ( int k = 0 ; k<6 ; k++ ){
    S += (Q1[k]*Q2[k]) + (QB1[k]*QB2[k]);
    A += (Q1[k]*QB2[k]) + (QB1[k]*Q2[k]);
  }
   
  //c   - compute seven combinations
  vector <double> H(7);
  H[0]	= G1*G2;
  H[1] = SumQ1*SumQ2 + SumQB1*SumQB2 - S;
  H[2] = S;
  H[3] = A;
  H[4] = SumQ1*SumQB2 + SumQB1*SumQ2 - A;
  H[5] = (SumQ1+SumQB1)*G2;
  H[6] = G1*(SumQ2+SumQB2);

  if (NSubproc == 6) {
    H[5] += H[6];
    H.resize(6);
  }
  return H;   

}


//______________________________________________________________________________


void FastNLOReader::SetExternalFuncForMuR( double (*Func)(double,double)  , bool Verbose ){
  if ( BBlocksSMCalc[0][0]->NScaleDep != 3 ) {
    printf("FastNLOReader::SetExternalFuncForMuR. WARNING! This is not a flexible-scale table.\n");
    printf("  SetFunctionalForm has no impact.\n");
    printf("  Please use a flexible-scale table, if you want to change your scale definition.\n");
    return;
  }

  Fct_MuR = Func;
  SetFunctionalForm( kExtern , kMuR );
  if ( Verbose ) {
    printf(" *  FastNLOReader::SetExternalFuncForMuR(). Test.\n");
    printf(" *    Scale1 = 1 ,      Scale2 = 1        ->  mu = func(1,1)             = %9.4f\n",(*Fct_MuR)(1,1));
    printf(" *    Scale1 = 91.1876, Scale2 = 91.1876  ->  mu = func(91.1876,91.1876) = %9.4f\n",(*Fct_MuR)(91.1876,91.1876));
    printf(" *    Scale1 = 1,       Scale2 = 91.1876  ->  mu = func(1,91.1876)       = %9.4f\n",(*Fct_MuR)(1,91.1876));
    printf(" *    Scale1 = 91.1876, Scale2 = 1        ->  mu = func(91.1876,1)       = %9.4f\n",(*Fct_MuR)(91.1876,1));
  }
  FillAlphasCache();
}


//______________________________________________________________________________


void FastNLOReader::SetExternalFuncForMuF( double (*Func)(double,double)  , bool ReFillCache , bool Verbose ){
  if ( BBlocksSMCalc[0][0]->NScaleDep != 3 ) {
    printf("FastNLOReader::SetExternalFuncForMuF. WARNING! This is not a flexible-scale table.\n");
    printf("  SetFunctionalForm has no impact.\n");
    printf("  Please use a flexible-scale table, if you want to change your scale definition.\n");
    return;
  }

  Fct_MuF = Func;
  SetFunctionalForm( kExtern , kMuF );
  if ( Verbose ) {
    printf(" *  FastNLOReader::SetExternalFuncForMuF(). Test.\n");
    printf(" *    Scale1 = 1 ,      Scale2 = 1        ->  mu = func(1,1)             = %9.4f\n",(*Fct_MuF)(1,1));
    printf(" *    Scale1 = 91.1876, Scale2 = 91.1876  ->  mu = func(91.1876,91.1876) = %9.4f\n",(*Fct_MuF)(91.1876,91.1876));
    printf(" *    Scale1 = 1,       Scale2 = 91.1876  ->  mu = func(1,91.1876)       = %9.4f\n",(*Fct_MuF)(1,91.1876));
    printf(" *    Scale1 = 91.1876, Scale2 = 1        ->  mu = func(91.1876,1)       = %9.4f\n",(*Fct_MuF)(91.1876,1));
  }
  if ( ReFillCache ){
    FillPDFCache();
  }
}



//______________________________________________________________________________


void FastNLOReader::StripWhitespace(string* s){
  string fastlast = &(*s)[s->size()-1];
  while ( !fastlast.compare(" ")){
    string::iterator it = s->end();
    s->erase(it-1);
    fastlast = &(*s)[s->size()-1];
  }
}

//______________________________________________________________________________


int FastNLOReader::ContrId( const ESMCalculation eCalc, const ESMOrder eOrder ) const {
  
  int Id = -1;
  
  if ( BBlocksSMCalc.empty() || bUseSMCalc[eCalc].empty() ){
    return Id;
  }

  // Requested order
  string requested = fOrdName[eCalc][eOrder];
  // Loop over all available orders of contribution type eCalc 
  for(unsigned int i=0; i<BBlocksSMCalc[eCalc].size(); i++){
    // Found order
    string available = fOrdName[BBlocksSMCalc[eCalc][i]->IContrFlag1-1][BBlocksSMCalc[eCalc][i]->IContrFlag2-1];
    if ( available == requested ) {Id = i;} 
  }
  return Id;
  
}

//______________________________________________________________________________


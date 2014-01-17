// Author: Daniel Britzger
// DESY, 29/07/2013
#include <string>
#include <algorithm>
#include "fastnlotk/fastNLOCreate.h"
#include "fastnlotk/fastNLOTools.h"
#include "fastnlotk/read_steer.h"

#include "fastnlotk/fastNLOCoeffAddFlex.h"
#include "fastnlotk/fastNLOCoeffAddFix.h"

using namespace std;


// fastNLOCreate* fastNLOCreate::fInstance = NULL;

int fastNLOCreate::nInst = 0;


// ___________________________________________________________________________________________________
fastNLOCreate::fastNLOCreate() 
{
   SetClassName("fastNLOCreate");
}



// ___________________________________________________________________________________________________
fastNLOCreate::fastNLOCreate(string steerfile) 
{
   //speaker::SetGlobalVerbosity(say::DEBUG);
   SetClassName("fastNLOCreate");
   nInst++;
   if ( nInst > 1 ) {
      error["fastNLOCreate"]<<"Only one instance of fastNLOCreate is currently possible in one program call."<<endl;
      error>>"This is due to lazy implementation of the reading of the parser-values from 'read_steer'."<<endl;
      error>>"To fix this limitation, read in every steerfile/fastNLOCreate into a dedicated read_steer-namespace and access the values from there."<<endl;
      exit(1);
   }
   ResetHeader();
   ReadSteering(steerfile);

   // Try to get warm-up values.
   // Otherwise a warm-up run will be initialized.
   GetWarmupValues();
   
   // init bin grid
   if ( fIsWarmup )   ReadBinning();			// if warmup, then always read binning from steering.
   else if ( BOOL(ReadBinningFromSteering) ) {
      ReadBinning();
      CheckWarmupConsistency();
   }
   else  UseBinGridFromWarmup();

   // now create one coefficient tasble.
   InitCoeffTable();
   SetOrderOfAlphasOfCalculation(fIOrd);

   // Init interpolation kernels
   if ( !fIsWarmup ) {
      InitInterpolationKernels();
      InitGrids();
   }
}


// ___________________________________________________________________________________________________
fastNLOCreate::~fastNLOCreate()
{
   // todo. cleanup arrays of kernels.
}


// ___________________________________________________________________________________________________
void fastNLOCreate::ReadSteering(string steerfile) 
{
   debug["ReadSteering"]<<"Steerfile = "<<steerfile<<endl;
   fSteerfile =  steerfile;
   READ(steerfile);
   PRINTALL();

   // header
   SetScenName(STRING(ScenarioName));
   SetItabversion(20000);

   // scenario specific things
   Ipublunits   = INT(PublicationUnits);
   ScDescript   = STRING_ARR(ScenarioDescription);
   Ecms         = DOUBLE(CenterOfMassEnergy);   // is often superseeded by generator-specific code.
   ILOord       = INT(LeadingOrder);            // is often superseeded by generator-specific code.
   INormFlag    = 0;
   fIOrd        = INT(OrderInAlphasOfCalculation);// is often superseeded by generator-specific code.
   SetFilename(STRING(OutputFilename));

   fIsFlexibleScale = BOOL(FlexibleScaleTable);
   fApplyPDFReweight = BOOL(ApplyPDFReweighting);
   SetOutputPrecision(INT(OutputPrecision));

   //if ( !fIsFlexibleScale )  ReadScaleFactors() ; // is called only when setting order of calculation

   //    //KR: Added possibility to store and read start of new rapidity bin in nobs
   //    //vector <int> RapIndex; ?
}


// ___________________________________________________________________________________________________
void fastNLOCreate::ReadScaleFactors()
{
   //! read scale factors from steering
   //! and init member fScaleFac

   if ( fIsFlexibleScale ) {warn["ReadScaleFactors"]<<"This function is only reasonable for fixed-scale tables!"<<endl;}
   vector<double> svar = DOUBLE_ARR(ScaleVariationFactors);
   fScaleFac.resize(svar.size());
   if ( svar.empty() ) { 
      // 'ScaleVariationFactors' not found -> using default
      warn["ReadScaleFactors"]<<"No list of scale-factors found in steering file. Using only scale-factor of '1'."<<endl;
      fScaleFac.push_back(1.0);
   }
   else if (  GetTheCoeffTable() && GetTheCoeffTable()->IsLO() ) {
      // scale-factors not needed -> using default of 1.0
      info["ReadScaleFactors"]<<"This is a leading-order run. There is no scale dependence in LO. Using only scalefactor of 1."<<endl;
      fScaleFac.resize(1);
      fScaleFac[0] = 1.0;
   }
   else if ( fIsWarmup ) {
      // scale-factors not needed -> using default of 1.0
      info["ReadScaleFactors"]<<"This is a warmup run. Using only scalefactor of 1."<<endl;
      fScaleFac.resize(1);
      fScaleFac[0] = 1.0;
   }
   else  {
      // sort list according to fastNLO conventions
      //  0:  1.0
      //  1-n:  nmin...nmax (without 1.0)
      vector<double> vtemp;
      bool foundunity = false;
      for ( unsigned int k = 0 ;k<svar.size(); k++ ){
	 if ( fabs(svar[k]-1.0)<1.e-8  && foundunity ) {
	    info["ReadScaleFactors"]<<"Found scale factor 1.0 two times in list ScaleVariationFactors. Ignoring second appearance."<<endl;
	    fScaleFac.resize(fScaleFac.size()-1);
	 }
	 else if ( fabs(svar[k]-1.0)<1.e-8 ) foundunity = true;
	 else vtemp.push_back(svar[k]);
      }
      if (!foundunity ) {
	 error["ReadScaleFactors"]<<"Could not found scalefactor of 1.0 in list ScaleVariationFactors. Exiting."<<endl;
	 exit(1);
      }
      sort(vtemp.begin(),vtemp.end());

      fScaleFac[0] = 1.0;
      int s = 0;
      for ( unsigned int k = 0 ;k<vtemp.size(); k++ ) {
	 fScaleFac[s+1] = vtemp[k];
	 if ( fScaleFac[s+1] == fScaleFac[s] ) {
	    info["ReadScaleFactors"]<<"Found scale factor '"<<fScaleFac[k+1]<<"' two times in list ScaleVariationFactors. Ignoring second appearance."<<endl;
	    fScaleFac.resize(fScaleFac.size()-1);
	    s--;
	 }
	 s++;
      }
   }
   info["ReadScaleFactors"]<<"Using following scale factors:"<<endl;
   for ( unsigned int k = 0 ;k<fScaleFac.size(); k++ )
      info>>"\tiScaleVar "<<k<<":\t"<< fScaleFac[k]<<endl;
}


// ___________________________________________________________________________________________________
void fastNLOCreate::InitCoeffTable()
{
   debug["InitCoeffTable"]<<endl;
   //! create a coeff table
   //CreateCoeffTable(0);
   CreateCoeffTable();

   //! set 'usual' variables for perturbative calculations
   InitVariablesInCoefficientTable();

   //! read in process specific variables
   ReadCoefficientSpecificVariables();
}


// ___________________________________________________________________________________________________
int fastNLOCreate::CreateCoeffTable()
{
   debug["CreateCoeffTable"]<<endl;
   if ( !fCoeff.empty() ){
      error["CreateCoeffAddFix"]<<"Vector of coefficients must be empty, since only one coefficient table is allowed."<<endl;
      exit(1);
   }
   if (fIsFlexibleScale )
      return fastNLOTable::CreateCoeffTable(fCoeff.size(), new fastNLOCoeffAddFlex(NObsBin,ILOord) );
   else
      return fastNLOTable::CreateCoeffTable(fCoeff.size(), new fastNLOCoeffAddFix(NObsBin) );
}


// ___________________________________________________________________________________________________
void fastNLOCreate::ReadBinning()
{
   // optimize read-in of bin grids
   // ToDo. Check sanity of bin-grid

   NDim         = INT(DifferentialDimension);
   //Scenario.SetNDim(NDim);
   if ( (int)STRING_ARR(DimensionLabels).size() < NDim ){
      error["ReadBinning"]<<"Each dimension needs a bin label. Exiting."<<endl;
      exit(1);
   }
   DimLabel     = STRING_ARR(DimensionLabels);
   if ( (int)INT_ARR(DimensionIsDifferential).size() < NDim ){
      error["ReadBinning"]<<"Each dimension need to specify if differential or not. Exiting."<<endl;
      exit(1);
   }
   IDiffBin     = INT_ARR(DimensionIsDifferential);
   DimLabel.resize(NDim); //safety
   IDiffBin.resize(NDim);

   bool AllDiff = true;
   bool AllBinInt = true;
   for ( unsigned int i = 0 ; i<IDiffBin.size() ; i++ ){
      AllDiff = AllDiff && (IDiffBin[i] == 1);
      AllBinInt = AllBinInt && (IDiffBin[i] != 1);
   }
   if ( !AllDiff && !AllBinInt ) {
      error["ReadBinning"]<<"All dimensions must be consistently either bin-integrated, or truly differential dimensions. Exiting."<<endl;
      exit(1);
   }

   if ( AllDiff && NDim == 3) { error["ReadBinning"]<<"Fully differential and triple-differential binning not yet implemented. exiting"<<endl;exit(1);}

   // read single-differential bin grid
   if ( NDim == 1 ) {
      vector<double> bgrid = DOUBLE_ARR(SingleDiffBinning);
      for ( unsigned int i = 0 ; i<bgrid.size()-1 ; i++ ){
         if ( bgrid[i] >= bgrid[i+1] ) {
            error["ReadBinning"]<<"The upper bin edge is below the lower one in bin "<<i+1<<". Exiting."<<endl;
            exit(1);
         }
      }
      if ( AllBinInt ) { // bin integrated bin grid
         NObsBin = bgrid.size()-1;
         Bin.resize(NObsBin);
         for ( unsigned int i = 0 ; i<bgrid.size()-1 ; i++ ){
            Bin[i].resize(1);
            Bin[i][0] = make_pair(bgrid[i],bgrid[i+1] );
         }
      }
      else { // truly differential bin grid
         NObsBin = bgrid.size();
         Bin.resize(NObsBin);
         for ( unsigned int i = 0 ; i<bgrid.size() ; i++ ){
            Bin[i].resize(1);
            Bin[i][0] = make_pair(bgrid[i],bgrid[i]) ;
         }
      }
   }

   // read double-differential bin grid
   else if ( NDim==2 ) {
      vector<vector<double> > in = DOUBLE_TAB(DoubleDifferentialBinning);
      NObsBin=0;
      Bin.clear();
      for ( unsigned int r = 0 ; r<in.size() ; r++ ){
         unsigned int nBin2Max = AllBinInt ? in[r].size()-1 : in[r].size();
         for ( unsigned int c = 2 ; c<nBin2Max ; c++ ){
            Bin.push_back(vector<pair<double,double> >(NDim));
            // sanity dim 1:
            if ( AllBinInt ) {
               if ( in[r][0]>=in[r][1] ){
                  error["ReadBinning"]<<"The upper bin edge ("<<in[r][1]<<") is below the lower one ("<<in[r][0]<<") in row "<<r+1<<". Exiting."<<endl;
                  exit(1);
               }
               if ( AllBinInt && r>0 && in[r][0]!=in[r-1][1] ) {
                  error["ReadBinning"]<<"The lower bin edge ("<<in[r][0]
                                      <<") is not identical to the upper bin edge to the previous bin ("<<in[r-1][1]<<") around row "<<r+2<<". Exiting."<<endl;
                  exit(1);
               }
               Bin[NObsBin][1] = make_pair(in[r][0],in[r][1]);
               // sanity dim 0:
               if ( in[r][c] >= in[r][c+1] ) {
                  error["ReadBinning"]<<"The upper bin edge ("<<in[r][c+1]<<") is below the lower one ("<<in[r][c]<<") in row "<<r+1<<" and column "<<c+1<<". Exiting."<<endl;
                  exit(1);
               }
               Bin[NObsBin][0] = make_pair(in[r][c],in[r][c+1]);
            }
            else {
               Bin[NObsBin][1] = make_pair(in[r][0],in[r][0]);
               Bin[NObsBin][0] = make_pair(in[r][c],in[r][c]);
            }
            NObsBin++; // count
         }
      }
   }

   // read in triple-differential binning
   else if ( NDim==3 ) {
      warn["ReadBinning"]<<"The code for reading of "<<NDim<<"-dimensional binnings was not fully tested. Please verify the code an remove this statement."<<endl;
      vector<vector<double> > in = DOUBLE_TAB(TripleDifferentialBinning);
      NObsBin=0;
      Bin.clear();
      for ( unsigned int r = 0 ; r<in.size() ; r++ ){
         if ( in[r].size() < 6 ) {
            warn["ReadBinning"]<<"At least six numbers are necessary to specify a 3-dimensional binning in row"<<r+1<<endl;
         }
         for ( unsigned int c = 4 ; c<in[r].size()-1 ; c++ ){
            NObsBin++;
            Bin.push_back(vector<pair<double,double> >(NDim));
            // sanity dim 2:
            if ( in[r][0]>=in[r][1] ){
               error["ReadBinning"]<<"The upper bin edge ("<<in[r][1]<<") is below the lower one ("<<in[r][0]<<") in row "<<r+1<<". Exiting."<<endl;
               exit(1);
            }
            Bin[NObsBin-1][2] = make_pair(in[r][0],in[r][1]);
            // sanity dim 1:
            if ( in[r][2]>=in[r][3] ){
               error["ReadBinning"]<<"The upper bin edge ("<<in[r][3]<<") is below the lower one ("<<in[r][2]<<") in row "<<r+1<<". Exiting."<<endl;
               exit(1);
            }
            Bin[NObsBin-1][1] = make_pair(in[r][2],in[r][3]);
            // sanity dim 0:
            if ( in[r][c] >= in[r][c+1] ) {
               error["ReadBinning"]<<"The upper bin edge ("<<in[r][c+1]<<") is below the lower one ("<<in[r][c]<<") in row "<<r+1<<" and column "<<c+1<<". Exiting."<<endl;
               exit(1);
            }
            Bin[NObsBin-1][0] = make_pair(in[r][c],in[r][c+1]);
         }
      }
   }
   else {
      error["ReadBinning"]<<"Reading of "<<NDim<<"-binnings from steering is not yet implemented. Exiting"<<endl;
      exit(1);
   }


   // ---------------------------------
   //       Bin width
   // ---------------------------------
   if ( BOOL(CalculateBinWidth) ){
      BinSize.resize(NObsBin);
      bool idi = false;
      for ( int i = 0 ; i<NObsBin ; i++ ){
         BinSize[i] = 1;
         for ( int d=0 ; d<NDim ; d++ ){
            if (IDiffBin[d]==0 ) {
               // nothing todo
            }
            else if ( IDiffBin[d]==1) {
               // nothing todo
               //warn["ReadBinning"]<<"Don't know how to handle truly differential bins for bin widths."<<endl;
            }
            else if (IDiffBin[d]==2 ) {
               BinSize[i] *= Bin[i][d].second-Bin[i][d].first;//UpBin[i][d]-LoBin[i][d];
               idi = true;
            }
         }
         // divide by binwidthfactor, but only if at least one dimension is differential
         if ( idi ) BinSize[i] *= DOUBLE(BinWidthFactor);
      }
      if (!idi) info["CalculateBinWidth"]<<"BinWidthFactor is not being used, since no observable is calculated differential."<<endl;
   }
   else {
      // read in bin width
      warn["ReadBinning"]<<"Reading of bindwidth only poorly  implemented! Improve it and remove this message."<<endl;
      if ( (int)DOUBLE_ARR(BinWidth).size()!=NObsBin) warn["ReadBinning"]<<"Number of bins of 'BinWidth' not consistent with bin grid."<<endl;
      BinSize=DOUBLE_ARR(BinWidth);
      BinSize.resize(NObsBin);
      for ( int i = 0 ; i<NObsBin ; i++ ){
	 if ( BinSize[i]==0 ) BinSize[i] = 1.0;
      }
   }

   info["ReadBinning"]<<"Read in "<<NDim<<"-dimensional bin grid with "<<NObsBin<<" bins succesfully."<<endl;
}


// ___________________________________________________________________________________________________
void fastNLOCreate::GetWarmupValues()
{
   //!
   //! GetWarmupValues.
   //! Checks if warmup-table exists and initialized
   //! member variable fIsWarmup
   //!
   debug["GetWarmupValues"]<<endl;

   // try to get warmup values
   vector<vector<double> > warmup = DOUBLE_TAB(WarmupValues);
   fIsWarmup = warmup.empty();

   // try again, with hard-coded convention:
   if ( fIsWarmup ) {
      info["GetWarmupValues"]<<"Could not get warmup table from steerfile. Now trying to read steerfile: "<<GetWarmupTableFilename()<<endl;
      READ(GetWarmupTableFilename());   
      warmup = DOUBLE_TAB(WarmupValues);
      fIsWarmup = warmup.empty();
      if ( !fIsWarmup ) info["GetWarmupValues"]<<"Reading of file "<<GetWarmupTableFilename()<<" contained warmup values!"<<endl;
   }
   
   // inform user about success
   info["GetWarmupValues"]<<"This will be "<<(fIsWarmup?"":"not")<<" a warmup run."<<endl;
}



// ___________________________________________________________________________________________________
void fastNLOCreate::UseBinGridFromWarmup()
{
   //! initialialize all binning related variables
   //! with values stored in the warmup file.
   vector<vector<double> > warmup =  DOUBLE_TAB(WarmupValues);
   NObsBin	= warmup.size();
   NDim		= INT(Warmup.DifferentialDimension);
   if ( (int)warmup[0].size() != (7+2*NDim) && (int)warmup[0].size() != (5+2*NDim) ) {
      error["UseBinGridFromWarmup"]<<"This warmup table has an unknown size of columns. Expecting "<<(7+2*NDim)<<" for flexible-scale, or "<<(5+2*NDim)<<" for fixed-scale tables. Exiting."<<endl;
      exit(1);
   }
   fIsFlexibleScale = ((int)warmup[0].size() == (7+2*NDim));
   IDiffBin     = INT_ARR(Warmup.DimensionIsDifferential);
   DimLabel	= STRING_ARR(Warmup.DimensionLabels);
   
   // make binning 
   const int i0 = fIsFlexibleScale ? 6 : 4;
   Bin.resize(NObsBin);
   BinSize.resize(NObsBin);
   for ( int i = 0 ; i < NObsBin ; i ++ ) {
      Bin[i].resize(NDim);
      if ( NDim==1 ) {
	 Bin[i][0] = make_pair(warmup[i][i0],warmup[i][i0+1]) ;
      }
      else if ( NDim==2 ) {
	 Bin[i][1] = make_pair(warmup[i][i0],warmup[i][i0+1]) ;
	 Bin[i][0] = make_pair(warmup[i][i0+2],warmup[i][i0+3]) ;
      }
      else if ( NDim==3 ) {
	 Bin[i][2] = make_pair(warmup[i][i0],warmup[i][i0+1]) ;
	 Bin[i][1] = make_pair(warmup[i][i0+2],warmup[i][i0+3]) ;
	 Bin[i][0] = make_pair(warmup[i][i0+4],warmup[i][i0+5]) ;
      }
      BinSize[i] = warmup[i][i0+NDim*2];
   }
}



// ___________________________________________________________________________________________________
bool fastNLOCreate::CheckWarmupConsistency()
{
   //! check if warmup values are consistent with steering card
   //! check if number of bins is consistent

   vector<vector<double> > warmup =  DOUBLE_TAB(WarmupValues);
   bool ret = true;
   cout<<endl;
   const string wrmuphelp = "Please remove warmup-file in order to calculate a new warmup-file which is compatible to your steering,\nor alternatively use 'ReadBinningFromSteering=false', then all binning-related information is taken from the warmup file.\n";
   
   if ( (int)warmup.size() != NObsBin ) {
      error["GetWarmupValues"]
	 <<"Table of warmup values is not compatible with steering file.\n"
	 <<"Different number of bins ("<<warmup.size()<<" instead of "<<NObsBin<<".\n"
	 <<wrmuphelp
	 <<"Exiting."<<endl;
      ret = false;
      exit(1);
   }
   if ( INT(Warmup.DifferentialDimension) != NDim  ) {
      error["GetWarmupValues"]
	 <<"Table of warmup values is not compatible with steering file.\n"
	 <<"Found different number of dimensions. NDim="<<NDim<<", Warmup.DifferentialDimension="<<INT(Warmup.DifferentialDimension)<<".\n"
	 <<wrmuphelp
	 <<"Exiting."<<endl;
      ret = false;
      exit(1);
   }

   // CoeffTable is not available during intialization
   //    if ( STRING(Warmup.ScaleDescriptionScale1) != GetTheCoeffTable()->ScaleDescript[0][0] ) {
   //       warn["GetWarmupValues"]
   // 	 <<"Table of warmup values is potentially incompatible with steering file.\n"
   // 	 <<"Found different scale description (ScaleDescriptionScale1). ScaleDescriptionScale1="<<GetTheCoeffTable()->ScaleDescript[0][0]
   // 	 <<", but Warmup.ScaleDescriptionScale1='"<<STRING(Warmup.ScaleDescriptionScale1)<<"'."<<endl<<endl;
   //       ret = false;
   //    }
   //    if ( fIsFlexibleScale && STRING(Warmup.ScaleDescriptionScale2) != GetTheCoeffTable()->ScaleDescript[0][1] ){
   //       warn["GetWarmupValues"]
   // 	 <<"Table of warmup values is potentially incompatible with steering file.\n"
   // 	 <<"Found different scale description (ScaleDescriptionScale2). ScaleDescriptionScale2='"<<GetTheCoeffTable()->ScaleDescript[0][0]<<"'"
   // 	 <<", but Warmup.ScaleDescriptionScale2='"<<STRING(Warmup.ScaleDescriptionScale1)<<"'."<<endl<<endl;
   //       ret = false;
   //    }

   if ( INT_ARR(Warmup.DimensionIsDifferential)[0] != IDiffBin[0]  ) {
      warn["GetWarmupValues"]
	 <<"Table of warmup values seems to be incompatible with steering file.\n"
	 <<"Found different diff-label for dimension 0  (IDiffBin). DimensionIsDifferential='"<<IDiffBin[0]<<"'"
	 <<", but Warmup.DimensionIsDifferential[0]='"<<DOUBLE_ARR(Warmup.DimensionIsDifferential)[0]<<"'. Exiting."<<endl;
      ret = false;
      exit(1);
   }
   if ( NDim > 1 && INT_ARR(Warmup.DimensionIsDifferential)[1] != IDiffBin[1]  ) {
      warn["GetWarmupValues"]
	 <<"Table of warmup values seems to be incompatible with steering file.\n"
	 <<"Found different diff-label for dimension 0  (IDiffBin). DimensionIsDifferential='"<<IDiffBin[1]<<"'"
	 <<", but Warmup.DimensionIsDifferential[0]='"<<DOUBLE_ARR(Warmup.DimensionIsDifferential)[1]<<"'. Exiting."<<endl;
      ret = false;
      exit(1);
   }
	 
   // check bining in detail
   for ( int i = 0 ; i < GetNObsBin() ; i ++ ) {
      int i0 = fIsFlexibleScale ? 6 : 4;
      if ( NDim == 1 ) {
	 if ( Bin[i][0].first != warmup[i][i0] || Bin[i][0].second != warmup[i][i0+1] ) {
	    error["GetWarmupValues"]
	       <<"Table of warmup values seems to be incompatible with steering file.\n"
	       <<"Found different binning for bin "<<i<<", steering: ["<<Bin[i][0].first<<","<<Bin[i][0].second
	       <<",], warmup: ["<<warmup[i][i0]<<","<<warmup[i][i0+1]<<"].\n"
	       <<wrmuphelp
	       <<"Exiting."<<endl;
	    ret = false;
	    exit(1);
	 }
      }
      else if ( NDim == 2 ) {
	 if ( Bin[i][1].first != warmup[i][i0] || Bin[i][1].second != warmup[i][i0+1] 
	      || Bin[i][0].first != warmup[i][i0+2] || Bin[i][0].second != warmup[i][i0+3] ) {
	    error["GetWarmupValues"]
	       <<"Table of warmup values seems to be incompatible with steering file.\n"
	       <<"Found different binning for bin "<<i<<", steering: ["<<Bin[i][1].first<<","<<Bin[i][1].second<<",] ["<<Bin[i][0].first<<","<<Bin[i][0].second
	       <<"], warmup: ["<<warmup[i][i0]<<","<<warmup[i][i0+1]<<"] ["<<warmup[i][i0+2]<<","<<warmup[i][i0+3]<<"].\n"
	       <<wrmuphelp
	       <<"Exiting."<<endl;
	    ret = false;
	    exit(1);
	 }
      }
      //check bin width
      double bwwrm = 0;
      if ( NDim == 1 ) bwwrm = warmup[i][i0+2];
      else if ( NDim == 2 ) bwwrm = warmup[i][i0+4];
      else if ( NDim == 3 ) bwwrm = warmup[i][i0+6];
      if ( fabs(BinSize[i] - bwwrm ) > 1.e-6 ) {
	 warn["GetWarmupValues"]
	    <<"Table of warmup values seems to be incompatible with steering file.\n"
	    <<"Found different bin size for bin "<<i<<". Steering: "<<BinSize[i]
	    <<", warmup: "<<bwwrm<<".\n"
	    <<wrmuphelp
	    <<"Please check consistency!"<<endl<<endl;
	 ret = false;
      }
   }
   return ret;
}



// ___________________________________________________________________________________________________
void fastNLOCreate::SetOrderOfAlphasOfCalculation(unsigned int ord)
{
   debug["SetOrderOfAlphasOfCalculation"]<<"ord="<<ord<<endl;
   //! set order of alpha_s of this calculation
   //! it must be: iLeadingOrder + iHigherOrder ;
   //! for instance: 3-jet-production in NLO = 4!
   fastNLOCoeffAddBase* c = GetTheCoeffTable();
   const int IOrdInitial = fIOrd;
   fIOrd = ord;
   c->Npow = ord;
   c->IContrFlag2 = ord-GetLoOrder()+1;
   c->CtrbDescript.resize(1);
   string ss[5] = {"LO","NLO","NNLO","NNNLO","unknown"};
   int iContrb = ord-GetLoOrder() <4 ? ord-GetLoOrder() : 4;
   c->CtrbDescript[0] = ss[iContrb];
   if ( !(GetLoOrder() == 1 || GetLoOrder() == 2 || GetLoOrder() == 3 || GetLoOrder() == 4) ) {
      warn["SetOrderOfAlphasOfCalculation"]<<"The leading order is of order of "<<GetLoOrder()<<" in alpha_s."<<endl;
      warn>>"\tThis may be unreasonable."<<endl;
   }
   if      ( ( ord - GetLoOrder()) == 0 ) {
      c->NSubproc               = INT(NSubProcessesLO);
      c->IPDFdef3               = INT(IPDFdef3LO);
   }
   else if ( ( ord - GetLoOrder()) == 1 ) {
      c->NSubproc               = INT(NSubProcessesNLO);
      c->IPDFdef3               = INT(IPDFdef3NLO);
   }
   else if ( ( ord - GetLoOrder()) == 2 ) {
      c->NSubproc               = INT(NSubProcessesNNLO);
      c->IPDFdef3               = INT(IPDFdef3NNLO);
   }
   else {
      error["SetOrderOfAlphasOfCalculation"]<<"Unknown order of pertubation theory: order="<<ord-GetLoOrder()<<" (ord="<<ord<<",ILOord="<<ILOord<<"). Exiting."<<endl;
      exit(1);
   }
   info["SetOrderOfAlphasOfCalculation"]<<"Using "<<c->NSubproc<<" subprocesses and IPDFdef3="<<c->IPDFdef3<<endl;


   // init array with counter processes (symmetric and asymmetric ones)
   if ( c->NPDFPDG.size() == 2 && c->NPDFDim == 1 ){
      fSymProc.resize(c->NSubproc);
      for ( int p = 0 ; p<c->NSubproc ; p++ ) fSymProc[p]=p;
      vector<vector<int> > asym = INT_TAB(AsymmetricProcesses);
      for ( unsigned int i = 0 ; i<asym.size() ; i ++ ) {
	 if ( asym[i][0]<(int)fSymProc.size() ) { // safety
	    if ( asym[i][1] >= GetNSubprocesses() || fSymProc[asym[i][0]] >= GetNSubprocesses() ) {
	       warn["AsymmetricProcesses"]<<"Subprocess "<<fSymProc[asym[i][0]]<<" is requested to be asymmetric with subprocess "<<asym[i][1]<<", but there are only "<<GetNSubprocesses()-1<<" subprocesses in this calculation. Ignoring call."<<endl;
	    }
	    else fSymProc[asym[i][0]] = asym[i][1];
	 }
      }
   }
   
   // Scale factors have to be updated, since this may be either a lo or nlo run.
   if ( !fIsFlexibleScale )  ReadScaleFactors();
   // NSubproc may have changed. We have to reinitialize the grids
   //const int nSFInitial = fScaleFac.size();
   if ( IOrdInitial!= fIOrd && !fIsWarmup ) {
      if ( !fIsFlexibleScale ) InitInterpolationKernels();
      InitGrids();
   }
}



// ___________________________________________________________________________________________________
void fastNLOCreate::SetLoOrder(int LOOrd)
{
   debug["SetLoOrder"]<<endl;
   fastNLOTable::SetLoOrder(LOOrd);
   if ( fIsFlexibleScale )
      ((fastNLOCoeffAddFlex*)GetTheCoeffTable())->fILOord = LOOrd;
}



// ___________________________________________________________________________________________________
void fastNLOCreate::InitVariablesInCoefficientTable()
{
   debug["InitVariablesInCoefficientTable"]<<endl;
   fastNLOCoeffAddBase* c = GetTheCoeffTable();
   c->IDataFlag = 0;            // No data, but theory
   c->IAddMultFlag = 0;         // additive contribution.
   c->IContrFlag1 = 1;          // fixed order: 1
   c->IContrFlag2 = 42;         // init with arbitrary number. to be specified later.
   c->IRef  = 0;                // it is not a reference calculation
   c->IScaleDep = 100;
   c->NScaleDep = 0;
   //c->NFragFunc               = 0;
   c->NFFDim            = 0;
   c->Nevt              = 0;
   c->SetIXsectUnits(12);       // it is often pb
}



// ___________________________________________________________________________________________________
void fastNLOCreate::ReadCoefficientSpecificVariables()
{
    debug["ReadCoefficientSpecificVariables"]<<endl;
   // todo: make it more user friendly
   // todo: include some sanity checks
   fastNLOCoeffAddBase* c = GetTheCoeffTable();
   c->NPDFPDG.resize(INT(NPDF));
   if ( c->NPDFPDG.size() >0 ) c->NPDFPDG[0] = INT(PDF1);
   if ( c->NPDFPDG.size() >1 ) c->NPDFPDG[1] = INT(PDF2);
   c->NPDFDim           = INT(NPDFDim);
   c->CodeDescript      = STRING_ARR(CodeDescription);
   c->IPDFdef1          = INT(IPDFdef1);
   c->IPDFdef2          = INT(IPDFdef2);
   c->IPDFdef3          = -1 ;          // safe initialization, is initialized in SetOrderOfAlphasOfCalculation(int); INT(IPDFdef3);
   c->NSubproc          = -1;           // safe initialization, is initialized in SetOrderOfAlphasOfCalculation(int);
   c->SetIXsectUnits(INT(UnitsOfCoefficients));
   c->NScaleDim = 1;  // NEVER SET NScaleDim TO ANY OTHER VALUE THAN 1 !!!

   //IPDFdef3 = NSubproc == 7 ? 2 : 1;
   //printf("         Set IPDFdef3 = %d, consistent with %d subprocesses.\n",IPDFdef3,NSubproc);
   const int NScales = 1;  
   c->Iscale.resize(NScales);
   c->Iscale[0] = 0;
   c->ScaleDescript.resize(NScales);

   if ( fIsFlexibleScale ){
      c->NScaleDep              = 3; // temporaily. Until known if generator runs in LO, NLO or NNLO.
      c->ScaleDescript[0].resize(2);
      if ( BOOL(ReadBinningFromSteering) ) {
	 c->ScaleDescript[0][0] = STRING(ScaleDescriptionScale1);
	 c->ScaleDescript[0][1] = STRING(ScaleDescriptionScale2);
      } else {
	 c->ScaleDescript[0][0] = STRING(Warmup.ScaleDescriptionScale1);
	 c->ScaleDescript[0][1] = STRING(Warmup.ScaleDescriptionScale2);
      }
      // ---- those numbers are partly ambigously defined in v2.1 ---- //
      // proper code would be, but needs also adjustment in reading, writing and getter functions!
      //       int NScales = 2;  
      //       c->Iscale.resize(NScales);	// Iscale is not used, but kept for consistency
      //       c->ScaleDescript.resize(NScales);
      //       c->Iscale[0] = 0;
      //       c->Iscale[1] = 1;
      //       c->ScaleDescript[0].resize(1);
      //       c->ScaleDescript[1].resize(1);
      //       c->ScaleDescript[0][0] = STRING(ScaleDescriptionScale1);
      //       c->ScaleDescript[1][0] = STRING(ScaleDescriptionScale2);
   } else {
      c->ScaleDescript[0].resize(1);
      if ( BOOL(ReadBinningFromSteering) )
	 c->ScaleDescript[0][0] = STRING(ScaleDescriptionScale1);
      else 
	 c->ScaleDescript[0][0] = STRING(Warmup.ScaleDescriptionScale1);
      c->SetNScaleDep(0);		// This is a fixed-scale table
   }
   
   //    if ( c->IsLO() ) c->IScaleDep = 0;
   //    else c->IScaleDep = 1;
   c->IScaleDep = 1;
}



// ___________________________________________________________________________________________________
int fastNLOCreate::GetBin(){
   //! get bin number, using
   //! observables from Scenario

   const int idiff = GetNumDiffBin();
   // -------------------------------
   // check cache and return if available
   if ( idiff == 1 ) {
      if ( fLastScen._o[0] == fScenario._o[0] ) return fObsBin;
   }
   else if ( idiff == 2 ) {
      if ( fLastScen._o[0] == fScenario._o[0] && fLastScen._o[1] == fScenario._o[1])  return fObsBin;
   }
   else if ( idiff == 3 ) {
      if ( fLastScen._o[0] == fScenario._o[0] && fLastScen._o[1] == fScenario._o[1] && fLastScen._o[2] == fScenario._o[2])  return fObsBin;
   }
   else {
      error["GetBin"]<<"Sorry. triple-differential binning not yet implemented. exiting."<<endl;
   }


   // -------------------------------
   // calc bin number and keep Observables
   if ( idiff == 1 ) fObsBin = GetBinNumber(fScenario._o[0]);
   else if ( idiff == 2 )  fObsBin = GetBinNumber(fScenario._o[1],fScenario._o[0]);
   //else if ( idiff == 3 )  fObsBin = GetBinNumber(fScenario._o[2],fScenario._o[1],fScenario._o[0]);
   else {
      error["GetBin"]<<"Sorry. triple-differential binning not yet implemented. exiting."<<endl;
   }
   fLastScen = fScenario;

   return fObsBin;
}



// ___________________________________________________________________________________________________
void fastNLOCreate::FillAllSubprocesses(const vector<vector<fnloEvent> >& events, const fnloScenario& scen)
{
   //! fill all subprocessess for all scale variations (into a fixed-scale table)
   //! events is expected to be of the form:
   //!   events[nscalevar][nsubproc]

   const bool bFasterCode = true; // experimental developement: try to make code faster
   if ( bFasterCode && !fIsWarmup && !fIsFlexibleScale) {
      // make filling code a little bit faster ... ~40%
      // if filling step "+=" is commented, then code is a factor ~6 faster
      fEvent = events[0][0];
      fScenario = scen;

      //if ( fEvent._w == 0 ) return; // nothing todo.

      const int ObsBin = (fScenario._iOB == -1) ? GetBin() : fScenario._iOB;
      if ( ObsBin < 0 ) return;
      if ( ObsBin >= GetNObsBin() ) return;
      fStats._nEvPS++;

   
      fastNLOCoeffAddFix* c = (fastNLOCoeffAddFix*)GetTheCoeffTable();
      // do interpolation
      double xmin = std::min(fEvent._x1,fEvent._x2);
      double xmax = std::max(fEvent._x1,fEvent._x2);
      vector<pair<int,double> > nxlo = fKernX[ObsBin]->GetNodeValues(xmin);
      vector<pair<int,double> > nxup = fKernX[ObsBin]->GetNodeValues(xmax);

      if ( fApplyPDFReweight ) {
	 fKernX[ObsBin]->CheckX(xmin);
	 fKernX[ObsBin]->CheckX(xmax);
	 ApplyPDFWeight(nxlo,xmin,fKernX[ObsBin]->GetGridPtr()); // changes node values
	 ApplyPDFWeight(nxup,xmax,fKernX[ObsBin]->GetGridPtr()); // changes node values
      }


      v4d& st = c->SigmaTilde[ObsBin];
      for ( unsigned int is = 0 ; is<events.size() ; is++ ) {
	 double mu = fScenario._m1 * fScaleFac[is];
	 const vector<pair<int,double> >& nmu  = fKernMuS[ObsBin][is]->GetNodeValues(mu);
	 for ( unsigned int m1 = 0 ; m1<nmu.size() ; m1++ ) {
	    v2d& stm1 = st[is][nmu[m1].first];
	    for ( unsigned int p = 0 ; p<events[is].size() ; p++ ) {
	       double wgt = events[is][p]._w * nmu[m1].second / BinSize[ObsBin];
 	       for ( unsigned int x1 = 0 ; x1<nxup.size() ; x1++ ) {
 		  for ( unsigned int x2 = 0 ; x2<nxlo.size() ; x2++ ) {
		     int xmaxbin = nxup[x1].first;
		     int xminbin = nxlo[x2].first;
		     int proc = events[is][p]._p;
  		     HalfMatrixCheck(xminbin,xmaxbin,proc);
		     int ixHM = GetXIndex(ObsBin,xminbin,xmaxbin);
// 		     double ww = nxup[x1].second * nxlo[x2].second * wgt + stm1[ixHM][proc];
// 		     double& sti = stm1[ixHM][proc];
// 		     sti = ww;	// this assignment is time consuming ?!
  		     stm1[ixHM][proc] += nxup[x1].second * nxlo[x2].second * wgt;//nmu[m1].second * wp[p];
//  		     c->SigmaTilde[ObsBin][is][nmu[m1].first][ixHM][proc] += nxup[x1].second * nxlo[x2].second * wgt;//nmu[m1].second * wp[p];
		  }
	       }
	    }
	 }
      }
   }
   else {
      for ( unsigned int is = 0 ; is<events.size() ; is++ ) {
	 FillAllSubprocesses( events[is], scen, is );
      }
   }
}



// ___________________________________________________________________________________________________
void fastNLOCreate::FillAllSubprocesses(const vector<fnloEvent>& events, const fnloScenario& scen, int scalevar)
{
   //! fill a list of subprocesses into the fastNLO table

   if ( (int)events.size() != GetNSubprocesses() ){
      error["FillAllSubprocess"]<<"This table expects "<<GetNSubprocesses()<<" subprocesses, but only "<<events.size()<<" are provided. Exiting."<<endl;
      exit(1);
   }
   for ( unsigned int p = 0 ; p<events.size() ; p++ ) {
      FillOneSubprocess(events[p],scen,scalevar);
   }
}



// ___________________________________________________________________________________________________
void fastNLOCreate::FillOneSubprocess(const fnloEvent& event, const fnloScenario& scen, int scalevar)
{
   fEvent = event;
   fScenario = scen;
   Fill(scalevar);
}



// ___________________________________________________________________________________________________
void fastNLOCreate::Fill(int scalevar)
{
   //!
   //! Fill values, which are stored in 'Event' and 'Scenario' into fastNLO table.
   //!
   //debug["Fill"]<<"Filling subprocess contributions into table."<<endl;

   //GetTheCoeffTable()->Nevt++; // todo: counting of events must be properly implemented
   fStats._nProc++; //keep statistics

   if ( fIsWarmup && scalevar==0 ) UpdateWarmupArrays();
   else FillContribution(scalevar);

   fEvent.ResetButX();
}



// ___________________________________________________________________________________________________
void fastNLOCreate::FillContribution(int scalevar)
{
   //! read information from 'Event' and 'Scenario'
   //! do the interpolation
   //! and fill into the tables.

   if ( fEvent._n > 0 ) SetNumberOfEvents(fEvent._n);

   const int ObsBin = (fScenario._iOB == -1) ? GetBin() : fScenario._iOB;
   if ( ObsBin < 0 ) return;
   if ( ObsBin >= GetNObsBin() ) return;
   fStats._nEvPS++;

   fastNLOCoeffAddBase* c = GetTheCoeffTable();

   // ---- DIS ---- //
   if ( c->GetNPDF() == 1 && fastNLOCoeffAddFlex::CheckCoeffConstants(c,true)) {
      // todo
      FillContributionFlexDIS((fastNLOCoeffAddFlex*)GetTheCoeffTable(),  ObsBin);
      //{error["FillContribution"]<<"Don't know how to fill this table. Exiting."<<endl; exit(1); }
   }
   else if ( c->GetNPDF() == 1 && fastNLOCoeffAddFix::CheckCoeffConstants(c,true)) {
      // todo
      {error["FillContribution"]<<"Don't know how to fill this table. Exiting."<<endl; exit(1); }
   }
   // ---- pp/ppbar ---- //
   else if ( c->GetNPDF() == 2 && fastNLOCoeffAddFlex::CheckCoeffConstants(c,true) )
      FillContributionFlexHHC((fastNLOCoeffAddFlex*)GetTheCoeffTable(),  ObsBin);
   else if ( c->GetNPDF() == 2 && fastNLOCoeffAddFix::CheckCoeffConstants(c,true) )
      FillContributionFixHHC((fastNLOCoeffAddFix*)GetTheCoeffTable(),  ObsBin, scalevar);
   else {
      error["FillContribution"]<<"Don't know how to fill this table. Exiting."<<endl; exit(1);
   }
}



// ___________________________________________________________________________________________________
void fastNLOCreate::FillContributionFixHHC(fastNLOCoeffAddFix* c, int ObsBin, int scalevar)
{
   //! read informatio from 'Event' and 'Scenario'
   //! do the interpolation
   //! and fill into the tables.
   debug["FillContributionFixHHC"]<<endl;

   if ( fEvent._w == 0 ) return; // nothing todo.

   // do interpolation
   double xmin = std::min(fEvent._x1,fEvent._x2);
   double xmax = std::max(fEvent._x1,fEvent._x2);
   vector<pair<int,double> > nxlo = fKernX[ObsBin]->GetNodeValues(xmin);
   vector<pair<int,double> > nxup = fKernX[ObsBin]->GetNodeValues(xmax);
   const double mu = fScenario._m1 * fScaleFac[scalevar];
   const vector<pair<int,double> >& nmu  = fKernMuS[ObsBin][scalevar]->GetNodeValues(mu);

   if ( fApplyPDFReweight ) {
      fKernX[ObsBin]->CheckX(xmin);
      fKernX[ObsBin]->CheckX(xmax);
      ApplyPDFWeight(nxlo,xmin,fKernX[ObsBin]->GetGridPtr());
      ApplyPDFWeight(nxup,xmax,fKernX[ObsBin]->GetGridPtr());
   }


   // fill grid
   if ( CheckWeightIsNan() ) return;
   double wgt = fEvent._w / BinSize[ObsBin];
   for ( unsigned int x1 = 0 ; x1<nxup.size() ; x1++ ) {
      for ( unsigned int x2 = 0 ; x2<nxlo.size() ; x2++ ) {
	 int xmaxbin = nxup[x1].first;
         int xminbin = nxlo[x2].first;
	 int p = fEvent._p;
         HalfMatrixCheck(xminbin,xmaxbin,p);
         int ixHM = GetXIndex(ObsBin,xminbin,xmaxbin);

         for ( unsigned int m1 = 0 ; m1<nmu.size() ; m1++ ) {
	    double w = wgt * nxup[x1].second * nxlo[x2].second * nmu[m1].second / BinSize[ObsBin];
//   	    cout<<"   Fill * : i="<<ObsBin<<" ix="<<ixHM<<", im1="<<nmu[m1].first<<", p="<<p<<", w="<<fEvent._w  * wfnlo
//   		<<",\twnlo="<<fEvent._w<<",\twx="<<nxup[x1].second * nxlo[x2].second<<",\tws="<<nmu[m1].second<<endl;
	    c->SigmaTilde[ObsBin][scalevar][nmu[m1].first][ixHM][p] += w;
         }
      }
   }
}



// ___________________________________________________________________________________________________
void fastNLOCreate::FillContributionFlexHHC(fastNLOCoeffAddFlex* c, int ObsBin)
{
   //! read informatio from 'Event' and 'Scenario'
   //! do the interpolation
   //! and fill into the tables.
   debug["FillContributionFlexHHC"]<<endl;

   if ( fEvent._w == 0 && fEvent._wf==0 && fEvent._wr==0 && fEvent._wrr==0 && fEvent._wff==0 && fEvent._wrf==0 ) return; // nothing todo.

   // do interpolation
   //cout<<"try to interpol. ObsBin="<<ObsBin<<" ,x1="<<fEvent._x1<<", x2="<<fEvent._x2<<", mu1="<<Scenario._m1<<", mu2="<<Scenario._m2<<endl;
   double xmin = std::min(fEvent._x1,fEvent._x2);
   double xmax = std::max(fEvent._x1,fEvent._x2);
   vector<pair<int,double> > nxlo = fKernX[ObsBin]->GetNodeValues(xmin);
   vector<pair<int,double> > nxup = fKernX[ObsBin]->GetNodeValues(xmax);
   vector<pair<int,double> > nmu1 = fKernMu1[ObsBin]->GetNodeValues(fScenario._m1);
   vector<pair<int,double> > nmu2 = fKernMu2[ObsBin]->GetNodeValues(fScenario._m2);


//       cout<<"neu: ObsBin = "<<ObsBin<<", Mu1="<<fScenario._m1<<", Mu2="<<fScenario._m2<<endl;
//       cout<<"     xmin="<<xmin<<"\txmax="<<xmax<<endl;
//       cout<<"     deltamin="<<fKernX[ObsBin]->GetDelta(xmin)<<"\tdeltamax="<<fKernX[ObsBin]->GetDelta(xmax)<<endl;
//       cout<<" ---- x-nodes ----- "<<endl;
//       cout<<"     xlo="<<nxlo[0].second<<"\txlo="<<nxlo[1].second<<"\txlo="<<nxlo[2].second<<"\txlo="<<nxlo[3].second<<endl;
//       cout<<"     xup="<<nxup[0].second<<"\txup="<<nxup[1].second<<"\txup="<<nxup[2].second<<"\txup="<<nxup[3].second<<endl;

//       cout<<" - - - - - x-grid - - - - "<<endl;
//       for ( unsigned int k = 0 ; k<fKernX[ObsBin]->fgrid.size() ; k ++ ) {
//       cout<<"k="<<k<<"\tXNode="<<fKernX[ObsBin]->fgrid[k]<<endl;
//       }
//       cout<<" ---- scalenodes -----"<<endl;
//       cout<<"     mu1="<<nmu1[0].second<<"\tmu1="<<nmu1[1].second<<"\tmu1="<<nmu1[2].second<<"\tmu1="<<nmu1[3].second<<endl;
//       cout<<"     mu2="<<nmu2[0].second<<"\tmu2="<<nmu2[1].second<<"\tmu2="<<nmu2[2].second<<"\tmu2="<<nmu2[3].second<<endl;

//       cout<<" ---- scale1 ----- mu1="<<fScenario._m1<<endl;
//       cout<<"       mu1="<<nmu1[0].second<<"\t  mu1="<<nmu1[1].second<<"\t  mu1="<<nmu1[2].second<<"\t  mu1="<<nmu1[3].second<<endl;
//       cout<<"       node2="<<fKernMu1[ObsBin]->FindLargestPossibleNode(fScenario._m1)<<endl;
//       cout<<"     NscalenodeScale1="<<fKernMu1[ObsBin]->fHgrid.size()<<endl;
//       cout<<"     delta="<<fKernMu1[ObsBin]->GetDelta(fScenario._m1)<<endl;
//       //cout<<"     nscale1="<<nscale1<<endl;
//       cout<<"     HScaleNode[0]="<<fKernMu1[ObsBin]->fHgrid[0]<<", HNode[1]="<<fKernMu1[ObsBin]->fHgrid[1]<<", HNode[2]="<<fKernMu1[ObsBin]->fHgrid[2]<<", HNode[3]="<<fKernMu1[ObsBin]->fHgrid[3]<<", HScaleNode[4]="<<fKernMu1[ObsBin]->fHgrid[4]<<endl;
//       cout<<"     ScaleNode[0]="<<fKernMu1[ObsBin]->fgrid[0]<<", Node[1]="<<fKernMu1[ObsBin]->fgrid[1]<<", Node[2]="<<fKernMu1[ObsBin]->fgrid[2]<<", Node[3]="<<fKernMu1[ObsBin]->fgrid[3]<<", ScaleNode[4]="<<fKernMu1[ObsBin]->fgrid[4]<<endl;

//       cout<<" ---- scale2 ----- mu2="<<fScenario._m2<<endl;
//       cout<<"       mu2="<<nmu2[0].second<<"\t  mu1="<<nmu2[1].second<<"\t  mu1="<<nmu2[2].second<<"\t  mu1="<<nmu2[3].second<<endl;
//       cout<<"       node2="<<fKernMu2[ObsBin]->FindLargestPossibleNode(fScenario._m2)<<endl;
//       cout<<"     NscalenodeScale2="<<fKernMu2[ObsBin]->fHgrid.size()<<endl;
//       cout<<"     delta="<<fKernMu2[ObsBin]->GetDelta(fScenario._m2)<<endl;
//       cout<<"     HScaleNode[0]="<<fKernMu2[ObsBin]->fHgrid[0]<<", HNode[1]="<<fKernMu2[ObsBin]->fHgrid[1]<<", HNode[2]="<<fKernMu2[ObsBin]->fHgrid[2]<<", HNode[3]="<<fKernMu2[ObsBin]->fHgrid[3]<<", HScaleNode[4]="<<fKernMu2[ObsBin]->fHgrid[4]<<endl;
//       cout<<"     ScaleNode[0]="<<fKernMu2[ObsBin]->fgrid[0]<<", Node[1]="<<fKernMu2[ObsBin]->fgrid[1]<<", Node[2]="<<fKernMu2[ObsBin]->fgrid[2]<<", Node[3]="<<fKernMu2[ObsBin]->fgrid[3]<<", ScaleNode[4]="<<fKernMu2[ObsBin]->fgrid[4]<<endl;

   if ( fApplyPDFReweight ) {
      fKernX[ObsBin]->CheckX(xmin);
      fKernX[ObsBin]->CheckX(xmax);
      ApplyPDFWeight(nxlo,xmin,fKernX[ObsBin]->GetGridPtr());
      ApplyPDFWeight(nxup,xmax,fKernX[ObsBin]->GetGridPtr());
   }
   //       cout<<" --  after reweight: --  "<<endl;
   //       cout<<"     n1min="<<nxlo[0].second<<"\ttn1min="<<nxlo[1].second<<"\txlo="<<nxlo[2].second<<"\txlo="<<nxlo[3].second<<endl;
   //       cout<<"     n1max="<<nxup[0].second<<"\ttn1max="<<nxup[1].second<<"\txup="<<nxup[2].second<<"\txup="<<nxup[3].second<<endl;
   //       cout<<"     mu1="<<nmu1[0].second<<"\tmu1="<<nmu1[1].second<<"\tmu1="<<nmu1[2].second<<"\tmu1="<<nmu1[3].second<<endl;
   //       cout<<"     mu2="<<nmu2[0].second<<"\tmu2="<<nmu2[1].second<<"\tmu2="<<nmu2[2].second<<"\tmu2="<<nmu2[3].second<<endl;
   //       cout<<"  0-nodes: mi1="<< nmu1[0].first<<",  mi2="<< nmu2[0].first<<", xup1="<<nxup[0].first<<", xdn1="<<nxlo[0].first<<endl;


   // fill grid
   if ( CheckWeightIsNan() ) return;
   for ( unsigned int x1 = 0 ; x1<nxup.size() ; x1++ ) {
      for ( unsigned int x2 = 0 ; x2<nxlo.size() ; x2++ ) {
         int xmaxbin = nxup[x1].first;
         int xminbin = nxlo[x2].first;
         int p = fEvent._p;
         HalfMatrixCheck(xminbin,xmaxbin,p);
         int ixHM = GetXIndex(ObsBin,xminbin,xmaxbin);

         for ( unsigned int m1 = 0 ; m1<nmu1.size() ; m1++ ) {
            for ( unsigned int mu2 = 0 ; mu2<nmu2.size() ; mu2++ ) {
               double wfnlo = nxup[x1].second * nxlo[x2].second * nmu1[m1].second * nmu2[mu2].second / BinSize[ObsBin];
               if ( isnan(wfnlo) ) {
                  error[""]<<"wfnlo is a nan."<<endl;
                  fKernX[ObsBin]->PrintGrid();
                  fKernMu1[ObsBin]->PrintGrid();
                  fKernMu2[ObsBin]->PrintGrid();
                  cout<<"ix1="<<x1<<", ix2="<<x2<<", im1="<<m1<<", im2="<<mu2<<endl;
                  cout<<"x1="<<nxlo[x1].second<<", x1="<<x1<<", xval="<<xmin<<endl;
                  cout<<"x2="<<nxup[x2].second<<", x2="<<x2<<", xval="<<xmax<<endl;
                  cout<<"m1="<< nmu1[m1].second<<", m1="<<m1<<", mu1val="<<fScenario._m1<<endl;
                  cout<<"m2="<<nmu2[mu2].second<<", m2="<<mu2<<", mu2val="<<fScenario._m2<<endl;
                  exit(1);
               }
	       // 	       cout<<"ObsBin="<<ObsBin<<", ixHM="<<ixHM<<", m1="<<nmu1[m1].first<<", m2="<< nmu2[mu2].first<<", p="<<p
	       // 		   <<" ,i(x="<<x1<<",x2="<<x2<<",m1="<<m1<<",m2="<<mu2<<") [xlo="<<xmin<<",xup="<<xmax<<",m1="<<fScenario._m1<<",m2="<<fScenario._m2<<"]"<<endl;
	       // 	       cout<<" ggg-n-  : O="<<ObsBin<<", ix="<<ixHM<<", im1="<<nmu1[m1].first<<", im2="<<nmu2[mu2].first<<", p="<<p<<", wfnlo="<<wfnlo<<", wxu="<<nxup[x1].second<<", wxd="<<nxlo[x2].second<<", wm1="<<nmu1[m1].second<<", wm2="<<nmu2[mu2].second<<endl;
               if ( fEvent._w  != 0 ) {
		  //cout<<"   Fill * : ix="<<ixHM<<", im1="<<nmu1[m1].first<<", im2="<<nmu2[mu2].first<<", p="<<p<<", w="<<fEvent._w  * wfnlo<<endl;
		  c->SigmaTildeMuIndep[ObsBin][ixHM][nmu1[m1].first][nmu2[mu2].first][p]  += fEvent._w  * wfnlo;
               }
               if ( fEvent._wf != 0 ){
		  //cout<<"   Fill F : ix="<<ixHM<<", im1="<<nmu1[m1].first<<", im2="<<nmu2[mu2].first<<", p="<<p<<", w="<<fEvent._wf  * wfnlo<<endl;
		  c->SigmaTildeMuFDep [ObsBin][ixHM][nmu1[m1].first][nmu2[mu2].first][p]  += fEvent._wf * wfnlo;
               }
               if ( fEvent._wr != 0 ) {
		  //cout<<"   Fill R : ix="<<ixHM<<", im1="<<nmu1[m1].first<<", im2="<<nmu2[mu2].first<<", p="<<p<<", w="<<fEvent._wr  * wfnlo<<endl;
                  c->SigmaTildeMuRDep [ObsBin][ixHM][nmu1[m1].first][nmu2[mu2].first][p]  += fEvent._wr * wfnlo;
               }
               if ( fEvent._wrr != 0 ) {
                  c->SigmaTildeMuRRDep [ObsBin][ixHM][nmu1[m1].first][nmu2[mu2].first][p]  += fEvent._wrr * wfnlo;
               }
               if ( fEvent._wff != 0 ) {
                  c->SigmaTildeMuFFDep [ObsBin][ixHM][nmu1[m1].first][nmu2[mu2].first][p]  += fEvent._wff * wfnlo;
               }
               if ( fEvent._wrf != 0 ) {
                  c->SigmaTildeMuRFDep [ObsBin][ixHM][nmu1[m1].first][nmu2[mu2].first][p]  += fEvent._wrf * wfnlo;
               }
            }
         }
      }
   }
}



// ___________________________________________________________________________________________________
void fastNLOCreate::FillContributionFlexDIS(fastNLOCoeffAddFlex* c, int ObsBin)
{
   //! read information from 'Event' and 'Scenario'
   //! do the interpolation
   //! and fill into the tables.
   debug["FillContributionFlexHHC"]<<endl;
   
   if ( fEvent._w == 0 && fEvent._wf==0 && fEvent._wr==0 ) return; // nothing todo.

   // do interpolation
   //cout<<"try to interpol. ObsBin="<<ObsBin<<" ,x1="<<fEvent._x1<<", x2="<<fEvent._x2<<", mu1="<<Scenario._m1<<", mu2="<<Scenario._m2<<endl;

   // todo, just: 'x'
   double x = fEvent._x1;
   vector<pair<int,double> > nx = fKernX[ObsBin]->GetNodeValues(x);
   vector<pair<int,double> > nmu1 = fKernMu1[ObsBin]->GetNodeValues(fScenario._m1);
   vector<pair<int,double> > nmu2 = fKernMu2[ObsBin]->GetNodeValues(fScenario._m2);


   if ( fApplyPDFReweight ) {
      fKernX[ObsBin]->CheckX(x);
      ApplyPDFWeight(nx,x,fKernX[ObsBin]->GetGridPtr());
   }

   // fill grid
   if ( CheckWeightIsNan() ) return;
   for ( unsigned int ix = 0 ; ix<nx.size() ; ix++ ) {
      int p = fEvent._p;
      int xIdx = nx[ix].first;
      //HalfMatrixCheck(xminbin,xmaxbin,p);
      //int ixHM = GetXIndex(ObsBin,xminbin,xmaxbin);

      for ( unsigned int m1 = 0 ; m1<nmu1.size() ; m1++ ) {
	 for ( unsigned int mu2 = 0 ; mu2<nmu2.size() ; mu2++ ) {
	    double wfnlo = nx[ix].second * nmu1[m1].second * nmu2[mu2].second / BinSize[ObsBin];
	    if ( isnan(wfnlo) ) {
	       error[""]<<"wfnlo is a nan."<<endl;
	       fKernX[ObsBin]->PrintGrid();
	       fKernMu1[ObsBin]->PrintGrid();
	       fKernMu2[ObsBin]->PrintGrid();
	       cout<<"ix1="<<ix<<", im1="<<m1<<", im2="<<mu2<<endl;
	       cout<<"x1="<<nx[ix].second<<", ix="<<ix<<", xval="<<x<<endl;
	       cout<<"m1="<< nmu1[m1].second<<", m1="<<m1<<", mu1val="<<fScenario._m1<<endl;
	       cout<<"m2="<<nmu2[mu2].second<<", m2="<<mu2<<", mu2val="<<fScenario._m2<<endl;
	       exit(1);
	    }
	    if ( fEvent._w  != 0 ) {
	       //                 cout<<"   Fill * : ix="<<ixHM<<", im1="<<nmu1[m1].first<<", im2="<<nmu2[mu2].first<<", p="<<p<<", w="<<fEvent._w  * wfnlo<<endl;
	       c->SigmaTildeMuIndep[ObsBin][xIdx][nmu1[m1].first][nmu2[mu2].first][p]  += fEvent._w  * wfnlo;
	    }
	    if ( fEvent._wf != 0 ){
	       //                 cout<<"   Fill F : ix="<<ixHM<<", im1="<<nmu1[m1].first<<", im2="<<nmu2[mu2].first<<", p="<<p<<", w="<<fEvent._wf  * wfnlo<<endl;
	       c->SigmaTildeMuFDep [ObsBin][xIdx][nmu1[m1].first][nmu2[mu2].first][p]  += fEvent._wf * wfnlo;
	    }
	    if ( fEvent._wr != 0 ) {
	       //                 cout<<"   Fill R : ix="<<ixHM<<", im1="<<nmu1[m1].first<<", im2="<<nmu2[mu2].first<<", p="<<p<<", w="<<fEvent._wr  * wfnlo<<endl;
	       c->SigmaTildeMuRDep [ObsBin][xIdx][nmu1[m1].first][nmu2[mu2].first][p]  += fEvent._wr * wfnlo;
	    }
	 }
      }
   }
}


// ___________________________________________________________________________________________________
inline void fastNLOCreate::HalfMatrixCheck(int& xminbin, int& xmaxbin, int& subproc) const {
   //! check if half-matrix notation
   //! if half-matrix notation, and xmin-node is larger than xmax-node
   //! exchange suprocesses according to fSymProc and adjust x-nodes.
   //!
   if ( xminbin > xmaxbin  ) {
      if ( GetTheCoeffTable()->GetNPDFDim() == 1 ) { // half-matrix notation (otherwise nothing todo)
//          if ( (int)fSymProc.size() != GetTheCoeffTable()->GetNSubproc() )
//             error["HalfMatrixCheck"]<<"Necessary array with symmetric processes for half-matrix notation not initialized."<<endl;

         //cout<<"exchange supbrpc. xminbin="<<xminbin<<", xmaxbin="<<xmaxbin<<", p="<<subproc<<", pAsym="<<fSymProc[subproc]<<endl;
         int di = xminbin - xmaxbin;
         xmaxbin += di;        // modify indicees
         xminbin -= di;
         subproc = fSymProc[subproc];           // exchange asymmetric process
      }
   }
}


// ___________________________________________________________________________________________________
bool fastNLOCreate::CheckWeightIsNan() {
   //! check if weights contain isnan
   if ( isnan(fEvent._w)) {
      error["CheckWeightIsNan"]<<"(Scale-independent) weight is 'nan'"<<endl;
      return true;
   }
   if ( isnan(fEvent._wf)) {
      error["CheckWeightIsNan"]<<"Factorization scale dependent weight is 'nan'"<<endl;
      return true;
   }
   if ( isnan(fEvent._wr)) {
      error["CheckWeightIsNan"]<<"Renormalization scale dependent weight is 'nan'"<<endl;
      return true;
   }
   return false;
}


// ___________________________________________________________________________________________________
inline int fastNLOCreate::GetXIndex(const int& ObsBin,const int& x1bin,const int& x2bin) const {
   //! get index if 1 or two hadrons are involved
   //switch (GetTheCoeffTable()->GetNPDFDim() ) {
   // faster code ?
   //       if ( ((fastNLOCoeffAddBase*)fCoeff[0])->NPDFDim == 1 ) return x1bin + (x2bin*(x2bin+1)/2);
   //       else if ( ((fastNLOCoeffAddBase*)fCoeff[0])->NPDFDim == 0 ) return x1bin; // linear
   //       else return x1bin + x2bin * GetTheCoeffTable()->GetNxtot1(ObsBin); // full matrix

   //    switch (GetTheCoeffTable()->GetNPDFDim() ) {
   switch (((fastNLOCoeffAddBase*)fCoeff[0])->NPDFDim) {
   case 1:
      return x1bin + (x2bin*(x2bin+1)/2);    // half matrix
   case 0:
      return x1bin; // linear
   case 2:
      return x1bin + x2bin * GetTheCoeffTable()->GetNxtot1(ObsBin); // full matrix
   default:
      return -1; // this will cause a crash :)
   }
};



// ___________________________________________________________________________________________________
int fastNLOCreate::GetNxmax(const vector<double>* xGrid1, const vector<double>* xGrid2 ){
   switch (GetTheCoeffTable()->GetNPDFDim() ) {
   case 0:
      return xGrid1->size();
   case 1:
      if ( !xGrid2 ) error["GetNxmax"]<<"Error. Second x-grid must be specified."<<endl;
      if (xGrid1->size() != xGrid2->size() )error["GetNxmax"]<<"Grid sizes in half-matrix notation must have equal size."<<endl;
      return ((int)pow((double)xGrid1->size(),2)+xGrid1->size())/2;
   case 2:
      if ( !xGrid2 ) error["GetNxmax"]<<"Error. Second x-grid must be specified."<<endl;
      return xGrid1->size()*xGrid2->size();
   default:
      return 0;
   }
};


// ___________________________________________________________________________________________________
inline void fastNLOCreate::ApplyPDFWeight(vector<pair<int,double> >& nodes, const double x, const vector<double>* grid ) const{
//    double pdfwgtmax = PDFwgt(xmax);
//    for( int i1 = 0; i1 < 4; i1++) {
//       if ((nxmaxf-1+i1) >= 0 && (nxmaxf-1+i1) < Nxtot1[ObsBin] ) {
//          cefmax[i1] *= pdfwgtmax/PDFwgt(XNode1[ObsBin][nxmaxf-1+i1]);
//       }
//    }
   double wgtx = CalcPDFReweight(x);
   for( unsigned int in = 0; in < nodes.size(); in++) {
      double wgtn = CalcPDFReweight(grid->at(nodes[in].first));
      if ( wgtn==0 ) {error["ApplyPDFWeight"]<<"Cannot divide by 0."<<endl; exit(1);}
      //cout<<"in="<<in<<" wgtx="<<wgtx<<", x="<<x<<", wgtn="<<wgtn<<", nod="<<grid->at(nodes[in].first)<<endl;
      nodes[in].second *= wgtx/wgtn;
   }
}


// ___________________________________________________________________________________________________
inline double fastNLOCreate::CalcPDFReweight(double x) const {
   if ( x<=0 ) { error["CalcPDFReweight"]<<"Cannot calculate sqrt of negative numbers or divide by zero. x="<<x<<endl; exit(1);}
   double w=(1.-0.99*x)/sqrt(x);
   return w*w*w;
}


// ___________________________________________________________________________________________________
void fastNLOCreate::MultiplyCoefficientsByBinWidth() {
   //! Multiply all coefficients by binwidth

   if ( fIsFlexibleScale ) {
      fastNLOCoeffAddFlex* c = (fastNLOCoeffAddFlex*)GetTheCoeffTable();
      for (int i=0; i<GetNObsBin(); i++) {
         int nxmax = c->GetNxmax(i);
         for (unsigned int jS1=0; jS1<c->GetNScaleNode1(i); jS1++) {
            for (unsigned int kS2=0; kS2<c->GetNScaleNode2(i); kS2++) {
               for (int x=0; x<nxmax; x++) {
                  for (int n=0; n<c->GetNSubproc(); n++) {
                     c->SigmaTildeMuIndep[i][x][jS1][kS2][n] *= BinSize[i];
                     if ( c->GetNScaleDep() >= 5 ) {
                        c->SigmaTildeMuFDep [i][x][jS1][kS2][n] *= BinSize[i];
                        c->SigmaTildeMuRDep [i][x][jS1][kS2][n] *= BinSize[i];
                        if ( c->GetNScaleDep() >= 6 ) {
                           c->SigmaTildeMuRRDep [i][x][jS1][kS2][n] *= BinSize[i];
                           c->SigmaTildeMuFFDep [i][x][jS1][kS2][n] *= BinSize[i];
                           c->SigmaTildeMuRFDep [i][x][jS1][kS2][n] *= BinSize[i];
                        }
                     }
                  }
               }
            }
         }
      }
   }
   else {
      fastNLOCoeffAddFix* c = (fastNLOCoeffAddFix*)GetTheCoeffTable();
      for (int i=0; i<GetNObsBin(); i++) {
         for ( unsigned int s=0 ; s<c->SigmaTilde[i].size() ; s++ ) {
            for ( unsigned int x=0 ; x<c->SigmaTilde[i][s].size() ; x++ ) {
               for ( unsigned int l=0 ; l<c->SigmaTilde[i][s][x].size() ; l++ ) {
                  for ( unsigned int m=0 ; m<c->SigmaTilde[i][s][x][m].size() ; m++ ) {
                     c->SigmaTilde[i][s][x][l][m] *= BinSize[i];
                  }
               }
            }
         }
      }
   }
}


// ___________________________________________________________________________________________________
void fastNLOCreate::DivideCoefficientsByBinWidth() {
//! Divide all coefficients by binwidth
   if ( fIsFlexibleScale ) {
      fastNLOCoeffAddFlex* c = (fastNLOCoeffAddFlex*)GetTheCoeffTable();
      for (int i=0; i<GetNObsBin(); i++) {
         int nxmax = c->GetNxmax(i);
         for (unsigned int jS1=0; jS1<c->GetNScaleNode1(i); jS1++) {
            for (unsigned int kS2=0; kS2<c->GetNScaleNode2(i); kS2++) {
               for (int x=0; x<nxmax; x++) {
                  for (int n=0; n<c->GetNSubproc(); n++) {
                     c->SigmaTildeMuIndep[i][x][jS1][kS2][n] /= BinSize[i];
                     if ( c->GetNScaleDep() >= 5 ) {
                        c->SigmaTildeMuFDep [i][x][jS1][kS2][n] /= BinSize[i];
                        c->SigmaTildeMuRDep [i][x][jS1][kS2][n] /= BinSize[i];
                        if ( c->GetNScaleDep() >= 6 ) {
                           c->SigmaTildeMuRRDep [i][x][jS1][kS2][n] /= BinSize[i];
                           c->SigmaTildeMuFFDep [i][x][jS1][kS2][n] /= BinSize[i];
                           c->SigmaTildeMuRFDep [i][x][jS1][kS2][n] /= BinSize[i];
                        }
                     }
                  }
               }
            }
         }
      }
   }
   else {
      fastNLOCoeffAddFix* c = (fastNLOCoeffAddFix*)GetTheCoeffTable();
      for (int i=0; i<GetNObsBin(); i++) {
         for ( unsigned int s=0 ; s<c->SigmaTilde[i].size() ; s++ ) {
            for ( unsigned int x=0 ; x<c->SigmaTilde[i][s].size() ; x++ ) {
               for ( unsigned int l=0 ; l<c->SigmaTilde[i][s][x].size() ; l++ ) {
                  for ( unsigned int m=0 ; m<c->SigmaTilde[i][s][x][m].size() ; m++ ) {
                     c->SigmaTilde[i][s][x][l][m] /= BinSize[i];
                  }
               }
            }
         }
      }
   }
}


// ___________________________________________________________________________________________________
void fastNLOCreate::MultiplyCoefficientsByConstant(double coef) {
//! Divide all coefficients by binwidth
   if ( fIsFlexibleScale ) {
      fastNLOCoeffAddFlex* c = (fastNLOCoeffAddFlex*)GetTheCoeffTable();
      for (int i=0; i<GetNObsBin(); i++) {
         int nxmax = c->GetNxmax(i);
         for (unsigned int jS1=0; jS1<c->GetNScaleNode1(i); jS1++) {
            for (unsigned int kS2=0; kS2<c->GetNScaleNode2(i); kS2++) {
               for (int x=0; x<nxmax; x++) {
                  for (int n=0; n<c->GetNSubproc(); n++) {
                     c->SigmaTildeMuIndep[i][x][jS1][kS2][n] *= coef;
                     if ( c->GetNScaleDep() >= 5 ) {
                        c->SigmaTildeMuFDep [i][x][jS1][kS2][n] *= coef;
                        c->SigmaTildeMuRDep [i][x][jS1][kS2][n] *= coef;
                        if ( c->GetNScaleDep() >= 6 ) {
                           c->SigmaTildeMuRRDep [i][x][jS1][kS2][n] *= coef;
                           c->SigmaTildeMuFFDep [i][x][jS1][kS2][n] *= coef;
                           c->SigmaTildeMuRFDep [i][x][jS1][kS2][n] *= coef;
                        }
                     }
                  }
               }
            }
         }
      }
   }
   else {
      fastNLOCoeffAddFix* c = (fastNLOCoeffAddFix*)GetTheCoeffTable();
      for (int i=0; i<GetNObsBin(); i++) {
         for ( unsigned int s=0 ; s<c->SigmaTilde[i].size() ; s++ ) {
            for ( unsigned int x=0 ; x<c->SigmaTilde[i][s].size() ; x++ ) {
               for ( unsigned int l=0 ; l<c->SigmaTilde[i][s][x].size() ; l++ ) {
                  for ( unsigned int m=0 ; m<c->SigmaTilde[i][s][x][m].size() ; m++ ) {
                     c->SigmaTilde[i][s][x][l][m] *= coef;
                  }
               }
            }
         }
      }
   }
}


// ___________________________________________________________________________________________________
void fastNLOCreate::UpdateWarmupArrays(){
   //! Update the warmup-arrays fWMu1, fWx und fWMu2
   if ( fWx.empty() ) InitWarmupArrays();

   const int ObsBin = GetBin();
   debug["UpdateWarmupArrays"]<<"ObsBin="<<ObsBin<<"\tmu1="<<fScenario._m1<<"\tmu2="<<fScenario._m2<<"\tx1="<<fEvent._x1<<"\tx2="<<fEvent._x2<<endl;
   if ( ObsBin >= 0 ) {
      fWMu1[ObsBin].first       = std::min(fScenario._m1,fWMu1[ObsBin].first) ;
      fWMu1[ObsBin].second      = std::max(fScenario._m1,fWMu1[ObsBin].second) ;
      if ( GetTheCoeffTable()->IPDFdef1 == 3 ) {
         fWx[ObsBin].first      = std::min(std::min(fEvent._x1,fEvent._x2),fWx[ObsBin].first) ;
         fWx[ObsBin].second     = std::max(std::max(fEvent._x1,fEvent._x2),fWx[ObsBin].second) ;
      } else if (GetTheCoeffTable()->IPDFdef1 == 2 ){
         fWx[ObsBin].first      = std::min(fEvent._x1,fWx[ObsBin].first) ;
         fWx[ObsBin].second     = std::max(fEvent._x1,fWx[ObsBin].second) ;
      } else
         error["UpdateWarmupArrays"]<<"nothing reasonable implemented yet."<<endl;
      if ( fIsFlexibleScale ) {
         fWMu2[ObsBin].first    = std::min(fScenario._m2,fWMu2[ObsBin].first) ;
         fWMu2[ObsBin].second   = std::max(fScenario._m2,fWMu2[ObsBin].second) ;
      }
   }
}


// ___________________________________________________________________________________________________
void fastNLOCreate::InitWarmupArrays(){
   debug["InitWarmupArrays"]<<endl;
   //! initialize arrays to store and determined warm-up values
   //! initialize with reasonable values
   fWMu1.resize(GetNObsBin());
   fWMu2.resize(GetNObsBin());
   fWx.resize(GetNObsBin());
   for ( int i = 0 ; i < GetNObsBin() ; i ++ ) {
      fWMu1[i].first    = 10e10;
      fWMu1[i].second   = -10e10;
      fWMu2[i].first    = 10e10;
      fWMu2[i].second   = -10e10;
      fWx[i].first      = 10e10;
      fWx[i].second     = -10e10;
   }
}


// ___________________________________________________________________________________________________
void fastNLOCreate::WriteTable() {
   //if ( GetTheCoeffTable()->GetNevt(0,0) <= 0 ) {
   if ( GetTheCoeffTable()->Nevt <= 0 ) {
      warn["WriteTable"]<<"Number of events seems to be not filled. Please use SetNumberOfEvents(int) before writing table."<<endl;
   }
   fStats.PrintStats();
   if ( fIsWarmup ) {
      info["WriteTable"]<<"Writing warmup table instead of coefficient table."<<endl;
      WriteWarmupTable();
   }
   else {
      if ( ffilename == "" ) { 
	 error["WriteTable"]<<"No filename given."<<endl; 
	 exit(1);
      }
      // Number of events must be counted correctly.
      // I.e. the counting should be performed by the generator.
      // ->Divide by BinWidth
      fastNLOTable::WriteTable();
      // ->Multiply by BinWidth
   }
}


// ___________________________________________________________________________________________________
void fastNLOCreate::WriteTable(string filename) {
   SetFilename(filename);
   WriteTable();
}


// ___________________________________________________________________________________________________
void fastNLOCreate::WriteWarmupTable(){
   string tempfn = ffilename;
   string warmupfile = GetWarmupTableFilename();
   info["WriteWarmupTable"]<<"Writing warmup table to: "<<warmupfile<<endl;
   SetFilename(warmupfile);
   // open stream;
   OpenFileRewrite();
   // write to disk
   OutWarmup(warmupfile);
   // close file
   CloseStream();
   SetFilename(tempfn);
}



// ___________________________________________________________________________________________________
void fastNLOCreate::PrintWarmupValues(){
   OutWarmup("");
}


// ___________________________________________________________________________________________________
void fastNLOCreate::OutWarmup(string file){
   if ( fWx.empty() ) {
      warn["OutWarmup"]<<"Warmup arrays not initialized. Did you forgot to fill values?"<<endl;
//       warn["OutWarmup"]<<"  Continuting, but writing unreasonalby large/small values as warmup values..."<<endl;
//       InitWarmupArrays();
      error["OutWarmup"]<<" Do not write out unreasonable warmup table. Exiting."<<endl; exit(1);
   }
   std::ostream& sout = file == "" ?
      std::cout : (*ofilestream);

   sout<<"! This is a automatically generated file by fastNLO and holds the values of the warmup run. "<<endl;
   sout<<"! The values are valid for the scenario "<<GetScenName() << endl;
   sout<<"! and if calculated with the steerfile: "<< fSteerfile <<endl;
   sout<<"! but only if no serious changes have been performed since its creation."<<endl;
   sout<<"! "<<endl;
   sout<<"! Delete this file, if you want fastNLO to calculate a new one."<<endl;
   sout<<"! "<<endl;
   sout<<"! This file has been calculated using "<<GetTheCoeffTable()->Nevt<<" contributions."<<endl;
   sout<<"!   ( Mind: contributions != events. And contributions are not necessarily in phase space region."<<endl;
   sout<<"! Please check by eye for reasonability of the values."<<endl;
   sout<<" " <<endl;

   // write variables of warmup run
   sout<<"Warmup.OrderInAlphasOfCalculation\t"<<  fIOrd <<endl;
   sout<<"Warmup.ScaleDescriptionScale1\t\t\""<< GetTheCoeffTable()->ScaleDescript[0][0]<<"\""<<endl;
   if ( fIsFlexibleScale ) 
      sout<<"Warmup.ScaleDescriptionScale2\t\t\""<< GetTheCoeffTable()->ScaleDescript[0][1]<<"\"" <<endl;
   sout<<"Warmup.DifferentialDimension\t\t"<< NDim <<endl;
   sout<<"Warmup.DimensionLabels {\n  ";
   for ( int i = 0 ; i < NDim; i ++ ) sout<<"\""<<DimLabel[i]<<"\"  ";
   sout<<"\n} "<<endl;
   
   sout<<"Warmup.DimensionIsDifferential {\n  ";
   for ( int i = 0 ; i < NDim; i ++ ) sout<<"\""<<IDiffBin[i]<<"\"  ";
   sout<<"\n} "<<endl;
   sout<<endl;

   // write readable table
   char buf[4000];
   sout<<"WarmupValues {{"<<endl;
   //sout<<"   xmin      xmax     mu1min    m1max     mu2min    mu2max "<<endl;
   if ( fIsFlexibleScale ) {
      // table header
      sprintf(buf,"   %9s  %9s  %14s  %14s  %14s  %14s",
              "x_min","x_max",
              GetWarmupHeader(0,"min").c_str(),
              GetWarmupHeader(0,"max").c_str(),
              GetWarmupHeader(1,"min").c_str(),
              GetWarmupHeader(1,"max").c_str());
      sout<<buf; 
      if ( NDim == 1 )  
	 sprintf(buf,"  %9s_Lo  %9s_Up",DimLabel[0].c_str() ,DimLabel[0].c_str()); 
      else if ( NDim == 2 )
	 sprintf(buf,"  %9s_Lo  %9s_Up  %9s_Lo  %9s_Up",DimLabel[1].c_str() ,DimLabel[1].c_str(), DimLabel[0].c_str() ,DimLabel[0].c_str()); 
      else if ( NDim == 3 )
	 sprintf(buf,"  %9s_Lo  %9s_Up  %9s_Lo  %9s_Up  %9s_Lo  %9s_Up",
		 DimLabel[2].c_str() ,DimLabel[2].c_str(), DimLabel[1].c_str() ,DimLabel[1].c_str(), DimLabel[0].c_str() ,DimLabel[0].c_str()); 
      sout<<buf; 
      sprintf(buf,"  %12s","BinWidth"); 
      sout<<buf<<endl;
      
      // table values
      for ( int i = 0 ; i < GetNObsBin() ; i ++ ) {
         if ( fWx[i].first < 1.e-4 ) {
            warn["OutWarmup"]<<"The xmin value in bin "<<i<<" seems to be unreasonably low (xmin="<<fWx[i].first<<"). Taking xmin=1.e-4 instead."<<endl;
            fWx[i].first=1.e-4;
         }
         sprintf(buf,"   %9.2e  %9.2e  %14.2f  %14.2f  %14.3f  %14.3f",
                 fWx[i].first,fWx[i].second,fWMu1[i].first,fWMu1[i].second,fWMu2[i].first,fWMu2[i].second);
	 sout<<buf; 
	 if ( NDim == 1 )  
	    sprintf(buf,"  %12.3f  %12.3f",Bin[i][0].first , Bin[i][0].second); 
	 else if ( NDim == 2 )
	    sprintf(buf,"  %12.3f  %12.3f  %12.3f  %12.3f",Bin[i][1].first , Bin[i][1].second, Bin[i][0].first , Bin[i][0].second); 
	 else if ( NDim == 3 )
	    sprintf(buf,"  %12.3f  %12.3f  %12.3f  %12.3f  %12.3f  %12.3f",
		    Bin[i][2].first , Bin[i][1].second,Bin[i][2].first , Bin[i][1].second, Bin[i][0].first , Bin[i][0].second);
	 sout<<buf; 
	 sprintf(buf,"  %12.3f",BinSize[i]);
         sout<<buf<<endl;
      }
   }
   else {
      // is ScaleDescript available?
      if ( GetTheCoeffTable()->ScaleDescript[0].empty() ){ error["OutWarmup"]<<"Scale description is empty. but needed. Probably this has to be implemented."<<endl; exit(1);};
      // table header
      sprintf(buf,"   %9s  %9s  %16s  %16s",
              "x_min","x_max",
              GetWarmupHeader(0,"min").c_str(),
              GetWarmupHeader(0,"max").c_str() );
      sout<<buf;
      if ( NDim == 1 )  
	 sprintf(buf,"  %9s_Lo  %9s_Up",DimLabel[0].c_str() ,DimLabel[0].c_str()); 
      else if ( NDim == 2 )
	 sprintf(buf,"  %9s_Lo  %9s_Up  %9s_Lo  %9s_Up",DimLabel[1].c_str() ,DimLabel[1].c_str(), DimLabel[0].c_str() ,DimLabel[0].c_str()); 
      else if ( NDim == 3 )
	 sprintf(buf,"  %9s_Lo  %9s_Up  %9s_Lo  %9s_Up  %9s_Lo  %9s_Up",
		 DimLabel[2].c_str() ,DimLabel[2].c_str(), DimLabel[1].c_str() ,DimLabel[1].c_str(), DimLabel[0].c_str() ,DimLabel[0].c_str()); 
      sout<<buf; 
      sprintf(buf,"  %12s","BinWidth"); 
      sout<<buf<<endl;

      // table values
      for ( int i = 0 ; i < GetNObsBin() ; i ++ ) {
         if ( fWx[i].first < 1.e-4 ) {
            warn["OutWarmup"]<<"The xmin value in bin "<<i<<" seems to be unreasonably low (xmin="<<fWx[i].first<<"). Taking xmin=1.e-4 instead."<<endl;
            fWx[i].first=1.e-4;
         }
         sprintf(buf,"   %9.2e  %9.2e  %16.2f  %16.2f",
                 fWx[i].first, fWx[i].second, fWMu1[i].first, fWMu1[i].second);
 	 sout<<buf; 
	 if ( NDim == 1 )  
	    sprintf(buf,"  %12.3f  %12.3f",Bin[i][0].first, Bin[i][0].second); 
	 else if ( NDim == 2 )
	    sprintf(buf,"  %12.3f  %12.3f  %12.3f  %12.3f", Bin[i][1].first, Bin[i][1].second, Bin[i][0].first, Bin[i][0].second); 
	 else if ( NDim == 3 )
	    sprintf(buf,"  %12.3f  %12.3f  %12.3f  %12.3f  %12.3f  %12.3f",
		    Bin[i][2].first , Bin[i][2].second, Bin[i][1].first, Bin[i][1].second, Bin[i][0].first, Bin[i][0].second);
	 sout<<buf; 
	 sprintf(buf,"  %12.3f",BinSize[i]);
	 sout<<buf<<endl;
      }
   }
   sout<<"}}"<<endl;
   //*(new ofstream(filename.Data()));
}


// ___________________________________________________________________________________________________
string fastNLOCreate::GetWarmupHeader(int iScale, string minmax ){
   string Descript = GetTheCoeffTable()->ScaleDescript[0][iScale];
   // replace all 'spaces' with _
   replace(Descript.begin(), Descript.end(), ' ', '_');
   // add 'minmax'
   string ret = "";
   ret += Descript;
   ret += "_";
   ret += minmax;
   return ret;
}



// ___________________________________________________________________________________________________
string fastNLOCreate::GetWarmupTableFilename(){
   string ret = fSteerfile;
   size_t pos = ret.find(".str");
   if(pos != std::string::npos) ret.erase(pos,4);
   pos = ret.find(".steer");
   if(pos != std::string::npos) ret.erase(pos,6);
   ret += "_";
   ret += GetScenName();
   ret += "_warmup.txt";
   return ret;
}


// ___________________________________________________________________________________________________
void  fastNLOCreate::InitGrids() {
   debug["InitGrids"]<<endl;
   if ( fKernX.empty() ) error["InitGrids"]<<"Interpolation kernels must be initialized before calling this function."<<endl;
   if ( fIsFlexibleScale ) {
      fastNLOCoeffAddFlex* c = (fastNLOCoeffAddFlex*)GetTheCoeffTable();
      if ( (c->GetNPDF()==2 && c->GetNPDFDim() == 1) || (c->GetNPDF()==1)   ) {;} // ok!
      else {
	 error["InitGrids"]<<"Only half-matrix or DIS implemented."<<endl; exit(1);
      }
      c->ScaleNode1.resize(GetNObsBin());
      c->ScaleNode2.resize(GetNObsBin());
      c->XNode1.resize(GetNObsBin());
      vector<vector<vector<vector<vector<double> > > > > stype(GetNObsBin());
      for ( int i = 0 ; i < GetNObsBin() ; i ++ ) {
         c->ScaleNode1[i] = fKernMu1[i]->GetGrid();
         c->ScaleNode2[i] = fKernMu2[i]->GetGrid();
         c->XNode1[i]     = fKernX[i]->GetGrid();

         // SigmaTilde [NObsBins] ['n' x-nodes] [n s1-Nodes] [n s2-Nodes] [nsubproc]
         int nxmax = GetNxmax (fKernX[i]->GetGridPtr() , fKernX[i]->GetGridPtr() );
         stype[i].resize(nxmax);
         for ( unsigned int x = 0 ; x<stype[i].size() ; x++ ) {
            stype[i][x].resize(c->ScaleNode1[i].size());
            for ( unsigned int m1 = 0 ; m1<stype[i][x].size() ; m1++ ) {
               stype[i][x][m1].resize(c->ScaleNode2[i].size());
               for ( unsigned int mu2 = 0 ; mu2<stype[i][x][m1].size() ; mu2++ ) {
                  stype[i][x][m1][mu2].resize(GetNSubprocesses());
               }
            }
         }
      }
      c->SigmaTildeMuIndep      = stype;
      c->SigmaTildeMuRDep       = stype;
      c->SigmaTildeMuFDep       = stype;

      c->SigmaTildeMuRRDep      = stype;
      c->SigmaTildeMuFFDep      = stype;
      c->SigmaTildeMuRFDep      = stype;
      //c->SigmaTildeMuIndep(GetNObsBin());///
   }

   else {
      fastNLOCoeffAddFix* c = (fastNLOCoeffAddFix*)GetTheCoeffTable();
      if ( (c->GetNPDF()==2 && c->GetNPDFDim() == 1)  ) {;} // ok!
      else {
	 error["InitGrids"]<<"Only half-matrix is implemented for grids for fixed-scale tables."<<endl; exit(1);
      }
      
      int nscalevar = fScaleFac.size();
      if (nscalevar==0 ) {
	 error["InitGrids"]<<"No scale factors found."<<endl;
      }
      c->Nscalevar.resize(1);			// 1 = NScaleDim
      c->Nscalevar[0]  = nscalevar;		// testing

      c->ScaleFac.resize(1);			// 1 = NScaleDim
      c->ScaleFac[0] = fScaleFac;
      c->XNode1.resize(GetNObsBin());
      for ( int i = 0 ; i < GetNObsBin() ; i ++ ) {
         c->XNode1[i]     = fKernX[i]->GetGrid();
      }

      int nscalenode = fKernMuS[0][0]->GetGrid().size();
      // scale nodes
      fastNLO::ResizeVector( c->ScaleNode, GetNObsBin(), 1, nscalevar, nscalenode );
      for ( int i = 0 ; i < GetNObsBin() ; i ++ ) {
         for(int k=0;k<nscalevar;k++){
	    c->ScaleNode[i][0][k] = fKernMuS[i][k]->GetGrid();
	 }
      }
      c->ResizeSigmaTilde();
   }

}


// ___________________________________________________________________________________________________
void  fastNLOCreate::InitInterpolationKernels() {
   debug["InitInterpolationKernels"]<<endl;
   if ( fIsWarmup ) {
      error["InitInterpolationKernels"]<<"Interpolation kernels can only be initialized in production runs. Warmup values must be known."<<endl;
   }

   fKernX.resize(GetNObsBin());
   if ( fIsFlexibleScale ) {
      fKernMu1.resize(GetNObsBin());
      fKernMu2.resize(GetNObsBin());
   }
   else {
      fKernMuS.resize(GetNObsBin());
      for ( int i = 0 ; i < GetNObsBin() ; i ++ ) {
	 fKernMuS[i].resize(fScaleFac.size());
      }
   }
   // todo. clean up memory
   vector<double> wrmX = DOUBLE_COL(WarmupValues,x_min);
   vector<double> wrmMu1Up, wrmMu1Dn;
   wrmMu1Dn = read_steer::getdoublecolumn("WarmupValues",GetWarmupHeader(0,"min"));
   wrmMu1Up = read_steer::getdoublecolumn("WarmupValues",GetWarmupHeader(0,"max"));
   if ( (int)wrmMu1Dn.size()!=GetNObsBin() || (int)wrmMu1Up.size()!= GetNObsBin() ){
      error["InitInterpolationKernels"]<<"Could not read warmup values for Mu1. Exiting."<<endl; exit(1);
   }
   vector<double> wrmMu2Up, wrmMu2Dn;
   if ( fIsFlexibleScale ){
      wrmMu2Dn = read_steer::getdoublecolumn("WarmupValues",GetWarmupHeader(1,"min"));
      wrmMu2Up = read_steer::getdoublecolumn("WarmupValues",GetWarmupHeader(1,"max"));
      if ( (int)wrmMu2Dn.size()!=GetNObsBin() || (int)wrmMu2Up.size()!= GetNObsBin() ){
         error["InitInterpolationKernels"]<<"Could not read warmup values for Mu2. Exiting."<<endl; exit(1);
      }
   }

   for ( int i = 0 ; i < GetNObsBin() ; i ++ ) {
      // ------------------------------------------------
      // init x-interpolation kernels
      // ------------------------------------------------
      debug["InitInterpolationKernels"]<<"Make x grid for obsbin="<<i<<endl;
      fKernX[i] = MakeInterpolationKernels(STRING(X_Kernel),wrmX[i],1);

      if ( BOOL(X_NoOfNodesPerMagnitude) )
         fKernX[i]->MakeGridsWithNNodesPerMagnitude( fastNLOInterpolBase::TranslateGridType(STRING(X_DistanceMeasure)), INT(X_NNodes)  );
      else
         fKernX[i]->MakeGrids( fastNLOInterpolBase::TranslateGridType(STRING(X_DistanceMeasure)), INT(X_NNodes));

      fKernX[i]->RemoveLastNode();

      // ------------------------------------------------
      // init scale1-interpolation kernels
      // ------------------------------------------------
      if ( fIsFlexibleScale ){
	 debug["InitInterpolationKernels"]<<"Make Mu1 grid for obsbin="<<i<<endl;
	 fKernMu1[i] = MakeInterpolationKernels(STRING(Mu1_Kernel),wrmMu1Dn[i],wrmMu1Up[i]);
	 fKernMu1[i]->MakeGrids( fastNLOInterpolBase::TranslateGridType(STRING(Mu1_DistanceMeasure)), INT(Mu1_NNodes));
	 //fKernMu1[i]->PrintGrid();
	 
	 // ------------------------------------------------
	 // init scale2-interpolation kernels
	 // ------------------------------------------------
	    debug["InitInterpolationKernels"]<<"Make Mu2 grid for obsbin="<<i<<endl;
	    fKernMu2[i] = MakeInterpolationKernels(STRING(Mu2_Kernel),wrmMu2Dn[i],wrmMu2Up[i]);
	    fKernMu2[i]->MakeGrids( fastNLOInterpolBase::TranslateGridType(STRING(Mu2_DistanceMeasure)), INT(Mu2_NNodes));
      }
      else {
	 for ( unsigned int k = 0 ; k<fScaleFac.size() ; k++ ) {
	    double muDn = fScaleFac[k]*wrmMu1Dn[i];
	    double muUp = fScaleFac[k]*wrmMu1Up[i];
	    fKernMuS[i][k] = MakeInterpolationKernels(STRING(Mu1_Kernel),muDn,muUp);
	    fKernMuS[i][k]->MakeGrids( fastNLOInterpolBase::TranslateGridType(STRING(Mu1_DistanceMeasure)), INT(Mu1_NNodes));
 	 }
      }
   }
}



// ___________________________________________________________________________________________________
fastNLOInterpolBase* fastNLOCreate::MakeInterpolationKernels(string KernelName, double xdn, double xup) {
   //! This function identifies the string-identifier
   //! and creates the corresponding fastNLO Interpolation kernel

   if ( KernelName == "CatmulRom" )
      return (fastNLOInterpolBase*)(new fastNLOInterpolCatmulRom(xdn,xup));
   // else if ( KernelName == "...") // todo implement other kernels here!
   //   return ...
   else {
      warn["MakeInterpolationKernels"]<<"Cannot find kernel routine:" <<KernelName<<" or kernel not (yet) implemented. Exiting."<<endl;
      exit(1);
   }
   return NULL; // default return
}

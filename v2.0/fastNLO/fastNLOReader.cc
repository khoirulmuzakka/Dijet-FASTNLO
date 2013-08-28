// Author: Daniel Britzger
// DESY, 14/08/2013
#include <string>
#include <algorithm>
#include <cfloat>
#include "fastNLOReader.h"
#include "read_steer.h"

#include "fastNLOCoeffAddFlex.h"
#include "fastNLOCoeffAddFix.h"

using namespace std;
using namespace fastNLO;
using namespace say;


fastNLOReader::fastNLOReader(): fastNLOTable() {
   SetClassName("fastNLOReader");
}

fastNLOReader::fastNLOReader(string filename) : fastNLOTable(filename) {
   SetClassName("fastNLOReader");
   debug["fastNLOReader"]<<"New fastNLOReader reading filename="<<filename<<endl;
   fCoeffData          = NULL;
   Coeff_LO_Ref        = NULL;
   Coeff_NLO_Ref       = NULL;
   fUnits               = fastNLO::kPublicationUnits;
   fMuRFunc             = fastNLO::kScale1;
   fMuFFunc             = fastNLO::kScale1;
   fPDFSuccess          = false;
   fAlphasCached        = 0.;
   fPDFCached           = 0.;
   SetFilename(filename);
 }

//______________________________________________________________________________



fastNLOReader::~fastNLOReader(void) {
}




void fastNLOReader::SetFilename(string filename) {
   debug["SetFilename"]<<"New filename="<<filename<<endl;
   ffilename    = filename;
   Init();
}


//______________________________________________________________________________



void fastNLOReader::Init() {
   debug["Init"]<<endl;

   // Initialize lists for BlockB's
   BBlocksSMCalc.resize(10);
   BBlocksNewPhys.resize(10);
   bUseSMCalc.resize(BBlocksSMCalc.size());
   bUseNewPhys.resize(BBlocksNewPhys.size());

   // Initialize Coeff's
   fastNLOCoeffBase* Coeff_LO   = NULL;
   fastNLOCoeffBase* Coeff_NLO  = NULL;
   fastNLOCoeffBase* Coeff_THC1 = NULL;
   fastNLOCoeffBase* Coeff_THC2 = NULL;
   fastNLOCoeffBase* Coeff_NPC1 = NULL;

   // run over all coefficient tables, identify and sort contributions.
   for (unsigned int i= 0; i<fCoeff.size() ; i++ ){
      fastNLOCoeffBase* c = GetCoeffTable(i);
      // give contribution a reasonable name
      //char nbuf[400];
      //       sprintf(nbuf,"Coeff. %s %s %s",
      // 	      _ContrName[c->GetIContrFlag1()-1].c_str(),_OrdName[c->GetIContrFlag1()-1][c->GetIContrFlag2()-1].c_str(),_fNSDep[c->GetNScaleDep()].c_str());
      //       c->SetName(nbuf);

      // data
      if ( fastNLOCoeffData::CheckCoeffConstants(c,true) ) {
	 debug["Init"]<<"Found data table."<<endl;
	 //c->SetName("Data");
	 if ( fCoeffData ) warn["Init"]<<"Already one data table present. Substituting."<<endl;
	 fCoeffData = (fastNLOCoeffData*)c;
      }
      // additive contributions
      else if ( fastNLOCoeffAddBase::CheckCoeffConstants(c,true) ) {
	 // Reference table, implemented only for LO or NLO
	 if ( ((fastNLOCoeffAddBase*)c)->IsReference() ) { 
	    debug["Init"]<<"Found reference table."<<endl;
	    if ( c->IsLO() ) {
// 	       if ( fastNLOCoeffAddFix::CheckCoeffConstants(&c,true) )   c->SetName("Coeff. LO Reference. v2.0.");
// 	       if ( fastNLOCoeffAddFlex::CheckCoeffConstants(&c,true) )   c->SetName("Coeff. LO Reference. v2.1."); // not forseen
	       Coeff_LO_Ref           = (fastNLOCoeffAddBase*)c;
	    } 
	    else if ( c->IsNLO() ) {
// 	       if ( fastNLOCoeffAddFix::CheckCoeffConstants(&c,true) )   c->SetName("Coeff. NLO Reference. v2.0.");
// 	       if ( fastNLOCoeffAddFlex::CheckCoeffConstants(&c,true) )   c->SetName("Coeff. NLO Reference. v2.1."); // not foreseen
	       Coeff_NLO_Ref  = (fastNLOCoeffAddBase*)c;
	    } else {
	       error["ReadTable"]<<"Reference tables are only implemented for fixed order (LO and NLO), stopped!\n";
	       exit(1);
	    }
	 }
	 // Additive fixed order (perturbative) contribution
	 else if ( c->GetIContrFlag1() == 1 ) {
	    if ( c->IsLO() )		Coeff_LO  = c;
	    else if ( c->IsNLO() )	Coeff_NLO = c;
	 }
	 // Threshold corrections
	 else if ( c->GetIContrFlag1() == 2 ) {
	    if ( c->GetIContrFlag2() == 1 )		Coeff_THC1 = c;
	    else if ( c->GetIContrFlag2() == 2 )	Coeff_THC2 = c;
	    else {
               error["Init"]<<"Threshold correction implemented only up to 2-loops, exiting!\n";
               exit(1);
            }
	 }
      }
      // multiplicative corrections
      else if ( fastNLOCoeffMult::CheckCoeffConstants(c,true) ) {
	 // Non-perturbative corrections
	 if ( c->GetIContrFlag1()==4 )	Coeff_NPC1 = c;
	 else {
	    error["ReadTable"]<<"Further multiplicative corrections not yet implemented, stopped!\n";
	    exit(1);
	 }
      }
   }

   // Assign NPC, switch off by default
   if (Coeff_NPC1) {
      BBlocksSMCalc[kNonPerturbativeCorrection].push_back(Coeff_NPC1);
      bUseSMCalc[kNonPerturbativeCorrection].push_back(false);
   }

   // Assign THC, switch off by default
   if (Coeff_THC1) {
      BBlocksSMCalc[kThresholdCorrection].push_back(Coeff_THC1);
      bUseSMCalc[kThresholdCorrection].push_back(false);
   }
   if (Coeff_THC2) {
      BBlocksSMCalc[kThresholdCorrection].push_back(Coeff_THC2);
      bUseSMCalc[kThresholdCorrection].push_back(false);
   }

   // Assign fixed order calculations (LO must be [0]), switch on by default
   if (Coeff_LO)  {
      BBlocksSMCalc[kFixedOrder].push_back(Coeff_LO);
      bUseSMCalc[kFixedOrder].push_back(true);
   } else {
      error["Init"]<<"Could not find any LO Calculation. Exiting!\n";
      exit(1);
   }
   if (Coeff_NLO) {
      BBlocksSMCalc[kFixedOrder].push_back(Coeff_NLO);
      bUseSMCalc[kFixedOrder].push_back(true);
   } else {
      error["Init"]<<"Could not find any NLO Calculation. Exiting!\n";
      exit(1);
   }

   //int iprint = 2;
   //PrintFastNLOTableConstants(iprint);
   InitScalevariation();
}


//______________________________________________________________________________


void fastNLOReader::InitScalevariation() {
   debug["InitScalevariation"]<<endl;
   fScaleFacMuR  = 1.;
   fScaleFacMuF  = 1.;
   fScalevar     = -1;

   if (!GetIsFlexibleScaleTable()) {
      fastNLOCoeffAddFix* cNLO = (fastNLOCoeffAddFix*)BBlocksSMCalc[kFixedOrder][kNextToLeading];
      for (int iscls=0; iscls< GetNScaleVariations(); iscls++) {
         const double muFac = cNLO->GetScaleFactor(iscls);
         if (fabs(muFac-1.0) < 1.e-7) {
            SetScaleVariation(iscls,true);
            break;
         }
      }
      if (fScalevar == -1) {
         error["InitScalevariation"]<<"Could not found scale variation with scale factor 1.0. Exiting.\n";
         exit(1);
      }
   } else {
      // this is a MuVar table. You can vary mu_f and mu_r independently by any factor
      // and you can choose the functional form of mu_f and mu_r as functions of
      // scale1 and scale1 (called partly scaleQ2 and scalePt).

      fastNLOCoeffAddFlex* cNLO = (fastNLOCoeffAddFlex*)BBlocksSMCalc[kFixedOrder][kNextToLeading];
      if (cNLO->GetScaleDescr()[0].size() <0) { // ???
         warn["InitScalevariation"]<<"No scaledescription available.\n"; 
         SetFunctionalForm(kScale1 , kMuR);
         SetFunctionalForm(kScale1 , kMuF);
         return;
      }

      // ---- DIS ---- //
      if (cNLO->GetNPDF() == 1) {
         SetFunctionalForm(kQuadraticMean , kMuR);
         SetFunctionalForm(kScale1 , kMuF);
      }
      // ---- HHC --- //
      else if (cNLO->GetNPDF() == 2) {
         SetFunctionalForm(kScale1 , kMuR);
         SetFunctionalForm(kScale1 , kMuF);
      } else {
         error<<"Unknown process.\n";
         exit(1);
      }
   }
}




//______________________________________________________________________________


double fastNLOReader::SetScaleVariation(int scalevar , bool FirstCall) {
   debug["SetScaleVariation"]<<"Setting to scalevar="<<scalevar<<endl;
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

   if (GetIsFlexibleScaleTable()) {
      info["SetScaleVariation"]<<"This is a flexible-scale table. No Scalevariation tables available!"<<endl;
      man<<"You can choose freely (within reason) a factorization scale factor. Your Scalevar has to be '0'.\n";
      man<<"Please use SetScaleFacMuR(double) and SetScaleFacMuF(double) to set scale factors.\n";
      return 0;
   }
   else {
      // Check for maximal scale variation of all rel. and active SM calcs
      int scalevarmax = GetNScaleVariations();
      fastNLOCoeffAddFix* cNLO = (fastNLOCoeffAddFix*)B_NLO();  
	 if (scalevar >= scalevarmax) {
	    warn["SetScaleVariation"]<<"This table has only "<<scalevarmax<<" scale variation(s) stored!"<<endl;
	    man<<"For the currently active contributions. You wanted to access the non-existing number "<<scalevar<<endl;
	    man<<"Using '0' instead."<<endl;;
	    fScalevar = 0;
	    return cNLO->GetScaleFactor(0);
	 }

      fScalevar     = scalevar;
      fScaleFacMuF  = cNLO->GetScaleFactor(fScalevar);
      info["SetScaleVariation"]
         <<"Selecting MuF table according to a multiplicative scale factor of the factorization scale of "
         <<fScaleFacMuF<<" times the nominal scale."<<endl;

      // check for threshold corrections.
      if (!BBlocksSMCalc[kThresholdCorrection].empty()) {
	 bool lkth = false;
	 for (unsigned int i = 0 ; i <BBlocksSMCalc[kThresholdCorrection].size() ; i++) {
	    if (bUseSMCalc[kThresholdCorrection][i]) {
	       lkth = true;
	    }
	 }
	 if (lkth && fabs(fScaleFacMuR-fScaleFacMuF) > DBL_MIN) {
	    fScaleFacMuR = fScaleFacMuF;
	    warn["SetScaleVariation."]<<"Threshold corrections do not allow variations of the renormalization scale!"<<endl;
	    man<<"The scale factor for MuR has been set equal to the one for MuF = "<<fScaleFacMuF<<endl;
	    man<<"Either select a different simultaneous scale variation i, if possible, via fastNLOReader::SetScaleVariation(i)"<<endl;
	    man<<"or deactivate first all threshold corrections using fastNLOReader::SetContributionON(kTresholdCorrections,Id,false)."<<endl;
	 }
      }
      return cNLO->GetScaleFactor(fScalevar);
   }
}



//______________________________________________________________________________

void fastNLOReader::SetContributionON(ESMCalculation eCalc , unsigned int Id , bool SetOn) {
   //
   // 'Activate' contribution to be considered in the calculation
   // of the cross section.
   //

   // sanity check
   if (bUseSMCalc[eCalc].empty() || BBlocksSMCalc.empty()) {
      warn["SetContributionON"]
	 <<"This contribution ("<<_ContrName[eCalc]<<") does not exist in this table. Cannot switch it On/Off. Ignoring call.\n";
      return;
   }
   if (bUseSMCalc[eCalc].size() < Id || BBlocksSMCalc[eCalc].size() < Id || !BBlocksSMCalc[eCalc][Id]) {
      warn["SetContributionON"]
            <<"This Id="<<Id<<" does not exist for this contribtion. Cannot switch it On/Off. Ignoring call.\n";
      return;
   }

   info<<(SetOn?"Activating":"Deactivating")
       <<" contribution '"<<_ContrName[eCalc]
       <<" with Id="<<Id<<endl;

   if (!bUseSMCalc[eCalc][Id] && SetOn) {
      fastNLOCoeffAddBase* c = (fastNLOCoeffAddBase*)BBlocksSMCalc[eCalc][Id];
      if (!c->GetIAddMultFlag()) {
         // Fill alpha_s cache
         debug["SetContributionON"]<<"Call FillAlphasCache for contribution eCalc="<<eCalc<<"\tId="<<Id<<endl;
         if (!GetIsFlexibleScaleTable(c)) {
            FillAlphasCacheInBlockBv20((fastNLOCoeffAddFix*)c);
         } else {
            FillAlphasCacheInBlockBv21((fastNLOCoeffAddFlex*)c);
         }
         // Fill PDF cache
         debug["SetContributionON"]<<"Call FillPDFCache for contribution eCalc="<<eCalc<<"\tId="<<Id<<endl;
         // linear: DIS-case
         // ---- DIS ---- //
         if (c->GetIPDFdef1() == 2) {
            if (c->GetNPDFDim() == 0) {
               if (!GetIsFlexibleScaleTable(c)) {
                  FillBlockBPDFLCsDISv20((fastNLOCoeffAddFix*)c);
               } else {
                  FillBlockBPDFLCsDISv21((fastNLOCoeffAddFlex*)c);
               }
            }
         }
         // ---- pp ---- //
         else if (c->GetIPDFdef1() == 3) {
            if (c->GetNPDFDim() == 1) {
               if (!GetIsFlexibleScaleTable()) {
                  FillBlockBPDFLCsHHCv20((fastNLOCoeffAddFix*)c);
               } else {
                  FillBlockBPDFLCsHHCv21((fastNLOCoeffAddFlex*)c);
               }
            } else {
               error<<"Only half matrices for hh is implemented.\n";
               exit(1);
            }
         } else {
            error<<"Tables not yet implemented.\n";
         }
      }
   }
   // set the new value
   bUseSMCalc[eCalc][Id] = SetOn;

}


//______________________________________________________________________________


int fastNLOReader::GetNScaleVariations() const {
   if (GetIsFlexibleScaleTable()) {
      info["GetNScaleVariations"]<<"This is a 'flexible-scale' table, therefore you can choose all desired scale variations."<<endl;
      return 0;
   }
   else {
      // Check for maximal scale variation of all rel. and active SM calcs
      // Assume a maximum of 10!
      unsigned int scalevarmax = 10;
      for (unsigned int j = 0 ; j<BBlocksSMCalc.size() ; j++) {
	 if (!BBlocksSMCalc.empty()) {
	    for (unsigned int i = 0 ; i<BBlocksSMCalc[j].size() ; i++) {
	       fastNLOCoeffAddFix* c = (fastNLOCoeffAddFix*)BBlocksSMCalc[j][i];
	       // Do not check pQCD LO or mult. corrections
	       if (bUseSMCalc[j][i] && !c->GetIAddMultFlag() &&
		   !(j==kFixedOrder && i==kLeading)) {
		  if (c->GetNScalevar() < (int)scalevarmax) {
		     scalevarmax = c->GetNScalevar();
		  }
	       }
	    }
	 }
      }
      debug["GetNScaleVariations"]<<"Found "<<scalevarmax<<" scale variations."<<endl;
      return scalevarmax;
      return 0;
   }
}

//______________________________________________________________________________


vector < double > fastNLOReader::GetScaleFactors() const {
   if (GetIsFlexibleScaleTable()) {
      info["GetScaleFactors"]<<"This is a 'flexible scale table', therefore you can choose all desired scale variations."<<endl;
      return vector<double>();
   }
   else 
      return ((fastNLOCoeffAddFix*)BBlocksSMCalc[kFixedOrder][kNextToLeading])->GetAvailableScaleFactors();
}


//______________________________________________________________________________


vector < double > fastNLOReader::GetCrossSection() {
   // Get fast calculated NLO cross section
   if (XSection.empty()) CalcCrossSection();
   return XSection;
}


//______________________________________________________________________________


vector < double > fastNLOReader::GetKFactors() {
   // Get ratio of fast calculated NLO to LO cross section
   if (XSection.empty()) CalcCrossSection();
   return kFactor;
}

//______________________________________________________________________________

vector < double > fastNLOReader::GetQScales(int irelord) {
   // Get XSection weighted Q scale in bin
   if (XSection.empty()) CalcCrossSection();
   if (irelord == 0) {
      return QScale_LO;
   } else {
      return QScale;
   }
}

//______________________________________________________________________________


vector < double > fastNLOReader::GetReferenceCrossSection() {
   // Get reference cross section from direct nlojet++ calculation
   if (XSectionRef.empty() && XSectionRef_s1.empty()) {
      CalcReferenceCrossSection();
   }
   if (GetIsFlexibleScaleTable()) {
      if (fMuFFunc == kScale1 && fMuRFunc == kScale1)                   return XSectionRef_s1;
      else if (fMuFFunc == kScale2 && fMuRFunc == kScale2)              return XSectionRef_s2;
      else if (fMuFFunc == kQuadraticMean && fMuRFunc == kQuadraticMean)return XSectionRefMixed;
      else return XSectionRefMixed;
   } else return XSectionRef; // XSectionRef from BlockB-Ref
   return XSectionRef;
}

//______________________________________________________________________________



void fastNLOReader::CalcReferenceCrossSection() {
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

   if (!GetIsFlexibleScaleTable()) {
      if (Coeff_LO_Ref && Coeff_NLO_Ref) {
	 for (int i=0; i<NObsBin; i++) {
	    double unit = fUnits==kAbsoluteUnits ? BinSize[i] : 1.;
	    for (int l=0; l<Coeff_LO_Ref->GetNSubproc(); l++) {
	       XSectionRef[i] +=  ((fastNLOCoeffAddFix*)Coeff_LO_Ref)->GetSigmaTilde(i,0,0,0,l) * unit; // no scalevariations in LO tables
	    }
	    for (int l=0; l<Coeff_NLO_Ref->GetNSubproc(); l++) {
	       XSectionRef[i] +=  ((fastNLOCoeffAddFix*)Coeff_NLO_Ref)->GetSigmaTilde(i,fScalevar,0,0,l) * unit;
	    }
	 }
      }
      else 
	 warn["CalcReferenceCrossSection"]<<"No reference cross sections for LO and NLO available.\n";
   }
   else {
      for (int i=0; i<NObsBin; i++) {
         double unit = fUnits==kAbsoluteUnits ? BinSize[i] : 1.;
	 fastNLOCoeffAddFlex* cLO = (fastNLOCoeffAddFlex*)BBlocksSMCalc[kFixedOrder][kLeading];
         for (int n=0; n<cLO->GetNSubproc(); n++) {
            XSectionRefMixed[i]             += cLO->SigmaRefMixed[i][n] * unit;
            XSectionRef_s1[i]               += cLO->SigmaRef_s1[i][n] * unit;
            XSectionRef_s2[i]               += cLO->SigmaRef_s2[i][n] * unit;
         }
	 fastNLOCoeffAddFlex* cNLO = (fastNLOCoeffAddFlex*)BBlocksSMCalc[kLeading][kNextToLeading];
	 for (int n=0; n<cNLO->GetNSubproc(); n++) {
            XSectionRefMixed[i]             += cNLO->SigmaRefMixed[i][n] * unit;
            XSectionRef_s1[i]               += cNLO->SigmaRef_s1[i][n] * unit;
            XSectionRef_s2[i]               += cNLO->SigmaRef_s2[i][n] * unit;
         }
      }
   }
}


//______________________________________________________________________________


bool fastNLOReader::PrepareCache() {
   // check pdf cache
   const double PDFcks = CalcNewPDFChecksum();
   if (fPDFCached==0. || (fPDFCached!=0. && fabs(PDFcks/fPDFCached -1.) > 1.e-7)) {
      debug["PrepareCache"]<<"Need to refill PDFCache, since PDFCecksum="<<PDFcks<<" and fPDFCached="<<fPDFCached<<endl;
      FillPDFCache(PDFcks);
   } else  debug["PrepareCache"]<<"No need to refill PDFCache."<<endl;

   // check pdf cache
   if (!fPDFSuccess) {
      error["PrepareCache"]<<"Cannot calculate cross sections. PDF has not been initalized successfully."<<endl;
      return false;
   }

   // check alpha_s cache
   const double asref = CalcReferenceAlphas();
   if (fAlphasCached == 0. || fAlphasCached != asref) {
      debug["PrepareCache"]<<"Need to refill AlphasCache, since fAlphasCached="<<fAlphasCached<<endl;
      FillAlphasCache();
   }
   // do we now have an alphas?
   if (fAlphasCached==0. || fAlphasCached != asref) {
      error["PrepareCache"]<<"Filling of alpha_s cache failed. fAlphasCached="<<fAlphasCached<<"\tasref="<<asref<<endl;
      return false;
   }
   return true;
}


//______________________________________________________________________________



void fastNLOReader::CalcCrossSection() {
   debug["CalcCrossSection"]<<endl;
   //
   //  Initialize the internal arrays with the NLO cross sections
   //  with the information from the FastNLO file, the pdf and
   //  the defined alpha_s
   //
   // xs = (sum(all active additive (perturbative) contr,) + sum(all active SM corrections) + sum(all active new physics contr.) ) * prod(all active multipl. contr.)

   XSection_LO.clear();
   XSection.clear();
   XSection_LO.resize(NObsBin);
   XSection.resize(NObsBin);
   kFactor.clear();
   kFactor.resize(NObsBin);
   QScale_LO.clear();
   QScale.clear();
   QScale_LO.resize(NObsBin);
   QScale.resize(NObsBin);

   // handle alpha_s and PDF Cache
   bool CacheOK = PrepareCache();
   if (!CacheOK) {
      error["CalcCrossSection"]<<"Caching failed. Cannot calculate cross sections."<<endl;
      return;
   }

   // perturbative (additive) contributions
   for (unsigned int j = 0 ; j<BBlocksSMCalc.size() ; j++) {
      for (unsigned int i = 0 ; i<BBlocksSMCalc[j].size() ; i++) {
	 if ( bUseSMCalc[j][i] ) {
	    if ( fastNLOCoeffAddFlex::CheckCoeffConstants(BBlocksSMCalc[j][i],true) )
	       CalcCrossSectionv21((fastNLOCoeffAddFlex*)BBlocksSMCalc[j][i]);
	    else if ( fastNLOCoeffAddFix::CheckCoeffConstants(BBlocksSMCalc[j][i],true) )
	       CalcCrossSectionv20((fastNLOCoeffAddFix*)BBlocksSMCalc[j][i]);
	 }
      }
   }


   // contributions from the a-posteriori scale variation
   if (!GetIsFlexibleScaleTable()) {
      fastNLOCoeffAddFix* cNLO = (fastNLOCoeffAddFix*)B_NLO();
      if ( fabs( fScaleFacMuR - cNLO->GetScaleFactor(fScalevar) ) > DBL_MIN ) {
	 CalcAposterioriScaleVariation();
      }
   }

   // calculate LO cross sections
   if (!GetIsFlexibleScaleTable())
      CalcCrossSectionv20((fastNLOCoeffAddFix*)B_LO(),true);
   else
      CalcCrossSectionv21((fastNLOCoeffAddFlex*)B_LO(),true);


   // non-perturbative corrections (multiplicative corrections)
   for (unsigned int j = 0 ; j<BBlocksSMCalc.size() ; j++) {
      for (unsigned int i = 0 ; i<BBlocksSMCalc[j].size() ; i++) {
   	 if (bUseSMCalc[j][i] ) {
	    if ( fastNLOCoeffMult::CheckCoeffConstants( BBlocksSMCalc[j][i] , true ) ) {
	       fastNLOCoeffMult* cMult = (fastNLOCoeffMult*) BBlocksSMCalc[j][i];
	       if ( cMult->GetIContrFlag1() == 4 && cMult->GetIContrFlag2() == 1) {
		  debug["CalcCrossSection"]<<"Adding multiplicative non-perturbative correction."<<endl;
		  for (int iB=0; iB<NObsBin; iB++) {
		     XSection[iB] *= cMult->GetMultFactor(iB);
		     //            XSection_LO[iB]     *= BBlocksSMCalc[j][i]->fact[iB];
		  }
	       } else {
		  error["CalcCrossSection"]<<"Found unknown multiplicative correction. Printing coeff table and exiting..."<<endl;
		  cMult->Print();
		  exit(1);
	       }
	    }
	 }
      }
   }

   // ---- k-factor calculation ---- //
   debug["CalcCrossSection"]<<"Calculate k-factors: xs/xs_LO"<<endl;
   for (int i=0; i<NObsBin; i++) {
      kFactor[i] = XSection[i] / XSection_LO[i];
   }

   // ---- Q-scale calculation ---- //
   debug["CalcCrossSection"]<<"Calculate Q-scales: xsQ/xs"<<endl;
   for (int i=0; i<NObsBin; i++) {
      QScale_LO[i] = QScale_LO[i]/XSection_LO[i];
      QScale[i]    = QScale[i]/XSection[i];
   }
}
//______________________________________________________________________________


void fastNLOReader::CalcAposterioriScaleVariation() {
   double scalefac       = fScaleFacMuR/fScaleFacMuF;
   debug["CalcAposterioriScaleVariation"]<<"scalefac="<<scalefac<<endl;
   if ( GetIsFlexibleScaleTable() ) { error["CalcAposterioriScaleVariation"]<<"This function is only reasonable for non-flexible scale tables."<<endl; exit(1);}
   fastNLOCoeffAddFix* cLO  = (fastNLOCoeffAddFix*) B_LO();
   vector<double>* XS    = &XSection;
   vector<double>* QS    = &QScale;
   const double n     = cLO->GetNpow();
   const double L     = log(scalefac);
   const double beta0 = (11.*3.-2.*5)/3.;
   for (int i=0; i<NObsBin; i++) {
      int nxmax = cLO->GetNxmax(i);
      double unit = fUnits==kAbsoluteUnits ? BinSize[i] : 1.;
      for (int j=0; j<cLO->GetTotalScalenodes(); j++) {
         double asnp1 = pow(cLO->AlphasTwoPi_v20[i][j],(n+1)/n);//as^n+1
         for (int k=0; k<nxmax; k++) {
            for (int l=0; l<cLO->GetNSubproc(); l++) {
               double clo  = cLO->GetSigmaTilde(i,0,j,k,l) *  cLO->PdfLc[i][j][k][l] * unit;
               double xsci = asnp1 * clo * n * L * beta0;
               double mur  = fScaleFacMuR * cLO->GetScaleNode(i,0,j);
               XS->at(i) +=  xsci;
               QS->at(i) +=  xsci*mur;
            }
         }
      }
   }
}



//______________________________________________________________________________


void fastNLOReader::CalcCrossSectionv21(fastNLOCoeffAddFlex* c , bool IsLO) {
   //debug["CalcCrossSectionv21"]<<"B->fname="<<B->GetName()<<"\tNpow="<<B->GetNpow()<<"\tIsLO="<<IsLO<<endl;
   //
   //  Cross section calculation for DIS and HHC tables in v2.1 format
   //

   vector<double>* XS = IsLO ? &XSection_LO : &XSection;
   vector<double>* QS = IsLO ? &QScale_LO : &QScale;
   //c->fact.resize(NObsBin);
   for (int i=0; i<NObsBin; i++) {
      //B->fact[i]=0;
      int nxmax = c->GetNxmax(i);
      double unit = (fUnits==kAbsoluteUnits) ? BinSize[i] : 1.;
      for (unsigned int jS1=0; jS1<c->GetNScaleNode1(i); jS1++) {
         for (unsigned int kS2=0; kS2<c->GetNScaleNode2(i); kS2++) {
            double Q2		= pow(c->GetScaleNode1(i,jS1),2);
            double mur		= CalcMu(kMuR , c->GetScaleNode1(i,jS1) ,  c->GetScaleNode2(i,kS2) , fScaleFacMuR);
            double muf		= CalcMu(kMuF , c->GetScaleNode1(i,jS1) ,  c->GetScaleNode2(i,kS2) , fScaleFacMuF);
            double mur2		= pow(mur,2);
            double muf2		= pow(muf,2);
            for (int x=0; x<nxmax; x++) {
               for (int n=0; n<c->GetNSubproc(); n++) {
		  double as		= c->AlphasTwoPi[i][jS1][kS2]; 
		  double pdflc		= c->PdfLcMuVar[i][x][jS1][kS2][n];
                  if (pdflc == 0.) continue;
                  double fac  = as * pdflc * unit;
                  double xsci =  c->SigmaTildeMuIndep[i][x][jS1][kS2][n] * fac;
		  if ( c->GetNScaleDep() >= 5 ) {
		     xsci             += c->SigmaTildeMuFDep [i][x][jS1][kS2][n] * log(muf2) * fac;
		     xsci             += c->SigmaTildeMuRDep [i][x][jS1][kS2][n] * log(mur2) * fac;
		     if ( c->GetIPDFdef1() == 2 ) {   // DIS tables use log(mu/Q2) instead of log(mu)
			xsci -= c->SigmaTildeMuFDep [i][x][jS1][kS2][n] * log(Q2) * fac;
			xsci -= c->SigmaTildeMuRDep [i][x][jS1][kS2][n] * log(Q2) * fac;
		     }
		  }
                  XS->at(i)   += xsci;
                  //B->fact[i]  += xsci;
                  QS->at(i)   += xsci*mur;
               }
            }
         }
      }
   }
}


//______________________________________________________________________________


void fastNLOReader::CalcCrossSectionv20(fastNLOCoeffAddFix* c , bool IsLO) {
   debug["CalcCrossSectionv20"]<<"Npow="<<c->GetNpow()<<"\tIsLO="<<IsLO<<endl;
   //
   //  Cross section calculation in v2.0 format
   //

   int scaleVar          = c->GetNpow() == ILOord ? 0 : fScalevar;
   vector<double>* XS    = IsLO ? &XSection_LO : &XSection;
   vector<double>* QS    = IsLO ? &QScale_LO : &QScale;
   //B->fact.resize(NObsBin);
   for (int i=0; i<NObsBin; i++) {
      //B->fact[i] = 0;
      int nxmax = c->GetNxmax(i);
      double unit = fUnits==kAbsoluteUnits ? BinSize[i] : 1.;
      for (int j=0; j<c->GetTotalScalenodes(); j++) {
	 double scalefac = fScaleFacMuR/c->GetScaleFactor(scaleVar);
	 double mur      = scalefac * c->GetScaleNode(i,scaleVar,j);
         for (int k=0; k<nxmax; k++) {
            for (int l=0; l<c->GetNSubproc(); l++) {
	       // hier hier todo
	       double xsci     = c->GetSigmaTilde(i,scaleVar,j,k,l) *  c->AlphasTwoPi_v20[i][j]  * c->PdfLc[i][j][k][l] * unit;
               XS->at(i)      +=  xsci;
               //B->fact[i]     +=  xsci;
               QS->at(i)      +=  xsci*mur;
            }
         }
      }
   }
}

//______________________________________________________________________________


void fastNLOReader::SetUnits(EUnits Unit) {
   if (fUnits != Unit) {
      fUnits  = Unit;
      //CalcCrossSection();
   } else {
      // nothing todo
   }
}


//______________________________________________________________________________


void fastNLOReader::FillAlphasCache() {
   debug["FillAlphasCache"]<<endl;
   //
   //  Fill the internal alpha_s cache.
   //  This is usally called automatically. Only if you
   //  make use of ReFillCache==false options, you have
   //  to take care of this filling by yourself.
   //

   // check if the alpha_s value is somehow reasonable
   debug["FillAlphasCache"]<<"Sanity check!"<<endl;
   TestAlphas();

   // is there a need for a recalclation?
   const double asNew = CalcReferenceAlphas();
   if (asNew == fAlphasCached) {
      debug["FillAlphasCache"]<<"No need for a refilling of AlphasCache. asNew==fAlphasCached="<<asNew<<endl;
   } else {
      fAlphasCached = asNew;
      for (unsigned int j = 0 ; j<BBlocksSMCalc.size() ; j++) {
	 for (unsigned int i = 0 ; i<BBlocksSMCalc[j].size() ; i++) {
	    // Check that this contribution type j and no. i should actually be used
	    // Otherwise deactivation of e.g. threshold corr. is not respected here
	    if (bUseSMCalc[j][i] ) {
	       fastNLOCoeffBase* c = BBlocksSMCalc[j][i];
	       if ( fastNLOCoeffAddFlex::CheckCoeffConstants(c,true) )
		  FillAlphasCacheInBlockBv21((fastNLOCoeffAddFlex*)c);
	       else if ( fastNLOCoeffAddFix::CheckCoeffConstants(c,true) )
		  FillAlphasCacheInBlockBv20((fastNLOCoeffAddFix*)c);
	       else {
		  error["FillAlphasCache"]<<"Could not identify contribution. Printing."<<endl;
		  c->Print();
	       }
	    }
	 }
      }
   }
}


//______________________________________________________________________________


void fastNLOReader::FillAlphasCacheInBlockBv20(fastNLOCoeffAddFix* c) {
   //
   //  Internal method for filling alpha_s cache
   //
   
   // todo: the flag IScaleDep should also indicate whether this contribution may contain scale variations
   int scaleVar          = c->GetNpow() == ILOord ? 0 : fScalevar;

   // Sanity check that scaleVar is in allowed range
   // For thresh. corr. can otherwise lead to inf and then segfault!
   if (scaleVar >= GetNScaleVariations()) {
      error<<"Trying to refresh  cache for non-existing scale variation no. "<<scaleVar<<" while only "<<GetNScaleVariations()<<" ex\
ist in total. Aborted."<<endl;
      exit(1);
   }
   double scalefac       = fScaleFacMuR/c->GetScaleFactor(scaleVar);
   debug["FillAlphasCacheInBlockBv20"]<<"scalefac="<<scalefac<<"\tscaleVar="<<scaleVar<<endl;

   for (int i=0; i<NObsBin; i++) {
      for (int j=0; j<c->GetTotalScalenodes(); j++) {
         double mur        = scalefac * c->GetScaleNode(i,scaleVar,j);
         double as         = CalcAlphas(mur);
         c->AlphasTwoPi_v20[i][j] = pow(as/TWOPI , c->GetNpow());
      }
   }
}


//______________________________________________________________________________


void fastNLOReader::FillAlphasCacheInBlockBv21(fastNLOCoeffAddFlex* c) {
   //
   //  Internal method for filling alpha_s cache
   //

   for (int i=0; i<NObsBin; i++) {
      for (unsigned int jS1=0; jS1<c->GetNScaleNode1(i); jS1++) {
         for (unsigned int kS2=0; kS2<c->GetNScaleNode2(i); kS2++) {
            double mur              = CalcMu(kMuR , c->GetScaleNode1(i,jS1) ,  c->GetScaleNode2(i,kS2) , fScaleFacMuR);
            double as               = CalcAlphas(mur);
            double alphastwopi      = pow(as/TWOPI, c->GetNpow() );
            c->AlphasTwoPi[i][jS1][kS2] = alphastwopi;
         }
      }
   }
}


//______________________________________________________________________________



double fastNLOReader::CalcAlphas(double Q) {
   //
   //  Internal method for calculating the alpha_s(mu)
   //
   return EvolveAlphas(Q);
}


//______________________________________________________________________________


double fastNLOReader::CalcReferenceAlphas() {
   double mu = 0;
   if (GetIsFlexibleScaleTable()) {
      if (fMuRFunc==kExtern) mu = (*Fct_MuR)(91.,1.)*(fScaleFacMuR+0.1);
      else mu = 91.1876111111+(fMuRFunc*0.1)+(fScaleFacMuR);
   } else mu = 91.187611111115*(fScaleFacMuR+0.1)+fScalevar*0.1;
   double as = CalcAlphas(mu);
   if (isnan(as)) {
      error["CalcReferenceAlphas"]<<"Reference alphas is a 'nan' for scale mu="<<mu<<endl;
      //exit(1);
   }
   return as;
}


//______________________________________________________________________________


double fastNLOReader::CalcNewPDFChecksum() {
   // calculate a PDF checksum to
   // decide, whether PDF cache has to be refilled

   // init PDF and check success
   debug["CalcNewPDFChecksum"]<<"Call InitPDF() in user module."<<endl;
   fPDFSuccess = InitPDF();
   debug["CalcNewPDFChecksum"]<<"Return value InitPDF() = "<<fPDFSuccess<<endl;
   if (!fPDFSuccess) {
      warn["CalcPDFChecksum"]<<"PDF initialization failed. Please check PDF interface in your FastNLO user module."<<endl;
      return 0.;
   }

   // calculate checksum for some scales and flavors
   double muf = 0;
   if (GetIsFlexibleScaleTable()) {
      if (fMuFFunc==kExtern) muf = (*Fct_MuF)(91.,1.)/91.*(fScaleFacMuF+0.1) ;
      else muf = 91.1+0.1*fMuFFunc+fScaleFacMuF;
   } else muf=(fScaleFacMuF+0.1)+fScalevar*0.1;
   double cks = CalcChecksum(muf);
   return cks;
}


//______________________________________________________________________________


double fastNLOReader::CalcChecksum(double mufac) {
   debug["CalcChecksum"]<<"Calculate checksum of 13 flavors, 3 mu_f values, and 3 x-values, for scalefac="<<mufac<<endl;
   double cks = 0;
   vector<double> xfx(13);
   const double mf[3] = { 3,10,91.18};
   const double x[3] = {1.e-1,1.e-2,1.e-3};
   for (int jf = 0 ; jf<3 ; jf++) {
      double mu = mf[jf]* mufac;//(fScaleFacMuF+0.1)+fScalevar*0.1;
      for (int ix = 0 ; ix<3 ; ix++) {
         xfx = GetXFX(x[ix],mu);
         for (unsigned int fl = 0 ; fl<xfx.size() ; fl++) {
            cks+=xfx[fl];
         }
      }
   }
   debug["CalcChecksum"]<<"Calculated checksum = "<<cks<<endl;
   return cks;
}


//______________________________________________________________________________



bool fastNLOReader::TestAlphas() {
   const double as = CalcAlphas(91.18);
   if (as < 0.01 || as > 0.5) {
      warn["TestAlphas"]<<"The alphas value, returned by the user class seems to be unreasonably small/large."<<endl;
      warn["TestAlphas"]<<"The evolution code calculated alphas(Mz~91.18GeV) = "<<as<<endl;
      return false;
   }
   debug["TestAlphas"]<<"Sanity check of alpha_s(MZ=91.18) = "<<as<<endl;
   return true;
}


//______________________________________________________________________________



bool fastNLOReader::TestXFX() {
   vector<double> pdftest = GetXFX(1.e-2,10);
   if (pdftest.size() != 13) {
      error["TestXFX"]<<"The pdf array must have the size of 13 flavors.\n";
      return false;
   }
   // if ( pdftest[6] == 0. )printf("fastNLOReader. Warning. There seems to be no gluon in the pdf.\n");
   // double sum = 0;
   // for ( int i = 0 ; i<13 ; i++ ) sum+=fabs(pdftest[i]);
   // if ( sum== 0. ) printf("fastNLOReader. Error. All 13 pdf probabilities are 0. There might be sth. wrong in the pdf interface. Please check FastNLOUser::GetXFX().\n");
   for (int i = 0 ; i<13 ; i++) {
      if (pdftest[i] > 1.e10 || (pdftest[i] < 1.e-10 && pdftest[i] > 1.e-15)) {
         warn["TestXFX"]<<"The pdf probability of the "<<i<<"'s flavor seeems to be unreasonably large/small (pdf="<<pdftest[i]<<").\\n";
      }
   }
   return true;
}



//______________________________________________________________________________


void fastNLOReader::FillPDFCache(double chksum) {
   debug["FillPDFCache"]<<"Passed chksum="<<chksum<<". Do not recalculate checksum (which calls InitPDF()) if chksum!=0."<<endl;
   //
   //  Fill the internal pdf cache.
   //  This function has to be called by the user, since the
   //  pdf parameters and evolutions are calculated externally.
   //

   // reset checknum
   // check if the alpha_s value is somehow reasonable
   double PDFnew = chksum;
   if (chksum == 0.) {
      debug["FillPDFCache"]<<"Calculate Checksum!"<<endl;
      PDFnew = CalcNewPDFChecksum();
      if (PDFnew==0.) {
         warn["FillPDFCache"]<<"PDF Checksum is zero."<<endl;
      }
      debug["FillPDFCache"]<<"PDF Checksum = "<<PDFnew<<endl;
   }

   // is there a need for a recalculation?
   if (fPDFCached != 0. && fabs(PDFnew/fPDFCached - 1.) < 1.e-7) {
      debug["FillPDFCache"]<<"No need for a refilling of PDFCache. fPDFCached=RefreshPDFChecksum()"<<PDFnew<<endl;
   } else {
      debug["FillPDFCache"]<<"Refilling PDF cache"<<endl;
      fPDFCached = PDFnew;

      // check (or not) if the pdf is somehow reasonable
      TestXFX();

      for (unsigned int j = 0 ; j<BBlocksSMCalc.size() ; j++) {
	 for (unsigned int i = 0 ; i<BBlocksSMCalc[j].size() ; i++) {
	    // Check that this contribution type j and no. i should actually be used
	    // Otherwise deactivation of e.g. threshold corr. is not respected here
	    if (bUseSMCalc[j][i] ) {
	       // linear: DIS-case
	       // ---- DIS ---- //
	       fastNLOCoeffAddBase* c = (fastNLOCoeffAddBase*)BBlocksSMCalc[j][i];
	       if (c->GetIPDFdef1() == 2) {
		  if (c->GetNPDFDim() == 0) {
		     if (!GetIsFlexibleScaleTable(c)) 
			FillBlockBPDFLCsDISv20((fastNLOCoeffAddFix*)c);
		     else 
			FillBlockBPDFLCsDISv21((fastNLOCoeffAddFlex*)c);
		  }
	       }
	       // ---- pp ---- //
	       else if (c->GetIPDFdef1() == 3) {
		  if (c->GetNPDFDim() == 1) {
		     if (!GetIsFlexibleScaleTable(c)) FillBlockBPDFLCsHHCv20((fastNLOCoeffAddFix*)c);
		     else FillBlockBPDFLCsHHCv21((fastNLOCoeffAddFlex*)c);
		  } else {
		     error<<"Only half matrices for hh is implemented.\n";
		     exit(1);
		  }
	       } else {
		  error<<"IPDFdef of tables must be 1 or 2.\n";
	       }
	    }
	 }
      }
   }
}


//______________________________________________________________________________


void fastNLOReader::FillBlockBPDFLCsDISv20(fastNLOCoeffAddFix* c) {
   debug["FillBlockBPDFLCsDISv20"]<<endl;
   // todo: flag IScaleDep should indicate whether scale variations may exist or not.
   int scaleVar          = c->GetNpow() == ILOord ? 0 : fScalevar;
   double scalefac       = (c->GetScaleFactor(scaleVar) == fScaleFacMuF) ? 1. : fScaleFacMuF;
   vector<double> xfx(13); // PDFs of all partons
   if (!GetIsFlexibleScaleTable(c)) {
      for (int i=0; i<NObsBin; i++) {
         int nxmax = c->GetNxmax(i);
         for (int j=0; j<c->GetNScaleNode(); j++) {
            for (int k=0; k<nxmax; k++) {
               double xp     = c->GetXNode1(i,k);
               double muf    = scalefac * c->GetScaleNode(i,scaleVar,j);
               xfx = GetXFX(xp,muf);
	       c->PdfLc[i][j][k] = CalcPDFLinearCombination(c,xfx);
	       //                vector < double > buffer = CalcPDFLinearCombDIS(xfx , c->GetNSubproc());
	       //                for (int l=0; l<c->GetNSubproc(); l++) {
	       //                   c->PdfLc[i][j][k][l] = buffer[l];
	       //                }
            }
         }
      }
   }
}


//______________________________________________________________________________


void fastNLOReader::FillBlockBPDFLCsDISv21(fastNLOCoeffAddFlex* c) {
   debug["FillBlockBPDFLCsDISv21"]<<endl;//<<"CoeffTable = "<<endl;

   if (c->PdfLcMuVar.empty()) {
      error<< "PdfLcMuVar is empty in CoeffTable. Printing and exiting."<<endl;
      c->Print();
      exit(1);
   }

   for (int i=0; i<NObsBin; i++) {
      // speed up! if mu_f is only dependent on one variable, we can safe the loop over the other one
      for (int x=0; x<c->GetNxmax(i); x++) {
         //double xp = c->GetXNode1(i,x);
         double xp = c->GetXNode1(i,x);
	 if (fMuFFunc != kScale1 &&  fMuFFunc != kScale2) {   // that't the standard case!
            for (unsigned int jS1=0; jS1<c->GetNScaleNode1(i); jS1++) {
	       for (unsigned int kS2=0; kS2<c->GetNScaleNode2(i); kS2++) {
		  double muf = CalcMu(kMuF , c->GetScaleNode1(i,jS1) ,  c->GetScaleNode2(i,kS2) , fScaleFacMuF);
		  c->PdfLcMuVar[i][x][jS1][kS2] = CalcPDFLinearCombination(c,GetXFX(xp,muf));
		  //c->PdfLcMuVar[i][x][jS1][kS2] = CalcPDFLinearCombDIS(GetXFX(xp,muf) , c->GetNSubproc() );
               }
            }
         } 
	 else if (fMuFFunc == kScale2) { // speed up
            for (unsigned int kS2=0; kS2<c->GetNScaleNode2(i); kS2++) {
               double muf = CalcMu(kMuF , 0 ,  c->GetScaleNode2(i,kS2) , fScaleFacMuF);
               //vector < double > buffer = CalcPDFLinearCombDIS(GetXFX(xp,muf) , c->GetNSubproc() );
	       vector<double > buffer = CalcPDFLinearCombination(c,GetXFX(xp,muf));
               for (unsigned int jS1=0; jS1<c->GetNScaleNode1(i); jS1++) {
                  c->PdfLcMuVar[i][x][jS1][kS2] = buffer;
               }
            }
         } 
	 else if (fMuFFunc == kScale1) { // speed up
	    for (unsigned int jS1=0; jS1<c->GetNScaleNode1(i); jS1++) {
               double muf = CalcMu(kMuF , c->GetScaleNode1(i,jS1) , 0 , fScaleFacMuF);
               //vector < double > buffer = CalcPDFLinearCombDIS(GetXFX(xp,muf) , c->GetNSubproc() );
	       vector<double > buffer = CalcPDFLinearCombination(c,GetXFX(xp,muf));
                for (unsigned int kS2=0; kS2<c->GetNScaleNode2(i); kS2++) {
		   c->PdfLcMuVar[i][x][jS1][kS2] = buffer;
               }
            }
         }
      }
   }
   debug["FillBlockBPDFLCsDISv21"]<<"done." <<endl;

}


//______________________________________________________________________________


void fastNLOReader::FillBlockBPDFLCsHHCv20(fastNLOCoeffAddFix* c) {
   int scaleVar          = c->GetNpow() == ILOord ? 0 : fScalevar; // Use IScaleDep
   double scalefac       = fScaleFacMuF/c->GetScaleFactor(scaleVar);
   debug["FillBlockBPDFLCsHHCv20"]<<"scalefac="<<scalefac<<endl;

   if ( c->GetNPDFDim() != 1 ) { error["FillBlockBPDFLCsHHCv20"]<<"full matrix notation not implemented."<<endl;exit(1);}
   
   vector < vector < double > > xfx; // PDFs of all partons
   for (int i=0; i<NObsBin; i++) {
      int nxmax = c->GetNxmax(i);
      int nxbins1 = c->GetNxtot1(i); // number of columns in half matrix
      xfx.resize(nxbins1);
      for (int j=0; j<c->GetNScaleNode(); j++) {
	 // determine all pdfs of hadron1
	 for (int k=0; k<nxbins1; k++) {
	    double xp     = c->GetXNode1(i,k);
	    double muf    = scalefac * c->GetScaleNode(i,scaleVar,j);
	    xfx[k]        = GetXFX(xp,muf);
	 }
	 int x1bin = 0;
	 int x2bin = 0;
	 for (int k=0; k<nxmax; k++) {
	    c->PdfLc[i][j][k] = CalcPDFLinearCombination(c,xfx[x2bin],xfx[x1bin]);
	    x1bin++;
	    if (x1bin>x2bin) {
	       x1bin = 0;
	       x2bin++;
	    }
	 }
      }
   }
}


//______________________________________________________________________________


void fastNLOReader::FillBlockBPDFLCsHHCv21(fastNLOCoeffAddFlex* c) {
   debug["FillBlockBPDFLCsHHCv21"]<<endl;
   if (c->PdfLcMuVar.empty()) {
      cout<< "PdfLcMuVar in CoeffTable is not accessible (resized)."<<endl;
      exit(1);
   }
   if ( c->GetNPDFDim() != 1 ) { error["FillBlockBPDFLCsHHCv20"]<<"'full matrix notation not implemented."<<endl;exit(1);}
  
   vector < vector < double > > xfx; // PDFs of all partons
   for (int i=0; i<NObsBin; i++) {
      int nxmax = c->GetNxmax(i);
      int nxbins1 = c->GetNxtot1(i); // number of columns in half matrix
      xfx.resize(nxbins1);

      if (fMuFFunc != kScale1 &&  fMuFFunc != kScale2)  {   // that't the standard case!
         for (unsigned int jS1=0; jS1<c->GetNScaleNode1(i); jS1++) {
            for (unsigned int kS2=0; kS2<c->GetNScaleNode2(i); kS2++) {
               // determine all pdfs of hadron1
               for (int k=0; k<nxbins1; k++) {
                  double muf = CalcMu(kMuF , c->GetScaleNode1(i,jS1) ,  c->GetScaleNode2(i,kS2) , fScaleFacMuF);
                  double xp   = c->GetXNode1(i,k);
                  xfx[k] = GetXFX(xp,muf);
               }
               int x1bin = 0;
               int x2bin = 0;
               for (int x=0; x<nxmax; x++) {
		  // CalcPDFLinearCombination calculats Anti-proton from proton
		  c->PdfLcMuVar[i][x][jS1][kS2] = CalcPDFLinearCombination(c,xfx[x2bin],xfx[x1bin]);
                  x1bin++;
                  if (x1bin>x2bin) {
                     x1bin = 0;
                     x2bin++;
                  }
               }
            }
         }
      } else if (fMuFFunc == kScale2) {   // speed up
	 for (unsigned int kS2=0; kS2<c->GetNScaleNode2(i); kS2++) {
            // determine all pdfs of hadron1
            for (int k=0; k<nxbins1; k++) {
               double muf = CalcMu(kMuF , 0 ,  c->GetScaleNode2(i,kS2) , fScaleFacMuF);
               double xp     = c->GetXNode1(i,k);
               xfx[k] = GetXFX(xp,muf);
            }
	    for (unsigned int jS1=0; jS1<c->GetNScaleNode1(i); jS1++) {
               int x1bin = 0;
               int x2bin = 0;
               for (int x=0; x<nxmax; x++) {
		  c->PdfLcMuVar[i][x][jS1][kS2] = CalcPDFLinearCombination(c,xfx[x2bin],xfx[x1bin]);
                  x1bin++;
                  if (x1bin>x2bin) {
                     x1bin = 0;
                     x2bin++;
                  }
               }
            }
         }
      } else if (fMuFFunc == kScale1) {   // speed up
	 for (unsigned int jS1=0; jS1<c->GetNScaleNode1(i); jS1++) {
            // determine all pdfs of hadron1
            for (int k=0; k<nxbins1; k++) {
               double muf = CalcMu(kMuF , c->GetScaleNode1(i,jS1) , 0 , fScaleFacMuF);
               double xp     = c->GetXNode1(i,k);
               xfx[k] = GetXFX(xp,muf);
            }
	    for (unsigned int kS2=0; kS2<c->GetNScaleNode2(i); kS2++) {
               int x1bin = 0;
               int x2bin = 0;
               for (int x=0; x<nxmax; x++) {
		  c->PdfLcMuVar[i][x][jS1][kS2] = CalcPDFLinearCombination(c,xfx[x2bin],xfx[x1bin]);
                  x1bin++;
                  if (x1bin>x2bin) {
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

void fastNLOReader::SetExternalFuncForMuR(double(*Func)(double,double)) {
   if (!GetIsFlexibleScaleTable()) {
      warn["SetExternalFuncForMuR"]<<"This is not a flexible-scale table and SetExternalFuncForMuR has no impact.\n";
      man<<"Please use a flexible-scale table, if you want to change your scale definition.\n";
      return;
   }

   Fct_MuR = Func;
   SetFunctionalForm(kExtern , kMuR);
   info["SetExternalFuncForMuR"]<<"Testing external function:"<<endl;
   info<<"Scale1 = 1 ,      Scale2 = 1        ->  mu = func(1,1)             = "<<(*Fct_MuR)(1,1)<<endl;
   info<<"Scale1 = 91.1876, Scale2 = 91.1876  ->  mu = func(91.1876,91.1876) = "<<(*Fct_MuR)(91.1876,91.1876)<<endl;
   info<<"Scale1 = 1,       Scale2 = 91.1876  ->  mu = func(1,91.1876)       = "<<(*Fct_MuR)(1,91.1876)<<endl;
   info<<"Scale1 = 91.1876, Scale2 = 1        ->  mu = func(91.1876,1)       = "<<(*Fct_MuR)(91.1876,1)<<endl;
}


//______________________________________________________________________________


void fastNLOReader::SetExternalFuncForMuF(double(*Func)(double,double)) {
   if (!GetIsFlexibleScaleTable()) {
      warn["SetExternalFuncForMuF"]<<"This is not a flexible-scale table and SetExternalFuncForMuF has no impact.\n";
      man<<"Please use a flexible-scale table, if you want to change your scale definition.\n";
      return;
   }

   Fct_MuF = Func;
   SetFunctionalForm(kExtern , kMuF);
   info["SetExternalFuncForMuF"]<<"Testing external function:"<<endl;
   info<<"Scale1 = 1 ,      Scale2 = 1        ->  mu = func(1,1)             = "<<(*Fct_MuF)(1,1)<<endl;
   info<<"Scale1 = 91.1876, Scale2 = 91.1876  ->  mu = func(91.1876,91.1876) = "<<(*Fct_MuF)(91.1876,91.1876)<<endl;
   info<<"Scale1 = 1,       Scale2 = 91.1876  ->  mu = func(1,91.1876)       = "<<(*Fct_MuF)(1,91.1876)<<endl;
   info<<"Scale1 = 91.1876, Scale2 = 1        ->  mu = func(91.1876,1)       = "<<(*Fct_MuF)(91.1876,1)<<endl;
}


//______________________________________________________________________________




void fastNLOReader::SetFunctionalForm(EScaleFunctionalForm func , fastNLO::EMuX MuX) {
   //
   //  For MuVar tables this method sets the functional form of
   //  the renormalization or the factorization scale.
   //     func:  Choose a pre-defined function
   //     kMuX:  is it for mu_r or for mu_f ?
   //

   if (!GetIsFlexibleScaleTable()) {
      warn<<"This is not a flexible-scale table. SetFunctionalForm cannot be used.\n";
      return;
   }

   // ---- setting scale ---- //
   if (MuX == kMuR) fMuRFunc = func;
   else fMuFFunc = func;


   // ---- cross check ---- //
   if (func == kScale2 || func == kQuadraticSum ||  func == kQuadraticMean || func == kQuadraticSumOver4
       || func == kLinearMean || func == kLinearSum  ||  func == kScaleMax|| func == kScaleMin) {
      int nnode = ((fastNLOCoeffAddFlex*)BBlocksSMCalc[0][1])->GetNScaleNode2(0);
      if (nnode <= 3) {
         error<<"There is no second scale variable available in this table. Using fastNLO::kScale1 only.\n";
         SetFunctionalForm(kScale1,MuX);
      }
      for (int i=0; i<NObsBin; i++) {
	 nnode = ((fastNLOCoeffAddFlex*)BBlocksSMCalc[0][1])->GetNScaleNode2(i);
	 if (nnode < 4) {
            warn<<"Scale2 has only very little nodes (n="<<nnode<<") in bin "<<i<<".\n";
         }
      }
   }
   //PrintScaleSettings(MuX);//not yet ported to v2.2
}


//______________________________________________________________________________


void fastNLOReader::SetMuRFunctionalForm(EScaleFunctionalForm func) {
   SetFunctionalForm(func,kMuR);
}
void fastNLOReader::SetMuFFunctionalForm(EScaleFunctionalForm func) {
   SetFunctionalForm(func,kMuF);
}


//______________________________________________________________________________



double fastNLOReader::CalcMu(fastNLO::EMuX kMuX , double scale1, double scale2, double scalefac) {
   //
   //  Calculate the scales with the defined function and the
   //  corresponding prefactor.
   //
   if (kMuX == kMuR && fScaleFacMuR != scalefac) error<<"Sth. went wrong with the scales.\n";
   if (kMuX == kMuF && fScaleFacMuF != scalefac) error<<"Sth. went wrong with the scales.\n";

   EScaleFunctionalForm Func = (kMuX == kMuR) ? fMuRFunc : fMuFFunc;
   double mu = 0;
   if (Func == fastNLO::kScale1)            mu      = scale1;
   else if (Func == fastNLO::kScale2)            mu      = scale2;
   else if (Func == fastNLO::kQuadraticSum)      mu      = FuncMixedOver1(scale1,scale2);
   else if (Func == fastNLO::kQuadraticMean)     mu      = FuncMixedOver2(scale1,scale2);
   else if (Func == fastNLO::kQuadraticSumOver4) mu      = FuncMixedOver4(scale1,scale2);
   else if (Func == fastNLO::kLinearMean)        mu      = FuncLinearMean(scale1,scale2);
   else if (Func == fastNLO::kLinearSum)         mu      = FuncLinearSum(scale1,scale2);
   else if (Func == fastNLO::kScaleMax)          mu      = FuncMax(scale1,scale2);
   else if (Func == fastNLO::kScaleMin)          mu      = FuncMin(scale1,scale2);
   else if (Func == fastNLO::kExpProd2)          mu      = FuncExpProd2(scale1,scale2);
   else if (Func == fastNLO::kExtern)           mu      = (kMuX==kMuR) ? (*Fct_MuR)(scale1,scale2) : (*Fct_MuF)(scale1,scale2);
   else error["CalcMu"]<<"Could not identify functional form for scales calculation.\n";

   return scalefac * mu;
}


//______________________________________________________________________________
double fastNLOReader::FuncMixedOver1(double scale1 , double scale2) {
   return (sqrt((pow(scale1,2) + pow(scale2,2))  / 1.));
}

//______________________________________________________________________________
double fastNLOReader::FuncMixedOver2(double scale1 , double scale2) {
   return (sqrt((pow(scale1,2) + pow(scale2,2))  / 2.));
}

//______________________________________________________________________________
double fastNLOReader::FuncMixedOver4(double scale1 , double scale2) {
   return (sqrt((pow(scale1,2) + pow(scale2,2))  / 4.));
}

//______________________________________________________________________________
double fastNLOReader::FuncLinearMean(double scale1 , double scale2) {
   return (scale1 + scale2) / 2.;
}

//______________________________________________________________________________
double fastNLOReader::FuncLinearSum(double scale1 , double scale2) {
   return scale1 + scale2;
}
//______________________________________________________________________________
double fastNLOReader::FuncMax(double scale1 , double scale2) {
   if (scale1 > scale2) return scale1;
   else return scale2;
}

//______________________________________________________________________________
double fastNLOReader::FuncMin(double scale1 , double scale2) {
   if (scale1 < scale2) return scale1;
   else return scale2;
}

//______________________________________________________________________________
double fastNLOReader::FuncExpProd2(double scale1 , double scale2) {
   return (scale1 * exp(0.3*scale2));
}



//______________________________________________________________________________



int fastNLOReader::ContrId(const ESMCalculation eCalc, const ESMOrder eOrder) const {
   int Id = -1;
   if (BBlocksSMCalc.empty() || bUseSMCalc[eCalc].empty()) {
      return Id;
   }

   // Requested order
   string requested = _OrdName[eCalc][eOrder];
   // Loop over all available orders of contribution type eCalc
   for (unsigned int i=0; i<BBlocksSMCalc[eCalc].size(); i++) {
      // Found order
      int iFlag1 = BBlocksSMCalc[eCalc][i]->GetIContrFlag1();
      int iFlag2 = BBlocksSMCalc[eCalc][i]->GetIContrFlag2();
      string available = _OrdName[iFlag1-1][iFlag2-1];
      if (available == requested) {
         Id = i;
      }
   }
   return Id;
}



//______________________________________________________________________________
//
//                      Print outs
//______________________________________________________________________________


void fastNLOReader::PrintTableInfo(const int iprint) const {
   // this function is inherited from fastNLOTable.
   fastNLOTable::PrintTableInfo(iprint);
}

//______________________________________________________________________________
void fastNLOReader::PrintFastNLOTableConstants(const int iprint) const {
   // this function is inherited from fastNLOTable.
   this->fastNLOTable::PrintFastNLOTableConstants(iprint);
}


//______________________________________________________________________________
void fastNLOReader::PrintCrossSections() const {
   //
   // Print Cross sections in NLO, k-factors and Reference table cross sections
   //

   //   if ( XSection.empty() )    CalcCrossSection();
   //   if ( XSectionRef.empty() && XSectionRef_s1.empty() )    CalcReferenceCrossSection();}

   vector < double > xs = XSection;

   printf(" #  \n");
   printf(" #  FastNLO Cross sections for\n");
   for (unsigned int i = 0 ; i < ScDescript.size() ; i++) {
      printf(" #     %s\n",ScDescript[i].c_str());
   }
   printf(" #  at sqrt(s) = %8.2f GeV\n", Ecms);
   printf(" #  \n");
   printf(" #  This is a %s-differential table in %s", ((NDim==1)?"single":"double"),GetDimLabel(0).c_str());
   if (NDim==2) printf(" and in %s",GetDimLabel(1).c_str());
   printf(".\n");
   printf(" #\n");

   string Aunits[16]  = { "[b] --   ","","","[mb] --  ","","","[mu b] --","","","[nb] --  ","","","[pb] --  ","","","[fb] --  "};
   string Nounits[16] = { " --      ","",""," --      ","",""," --      ","",""," --      ","",""," --      ","",""," --      "};
   string* unit = fUnits==kAbsoluteUnits ? Aunits : Nounits;


   if (NDim == 2) {
      double lobindim2 = -42;
      printf(" #  - Bin - |   ---  %5s  ---        -- XS-FNLO %s -- k-factor -- |\n",GetDimLabel(0).c_str(),unit[Ipublunits].c_str());
      printf(" #  --------------------------------------------------------------------\n");
      for (unsigned int i=0; i<xs.size(); i++) {
         if (GetLoBin(i,1) != lobindim2) {
            printf(" #                  ---->  from %9.3f to %9.3f in %s  <----\n",GetLoBin(i,1),GetUpBin(i,1),GetDimLabel(1).c_str());
            lobindim2 = GetLoBin(i,1);
         }
         printf(" #   %4.0f   | %9.3f - %9.3f       % 9.4e           % 5.2f      |\n",i*1.,GetLoBin(i,0),GetUpBin(i,0),xs[i],kFactor[i]);
      }
   }

   else {
      printf("   ---  %5s  ---        - Bin -       -- XS-FNLO --  \n",GetDimLabel(0).c_str());
      for (unsigned int i=0; i<xs.size(); i++) {
         printf("  %9.3f - %9.3f   %3.0f         % 9.4e\n",GetLoBin(i,0),GetUpBin(i,0),i*1.,xs[i]);
      }
   }
   printf(" #  --------------------------------------------------------------------\n");
}


//______________________________________________________________________________


void fastNLOReader::PrintCrossSectionsWithReference() {
   //
   //  Print Cross sections in NLO, k-factors and Reference table cross sections
   //
   //  Please mention, that the reference cross section can be easily deviating
   //  more than 20% (scales, pdfs, alpha_s, etc...). This does not mean that
   //  the table is wrong!
   //

   vector < double > xs = XSection;
   vector < double > xsref;
   if (XSection.empty())      CalcCrossSection();
   if (XSectionRef.empty() && XSectionRef_s1.empty())      CalcReferenceCrossSection();

   if (GetIsFlexibleScaleTable()) {
      if (fMuFFunc == kScale1 && fMuRFunc == kScale1)   {
         printf(" #  FastNLOReader::PrintCrossSectionsWithReference. Info. Taking reference cross sections 's1'\n");
         xsref = XSectionRef_s1;
      } else if (fMuFFunc == kScale2 && fMuRFunc == kScale2) {
         printf(" #  FastNLOReader::PrintCrossSectionsWithReference. Info. Taking reference cross sections 's2'\n");
         xsref = XSectionRef_s2;
      } else if (fMuFFunc == kQuadraticMean && fMuRFunc == kQuadraticMean) {
         printf(" #  FastNLOReader::PrintCrossSectionsWithReference. Info. Taking reference cross sections 'mixed'\n");
         xsref = XSectionRefMixed;
      } else {
         xsref = XSectionRefMixed;
         printf(" #  FastNLOReader::PrintCrossSectionsWithReference. Info. Taking reference cross sections 'mixed'\n");
      }
   } else xsref = XSectionRef;


   printf(" #  \n");
   printf(" #  FastNLO Cross sections for\n");
   for (unsigned int i = 0 ; i < ScDescript.size() ; i++) {
      printf(" #     %s\n",ScDescript[i].c_str());
   }
   printf(" #  at sqrt(s) = %8.2f GeV\n", Ecms);
   printf(" #  \n");
   printf(" #  This is a %s-differential table in %s", ((NDim==1)?"single":"double"),DimLabel[0].c_str());
   if (NDim==2) printf(" and %s",DimLabel[1].c_str());
   printf(" #  \n");
   printf(" #  Please mention, that the reference cross section can easily deviating up to more\n *  than 20%% due to different scale choices, alhpa_s value/evolution, PDFs, etc.");
   printf(" #  This does not mean, that this FastNLO table is wrong!\n\n");
   printf(" #  There are three reference cross sections stored for different scale choices.\n");
   printf(" #  If you have choosen mu_r=mu_f=%s, or mu_r=mu_f=%s or mu_r=mu_f=sqrt((%s^2+%s^2)/2), then you access automatically the corresponding reference cross section.\n",
          B_NLO()->GetScaleDescription(0).c_str(),B_NLO()->GetScaleDescription(1).c_str(),
	  B_NLO()->GetScaleDescription(0).c_str(),B_NLO()->GetScaleDescription(1).c_str());
   printf(" #  In any other case your reference cross section is calculated using mu_r=mu_f=sqrt((%s^2+%s^2)/2).\n",
	  B_NLO()->GetScaleDescription(0).c_str(),B_NLO()->GetScaleDescription(1).c_str());
   printf(" #  To be fully consistent with the nlojet++ reference cross section, you also have to adjust alpha_s and the alpha_s evolution accordingly.\n\n");
   
   printf("\n");
   printf(" #\n");

   string Aunits[16]  = { "[b] --   ","","","[mb] --  ","","","[mu b] --","","","[nb] --  ","","","[pb] --  ","","","[fb] --  "};
   string Nounits[16] = { " --      ","",""," --      ","",""," --      ","",""," --      ","",""," --      ","",""," --      "};
   string* unit = fUnits==kAbsoluteUnits ? Aunits : Nounits;

   if (NDim == 2) {
      double lobindim2 = -321312;
      printf(" #  - Bin - |   ---  %5s  ---        -- XS-FNLO %s -- k-factor -- |  -- XS-ref (NLOJET++) --    Diff [%%]\n",DimLabel[0].c_str(),unit[Ipublunits].c_str());
      printf(" #  -----------------------------------------------------------------------------------------------------------\n");
      for (unsigned int i=0; i<xs.size(); i++) {
         if (GetLoBin(i,1) != lobindim2) {
            printf(" #                    ---->  from %9.3f to %9.3f in %s  <----\n",GetLoBin(i,1),GetUpBin(i,1),GetDimLabel(1).c_str());
            lobindim2 = GetLoBin(i,1);
         }
         printf(" #   %4.0f   | %9.3f - %9.3f      % 9.4e           % 5.3f      |     % 9.4e            % 5.4f\n",
		i*1.,GetLoBin(i,0),GetUpBin(i,0),xs[i],kFactor[i],xsref[i],(xs[i]-xsref[i])/xsref[i]*100.);
      }
   }

   else {
      printf("FastNLOReader::PrintCrossSections( ). Info. Single differential printing of cross sections not yet nicely implemented.\n");
      printf("   ---  %s  ---        - Bin -    -- XS-FNLO  --       -- XS-ref (NLOJET++) --    Diff [%%]\n",DimLabel[0].c_str());
      for (unsigned int i=0; i<xs.size(); i++) {
         printf("  %9.3f - %9.3f   %3.0f         % 9.4e           % 9.4e          % 5.4f\n",GetLoBin(i,0),GetUpBin(i,0),i*1.,xs[i],xsref[i],(xs[i]-xsref[i])/xsref[i]*100.);
      }
   }
   printf(" #  ------------------------------------------------------------------------------------------------------------\n");
}


//______________________________________________________________________________
void fastNLOReader::PrintCrossSectionsDefault(const vector <double> kthc) const {
   //
   // Print observable binnning and cross sections at
   // LO, NLO and K factors like in Fortran Reader for comparison
   //
   // This function can only print double-differential crosss section tables.
   //

   // Check on existence of 2-loop threshold corrections
   //const int ithc2 = kthc.empty() ? -1 : ContrId( fastNLO::kThresholdCorrection, fastNLO::kNextToLeading);
   const int ithc2 = kthc.empty() ? -1 : ContrId(kThresholdCorrection,kNextToLeading);

   cout << _DSEP << endl;
   printf(" Cross Sections\n");
   if (!GetIsFlexibleScaleTable())
      printf(" The scale chosen here are: mu_f = % #6.3f * %s, and mu_r = % #6.3f * %s \n",
	     fScaleFacMuF, B_LO()->GetScaleDescription().c_str(), fScaleFacMuR, B_LO()->GetScaleDescription().c_str());
   cout << _SSEP << endl;

   if (NDim == 2) {
      // non-perturbative corrections (just first np correction)
      const int inpc1 = ContrId(kNonPerturbativeCorrection,kLeading);
      const vector < double > knpc = inpc1>-1 ? ((fastNLOCoeffMult*)BBlocksSMCalc[kNonPerturbativeCorrection][kLeading])->GetMultFactor() : vector<double>(NObsBin);

      string header0 = "  IObs  Bin Size IODim1 ";
      string header1 = "   IODim2 ";
      string header2 = " LO cross section   NLO cross section   K NLO";
      if (ithc2>-1)header2 += "     K THC";
      if (inpc1>-1)header2 += "     K NPC";
      unsigned int NDimBins[NDim];
      printf("%s [ %-12s ] %s [  %-12s  ] %s\n",
             header0.c_str(),DimLabel[0].c_str(),header1.c_str(),DimLabel[1].c_str(),header2.c_str());
      cout << _SSEP << endl;
      for (int i=0; i<NObsBin; i++) {
         for (int j=0; j<NDim; j++) {
            if (i==0)                                  NDimBins[j] = 1;
            else if (GetLoBin(i-1,j) < GetLoBin(i,j))       NDimBins[j]++;
            else if (GetLoBin(i,j) < GetLoBin(i-1,j))       NDimBins[j] = 1;
         }
         if (ithc2<0 && inpc1<0) {
            printf(" %5.i % -#10.4g %5.i % -#10.4g % -#10.4g %5.i  %-#8.2E  %-#8.2E %#18.11E %#18.11E %#9.5F",
                   i+1,BinSize[i],NDimBins[0],GetLoBin(i,0),GetUpBin(i,0),
                   NDimBins[1],GetLoBin(i,1),GetUpBin(i,1),XSection_LO[i],XSection[i],kFactor[i]);
         } else if (inpc1<0 && ithc2 != -1) {
            printf(" %5.i % -#10.4g %5.i % -#10.4g % -#10.4g %5.i  %-#8.2E  %-#8.2E %#18.11E %#18.11E %#9.5F %#9.5F",
                   i+1,BinSize[i],NDimBins[0],GetLoBin(i,0),GetUpBin(i,0),
                   NDimBins[1],GetLoBin(i,1),GetUpBin(i,1),XSection_LO[i],XSection[i],kFactor[i],kthc[i]);
         } else if (inpc1>-1 && ithc2 == -1) {
            printf(" %5.i % -#10.4g %5.i % -#10.4g % -#10.4g %5.i  %-#8.2E  %-#8.2E %#18.11E %#18.11E %#9.5F %#9.5F",
                   i+1,BinSize[i],NDimBins[0],GetLoBin(i,0),GetUpBin(i,0),
                   NDimBins[1],GetLoBin(i,1),GetUpBin(i,1),XSection_LO[i],XSection[i],kFactor[i],knpc[i]);
         } else {
            printf(" %5.i % -#10.4g %5.i % -#10.4g % -#10.4g %5.i  %-#8.2E  %-#8.2E %#18.11E %#18.11E %#9.5F %#9.5F %#9.5F",
                   i+1,BinSize[i],NDimBins[0],GetLoBin(i,0),GetUpBin(i,0),
                   NDimBins[1],GetLoBin(i,1),GetUpBin(i,1),XSection_LO[i],XSection[i],kFactor[i],kthc[i],knpc[i]);
         }
         printf("\n");
      }
   } else {
      warn["PrintCrossSectionsDefault"]<<"Print out optimized for two dimensions. No output for "<<NDim<<" dimensions."<<endl;
   }

}


//______________________________________________________________________________
void fastNLOReader::RunFastNLODemo() {
   //
   // This method prints out cross sections for different scale
   // variation tables. Though it also changes the currently stored
   // settings of this instance!
   //
   // PrintFastNLODemo is changing settings (like scale choices) of this reader.
   // 

   info["PrintFastNLODemo"]<<"PrintFastNLODemo is changing settings (like scale choices) of this reader."<<endl;

   // If flexible-scale table, set MuR and MuF functional forms
   if (GetIsFlexibleScaleTable()) {
      SetMuRFunctionalForm(fastNLO::kScale1);
      SetMuFFunctionalForm(fastNLO::kScale1);
      //SetMuRFunctionalForm(fastNLO::kExpProd2);
      //SetMuRFunctionalForm(fastNLO::kExpProd2);
   }

   // Check on existence of LO and NLO (Id = -1 if not existing)
   int ilo   = ContrId(fastNLO::kFixedOrder, fastNLO::kLeading);
   int inlo  = ContrId(fastNLO::kFixedOrder, fastNLO::kNextToLeading);
   if (ilo < 0 || inlo < 0) {
      error["PrintFastNLODemo"]<<"LO and/or NLO not found, nothing to be done!\n";
      return;
   }
   // Check on existence of 2-loop threshold corrections
   int ithc2 = ContrId(fastNLO::kThresholdCorrection, fastNLO::kNextToLeading);
   // Switched off by default. Don't do scale variations. Not available for the moment.
   //  if ( ithc2 > -1 ) {
   //    SetContributionON( fastNLO::kThresholdCorrection, ithc2, false, false );
   //  }

   // Pre-define desired order of scale variations
   const int nxmu = 4;
   double xmu[nxmu] = {1.0, 0.25, 0.5, 2.0};
   int   ixmu[nxmu] = { -1,   -1,  -1,  -1};
   // Get number of available scale variations and check on available scale factors,
   // in particular for MuF; set pointers
   int nscls = GetNScaleVariations();
   // With threshold corrections, allow only default scale (0)
   if (ithc2 > -1) {
      nscls = 1;
   }
   for (int iscls=0; iscls<nscls; iscls++) {
      SetScaleVariation(iscls);
      double fxmu = fScaleFacMuF;
      for (int i=0; i<nxmu; i++) {
         if (fabs(xmu[i]-fxmu) < 0.000001) {
            ixmu[i] = iscls;
         }
      }
   }

   // Loop over scales
   for (int iscls=0; iscls<nxmu; iscls++) {
      // First result is with NLO, LO result via division by K factor
      if (ixmu[iscls] > -1) {
         SetScaleVariation(ixmu[iscls]);
         CalcCrossSection();

         // Second result: Include threshold corrections for NLO if available
         vector < double > kthc;
         if (ithc2 > -1) {
            vector < double > stdk = kFactor;
            SetContributionON(fastNLO::kThresholdCorrection, ithc2, true);
            CalcCrossSection();
            kthc = kFactor;
            // Threshold K factor is NLO including 2-loop vs. NLO
            for (unsigned int i=0; i<kthc.size(); i++) {
               if (fabs(kFactor[i]) > DBL_MIN) {
                  kthc[i] = kFactor[i]/stdk[i];
               } else {
                  kthc[i] = -1.;
               }
            }
            SetContributionON(fastNLO::kThresholdCorrection, ithc2, false);
         }

         PrintCrossSectionsDefault(kthc);
      }
   }
}


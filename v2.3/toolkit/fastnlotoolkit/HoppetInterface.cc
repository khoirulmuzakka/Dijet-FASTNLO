#include <cstdlib>
#include <LHAPDF/LHAPDF.h>
#include "fastnlotk/fastNLOReader.h"
#include "fastnlotk/speaker.h"
#include "fastnlotk/HoppetInterface.h"
#include "hoppet_v1.h"

// Initial values
bool HoppetInterface::IsInitialized = false;
fastNLOReader *HoppetInterface::fnlo = NULL;
// PDG
double HoppetInterface::QMass[6] = {PDG_MD, PDG_MU, PDG_MS, PDG_MC, PDG_MB, PDG_MT};
double HoppetInterface::fMz = PDG_MZ;
double HoppetInterface::fAlphasMz = PDG_ASMZ;
// VFNS & NLO
std::string HoppetInterface:: fnScheme = "variable";
int HoppetInterface::fnFlavor = 0;
int HoppetInterface::fnLoop = 2;



void HoppetInterface::InitHoppet(fastNLOReader& lfnlo) {
   // Reinitialise also if fnLoop has been changed from initial value
   if ( ! IsInitialized || HoppetInterface::fnLoop != 2 ) {
      StartHoppet();
      fnlo = &lfnlo;
   }

   // If fnFlavor < 3 use VFNS up to Nf_max = 5 as default)
   if ( HoppetInterface::fnFlavor < 3 ) {
      hoppetsetpolemassvfn_(QMass[3], QMass[4], 1.E10);
      say::info["InitHoppet"] << "Using variable-flavour number scheme with Nfmax = 5. "
                              << "M_c, M_b are set to PDG values, and M_t to 10^10 GeV." <<  std::endl;
   } else {
      hoppetsetffn_(HoppetInterface::fnFlavor);
   }
   say::info["InitHoppet"] << "Using alpha_S(Q) = " << fAlphasMz << " at Q = " << fMz << " for " <<
      HoppetInterface::fnLoop << " loops in alpha_S evolution and Nf = " <<
      HoppetInterface::fnFlavor << " flavors." << std::endl;
   // Carry out evolution defining alpha_S(Q), Q, nLoops, muR/muF, LHAsub, starting scale Q_PDF
   hoppetEvolve(HoppetInterface::fAlphasMz, HoppetInterface::fMz, HoppetInterface::fnLoop, 1.0, &HoppetInterface::LHAsub, 2.00001);
   // Fills the HOPPET PDF represenation using PDF provided by LHAPDF.
   //hoppetAssign(&evolvepdf_);
}

void HoppetInterface::StartHoppet() {
   //   std::cout << "StartHoppet: Start. IsInitialized = " << IsInitialized << std::endl;
   // Reasonable defaults for initialisation of PDF evolution
   double ymax   = 12.0;
   double dy     = 0.1;
   double Qmin   = 1.0;
   double Qmax   = 28000.;
   double dlnlnQ = dy/4.0;
   int numorder  = -6; // Note: Has nothing to do with perturbative order!
   int fnLoop    = 2;  // NLO for a start; reinitialise with HoppetInterface::fnLoop if called again
   if ( IsInitialized ) {
      fnLoop = HoppetInterface::fnLoop;
   }
   //   std::cout << "StartHoppet: fnLoop = " << fnLoop << std::endl;
   //   std::cout << "StartHoppet: HI::fnLoop = " << HoppetInterface::fnLoop << std::endl;
   hoppetStartExtended( ymax, dy, Qmin, Qmax, dlnlnQ, fnLoop, numorder, factscheme_MSbar);
   IsInitialized = true;
   //   std::cout << "StartHoppet: Ende. IsInitialized = " << IsInitialized << std::endl;
}

void HoppetInterface::LHAsub(const double & x, const double & Q, double * pdf) {
   //
   //Provides PDF for Hoppet
   for (int i=0; i<13; i++)
   {
   pdf[i] = HoppetInterface::fnlo->GetXFX(x,Q)[i];
   }
}


std::vector<double> HoppetInterface::GetSpl(double x, double Q){
   if (!IsInitialized) {
      say::error["GetSpl"] << "Hoppet not correctly initialized!" <<  std::endl;
      std::exit(1);
   }

   //! Returns the splitting functions
   static std::vector<double> xfx(13);
   hoppetEvalSplit(x, Q, 1, 5, &xfx[0]);
   return xfx;
}

std::vector<double> HoppetInterface::GetXFX(double x, double Q){
   //! Returns PDF
   if (!IsInitialized) {
      say::error["GetSpl"] << "Hoppet not correctly initialized!" <<  std::endl;
      std::exit(1);
   }
   static std::vector<double> xfx(13);
   hoppetEval(x, Q, &xfx[0]);
   return xfx;
}

double HoppetInterface::EvolveAlphas(double Q) {
   if (!IsInitialized) {
      say::error["EvolveAlphas"] << "Hoppet not correctly initialized!" <<  std::endl;
      std::exit(1);
   }
   return hoppetAlphaS(Q);
}

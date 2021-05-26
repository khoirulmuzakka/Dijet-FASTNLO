//!
//! UsefulNlojetTools
//!
//! Collection of useful functions for the fastNLO interface
//! to NLOJET++.
//!
// KR: Drop differentiation in 2jet and 3jet ProcConsts.
//     Default is 2jet, set differences for 3jet via steering!

#include <cstdio>
#include <algorithm> //c++98
#include <utility>   //c++11
#include "fastnlotk/fastNLOEvent.h"
#include "fastnlotk/fastNLOTableConstants.h"
#include "fastnlotk/speaker.h"
#include "fnlo_int_nlojet/pdf-hhc-dummy.h"
#include <bits/hhc-phasespace.h>
#include <bits/hhc-jetfunc.h>

namespace UsefulNlojetTools {
/**
   namespace UsefulNlojetTools

   Collection of useful functions and constant for the interface
   between nlojet++ and fastNLO, if nlojet++ is run in hhc-mode
   (i.e. for pp and ppbar collisions).
*/

//_______________________________________________________________________
fastNLO::GeneratorConstants GenConsts();

//_______________________________________________________________________
fastNLO::ProcessConstants ProcConsts_HHC();

//_______________________________________________________________________
unsigned int GetNj();

// get center of mass energy
double GetEcms();

//_______________________________________________________________________
unsigned int GetLoOrder();

//_______________________________________________________________________
unsigned int GetOrderOfRun(const std::basic_string<char> &__file_name);

//_______________________________________________________________________
int NlojetToFastnloIDHHC(int id);
int FastnloIdToNlojetIdHHC(int id);

//_______________________________________________________________________
std::vector<fnloEvent>
GetFlexibleScaleNlojetContribHHC(const nlo::event_hhc &p,
                                 const nlo::amplitude_hhc &amp);

// KR: Is it possible to get NSubproc from the Toolkit instead?
//     This would allow to read the setting as defined in the steering in
//     contrast to this duplication ...
//     In addition, this method had set the # of subprocesses for LO 3-jet to 6
//     which was never tried before ...
//     ==> Fixed back to 7 now.
//     New default: NSubProcesses always 7!
//_______________________________________________________________________
unsigned int GetNSubproc();

//_______________________________________________________________________
std::vector<fnloEvent>
GetFixedScaleNlojetContribHHC(const nlo::event_hhc &p,
                              const nlo::amplitude_hhc &amp, double mu);

//_______________________________________________________________________
std::vector<std::vector<fnloEvent> >
GetFixedScaleNlojetContribHHC(const nlo::event_hhc &p,
                              const nlo::amplitude_hhc &amp, double mu,
                              const std::vector<double> &scalevar);
}

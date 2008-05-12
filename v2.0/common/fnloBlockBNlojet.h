#include "fnloBlockB.h"
#include <bits/photo-process.h>

class fnloBlockBNlojet : public fnloBlockB {
 public:
   fnloBlockBNlojet(fnloBlockA1 *blocka1, fnloBlockA2 *blocka2) :fnloBlockB(blocka1,blocka2){;}
   void FillEventPhoto(int ObsBin,double x, double scale1, const nlo::amplitude_photo& amp, nlo::pdf_and_coupling_photo& pdf);

};


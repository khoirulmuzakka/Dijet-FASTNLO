#ifndef __pdf_hhc_h__
#define __pdf_hhc_h__ 1

#include <algorithm>
#include <bits/hhc-process.h>
#include <stdio.h>

class pdf_hhc_dummy : public nlo::pdf_and_coupling_hhc {
public:
  //   constructor
  pdf_hhc_dummy(unsigned int mem = 0) { ; }

  //   strong coupling
  double alpha_qcd(unsigned int nf, double mr2) { return 1.0; }

  //   parton in lepton distribution function
  void hadronA(double x, double Q2, unsigned int, unsigned int, double *f) {
    printf("pdf_hhc_dummy::hadronA : PDF in dummy should never be called.\n");
    exit(1);
  }

  //   the parton in the hadron distribution function
  void hadronB(double x, double Q2, unsigned int, unsigned int, double *f) {
    printf("pdf_hhc_dummy::hadronB : PDF in dummy should never be called.\n");
    exit(1);
  }

  nlo::weight_hhc pdf(double x1, double x2, double mf2, unsigned int nu,
                      unsigned int nd);
};

#endif

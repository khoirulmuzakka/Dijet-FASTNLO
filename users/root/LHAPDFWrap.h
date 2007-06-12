#ifndef _LHAPDFWRAP_H
#define _LHAPDFWRAP_H

#include <vector>

/*
  This class is a wrapper around the LHAPDFv2/3 package for parton
  distribution functions of the proton.  
  Adapted for LHAPDFv4 by Mike Whalley
 */

extern "C" {
  void finitpdfset_(char&);
  void finitpdf_(int&);
  void fevolvepdf_(double&, double &, double*);
  void fevolvepdfp_(double&, double &, double&, int&, double*);
  void fnumberpdf_(int&);
  void falphaspdf_(double&, double &);
  void fgetorderpdf_(int&);
  void fgetorderas_(int&);
  void fgetdesc_();
  void fgetqmass_(int&, double&);
  void fgetthreshold_(int&, double&);
  void fgetnf_(int&);
  void fgetlam4_(int&, double&);
  void fgetlam5_(int&, double&);
}

class LHAPDFWrap {

public:
  LHAPDFWrap();
  ~LHAPDFWrap();
  // std c'tor d'tor
  LHAPDFWrap(char *name, int memset);
  // typical constructor with pdfset 'name' and subset 'memset' see
  // LHAPDFv2/PDFset for available sets.  'name' is the full path to
  // the grid or data file of the desired set.

  std::vector< double > xfx(const double &x, const double &Q);
  // returns a vector xf(x, Q) with index 0 < i < 12.
  // 0..5 = tbar, ..., ubar, dbar; 
  // 6 = g; 
  // 7..12 = d, u, ..., t
  double xfx(const double &x, const double &Q, int fl);
  // returns xf(x, Q) for flavour fl - this time the flavour encoding
  // is as in the LHAPDF manual...
  // -6..-1 = tbar,...,ubar, dbar
  // 1..6 = duscbt
  // 0 = g
  std::vector< double > xfxp(const double &x, const double &Q2, const double &P2, int ip);
  double xfxp(const double &x, const double &Q2, const double &P2, int ip, int fl);

  inline void initPDFSet(char *name) {finitpdfset_(*name);}
  // the pdfset by name, see subdir 'PDFset' of LHAPDFv2 for choices
  inline void initPDF(int memset) {finitpdf_(memset);}
  // the choice of pdf subset out of one distribution

  inline void getDescription() {fgetdesc_();}
  // prints a brief description of the current pdf set to stdout
  int numberPDF();
  // number of subsets available in the current distribution.
  double alphasPDF(double Q);
  // alphas used by the current pdf.
  int getOrderPDF();
  int getOrderAlphaS();
  // perturbative order of parton evolution and alpha_S resp.
  double getQMass(int f);
  // quark mass used for flavour f.
  double getThreshold(int f);
  // threshold for flavour f.
  int getNf();
  // number of flavours used in the current pdf set.
  double getLam4(int m);
  // value of qcd lambda4 for member m 
  double getLam5(int m);
  // value of qcd lambda5 for member m 
};

#endif

#include "LHAPDFWrap.h"
// adapted for LHAPDFv4 by Mike Whalley

LHAPDFWrap::LHAPDFWrap() {}

LHAPDFWrap::~LHAPDFWrap() {}

LHAPDFWrap::LHAPDFWrap(char *name, int memset) {
  finitpdfset_(*name);
  finitpdf_(memset);
}

std::vector<double> LHAPDFWrap::xfx(const double &x, const double &Q) {  
  double f[13], mx = x, mQ = Q;
  fevolvepdf_(mx, mQ, f);
  std::vector<double> r;
  for (int i=0; i<13; i++) r.push_back(f[i]);
  return r;
}

double LHAPDFWrap::xfx(const double &x, const double &Q, int fl) {  
  double f[13], mx = x, mQ = Q;
  fevolvepdf_(mx, mQ, f);
  return f[fl+6];
}

std::vector<double> LHAPDFWrap::xfxp(const double &x, const double &Q2, const double &P2, int ip) {  
  double f[13], mx = x, mQ2 = Q2, mP2 = P2;
  int  mip = ip;
  fevolvepdfp_(mx, mQ2, mP2, mip, f);
  std::vector<double> r;
  for (int i=0; i<13; i++) r.push_back(f[i]);
  return r;
}

double LHAPDFWrap::xfxp(const double &x, const double &Q2, const double &P2, int ip, int fl) {  
  double f[13], mx = x, mQ2 = Q2, mP2 = P2;
  int mip = ip;
  fevolvepdfp_(mx, mQ2, mP2, mip, f);
  return f[fl+6];
}

int LHAPDFWrap::numberPDF() {
  int N; 
  fnumberpdf_(N); 
  return N;
}

double LHAPDFWrap::alphasPDF(double Q) {
  double a;
  falphaspdf_(Q, a);
  return a;
}

int LHAPDFWrap::getOrderPDF() {
  int N;
  fgetorderpdf_(N);
  return N;
}

int LHAPDFWrap::getOrderAlphaS() {
  int N;
  fgetorderas_(N);
  return N;
}

double LHAPDFWrap::getQMass(int f) {
  double m;
  fgetqmass_(f, m);
  return m;
}

double LHAPDFWrap::getThreshold(int f) {
  double m;
  fgetthreshold_(f, m);
  return m;
}

int LHAPDFWrap::getNf() {
  int N; 
  fgetnf_(N); 
  return N;
}

double LHAPDFWrap::getLam4(int m){
  double l;
  fgetlam4_(m, l);
  return l;
}

double LHAPDFWrap::getLam5(int m){
  double l;
  fgetlam5_(m, l);
  return l;
}

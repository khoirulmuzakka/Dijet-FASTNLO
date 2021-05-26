//  Copyright (C) 2002 Zoltan Nagy
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#include "hhc2jet.h"
#include "ampg4.h"
#include "ampg5.h"
#include "ampq2g2.h"
#include "ampq2g3.h"
#include "ampq4.h"
#include "ampq4g1.h"
#include "bits/nlo-color.h"


//    PI_FACn = (2*pi)^(n-2) * 2^(n-2)
#define PI_FAC4 6234.18182617615598311119
#define PI_FAC5 492231.26711055587175698733



namespace nlo {

  hhc2jet::dipole_func hhc2jet::_S_dipole[9] = 
  { &hhc2jet::_M_di1, &hhc2jet::_M_di2, &hhc2jet::_M_di3, 
    &hhc2jet::_M_d01, &hhc2jet::_M_d02, &hhc2jet::_M_d03,
    &hhc2jet::_M_d12, &hhc2jet::_M_d13, &hhc2jet::_M_d23};

 
  hhc2jet::hhc2jet(const random_generator& rng, bool mchel, unsigned int nu, unsigned int nd, double al)
    : process_hhc(2, 2, nu, nd, al), _hhc_jet_base(nu+nd), _M_mchel(mchel) 
  {
    _M_g4   = new ampg4  (_M_ip, rng);
    _M_q2g2 = new ampq2g2(_M_ip, rng);
    _M_q4   = new ampq4  (_M_ip, rng);
			  
    _M_g5   = new ampg5  (_M_ip, rng);
    _M_q2g3 = new ampq2g3(_M_ip, rng);
    _M_q4g1 = new ampq4g1(_M_ip, rng);
  }
  
  hhc2jet::~hhc2jet() 
  {
    if(_M_g5)   delete _M_g5;
    if(_M_q2g3) delete _M_q2g3;
    if(_M_q4g1) delete _M_q4g1;
    if(_M_g4)   delete _M_g4;
    if(_M_q2g2) delete _M_q2g2;
    if(_M_q4)   delete _M_q4;
  }
  
  void hhc2jet::born_term(const event_type& p, weight_type& res) 
  {
    _M_ip.calculate(p);
    amp_tree(_M_g4, _M_q2g2, _M_q4, 0, res.begin());
    res *= PI_FAC4;

    for (int ikr = 0; ikr<7; ikr++) {
       if (isnan(res[ikr]) || isinf(res[ikr])) {
          std::cout << "NLOJet++_hhc2jet_born_term: NaN or Inf ERROR for ikr = " << ikr << ", res[ikr] = " << res[ikr] << std::endl;
       }
    }

  }
  
  void hhc2jet::real_term(const event_type& p, weight_type& res) 
  {
    _M_ip.calculate(p);
    amp_tree(_M_g5, _M_q2g3, _M_q4g1, 0, res.begin());
    res *= PI_FAC5;

    for (int ikr = 0; ikr<7; ikr++) {
       if (isnan(res[ikr]) || isinf(res[ikr])) {
          std::cout << "NLOJet++_hhc2jet_real_term: NaN or Inf ERROR for ikr = " << ikr << ", res[ikr] = " << res[ikr] << std::endl;
       }
    }

  }


  void hhc2jet::fini_term(double x1, double xjac1, double x2, double xjac2, 
			  const event_type& p, weight_type *res)
  {
    double al = alpha();
    static double loop[7];
    static su3_kp_i2 kp[7];

    _M_ip.calculate(p);
    amp_kp(al, _M_g4, _M_q2g2, _M_q4, kp);

    if(_M_mchel) amp_1loop_mch(_M_g4, _M_q2g2, _M_q4, loop);
    else amp_1loop(_M_g4, _M_q2g2, _M_q4, loop);

    double s = p[hadron(-1)]*p[hadron(0)];
    double e1 = (p[-1]*p[hadron( 0)])/s;
    double e2 = (p[ 0]*p[hadron(-1)])/s;

    //---- convolutions ----
    convolutions(e1, x1, xjac1, e2, x2, xjac2, al, kp, res);
    
    for(unsigned int i=0; i < 7; i++) {
      //---- 1-loop contributions ----
      res[2][i] += kp[i].loop + loop[i];
      
      //---- renormalization scale dependent term ----
      res[6][i] = 2.0*kp[i].tree*Gg(Nf);
    }

    //---- overall Pi factors ----
    for(unsigned int i=0; i < 7; i++)
      res[i] *= PI_FAC4;

    for (int ikr = 0; ikr<7; ikr++) {
       for (int jkr = 0; jkr<7; jkr++) {
          if (isnan(res[ikr][jkr]) || isinf(res[ikr][jkr])) {
             std::cout << "NLOJet++_hhc2jet_fini_term: NaN or Inf ERROR for ikr = " << ikr << ", jkr =" << jkr << ", res[ikr][jkr] = " << res[ikr][jkr] << std::endl;
          }
       }
    }

  }
  
  void hhc2jet::dipole_term(const event_type& p, const event_type& dp,
			    int i, int j, int k, weight_type& res) 
  {
    typedef split_fin<lorentzvector<double> > _SplitF;
    typedef split_ini<lorentzvector<double> > _SplitI;
    
    if(i <= 0) {      
      _M_sini = (k <= 0 ? (_SplitI *) &_M_sifi : (_SplitI *) &_M_siff);
      _M_sini -> set(p[i], p[j], p[k]);
    } else {
      _M_sfin = (k <= 0 ? (_SplitF *) &_M_sffi : (_SplitF *) &_M_sfff);
      _M_sfin -> set(p[i], p[j], p[k]);
    }
    
    int kt = (k == 3 ? j : k);
    int idx = (i==-1 ? j-1 : 2*i-(i*i-i)/2 + j+2);
    
    _M_ip.calculate(dp); 
    (this ->* _S_dipole[idx])(kt, i, res);
    res *= PI_FAC5;

    for (int ikr = 0; ikr<7; ikr++) {
       if (isnan(res[ikr]) || isinf(res[ikr])) {
          std::cout << "NLOJet++_hhc2jet_dipole_term: NaN or Inf ERROR for ikr = " << ikr << ", res[ikr] = " << res[ikr] << std::endl;
          std::cout << "NLOJet++_hhc2jet_dipole_term: kt = " << kt << ", idx = " << idx << std::endl;
       }
    }

  }

  //
  //             dipole functions
  //
#define HHC_FQG Vqg = (_M_sfin -> Vqg()).first
#define HHC_FQA Vqa = (_M_sfin -> Vqa()).first
#define HHC_FGG Vgg = (_M_sfin -> Vgg()).first
                   
#define HHC_IQG Vqg = (_M_sini -> Vqg()).first
#define HHC_IQQ Vqq = (_M_sini -> Vqq()).first
#define HHC_IGA Vga = (_M_sini -> Vga()).first
#define HHC_IGG Vgg = (_M_sini -> Vgg()).first

#define HHC_CCGG(amp,p1,p2) amp_ccgg(amp, kt, i, p1, p2, cc) 
#define HHC_CCQG(amp,p1,p2) amp_ccqg(amp, kt, i, p1, p2, cc) 
#define HHC_CCGQ(amp,p1,p2) amp_ccgq(amp, kt, i, p1, p2, cc) 
#define HHC_CCAG(amp,p1,p2) amp_ccag(amp, kt, i, p1, p2, cc) 
#define HHC_CCGA(amp,p1,p2) amp_ccga(amp, kt, i, p1, p2, cc) 
#define HHC_CCQA(amp,p1,p2) amp_ccqa(amp, kt, i, p1, p2, cc) 
#define HHC_CCAQ(amp,p1,p2) amp_ccaq(amp, kt, i, p1, p2, cc) 

#define HHC_CCQ4(p1,p2,p3,p4) amp_cc(_M_q4, kt, i, p1, p2, p3, p4, q4) 



  void hhc2jet::_M_di1(int kt, int i, weight_type& d) 
  {
    double cc[7], q4[2], HHC_IGG, HHC_IGA, HHC_IQQ, HHC_IQG;
    
    HHC_CCGG(_M_g4,1,2);
    HHC_CCAG(_M_q2g2,2,1);
    HHC_CCGA(_M_q2g2,2,1);
    HHC_CCQA(_M_q2g2,1,2);
    HHC_CCQ4(1,0,-1,2);
    d[0] = Vgg*cc[0]/6.0 + Nf*Vga*cc[1];
    d[1] = Vqq*cc[0]/2.0;
    d[2] = (Nf-1)*Vga*q4[0] + Vga*q4[1]/2.0;
    d[3] = Vqq*cc[2];
    d[4] = Vqq*cc[2]/2.0;
    d[5] = Vqg*cc[5]/6.0 + Vqq*cc[2];
    d[6] = Vqq*cc[2];

    HHC_CCGG(_M_q2g2,2,1);
    HHC_CCAQ(_M_q2g2,1,2);
    d[1] += Vqq*cc[0]/2.0;
    d[2] += Vga*cc[5]/2.0; 
  }
  
  void hhc2jet::_M_di2(int kt, int i, weight_type& d) 
  {
    double cc[7], q4[2], HHC_IGG, HHC_IGA, HHC_IQG, HHC_IQQ;
    
    HHC_CCGG(_M_g4,1,2);
    HHC_CCQG(_M_q2g2,1,2);
    HHC_CCGQ(_M_q2g2,1,2);
    HHC_CCQA(_M_q2g2,1,2);
    HHC_CCQ4(2,0,1,-1);
    d[0] = Vgg*cc[0]/6.0 + Nf*Vga*cc[1];
    d[1] = Vqg*cc[1]/2.0;
    d[2] = Vgg*cc[2]/2.0 + (Nf-1)*Vga*q4[0] + Vga*q4[1]/2.0;
    d[3] = 0.0;
    d[4] = Vqq*cc[2]/2.0;
    d[5] = Vqg*cc[5]/6.0;
    d[6] = 0.0;
  }
  
  void hhc2jet::_M_di3(int kt, int i, weight_type& d) 
  {
    double cc[7], q4[2], HHC_IGG, HHC_IQG, HHC_IQQ, HHC_IGA;
    
    HHC_CCGG(_M_g4,1,2);
    HHC_CCQA(_M_q2g2,1,2);
    HHC_CCQG(_M_q2g2,1,2);
    HHC_CCGQ(_M_q2g2,1,2);
    HHC_CCQ4(-1,0,1,2);
    d[0] = Vgg*cc[0]/6.0;
    d[1] = Vqg*cc[1]/2.0;
    d[2] = Vgg*cc[2]/2.0 + (Nf-1)*Vga*q4[0] + Vga*q4[1]/2.0;
    d[5] = Vqg*cc[5]/6.0;

    HHC_CCGG(_M_q2g2,1,2);
    HHC_CCQ4(0,-1,1,2);
    d[0] += Nf*Vgg*cc[0];
    d[1] += (Nf-0.5)*Vqq*cc[0];
    d[5] += (Nf-1)*Vqg*q4[0] + Vqg*q4[1];

    HHC_CCQ4(1,-1,2,0);
    d[3] = Vqg*q4[0];
    d[4] = Vqg*q4[1]/2.0;

    HHC_CCQ4(1,-1,0,2);
    d[6] = Vqg*q4[0];
  }
  
  void hhc2jet::_M_d01(int kt, int i, weight_type& d) 
  {
    double cc[7], q4[2], HHC_IGG, HHC_IGA, HHC_IQG, HHC_IQQ;
    
    HHC_CCGG(_M_g4,1,2);
    HHC_CCQG(_M_q2g2,2,1);
    HHC_CCGA(_M_q2g2,2,1);
    HHC_CCQA(_M_q2g2,2,1);
    HHC_CCQ4(1,-1,0,2);
    d[0] = Vgg*cc[0]/6.0 + Nf*Vga*cc[2];
    d[1] = Vga*(cc[5]/2.0 + (Nf-1)*q4[0] + q4[1]/2.0);
    d[2] = Vqq*cc[0]/2.0;
    d[3] = 0.0;
    d[4] = Vqq*cc[1]/2.0;
    d[5] = Vqg*cc[5]/6.0;
    d[6] = 0.0;

    HHC_CCGG(_M_q2g2,1,2);
    d[2] += Vqq*cc[0]/2.0;
  }
  
  void hhc2jet::_M_d02(int kt, int i, weight_type& d) 
  {
    double cc[7], q4[2], HHC_IGG, HHC_IGA, HHC_IQG, HHC_IQQ;
    
    HHC_CCGG(_M_g4,1,2);
    HHC_CCQG(_M_q2g2,1,2);
    HHC_CCGQ(_M_q2g2,1,2);
    HHC_CCQA(_M_q2g2,1,2);
    HHC_CCQ4(2,-1,1,0);
    d[0] = Vgg*cc[0]/6.0 + Nf*Vga*cc[2];
    d[1] = Vgg*cc[1]/2.0 + Vga*((Nf-1)*q4[0] + q4[1]/2.0);
    d[2] = Vqg*cc[2]/2.0;
    d[3] = Vqq*cc[1];
    d[4] = Vqq*cc[1]/2.0;
    d[5] = Vqg*cc[5]/6.0 + Vqq*cc[1];
    d[6] = Vqq*cc[1];
  }
  
  void hhc2jet::_M_d03(int kt, int i, weight_type& d) 
  {
    double cc[7], q4[2], HHC_IGG, HHC_IGA, HHC_IQG, HHC_IQQ;
    
    HHC_CCGG(_M_g4,1,2);
    HHC_CCQA(_M_q2g2,1,2);
    HHC_CCQG(_M_q2g2,1,2);
    HHC_CCGQ(_M_q2g2,1,2);
    HHC_CCQ4(0,-1,1,2);
    d[0] = Vgg*cc[0]/6.0;
    d[1] = Vgg*cc[1]/2.0 + Vga*((Nf-1)*q4[0] + q4[1]/2.0);
    d[2] = Vqg*cc[2]/2.0;
    d[5] = Vqg*cc[5]/6.0 + (Nf-1)*Vqg*q4[0] + Vqg*q4[1];
   
    HHC_CCGG(_M_q2g2,1,2);
    d[0] += Nf*Vgg*cc[0];
    d[2] += (Nf-0.5)*Vqq*cc[0];

    HHC_CCQ4(1,-1,2,0);
    d[3] = Vqg*q4[0];
    d[4] = Vqg*q4[1]/2.0;

    HHC_CCQ4(1,-1,0,2);
    d[6] = Vqg*q4[0];
  }
  
  void hhc2jet::_M_d12(int kt, int i, weight_type& d) 
  {
    double cc[7], HHC_FGG, HHC_FQA, HHC_FQG;
    
    HHC_CCGG(_M_g4,1,2);
    HHC_CCQG(_M_q2g2,1,2);
    HHC_CCGQ(_M_q2g2,1,2);
    HHC_CCQA(_M_q2g2,1,2);
    d[0] = Vgg*cc[0]/6.0 + Nf*Vqa*cc[0];
    d[1] = Vqg*cc[1]/2.0;
    d[2] = Vqg*cc[2]/2.0;
    d[3] = 0.0;
    d[4] = 0.0;
    d[5] = Vgg*cc[5]/6.0 + Nf*Vqa*cc[5];
    d[6] = 0.0;

    HHC_CCQG(_M_q2g2,2,1);
    HHC_CCGQ(_M_q2g2,2,1);
    d[1] += Vqa*(Nf-0.5)*cc[1];
    d[2] += Vqa*(Nf-0.5)*cc[2];
  }
  
  void hhc2jet::_M_d13(int kt, int i, weight_type& d) 
  {
    double cc[7], q4[2], HHC_FGG, HHC_FQG;
    
    HHC_CCGG(_M_g4,1,2);
    HHC_CCQG(_M_q2g2,1,2);
    HHC_CCGQ(_M_q2g2,1,2);
    HHC_CCQA(_M_q2g2,1,2);
    HHC_CCQ4(0,-1,1,2);
    d[0] = Vgg*cc[0]/6.0;
    d[1] = Vqg*cc[1]/2.0;
    d[2] = Vqg*cc[2]/2.0;
    d[5] = Vgg*cc[5]/6.0 + (Nf-1)*Vqg*q4[0] + Vqg*q4[1];

    HHC_CCGG(_M_q2g2,1,2);
    d[0] += Nf*Vqg*cc[0];

    HHC_CCQ4(1,-1,2,0);
    d[3] = Vqg*q4[0];
    d[4] = Vqg*q4[1]/2.0;

    HHC_CCQ4(1,-1,0,2);
    d[6] = Vqg*q4[0];
  }
  
  void hhc2jet::_M_d23(int kt, int i, weight_type& d) 
  {
    double cc[7], q4[2], HHC_FGG, HHC_FQG, HHC_FQA;
    
    HHC_CCGG(_M_g4,1,2);
    HHC_CCQG(_M_q2g2,1,2);
    HHC_CCGQ(_M_q2g2,1,2);
    HHC_CCQA(_M_q2g2,1,2);
    HHC_CCQ4(0,-1,1,2);
    d[0] = Vgg*cc[0]/6.0;
    d[1] = (Vgg + Vqa)*cc[1]/2.0;
    d[2] = (Vgg + Vqa)*cc[2]/2.0;
    d[5] = Vgg*cc[5]/6.0 + (Nf-1)*Vqg*q4[0] + Vqg*q4[1];

    HHC_CCGG(_M_q2g2,1,2);
    d[0] += Nf*Vqg*cc[0];

    HHC_CCQ4(1,-1,2,0);
    d[3] = Vqg*q4[0];
    d[4] = Vqg*q4[1]/2.0;

    HHC_CCQ4(1,-1,0,2);
    d[6] = Vqg*q4[0];
  }
}  //  namespace nlo

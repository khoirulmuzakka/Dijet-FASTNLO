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

//   nlo includes
#include "bits/nlo-color.h"
#include "bits/nlo-split.h"


#define mp(a,b) scalar_mp(a,b)
#define pm(a,b) scalar_pm(a,b)


namespace nlo {

  typedef std::pair<double, std::complex<double> > pair_type;
  typedef const lorentzvector<double>& clvr_type;


  void splitfff::set(clvr_type pi, clvr_type pj, clvr_type pk)
  {
    double sik = pi*pk, sjk = pk*pj;
    
    sij  = pi*pj;
    yijk = sij/(sij + sik + sjk);
    zi   = sik/(sik + sjk);
    
    lorentzvector<double> qij = pi + pj - yijk/(1.0 - yijk)*pk;
    
    hij = pm(pi,pj)*mp(qij,pi)*mp(qij,pj)/(mp(pi,pj)*pm(qij,pi)*pm(qij,pj));

//     static int hijfffcnt = 0;
//     hijfffcnt++;
//     if (hijfffcnt > 168000 ) {
//        std::cout << "Nlojet++:splitfff: hijfffcnt = " << hijfffcnt << std::endl;
//        std::cout << "Nlojet++:splitfff: Numerator: pm(pi,pj) = " << pm(pi,pj) << ", mp(qij,pi) = " << mp(qij,pi) << ", mp(qij,pj) = " << mp(qij,pj) << ", numerator = " << (pm(pi,pj)*mp(qij,pi)*mp(qij,pj)) << std::endl;
//        std::cout << "Nlojet++:splitfff: Denominator: mp(pi,pj) = " << mp(pi,pj) << ", pm(qij,pi) = " << pm(qij,pi) << ", pm(qij,pj) = " << pm(qij,pj) << ", denominator = " << (mp(pi,pj)*pm(qij,pi)*pm(qij,pj)) << std::endl;
//     }

    static int hijfffnaninfcnt = 0;
    if (isnan(real(hij)) || isnan(imag(hij)) || isinf(real(hij)) || isinf(imag(hij)) ) {
       hijfffnaninfcnt++;
       std::cout << "NloJet++_splitfff: ERROR! NaN or Inf in hij, hij = " << hij << std::endl;
       std::cout << "NloJet++_splitfff: NaN or Inf in hij occurrence # = " << hijfffnaninfcnt << std::endl;
       std::cout << "Nlojet++:splitfff: Numerator: pm(pi,pj) = " << pm(pi,pj) << ", mp(qij,pi) = " << mp(qij,pi) << ", mp(qij,pj) = " << mp(qij,pj) << ", numerator = " << (pm(pi,pj)*mp(qij,pi)*mp(qij,pj)) << std::endl;
       std::cout << "Nlojet++:splitfff: Denominator: mp(pi,pj) = " << mp(pi,pj) << ", pm(qij,pi) = " << pm(qij,pi) << ", pm(qij,pj) = " << pm(qij,pj) << ", denominator = " << (mp(pi,pj)*pm(qij,pi)*pm(qij,pj)) << std::endl;
    }
  }

  pair_type splitfff::Vqg() const {
    double vd = -(2.0/(1.0-zi*(1.0-yijk))-1.0-zi)/sij;
    return pair_type(vd, 0.0);
  }

  pair_type splitfff::Vqa() const 
  {
    double color = 0.5/Nc;
    double vd = -color*(1.0 - 2.0*zi*(1.0 - zi))/sij;
    std::complex<double> vc = -2.0*color*zi*(1.0 - zi)*hij/sij;
    if (isnan(real(hij)) || isnan(imag(hij)) || isinf(real(hij)) || isinf(imag(hij)) ) {
       std::cout << "NloJet++_splitfff_Vqa: ERROR! NaN or Inf in hij, hij = " << hij << std::endl;
    }

    return pair_type(vd, vc);
  }

  pair_type splitfff::Vgg() const 
  {
    double oy = 1.0-yijk, zj = 1.0-zi;
    double vd = -(2.0/(1.0-zi*oy) + 2.0/(1.0-zj*oy) - 4.0 + 2.0*zi*zj)/sij;
    std::complex<double> vc = 2.0*zi*zj*hij/sij;
    if (isnan(real(hij)) || isnan(imag(hij)) || isinf(real(hij)) || isinf(imag(hij)) ) {
       std::cout << "NloJet++_splitfff_Vgg: ERROR! NaN or Inf in hij, hij = " << hij << std::endl;
    }

    return pair_type(vd, vc);
  }

  void splitffi::set(clvr_type pi, clvr_type pj, clvr_type pa)
  {
    double sia = pi*pa, sja = pj*pa;

    sij  = pi*pj;
    xija = 1.0-sij/(sia+sja);
    zi   = sia/(sia+sja);
 
    lorentzvector<double> qij = pi + pj - (1.0 - xija)*pa;

    hij = pm(pi,pj)*mp(qij,pi)*mp(qij,pj)/(mp(pi,pj)*pm(qij,pi)*pm(qij,pj));

//     static int hijfficnt = 0;
//     hijfficnt++;
//     if (hijfficnt > 168000 ) {
//        std::cout << "Nlojet++:splitffi: hijfficnt = " << hijfficnt << std::endl;
//        std::cout << "Nlojet++:splitffi: Numerator: pm(pi,pj) = " << pm(pi,pj) << ", mp(qij,pi) = " << mp(qij,pi) << ", mp(qij,pj) = " << mp(qij,pj) << ", numerator = " << (pm(pi,pj)*mp(qij,pi)*mp(qij,pj)) << std::endl;
//        std::cout << "Nlojet++:splitffi: Denominator: mp(pi,pj) = " << mp(pi,pj) << ", pm(qij,pi) = " << pm(qij,pi) << ", pm(qij,pj) = " << pm(qij,pj) << ", denominator = " << (mp(pi,pj)*pm(qij,pi)*pm(qij,pj)) << std::endl;
//     }

    static int hijffinaninfcnt = 0;
    if (isnan(real(hij)) || isnan(imag(hij)) || isinf(real(hij)) || isinf(imag(hij)) ) {
       hijffinaninfcnt++;
       std::cout << "NloJet++_splitffi: ERROR! NaN or Inf in hij, hij = " << hij << std::endl;
       std::cout << "NloJet++_splitffi: NaN or Inf in hij occurrence # = " << hijffinaninfcnt << std::endl;
       std::cout << "Nlojet++:splitffi: Numerator: pm(pi,pj) = " << pm(pi,pj) << ", mp(qij,pi) = " << mp(qij,pi) << ", mp(qij,pj) = " << mp(qij,pj) << ", numerator = " << (pm(pi,pj)*mp(qij,pi)*mp(qij,pj)) << std::endl;
       std::cout << "Nlojet++:splitffi: Denominator: mp(pi,pj) = " << mp(pi,pj) << ", pm(qij,pi) = " << pm(qij,pi) << ", pm(qij,pj) = " << pm(qij,pj) << ", denominator = " << (mp(pi,pj)*pm(qij,pi)*pm(qij,pj)) << std::endl;
    }

  }

  pair_type splitffi::Vqg() const {
    double vd = -(2.0/(2.0 - zi - xija) - 1.0 - zi)/(sij*xija);
    return pair_type(vd, 0.0);
  }

  pair_type splitffi::Vqa() const 
  {
    double color = 0.5/Nc, zj = 1.0 - zi;
    double vd = -color*(1.0 - 2.0*zi*zj)/(sij*xija);
    std::complex<double> vc = -2.0*color*zi*zj*hij/(sij*xija);
    if (isnan(real(hij)) || isnan(imag(hij)) || isinf(real(hij)) || isinf(imag(hij)) ) {
       std::cout << "NloJet++_splitffi_Vqa: ERROR! NaN or Inf in hij, hij = " << hij << std::endl;
    }

    return pair_type(vd, vc);
  }

  pair_type splitffi::Vgg() const {
    double ox = 1.0 - xija, zj = 1.0 - zi;
    double vd = -(2.0/(zj + ox) + 2.0/(zi + ox) - 4.0 + 2.0*zi*zj)/(sij*xija);
    std::complex<double> vc = 2.0*zi*zj*hij/(sij*xija);
    if (isnan(real(hij)) || isnan(imag(hij)) || isinf(real(hij)) || isinf(imag(hij)) ) {
       std::cout << "NloJet++_splitffi_Vgg: ERROR! NaN or Inf in hij, hij = " << hij << std::endl;
    }

    return pair_type(vd, vc);
  }

  void splitiff::set(clvr_type pa, clvr_type pi, clvr_type pk)
  {
    double sik = pi*pk, sak = pk*pa;

    sai  = pa*pi;
    xika = 1.0-sik/(sai+sak);
    ui   = sai/(sai+sak);

    hai = mp(pa,pi)*pm(pi,pk)*mp(pk,pa)/(pm(pa,pi)*mp(pi,pk)*pm(pk,pa));
  }

  pair_type  splitiff::Vqg() const {
    double vd = -(2.0/(1.0 - xika + ui) - 1.0 - xika)/(sai*xika);
    return pair_type(vd, 0.0);
  }

  pair_type splitiff::Vga() const {
    double vd = -Nc*(1.0 - 2.0*xika*(1.0 - xika))/(Na*sai*xika);
    return pair_type(vd, 0.0);
  }

  pair_type splitiff::Vqq() const 
  {
    double color = Cf/Nc, ox = 1.0 - xika;
    double vd = -color*(xika + 2.0*ox/xika)/(sai*xika);
    std::complex<double> vc = 2.0*color*ox/xika*hai/(sai*xika);
    
    return pair_type(vd, vc);
  }

  pair_type splitiff::Vgg() const 
  {
    double ox = 1.0 - xika;
    double vd = -2.0*(1.0/(ox + ui) - 1.0 + (xika + 1.0/xika)*ox)/(sai*xika);
    std::complex<double> vc = 2.0*ox/xika*hai/(sai*xika);

    return pair_type(vd, vc);
  }

  void splitifi::set(clvr_type pa, clvr_type pi, clvr_type pb)
  {
    double sbi = pb*pi, sab = pa*pb;
    sai  = pa*pi;
    xiab = 1.0 - (sai+sbi)/sab;
    hai  = mp(pa,pi)*pm(pi,pb)*mp(pb,pa)/(pm(pa,pi)*mp(pi,pb)*pm(pb,pa));
  }

  pair_type splitifi::Vqg() const {
    double vd = -(2.0/(1.0 - xiab) - 1.0 - xiab)/(sai*xiab);
    return pair_type(vd, 0.0);
  }

  pair_type splitifi::Vga() const {
    double vd = -Nc*(1.0 - 2.0*xiab*(1.0 - xiab))/(Na*sai*xiab);
    return pair_type(vd, 0.0);
  }

  pair_type splitifi::Vqq() const 
  {
    double color = Cf/Nc, ox = 1.0 - xiab;
    double vd = -color*(xiab + 2.0*ox/xiab)/(sai*xiab);
    std::complex<double> vc = 2.0*color*ox/xiab*hai/(sai*xiab);

    return pair_type(vd, vc);
  }

  pair_type splitifi::Vgg() const 
  {
    double ox = 1.0 - xiab;
    double vd = -2.0*(xiab/ox + (xiab + 1.0/xiab)*ox)/(sai*xiab);
    std::complex<double> vc = 2.0*ox/xiab*hai/(sai*xiab);
    
    return pair_type(vd, vc);
  }
}  // namespace nlo

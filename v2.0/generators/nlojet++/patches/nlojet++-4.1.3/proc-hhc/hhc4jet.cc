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
#include "hhc4jet.h"
#include "ampg6.h"
#include "ampq2g4.h"
#include "ampq4g2.h"
#include "ampq6.h"


//    PI_FACn = (2*pi)^(n-2) * 2^(n-2)
#define PI_FAC6 38865023.041825071132937



namespace nlo {

  
  hhc4jet::hhc4jet(const random_generator& rng, bool mchel, unsigned int nu, unsigned int nd, double al)
    : process_hhc(4, 4, nu, nd, al), _hhc_jet_base(nu+nd), _M_mchel(mchel) 
  {
    _M_g6   = new ampg6  (_M_ip, rng);
    _M_q2g4 = new ampq2g4(_M_ip, rng);
    _M_q4g2 = new ampq4g2(_M_ip, rng);
    _M_q6   = new ampq6  (_M_ip, rng);
  }
  
  hhc4jet::~hhc4jet() 
  {
    if(_M_g6)   delete _M_g6;
    if(_M_q2g4) delete _M_q2g4;
    if(_M_q4g2) delete _M_q4g2;
    if(_M_q6)   delete _M_q6;
  }
  
  void hhc4jet::born_term(const event_type& p, weight_type& res) 
  {
    _M_ip.calculate(p);
    amp_tree(_M_g6, _M_q2g4, _M_q4g2, _M_q6, 0, res.begin());
    res *= PI_FAC6;

    for (int ikr = 0; ikr<7; ikr++) {
       if (isnan(res[ikr]) || isinf(res[ikr])) {
          std::cout << "NLOJet++_hhc4jet_born_term: NaN or Inf ERROR for ikr = " << ikr << ", res[ikr] = " << res[ikr] << std::endl;
       }
    }

  }
  
  void hhc4jet::real_term(const event_type&, weight_type&) 
  {}

  void hhc4jet::fini_term(double, double, double, double, const event_type&, weight_type *)
  {}
  
  void hhc4jet::dipole_term(const event_type&, const event_type&,
			    int, int, int, weight_type&) {}
}  //  namespace nlo

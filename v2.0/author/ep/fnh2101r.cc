//------ DON'T TOUCH THIS PART! ------
#include <bits/hhc-phasespace.h>
#include <bits/photo-phasespace.h>
#include <bits/hhc-process.h>
#include <bits/hhc-jetfunc.h>
#include <nlo++-module_add.h>

//----- used namespaces -----
using namespace nlo;
using namespace std;


//----- declaration of the user defined functons -----
void inputfunc(unsigned int&, unsigned int&, unsigned int&);
void psinput(phasespace_photo *, double&, double&);
user_base_hhc * userfunc();

typedef unsigned long int (*module_add_type)(bool, const list<basic_string<char> >&, const basic_string<char>&);
extern  module_add_type module_add;


//----- array of the symbols symbols -----
extern "C"{
  
  struct { 
	const char *name;
	void *address;
  } user_defined_functions[] = 
  {
	//   process index: hhc for gamma + p --> jets (photoproduction resolved photon)
	{"procindex", (void *) "photores"},
  
	//   input function 
	{"inputfunc", (void *) inputfunc},
  
	//   phase space input function 
	{"psinput", (void *) psinput},
  
	//   user defined functions
	{"userfunc",  (void *) userfunc},
  
	//   module to generate the readable result
	{"main_module_add", (void *) module_add},
  
	//  end of the list
	{0, 0}
  };
}
//------ END OF THE DO-NOT-TOUCH-PART ------



//------ USER DEFINED PART STARTS HERE ------
#include <algorithm>
#include <kT_clus.h>
#include "pdf-grv-cteq6.h"


//    Defines new sample type to do felavor decomposation in the initial state
const char *mySample_label[8] = {"gg", "qg", "gq", "qr", "qq", "qqb", "qrb","total"};
typedef weight<8U, mySample_label> mySample;

namespace nlo {
  template<>
  struct weight_conversion<weight_hhc, mySample>
    : public std::unary_function<weight_hhc, mySample>
  {
    mySample operator()(const weight_hhc& x)
    {
      mySample res;
      
      for(unsigned int i = 0; i < 7; i++)
		res[7] += (res[i] = x[i]);

      return res;
    }
  };
}

//   Here we specify some base user objects with the new sample type
//     (see photo-jetfunc.h and included header files therein) 
typedef basic_user<jetfunc_hhc, void,        mySample, weight_conversion> myUser0d;
typedef basic_user<jetfunc_hhc, double,      mySample, weight_conversion> myUser1d;
typedef basic_user<jetfunc_hhc, histpoint1d, mySample, weight_conversion> myUser1h;
typedef basic_user<jetfunc_hhc, histpoint2d, mySample, weight_conversion> myUser2h;


class UserPhoto : public basic_user_set<myUser0d, myUser1h, myUser2h> 
{
public:
  //   init and user function
  void initfunc(unsigned int);
  void userfunc(const event_hhc&, const amplitude_hhc&);
  
private:
  //   pdf
  pdf_grv_cteq6 pdf;

  // algorithms
  kT_clus_long jetclus;
   
  //  private types
  typedef lorentzvector<double> _Lv;

  // the jet structore
  bounded_vector<_Lv> cj, pj_lab, pj_3jet; 
  bounded_vector<unsigned int> jet;
  
  struct pT_sort {
    bool operator()(const _Lv& p1, const _Lv& p2) const {
      return p1.perp2() > p2.perp2();
    }
  };

  struct E_sort {
    bool operator()(const _Lv& p1, const _Lv& p2) const {
      return p1.T() > p2.T();
    }
  };

  //   making simpler the life
  double hard_scale_renorm, hard_scale_fact;
  const amplitude_hhc *amplitude_of_event;
  
};
  

//----- defines the module to sum up the results of the different runs -----
module_add_type module_add = 
main_module_add<basic_user_result<myUser0d::distbook_type, myUser1h::distbook_type, myUser2h::distbook_type> >; 


user_base_hhc * userfunc() {
  return new UserPhoto;
}

void inputfunc(unsigned int& nj, unsigned int& nu, unsigned int& nd)
{
  //  number of jets
  nj = 2U;
  
  //  number of the up and down type flavours
  nu = 2U;
  nd = 3U;
} 

void psinput(phasespace_photo *ps, double& el, double& eh)
{
  //  energy of incomings
  el = 27.6;
  eh = 920.0;
  
  //   You can use your own phase generator. 
  //   Here we use the default.
  ps = 0;
} 


void UserPhoto::initfunc(unsigned int)
{
  myUser0d::phys(0, "sigma(pT(1) > 25GeV, pT(2) > 15GeV)");
  myUser1h::phys(1, "Et,max eta1,2<1", 10, 25., 80.);
  myUser1h::phys(2, "Et,max etai<1,etaj>1", 10, 25., 80.);
  myUser1h::phys(3, "Et,max eta1,2>1", 10, 25., 80.);
  
}


void UserPhoto::userfunc(const event_hhc& p, const amplitude_hhc& amp)
{
  //----- initialize the user object -----
  amplitude_of_event = &amp;

  pdf.photon_momentum_fraction(0.1, 0.9);

 
  //----- do the cluster analysis-----
  jetclus.set_up(1,false,2);
  jetclus(p); jetclus.incl(cj, jet);
  
  unsigned int nj = cj.upper();
  if(nj < 2) return;   //  unable to resolve two jets

  //----- H1 cuts -----
  double pT = 0.0;
  double pTmin1 = 25.0, pTmin2 = 15.0, etamin = -0.5, etamax = 2.75;
  
  pj_lab.resize(1,0);
  for(unsigned int i = 1; i <= nj; i++)
     if((pT = cj[i].perp()) > pTmin2 && -cj[i].prapidity() < etamax && -cj[i].prapidity() > etamin) { // pz in wrong direction
      pj_lab.push_back(cj[i]);
    }
  
  nj = pj_lab.upper();

  if(nj < 2) return;   //  unable to resolve two jets
  
  //----- 2-jet mass -----
//   lorentzvector<double> q = pj_lab[1]+pj_lab[2];
//   double MJ = q.mag();
//   if(MJ < 65.0) return;
  

  //----- sorting the jets by transverse energy -----
  std::sort(pj_lab.begin(), pj_lab.end(), pT_sort()); 
  double pt1 = pj_lab[1].perp();
  double pt2 = pj_lab[2].perp();
  
  if(pt1 < pTmin1) return;

  double eta1 =  -pj_lab[1].prapidity();// pz in wrong direction 
  double eta2 =  -pj_lab[2].prapidity();// pz in wrong direction 

  // calculate cos(theta)
  //  double costheta = fabs( tanh((eta1-eta2)/2.) ); 
  
  //----- pdf -----
  amp.pdf_and_qcd_coupling(pdf, 389385730.);
  
  //----- hard scale -----
  //----- hard scale -----
  double meanpt = (pt1+pt2)/2.;
  hard_scale_renorm = hard_scale_fact = pow(meanpt,2.);

  //----- calculate the weight -----
  amplitude_hhc::contrib_type itype = amp.contrib();
  const weight_hhc& wt = amp(hard_scale_renorm, hard_scale_fact);

  //----- total 2-jet cross section -----
  myUser0d::physfill(0, wt);

  // Cut on x_gamma applied by restricting the y range
  double ymin = 1./(2.*0.8*27.6)*(
                                  pj_lab[1].perp()*exp(+pj_lab[1].prapidity())
                                + pj_lab[2].perp()*exp(+pj_lab[2].prapidity())
                                  );
  // x_gamma>0.8
  pdf.photon_momentum_fraction(0.1,ymin<0.9?ymin:0.9);
  //----- pdf -----
  amp.pdf_and_qcd_coupling(pdf, 389385730.);
  //----- calculate the weight -----
  const weight_hhc& wt2 = amp(hard_scale_renorm, hard_scale_fact);

  int id =0;
  if(eta1<1 && eta2<1){
     id = 1;
  }else{
     if(eta1>1 && eta2>1){
        id = 3;
     }else{
        id=2;
     }
  }

  switch(itype) {
  case amplitude_hhc::real: sphysfilld(id, pt1, sphysreal(wt2)); break;
  case amplitude_hhc::sub:  sphysfilld(id, pt1, sphyssub(wt2));  break;
  default: physfilld(id, pt1, wt2);
  }

  // Cut in xP
//   double Q2 = 2.*(p[-2]*p[-1]);
//   double xbj = Q2/(2.*p[hadron(0)]*(p[-1]-p[-2]));
//   if (xbj>0.1) return;



}




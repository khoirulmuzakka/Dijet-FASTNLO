//------ DON'T TOUCH THIS PART! ------
#include <bits/photo-phasespace.h>
#include <bits/photo-process.h>
#include <bits/photo-jetfunc.h>
#include <nlo++-module_add.h>

//----- used namespaces -----
using namespace nlo;
using namespace std;


//----- declaration of the user defined functons -----
void inputfunc(unsigned int&, unsigned int&, unsigned int&);
void psinput(phasespace_photo *, double&, double&);
user_base_photo * userfunc();

typedef unsigned long int (*module_add_type)(bool, const list<basic_string<char> >&, const basic_string<char>&);
extern  module_add_type module_add;


//----- array of the symbols symbols -----
extern "C"{
  
  struct { 
	const char *name;
	void *address;
  } user_defined_functions[] = 
  {
	//   process index: hhc for gamma + p --> jets (photoproduction direct photon)
	{"procindex", (void *) "photodir"},
  
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
#include "pdf-cteq6.h"
#include "fnloTable.h"


//    Defines new sample type to do felavor decomposation in the initial state
const char *mySample_label[4] = {"g", "u", "d", "total"};
typedef weight<4U, mySample_label> mySample;

namespace nlo {
  template<>
  struct weight_conversion<weight_photo, mySample>
    : public std::unary_function<weight_photo, mySample>
  {
    mySample operator()(const weight_photo& x)
    {
      mySample res;
      
      for(unsigned int i = 0; i < 3; i++)
		res[3] += (res[i] = x[i]);

      return res;
    }
  };
}

//   Here we specify some base user objects with the new sample type
//     (see photo-jetfunc.h and included header files therein) 
typedef basic_user<jetfunc_photo, void,        mySample, weight_conversion> myUser0d;
typedef basic_user<jetfunc_photo, double,      mySample, weight_conversion> myUser1d;
typedef basic_user<jetfunc_photo, histpoint1d, mySample, weight_conversion> myUser1h;
typedef basic_user<jetfunc_photo, histpoint2d, mySample, weight_conversion> myUser2h;


class UserPhoto : public basic_user_set<myUser0d, myUser1h, myUser2h> 
{
public:
  //   init and user function
   void initfunc(unsigned int);
   void userfunc(const event_photo&, const amplitude_photo&);

   void end_of_event();  
   virtual void phys_output(const std::basic_string<char>& fname, unsigned long nsave = 10000UL, bool txt = false);
  
private:
  //   pdf
  pdf_cteq6 pdf;

  // algorithms
  kT_clus_long jetclus;
   
  //  private types
  typedef lorentzvector<double> _Lv;

  // the jet structore
  bounded_vector<_Lv> cj, pj_lab; 
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
  const amplitude_photo *amplitude_of_event;

   //fastNLO starts here

   fnloTable *table;
   
   // ===== variables for the b-cubic interpolation =====
   // - the relative distances to the four nearest bins
   vector<double> cm ; 
   // - the weights for the cubic eigenfunctions (1-dim)
   vector<double> cefm; 

   double nevents;        // No of events calculated so far
   unsigned long nwrite;  // No of events after to write out the table

   string tablefilename; // The table file to write to
   bool textoutput; // If true, the table is written in plain ASCII instead of BASE64 encoded doubles (later for XML)
   amplitude_photo::integral_type itype; // Born, NLO etc.  
   time_t start_time;
   
   bool nlo;
   void writetable();
  
};
  

//----- defines the module to sum up the results of the different runs -----
module_add_type module_add = 
main_module_add<basic_user_result<myUser0d::distbook_type, myUser1h::distbook_type, myUser2h::distbook_type> >; 

user_base_photo * userfunc() {
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
  myUser1h::phys(1, "Et,max eta1,2<1", 10, 25., 80.);
  myUser1h::phys(2, "Et,max etai<1,etaj>1", 10, 25., 80.);
  myUser1h::phys(3, "Et,max eta1,2>1", 10, 25., 80.);


  //Set up fastNLO
  table = new fnloTable("notknownyet.tab");
  table->GetBlockA1()->SetScenName("fnh2101d");
  table->GetBlockA1()->SetNcontrib(1);
  table->GetBlockA1()->SetNmult(0);
  table->GetBlockA1()->SetNdata(0);
  table->GetBlockA2()->SetIpublunits(9);
  table->CreateBlockB(0);



  // ---- Initialize event counters
  nevents = 0;
  // Set some defaults
  if (nwrite==0) nwrite = 5000000;
  
  start_time = std::time(0);
}


void UserPhoto::userfunc(const event_photo& p, const amplitude_photo& amp)
{
  //----- initialize the user object -----
  amplitude_of_event = &amp;

  //----- photon momentum fraction -----
  double y = (p[-1]*p[hadron(0)])/(p[hadron(-1)]*p[hadron(0)]);
  if(y < 0.1 || y > 0.9) return; 

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

  //----- sorting the jets by transverse energy -----
  std::sort(pj_lab.begin(), pj_lab.end(), pT_sort()); 
  double pt1 = pj_lab[1].perp();
  double pt2 = pj_lab[2].perp();
  
  if(pt1 < pTmin1) return;

  double eta1 =  -pj_lab[1].prapidity();// pz in wrong direction 
  double eta2 =  -pj_lab[2].prapidity();// pz in wrong direction 

  // calculate x_gamma
  double xgamma = 1./(2.*y*27.6)*(pt1*exp(-eta1) + pt2*exp(-eta2) );
  
  //----- pdf -----
  amp.pdf_and_qcd_coupling(pdf, 389385730./137.0);
  
  //----- hard scale -----
  double meanpt = (pt1+pt2)/2.;
  hard_scale_renorm = hard_scale_fact = pow(meanpt,2.);
  
  //----- calculate the weight -----
  amplitude_photo::contrib_type itype = amp.contrib();
  const weight_photo& wt = amp(hard_scale_renorm, hard_scale_fact);

  //cut on xgamma
  if(xgamma<0.8) return;

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
  case amplitude_photo::real: sphysfilld(id, pt1, sphysreal(wt)); break;
  case amplitude_photo::sub:  sphysfilld(id, pt1, sphyssub(wt));  break;
  default: physfilld(id, pt1, wt);
  }

}

void UserPhoto::end_of_event(){
   // let NLOJET++ store its results
   basic_user_set<myUser0d, myUser1h, myUser2h>::end_of_event();
   
   nevents += 1;
   //-------- store table
   if (( (unsigned long)nevents % nwrite)==0){
      time_t hour, min, time = std::time(0) - start_time;
      
      hour = time/3600L;
      time -= hour*3600L;
      min  = time/60L;
      time -= min*60L;
      
      std::cout<<"--->     "
	  <<(hour < 10 ? "0" : "")<<hour
	  <<(min < 10 ? ":0" : ":")<<min
	  <<(time < 10 ? ":0" : ":")<<time<<std::endl;
      printf ("No. events: %.2G writing table....",nevents);
      cout.flush();
      writetable();
      printf("done.\n");
   }
}

void UserPhoto::phys_output(const std::basic_string<char>& __file_name, 
                          unsigned long __save, bool __txt) 
{
   // Suppress output of NLOJET++ files
   basic_user_set<myUser0d, myUser1h, myUser2h>::phys_output("",2000000000,false);
   tablefilename.assign(__file_name.c_str());
   tablefilename += ".tab";
   table->SetFilename(tablefilename);

   //Determine whether we are running LO or NLO
   const char* const file = __file_name.c_str(); 
   if(strstr(file,"born")!=NULL){
      nlo = false;
   }else{
      if(strstr(file,"nlo")!=NULL){
         nlo = true;
      }else{
         printf("This module can only be run at Born level or at NLO.\n");
         exit(1);
      }
   }
   textoutput = __txt;
   nwrite = __save;
}

void UserPhoto::writetable(){
   table->OpenFileWrite();
   table->WriteBlockA1();
   table->WriteBlockA2();
   table->WriteBlockB(0);
   table->CloseFileWrite();
}

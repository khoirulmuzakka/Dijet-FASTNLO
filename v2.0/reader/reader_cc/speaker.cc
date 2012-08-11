// Author: Daniel Britzger
// DESY, 06/08/2012

#include "speaker.h"
#include <string>
#include <iostream>
#include <map>


using namespace std;


std::map<unsigned long,speaker*>* speaker::list = NULL;
std::ostream* speaker::weg = NULL;
say::Verbosity speaker::fverb = say::INFO;
unsigned long speaker::ct = 0;
bool speaker::fe2cerr = true;

speaker::speaker(std::string prefix,say::Verbosity volume,bool err,bool quiet){
   if ( list==NULL ) list = new map<unsigned long,speaker*>();
   if ( weg==NULL ) {
      weg = new ostream(0);
      weg->clear(std::ios::badbit);
   }
   pref=prefix;
   fii=ct;
   (*list)[ct++] = this;
   fvol=volume;
   errs=err;
   fquiet= ( quiet || fvol<fverb );
   //std::cout<<"adding new speaker to list:"<<fii<<"\tsize="<<list->size()<<"\tfquiet="<<fquiet<<"\tvol="<<fvol<<"\tprefix="<<prefix<<std::endl;
}

speaker::~speaker(){
   list->erase(fii);
   if(list->empty()){
      delete list; list=NULL;
      delete weg; weg=NULL;
   }
}

std::ostream& speaker::operator() (std::string fct) const {
   if (fquiet) return *weg;
   //       *this<<"In "<<fct<<". ";
   if(errs && fe2cerr) return std::cerr<<fct;
   else return std::cout<<fct;
}

std::ostream& speaker::operator>> (std::string arg) const {
   return print(arg);
}

std::ostream& speaker::print(std::string mes) const {
   if (fquiet) return *weg;
   else {	
      if (errs&&fe2cerr) return std::cerr<<mes;
      else  return std::cout<<mes;
   }
}

std::ostream& speaker::operator[] (std::string fct) const {
   if ( fquiet ) return *weg;
   if ( !cn.empty()) return *this<<"["<<cn<<"::"<<fct<<"] ";
   else return *this<<"["<<fct<<"] ";
}

const speaker& speaker::prefix(std::string fct) const {
   if ( !fquiet ) {
      if (errs&&fe2cerr) std::cerr<<fct;
      else std::cout<<fct;
   }
   return *this; 
}


int speaker::SetGlobalVerbosity(say::Verbosity volume){
   fverb=volume;
   int c=0;
   for( map<unsigned long, speaker*>::const_iterator ii=(*list).begin(); ii!=(*list).end(); ++ii){
      (*ii).second->DoSpeak( (*ii).second->GetVolume()>=volume );
      c++;
   }
   return c;
}


PrimalScream::PrimalScream(std::string classname){//,std::string prefix=""){
   cn=classname;
   debug = speaker("Debug. ",say::DEBUG);
   debug.SetClassName(cn);
   man   = speaker("",say::MANUAL);
   man.SetClassName(cn);
   info  = speaker("Info. ",say::INFO);
   info.SetClassName(cn);
   warn  = speaker("Warning. ",say::WARNING);
   warn.SetClassName(cn);
   error = speaker("Error! ",say::ERROR,true);
   error.SetClassName(cn);
   debug["PrimalScream"]<<"Primal Scream initialized."<<std::endl;
}

void PrimalScream::SetVerbosity(say::Verbosity volume){
   debug.DoSpeak( debug.GetVolume() >= volume );
   man.DoSpeak( man.GetVolume() >= volume );
   info.DoSpeak( info.GetVolume() >= volume );
   warn.DoSpeak( warn.GetVolume() >= volume );
   error.DoSpeak( error.GetVolume() >= volume );
}

namespace say {
   speaker debug("Debug. ",say::DEBUG);
   speaker man("",say::MANUAL);
   speaker info("Info. ",say::INFO);
   speaker warn("Warning. ",say::WARNING);
   speaker error("Error! ",say::ERROR,true);
   //debug["namespace say"]<<"speakers initialized."<<std::endl;
   int SetGlobalVerbosity(Verbosity verbosity){ return speaker::SetGlobalVerbosity(verbosity);};
}



// namespace say {
//    enum Verbosity {DEBUG=-1000, MANUAL=-1, INFO=0, WARNING=1,ERROR=2,SILENT=1000};
// }
/*
class speaker {
public:
//    speaker(std::string prefix="",bool quiet=false,bool err=false) : 
//       weg(0) , quiet(quiet){
//       weg.clear(std::ios::badbit);
//       pref=prefix;
//       errs=err; };
//    speaker(const speaker& spk):weg(0){;}; 
//    std::ostream& operator() (std::string fct) {
//       if (quiet) return weg;
//       *this<<"In "<<fct<<". ";
//       if(errs) return std::cerr;
//       else return std::cout;
//    }
//    std::ostream& operator[] (std::string fct) {
//       if (quiet) return weg;
//       return *this<<"["<<fct<<"] ";
//       if(errs) return std::cerr;
//       else return std::cout;
//       //return *this;
//    }
//    speaker& operator+ (std::string fct){return this->prefix(fct);}
//    speaker& prefix(std::string fct) {
//       if (!quiet){
// 	 if (errs) std::cerr<<fct;
// 	 else std::cout<<fct;;
//       }
//       return *this;
//    }
//    template<typename T> std::ostream& operator<< (T arg) {
//       if (quiet) return weg;
//       else {	
// 	 if (errs) return std::cerr<<pref<<arg;
// 	 else return std::cout<<pref<<arg;
//       }
//    };
//    template<typename T> std::ostream& operator>> (T arg) {
//       return print(arg);
//    };
//    std::ostream& print(string mes) { 
//       if (quiet) return weg;
//       else {	
// 	 if (errs) return std::cerr<<mes;
// 	 else  return std::cout<<mes;
//       }
//    }
//    void DoSpeak(bool loud){quiet=!loud;};
//    bool GetSpeak() const {return !quiet;};
//    void SetPrefix(std::string prefix){pref=prefix;};
//    std::string GetPrefix(std::string prefix) const {return pref;};
   
   //speaker(const speaker& spk) : weg(0) {;}; 
   std::ostream& operator[] (std::string fct);
   speaker& operator+ (std::string fct){return this->prefix(fct);}
   speaker& prefix(std::string fct);
   
   template<typename T> std::ostream& operator<< (T arg) {
      if (fquiet) return *weg;
      else {	
	 if (errs && fe2cerr) return std::cerr<<pref<<arg;
	 else return std::cout<<pref<<arg;
      }
   }
   void DoSpeak(bool loud){fquiet=!loud;};
   bool GetSpeak() const {return !fquiet;};
   void SetPrefix(std::string prefix){pref=prefix;};
   std::string GetPrefix(std::string prefix) const {return pref;};
   void SetClassName(std::string classname){cn=classname;};
   std::string GetClassName(void) const {return cn;};
   say::Verbosity GetVolume(void) const { return vol;};
   void SetVolume(say::Verbosity volume) {vol=volume;};
   static void ErrorToErrStream(bool ToCerr){fe2cerr=ToCerr;};
   
protected:
protected:
   //std::ostream weg;
   static std::ostream* weg = NULL;
   bool fquiet;
   std::string pref;
   bool errs;
   say::Verbosity fvol;
   unsigned long fii;
   static unsigned long ct = 0;
   static bool fe2cerr = true;
   static say::Verbosity fverb;
   static std::map<unsigned long,speaker*>* list = NULL;
   std::string cn;

};

namespace say {
   extern speaker debug;
   extern speaker man;
   extern speaker info;
   extern speaker warn;
   extern speaker error;
   extern int SetGlobalVerbosity(Verbosity verbosity);
}


class PrimalScream {
public: 
   PrimalScream(std::string classname,std::string prefix="");
   void SetVerbosity(say::Verbosity volume);
   std::string cn;
   speaker debug;
   speaker man;
   speaker info;
   speaker warn;
   speaker error;
};

#endif //SPEAKER_H_
*/

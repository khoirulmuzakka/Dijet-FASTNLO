// Author: Daniel Britzger
// DESY, 06/08/2012

#ifndef SPEAKER_H_
#define SPEAKER_H_

#include <string>
#include <iostream>
#include <map>

namespace say {
enum Verbosity {DEBUG=-1000, MANUAL=-1, INFO=0, WARNING=1, ERROR=2, SILENT=1000};
}

class speaker {
public:
   speaker(std::string prefix="",say::Verbosity volume=say::INFO,bool err=false,bool quiet=false);
   speaker(const speaker& spk) {
      ;
   };
   ~speaker();
   //speaker(const speaker& spk) : weg(0) {;};
   std::ostream& operator[](std::string fct) const ;
   const speaker& operator+ (std::string fct) const {
      return this->prefix(fct);
   }
   const speaker& prefix(std::string fct) const ;
   std::ostream& operator()(std::string fct) const ;

   template<typename T> std::ostream& operator<< (T arg) const {
      if (fquiet) return *weg;
      else {
         if (errs && fe2cerr) return std::cerr<<pref<<arg;
         else return std::cout<<pref<<arg;
      }
   }
   std::ostream& operator>> (std::string arg) const ;
   std::ostream& print(std::string mes) const ;
   void DoSpeak(bool loud) {
      fquiet=!loud;
   };
   bool GetSpeak() const {
      return !fquiet;
   };
   void SetPrefix(std::string prefix) {
      pref=prefix;
   };
   std::string GetPrefix(std::string prefix) const {
      return pref;
   };
   void SetClassName(std::string classname) {
      cn=classname;
   };
   std::string GetClassName(void) const {
      return cn;
   };
   say::Verbosity GetVolume(void) const {
      return fvol;
   };
   void SetVolume(say::Verbosity volume) {
      fvol=volume;
   };
   static int SetGlobalVerbosity(say::Verbosity volume);
   static void ErrorToErrStream(bool ToCerr) {
      fe2cerr=ToCerr;
   };

protected:
   //std::ostream weg;
   static std::ostream* weg;
   bool fquiet;
   std::string pref;
   bool errs;
   say::Verbosity fvol;
   unsigned long fii;
   static unsigned long ct;
   static bool fe2cerr;
   static say::Verbosity fverb;
   static std::map<unsigned long,speaker*>* list;
   std::string cn;
};

namespace say {
extern speaker debug;
extern speaker man;
extern speaker info;
extern speaker warn;
extern speaker error;
extern speaker shout; // same as error but streamed to cout
extern int SetGlobalVerbosity(Verbosity verbosity);
}


class PrimalScream {
public:
   PrimalScream(std::string classname);//,std::string prefix="");
   void SetVerbosity(say::Verbosity volume);
   std::string cn;
protected:
   speaker debug;
   speaker man;
   speaker info;
   speaker warn;
   speaker error;
   speaker shout;
};

#endif //SPEAKER_H_

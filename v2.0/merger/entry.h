#ifndef __entry__
#define __entry__

#include <vector>
#include "fnloTable.h"

class Contrib;

class Contrib{
 public:
   bool operator==(Contrib val) const{
      return(
           IRef==val.IRef &&
           IDataFlag==val.IDataFlag && 
           IAddMultFlag==val.IAddMultFlag &&
           IContrFlag1==val.IContrFlag1 &&
           IContrFlag2==val.IContrFlag2 &&
           IContrFlag3==val.IContrFlag3  &&
             Npow==val.Npow);}
  string GetName1(){
     if(IDataFlag>0){
        return "data";
     }
     if(IAddMultFlag>0){
        return "correction factor";
     }
     switch(IContrFlag1){
     case 0: return "unknown";break;
     case 1: return "fixed order";break;
     case 2: return "SM correction";break;
     case 3: return "new physics";break;
     default:;
     }

  }
  string GetName2(int ILOord){
     string Refstring = "";
     if(IRef>0){
        Refstring = " (reference)";
     }
     if(IDataFlag>0){
        return "";
     }
     if(IAddMultFlag>0){
        return "";
     }
     switch(IContrFlag1){
     case 0: return "unknown"; break;
     case 1:
        switch(Npow-ILOord){
        case 0: return "LO"+Refstring; break;
        case 1: return "NLO"+Refstring; break;
        case 2: return "NNLO"+Refstring; break;
        case 3: return "NNLO"+Refstring; break;
        default:;
        }
     default:;
     }

   }
 public:
   int IRef;
   int IDataFlag;
   int IAddMultFlag;
   int IContrFlag1;
   int IContrFlag2;
   int IContrFlag3;
   int Npow;
};

class tableptr{
 public:
   vector <fnloTable*>::iterator table;
   int blockbno;
};

class Entry{
 public:
   Contrib contribution;
   vector <tableptr> tables;
};


#endif

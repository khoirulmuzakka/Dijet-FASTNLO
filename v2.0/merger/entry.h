#ifndef __entry__
#define __entry__

#include <vector>
#include "fnloTable.h"

class TContrib;

class TContrib{
 public:
   bool operator==(TContrib val) const{
      return(IDataFlag==val.IDataFlag && 
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
        case 0: return "LO"; break;
        case 1: return "NLO"; break;
        case 2: return "NNLO"; break;
        case 3: return "NNLO"; break;
        default:;
        }
     default:;
     }

   }
 public:
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
   TContrib contribution;
   vector <tableptr> tables;
};


#endif

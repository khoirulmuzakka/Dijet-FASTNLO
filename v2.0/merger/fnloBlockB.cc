#include "fnloBlockB.h"


int fnloBlockB::Read(istream *table){
   table->peek();
   if (table->eof()){
      printf("fnloBlockB::Read: Cannot read from file.\n");
      return(2);
   }

   int key = 0;
   *table >> key;
   if(key != tablemagicno){
      printf("fnloBlockB::Read: At beginning of block found %d instead of %d.\n",key,tablemagicno);
      return 1;
   };

   *table >> IXsectUnits;
   *table >> IDataFlag;
   *table >> IAddMultFlag;
   *table >> IContrFlag1;
   *table >> IContrFlag2;
   *table >> IContrFlag3;
   *table >> NContrDescr;
   CtrbDescript.resize(NContrDescr);
   char buffer[257];
   table->getline(buffer,256);
   for(int i=0;i<NContrDescr;i++){
      table->getline(buffer,256);
      CtrbDescript[i] = buffer;
      StripWhitespace(CtrbDescript[i]);
   }
   *table >> NCodeDescr;   
   CodeDescript.resize(NCodeDescr);
   table->getline(buffer,256);
   for(int i=0;i<NCodeDescr;i++){
      table->getline(buffer,256);
      CodeDescript[i] = buffer;
      StripWhitespace(CodeDescript[i]);
   }
   if(IDataFlag==1){
      *table >> Nuncorrel;
      UncDescr.resize(Nuncorrel);
      table->getline(buffer,256);
      for(int i=0;i<Nuncorrel;i++){
         table->getline(buffer,256);
         UncDescr[i] = buffer;
         StripWhitespace(UncDescr[i]);
      }
      *table >> Ncorrel;
      CorDescr.resize(Ncorrel);
      table->getline(buffer,256);
      for(int i=0;i<Ncorrel;i++){
         table->getline(buffer,256);
         CorDescr[i] = buffer;
         StripWhitespace(CorDescr[i]);
      }
      Xcenter.resize(BlockA2->GetNObsBin());
      Value.resize(BlockA2->GetNObsBin());
      UncorLo.resize(BlockA2->GetNObsBin());
      UncorHi.resize(BlockA2->GetNObsBin());
      CorrLo.resize(BlockA2->GetNObsBin());
      CorrHi.resize(BlockA2->GetNObsBin());
      for(int i=0;i<BlockA2->GetNObsBin();i++){
         *table >> Xcenter[i];
         *table >> Value[i];
         UncorLo[i].resize(Nuncorrel);
         UncorHi[i].resize(Nuncorrel);
         for(int j=0;j<Nuncorrel;j++){
            *table >> UncorLo[i][j];
            *table >> UncorHi[i][j];
         }
         CorrLo[i].resize(Ncorrel);
         CorrHi[i].resize(Ncorrel);
         for(int j=0;j<Ncorrel;j++){
            *table >> CorrLo[i][j];
            *table >> CorrHi[i][j];
         }
      }
      *table >> NErrMatrix;
      matrixelement.resize(NErrMatrix);
      for(int i=0;i<NErrMatrix;i++){
         matrixelement[i].resize((int)pow((double)BlockA2->GetNObsBin(),2));
         for(int j=0;j<(int)pow((double)BlockA2->GetNObsBin(),2);j++){
            *table >> matrixelement[i][j];
         }
      }
   }// end of IDataFlag==1
   if(IAddMultFlag==1){
      *table >> Nuncorrel;
      UncDescr.resize(Nuncorrel);
      table->getline(buffer,256);
      for(int i=0;i<Nuncorrel;i++){
         table->getline(buffer,256);
         UncDescr[i] = buffer;
         StripWhitespace(UncDescr[i]);
      }
      *table >> Ncorrel;
      CorDescr.resize(Ncorrel);
      table->getline(buffer,256);
      for(int i=0;i<Ncorrel;i++){
         table->getline(buffer,256);
         CorDescr[i] = buffer;
         StripWhitespace(CorDescr[i]);
      }
      fact.resize(BlockA2->GetNObsBin());
      UncorLo.resize(BlockA2->GetNObsBin());
      UncorHi.resize(BlockA2->GetNObsBin());
      CorrLo.resize(BlockA2->GetNObsBin());
      CorrHi.resize(BlockA2->GetNObsBin());
      for(int i=0;i<BlockA2->GetNObsBin();i++){
         *table >> fact[i];
         UncorLo[i].resize(Nuncorrel);
         UncorHi[i].resize(Nuncorrel);
         for(int j=0;j<Nuncorrel;j++){
            *table >> UncorLo[i][j];
            *table >> UncorHi[i][j];
         }
         CorrLo[i].resize(Ncorrel);
         CorrHi[i].resize(Ncorrel);
         for(int j=0;j<Ncorrel;j++){
            *table >> CorrLo[i][j];
            *table >> CorrHi[i][j];
         }
      }
   }// end of IAddMultFlag==1
   if(!(IDataFlag==1) && !(IAddMultFlag==1)){
      *table >> IRef;
      *table >> IScaleDep;
      *table >> Nevt;
      *table >> Npow;
      *table >> NPDF;
      if(NPDF>0){
         NPDFPDG.resize(NPDF);
         for(int i=0;i<NPDF;i++){
            *table >>  NPDFPDG[i];
         }      
      }
      *table >> NPDFDim;
      *table >> NFragFunc;
      if(NFragFunc>0){
         NFFPDG.resize(NFragFunc);
         for(int i=0;i<NFragFunc;i++){
            *table >>  NFFPDG[i];
         }      
      }
      *table >> NFFDim;
      *table >> NSubproc;
      *table >> IPDFdef1;
      *table >> IPDFdef2;
      *table >> IPDFdef3;
      if(IPDFdef1==0){
         for(int i=0;i<NSubproc;i++){
            // Missing: linear PDF combinations for IPDFdef1=0
            if(NPDF==1){
            }else{
               if(NPDF==2){
               }
            }
         }
      }
      Nxtot.resize(BlockA2->GetNObsBin());
      XNode1.resize(BlockA2->GetNObsBin());
      for(int i=0;i<BlockA2->GetNObsBin();i++){
         *table >> Nxtot[i];
         XNode1[i].resize(Nxtot[i]);
         for(int j=0;j<Nxtot[i];j++){
            *table >> XNode1[i][j];
         }         
      }      
      if(NPDFDim==2){
         Nxtot2.resize(BlockA2->GetNObsBin());
         XNode2.resize(BlockA2->GetNObsBin());
         for(int i=0;i<BlockA2->GetNObsBin();i++){
            *table >> Nxtot2[i];
            XNode2[i].resize(Nxtot2[i]);
            for(int j=0;j<Nxtot2[i];j++){
               *table >> XNode2[i][j];
            }         
         }      
      }
      if(NFragFunc>0){
         Nztot.resize(BlockA2->GetNObsBin());
         ZNode.resize(BlockA2->GetNObsBin());
         for(int i=0;i<BlockA2->GetNObsBin();i++){
            *table >> Nztot[i];
            ZNode[i].resize(Nztot[i]);
            for(int j=0;j<Nztot[i];j++){
               *table >> ZNode[i][j];
            }         
         }               
      }

      *table >> NScales;
      *table >> NScaleDim;
      Iscale.resize(NScales);
      for(int i=0;i<NScales;i++){
         *table >> Iscale[i];
      }      
      NscaleDescript.resize(NScaleDim);
      ScaleDescript.resize(NScaleDim);
      for(int i=0;i<NScaleDim;i++){
         *table >> NscaleDescript[i];
         ScaleDescript[i].resize(NscaleDescript[i]);
         table->getline(buffer,256);
         for(int j=0;j<NscaleDescript[i];j++){
            table->getline(buffer,256);
            ScaleDescript[i][j] = buffer;
            StripWhitespace(ScaleDescript[i][j]);
         }
      }
      Nscalevar.resize(NScaleDim);
      Nscalenode.resize(NScaleDim);
      for(int i=0;i<NScaleDim;i++){
         *table >> Nscalevar[i];
         *table >> Nscalenode[i];
      }
      ScaleFac.resize(NScaleDim);
      for(int i=0;i<NScaleDim;i++){
         ScaleFac[i].resize(Nscalevar[i]);
         for(int j=0;j<Nscalevar[i];j++){
            *table >> ScaleFac[i][j];
         }
      }
      ScaleNode.resize(BlockA2->GetNObsBin());
      for(int i=0;i<BlockA2->GetNObsBin();i++){
         ScaleNode[i].resize(NScaleDim);
         for(int j=0;j<NScaleDim;j++){
            ScaleNode[i][j].resize(Nscalevar[j]);
            for(int k=0;k<Nscalevar[j];k++){
               ScaleNode[i][j][k].resize(Nscalenode[j]);
               for(int l=0;l<Nscalenode[j];l++){
                  *table >> ScaleNode[i][j][k][l];
               }
            }
         }
      }
      
      SigmaTilde.resize(BlockA2->GetNObsBin());
      for(int i=0;i<BlockA2->GetNObsBin();i++){
         int nxmax =0;
         switch (NPDFDim) {
         case 0: nxmax = Nxtot[i];
            break;
         case 1: nxmax = ((int)pow((double)Nxtot[i],2)+Nxtot[i])/2;
            break;
         case 2: nxmax = Nxtot[i]*Nxtot2[i];
            break;
         default: ;
         }
         SigmaTilde[i].resize(NScaleDim);
         for(int j=0;j<NScaleDim;j++){
            SigmaTilde[i][j].resize(Nscalevar[j]);
            for(int k=0;k<Nscalevar[j];k++){
               SigmaTilde[i][j][k].resize(Nscalenode[j]);
               for(int l=0;l<Nscalenode[j];l++){
                  SigmaTilde[i][j][k][l].resize(nxmax);
                  for(int m=0;m<nxmax;m++){
                     SigmaTilde[i][j][k][l][m].resize(NSubproc);
                     for(int n=0;n<NSubproc;n++){
                        *table >> SigmaTilde[i][j][k][l][m][n];
                     }
                  }
               }
            }
         }
      }
   }// end of not data and not corrections
   key = 0;
   *table >> key;
   if(key != tablemagicno){
      printf("fnloBlockB::Read: At end of block found %d instead of %d.\n",key,tablemagicno);
      return 1;
   };
   // Put magic number back
   for(int i=0;i<(int)(log10((double)key)+1);i++){
      table->unget();
   }
   return 0;
}

int fnloBlockB::Write(ostream *table){

   *table << tablemagicno << endl;
   *table << IXsectUnits << endl;
   *table << IDataFlag << endl;
   *table << IAddMultFlag << endl;
   *table << IContrFlag1 << endl;
   *table << IContrFlag2 << endl;
   *table << IContrFlag3 << endl;
   *table << NContrDescr << endl;
   for(int i=0;i<NContrDescr;i++){
      *table << CtrbDescript[i] << endl;
   }
   *table << NCodeDescr << endl;   
   for(int i=0;i<NCodeDescr;i++){
      *table << CodeDescript[i] << endl;
   }
   if(IDataFlag==1){
      *table << Nuncorrel << endl;
      for(int i=0;i<Nuncorrel;i++){
         *table << UncDescr[i] << endl;
      }
      *table << Ncorrel << endl;
      for(int i=0;i<Ncorrel;i++){
         *table << CorDescr[i]  << endl;
      }
      for(int i=0;i<BlockA2->GetNObsBin();i++){
         *table << Xcenter[i] << endl;
         *table << Value[i] << endl;
         for(int j=0;j<Nuncorrel;j++){
            *table << UncorLo[i][j] << endl;
            *table << UncorHi[i][j] << endl;
         }
         for(int j=0;j<Ncorrel;j++){
            *table << CorrLo[i][j] << endl;
            *table << CorrHi[i][j] << endl;
         }
      }
      *table << NErrMatrix << endl;
      for(int i=0;i<NErrMatrix;i++){
         for(int j=0;j<(int)pow((double)BlockA2->GetNObsBin(),2);j++){
            *table << matrixelement[i][j] << endl;
         }
      }
   }// end of IDataFlag==1
   if(IAddMultFlag==1){
      *table << Nuncorrel << endl;
      for(int i=0;i<Nuncorrel;i++){
         *table << UncDescr[i]  << endl;
      }
      *table << Ncorrel << endl;
      for(int i=0;i<Ncorrel;i++){
         *table << CorDescr[i]  << endl;
      }
      for(int i=0;i<BlockA2->GetNObsBin();i++){
         *table << fact[i] << endl;
         for(int j=0;j<Nuncorrel;j++){
            *table << UncorLo[i][j] << endl;
            *table << UncorHi[i][j] << endl;
         }
         for(int j=0;j<Ncorrel;j++){
            *table << CorrLo[i][j] << endl;
            *table << CorrHi[i][j] << endl;
         }
      }
   }// end of IAddMultFlag==1
   if(!(IDataFlag==1) && !(IAddMultFlag==1)){
      *table << IRef << endl;
      *table << IScaleDep << endl;
      *table << Nevt << endl;
      *table << Npow << endl;
      *table << NPDF << endl;
      if(NPDF>0){
         for(int i=0;i<NPDF;i++){
            *table <<  NPDFPDG[i] << endl;
         }      
      }
      *table << NPDFDim << endl;
      *table << NFragFunc << endl;
      if(NFragFunc>0){
         for(int i=0;i<NFragFunc;i++){
            *table <<  NFFPDG[i] << endl;
         }      
      }
      *table << NFFDim << endl;
      *table << NSubproc << endl;
      *table << IPDFdef1 << endl;
      *table << IPDFdef2 << endl;
      *table << IPDFdef3 << endl;
      if(IPDFdef1==0){
         for(int i=0;i<NSubproc;i++){
            // Missing: linear PDF combinations for IPDFdef1=0
            if(NPDF==1){
            }else{
               if(NPDF==2){
               }
            }
         }
      }
      for(int i=0;i<BlockA2->GetNObsBin();i++){
         *table << Nxtot[i] << endl;
         for(int j=0;j<Nxtot[i];j++){
            *table << XNode1[i][j] << endl;
         }         
      }      
      if(NPDFDim==2){
         for(int i=0;i<BlockA2->GetNObsBin();i++){
            *table << Nxtot2[i] << endl;
            for(int j=0;j<Nxtot2[i];j++){
               *table << XNode2[i][j] << endl;
            }         
         }      
      }
      if(NFragFunc>0){
         for(int i=0;i<BlockA2->GetNObsBin();i++){
            *table << Nztot[i] << endl;
            for(int j=0;j<Nztot[i];j++){
               *table << ZNode[i][j] << endl;
            }         
         }               
      }

      *table << NScales << endl;
      *table << NScaleDim << endl;
      for(int i=0;i<NScales;i++){
         *table << Iscale[i] << endl;
      }      
      for(int i=0;i<NScaleDim;i++){
         *table << NscaleDescript[i] << endl;
         for(int j=0;j<NscaleDescript[i];j++){
            *table << ScaleDescript[i][j] << endl;
         }
      }
      for(int i=0;i<NScaleDim;i++){
         *table << Nscalevar[i] << endl;
         *table << Nscalenode[i] << endl;
      }
      for(int i=0;i<NScaleDim;i++){
         for(int j=0;j<Nscalevar[i];j++){
            *table << ScaleFac[i][j] << endl;
         }
      }
      for(int i=0;i<BlockA2->GetNObsBin();i++){
         for(int j=0;j<NScaleDim;j++){
            for(int k=0;k<Nscalevar[j];k++){
               for(int l=0;l<Nscalenode[j];l++){
                  *table << ScaleNode[i][j][k][l] << endl;
               }
            }
         }
      }
      
      for(int i=0;i<BlockA2->GetNObsBin();i++){
         int nxmax =0;
         switch (NPDFDim) {
         case 0: nxmax = Nxtot[i];
            break;
         case 1: nxmax = ((int)pow((double)Nxtot[i],2)+Nxtot[i])/2;
            break;
         case 2: nxmax = Nxtot[i]*Nxtot2[i];
            break;
         default: ;
         }
         for(int j=0;j<NScaleDim;j++){
            for(int k=0;k<Nscalevar[j];k++){
               for(int l=0;l<Nscalenode[j];l++){
                  for(int m=0;m<nxmax;m++){
                     for(int n=0;n<NSubproc;n++){
                        *table << SigmaTilde[i][j][k][l][m][n] << endl;
                     }
                  }
               }
            }
         }
      }
   }// end of not data and not corrections

   return 0;
}

void fnloBlockB::Add(fnloBlockB* other){
   double w1 = (double)Nevt / (Nevt+other->Nevt);
   double w2 = (double)other->Nevt / (Nevt+other->Nevt);
   for(int i=0;i<BlockA2->GetNObsBin();i++){
      int nxmax =0;
      switch (NPDFDim) {
      case 0: nxmax = Nxtot[i];
         break;
      case 1: nxmax = ((int)pow((double)Nxtot[i],2)+Nxtot[i])/2;
         break;
      case 2: nxmax = Nxtot[i]*Nxtot2[i];
         break;
      default: ;
      }
      for(int j=0;j<NScaleDim;j++){
         for(int k=0;k<Nscalevar[j];k++){
            for(int l=0;l<Nscalenode[j];l++){
               for(int m=0;m<nxmax;m++){
                  for(int n=0;n<NSubproc;n++){
                     SigmaTilde[i][j][k][l][m][n] = 
                        w1*SigmaTilde[i][j][k][l][m][n] + w2*other->SigmaTilde[i][j][k][l][m][n];
                  }
               }
            }
         }
      }
   }
   Nevt += other->Nevt;
;
}



bool fnloBlockB::IsCompatible(fnloBlockB* other){
   return true;
};

void fnloBlockB::StripWhitespace(string &str){
   for(string::iterator achar = str.end(); achar>str.begin();achar--) {
      if (*achar==0x20 || *achar==0x00){
         str.erase(achar);
      }else{
         break;
      }
   }
}

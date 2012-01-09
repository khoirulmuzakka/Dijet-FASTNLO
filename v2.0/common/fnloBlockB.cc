// KR: Add include because of header clean-up in gcc-4.3
#include <cstdlib>
#include <iostream>
#include <cmath>

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
   // *table >> IContrFlag3;	// IContrFlag3 is written here in v2.0 and v2.1 but  not in v2.0+
   // in v2.1. IContrFlag3 will be reintroduces again, and NScaleDep will be stored later in the table
   IContrFlag3 = 0;
   *table >> NScaleDep;
   int NContrDescr;
   *table >> NContrDescr;
   //   printf("  *  infnloBlockB::Read().  IDataFlag: %d, IAddMultFlag: %d, IContrFlag1: %d, IContrFlag2: %d, IContrFlag3: %d, NScaleDep: %d\n",IDataFlag,IAddMultFlag,IContrFlag1,IContrFlag2,IContrFlag3,NScaleDep );
   CtrbDescript.resize(NContrDescr);
   char buffer[257];
   table->getline(buffer,256);
   for(int i=0;i<NContrDescr;i++){
      table->getline(buffer,256);
      CtrbDescript[i] = buffer;
      //      StripWhitespace(CtrbDescript[i]);
   }
   int NCodeDescr;
   *table >> NCodeDescr;   
   CodeDescript.resize(NCodeDescr);
   table->getline(buffer,256);
   for(int i=0;i<NCodeDescr;i++){
      table->getline(buffer,256);
      CodeDescript[i] = buffer;
      //      StripWhitespace(CodeDescript[i]);
   }
   if(IDataFlag==1){
      *table >> Nuncorrel;
      UncDescr.resize(Nuncorrel);
      table->getline(buffer,256);
      for(int i=0;i<Nuncorrel;i++){
         table->getline(buffer,256);
         UncDescr[i] = buffer;
         //         StripWhitespace(UncDescr[i]);
      }
      *table >> Ncorrel;
      CorDescr.resize(Ncorrel);
      table->getline(buffer,256);
      for(int i=0;i<Ncorrel;i++){
         table->getline(buffer,256);
         CorDescr[i] = buffer;
         //         StripWhitespace(CorDescr[i]);
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
         //         StripWhitespace(UncDescr[i]);
      }
      *table >> Ncorrel;
      CorDescr.resize(Ncorrel);
      table->getline(buffer,256);
      for(int i=0;i<Ncorrel;i++){
         table->getline(buffer,256);
         CorDescr[i] = buffer;
         //         StripWhitespace(CorDescr[i]);
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
      //printf("  *  infnloBlockB::Read(). IRef : %d, IScaleDep: %d, Nevt: %d, Npow: %d, NPDF: %d, NPDFDim: %d\n", IRef ,IScaleDep  ,Nevt  , Npow ,NPDF , NPDFDim  );

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
      Nxtot1.resize(BlockA2->GetNObsBin());
      XNode1.resize(BlockA2->GetNObsBin());
      for(int i=0;i<BlockA2->GetNObsBin();i++){
         *table >> Nxtot1[i];
         XNode1[i].resize(Nxtot1[i]);
         for(int j=0;j<Nxtot1[i];j++){
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
      int NscaleDescript;
      ScaleDescript.resize(NScaleDim);
      for(int i=0;i<NScaleDim;i++){
         *table >> NscaleDescript;
         ScaleDescript[i].resize(NscaleDescript);
         table->getline(buffer,256);
         for(int j=0;j<NscaleDescript;j++){
            table->getline(buffer,256);
            ScaleDescript[i][j] = buffer;
            //            StripWhitespace(ScaleDescript[i][j]);
         }
      }


      //! v2.1 store NScaleDep here.
      //! v2.1 *table >> NScaleDep;
      
      if ( NScaleDep != 3 ) {
	 Nscalevar.resize(NScaleDim);
	 Nscalenode.resize(NScaleDim);
	 for(int i=0;i<NScaleDim;i++){
	    *table >> Nscalevar[i];
	    *table >> Nscalenode[i];
	 }
	 
	 // 	 printf("  *  infnloBlockB::Read().bins %d, NScalevar[0] %d, Nscalenode[0] %d,  NScaleDim %d  \n",
	 // 	 BlockA2->GetNObsBin(), Nscalevar[0] , Nscalenode[0] , NScaleDim );

	 
	 ScaleFac.resize(NScaleDim);
	 for(int i=0;i<NScaleDim;i++){
	    ScaleFac[i].resize(Nscalevar[i]);
	    for(int j=0;j<Nscalevar[i];j++){
	       *table >> ScaleFac[i][j];
	    }
	 }
	 
	 //printf("  *  infnloBlockB::Read().bins %d, NScalevar[0] %d, Nscalenode[0] %d, ScaleFac[0][0] %d,  NScaleDim %d  \n",
	 //BlockA2->GetNObsBin(), Nscalevar[0] , Nscalenode[0] , ScaleFac[0][0], NScaleDim );

	 //! DB: This for-loop-mess was replace by ResizeTable)( and ReadTable()      
	 //       ScaleNode.resize(BlockA2->GetNObsBin());
	 //       for(int i=0;i<BlockA2->GetNObsBin();i++){
	 //          ScaleNode[i].resize(NScaleDim);
	 //          for(int j=0;j<NScaleDim;j++){
	 //             ScaleNode[i][j].resize(Nscalevar[j]);
	 //             for(int k=0;k<Nscalevar[j];k++){
	 //                ScaleNode[i][j][k].resize(Nscalenode[j]);
	 //                for(int l=0;l<Nscalenode[j];l++){
	 //                   *table >> ScaleNode[i][j][k][l];
	 //                }
	 //             }
	 //          }
	 //       }

	 // @MW, @KR: I use this shorthand notation for reading vectors. Will it work for you?
	 //    I only see problems, if you make use of NScaleDim>1
	 ResizeTable( &ScaleNode , BlockA2->GetNObsBin(), 1 , Nscalevar[0] , Nscalenode[0] ); // should work, since NScaleDim==1, but is not yet tested for 100%
	 int nsn = ReadTable  ( &ScaleNode , table );
	 //printf("  *  infnloBlockB::Read(). Read %d lines of ScaleNode.\n",nsn);
	 
	 int XmaxFromI[1] = {0};
	 //printf(" &SigmaTilde  %i  %i  %i  *%i  %i\n", 
	 //BlockA2->GetNObsBin(), GetTotalScalevars(), GetTotalScalenodes(), XmaxFromI[0], NSubproc);
	 ResizeTable( &SigmaTilde , BlockA2->GetNObsBin(), GetTotalScalevars(), GetTotalScalenodes(), XmaxFromI, NSubproc );
	 int nst = ReadTable  ( &SigmaTilde , table );
	 //printf("  *  infnloBlockB::Read(). Read %d lines of SigmaTilde.\n",nst);
	 printf("  *  infnloBlockB::Read(). Read %d lines of FNLO v2.0 tables.\n",nst+nsn);
	 
	 //printf(" &PdfLc  %i  %i  #%i  %i\n", BlockA2->GetNObsBin(), GetTotalScalenodes(), XmaxFromI[0], NSubproc);
	 ResizeTable( &PdfLc , BlockA2->GetNObsBin(), GetTotalScalenodes(), XmaxFromI, NSubproc );
	 
	 //! DB: redundant code: was replaced by ResizeTable() and ReadTable() methods.
	 //       SigmaTilde.resize(BlockA2->GetNObsBin());
	 //       PdfLc.resize(BlockA2->GetNObsBin());
	 //       for(int i=0;i<BlockA2->GetNObsBin();i++){
	 //          int nxmax = GetNxmax(i);
	 //          int totalscalevars = GetTotalScalevars();
	 //          SigmaTilde[i].resize(totalscalevars);
	 //          for(int k=0;k<totalscalevars;k++){
	 //             int totalscalenodes =  GetTotalScalenodes();
	 //             SigmaTilde[i][k].resize(totalscalenodes);
	 //             PdfLc[i].resize(totalscalenodes);
	 //             for(int l=0;l<totalscalenodes;l++){
	 //                SigmaTilde[i][k][l].resize(nxmax);
	 //                PdfLc[i][l].resize(nxmax);
	 //                for(int m=0;m<nxmax;m++){
	 //                   SigmaTilde[i][k][l][m].resize(NSubproc);
	 //                   PdfLc[i][l][m].resize(NSubproc);
	 //                   for(int n=0;n<NSubproc;n++){
	 //                      *table >> SigmaTilde[i][k][l][m][n];
	 //                      PdfLc[i][l][m][n] = 0.;
	 //                   }
	 //                }
	 //             }
	 //          }
	 //       }

      }


      if ( NScaleDep == 3 ) {

	 //  ---- order of reading... ---- //
	 //    - nscalenode q2
	 //    - scalenode Q
	 //    - nscalenode pt
	 //    - scalenode pt
	 //    - simgatilde mu indep
	 //    - simgatilde mu_f dep
	 //    - simgatilde mu_r dep
	 //    - sigmarefmixed
	 //    - sigmaref scale 1
	 //    - sigmaref scale 2
	 // ------------------------------ //
	 int nn3 = 0;
	  
	 nn3 += ReadFlexibleVector  ( &ScaleNode1 , table );
	 nn3 += ReadFlexibleVector  ( &ScaleNode2 , table );
 	 NscalenodeScale1 = ScaleNode1[0].size();
 	 NscalenodeScale2 = ScaleNode2[0].size();

	 nn3 += ReadFlexibleVector  ( &SigmaTildeMuIndep , table , true );
	 nn3 += ReadFlexibleVector  ( &SigmaTildeMuFDep , table , true );
	 nn3 += ReadFlexibleVector  ( &SigmaTildeMuRDep , table , true );

	 nn3 += ReadFlexibleVector  ( &SigmaRefMixed , table , true );
	 nn3 += ReadFlexibleVector  ( &SigmaRef_s1 , table , true );
	 nn3 += ReadFlexibleVector  ( &SigmaRef_s2 , table , true );

// 	 *table >> NscalenodeScale1 ;
// 	 ResizeTable( &ScaleNode1 , BlockA2->GetNObsBin() , NscalenodeScale1 );
// 	 nn3 += ReadTable  ( &ScaleNode1 , table );

// 	 *table >> NscalenodeScale2 ;
// 	 ResizeTable( &ScaleNode2 , BlockA2->GetNObsBin() , NscalenodeScale2 );
// 	 nn3 += ReadTable  ( &ScaleNode2 , table );

// 	 int XMaxFromFromDim[1] = { 0 };
// 	 //ResizeTable( &PdfLcMuVar , BlockA2->GetNObsBin() , XMaxFromFromDim , NscalenodeScale1 , NscalenodeScale2 , NSubproc );
// 	 //ResizeTable( &AlphasTwoPi , BlockA2->GetNObsBin() , NscalenodeScale1 , NscalenodeScale2 );

// 	 ResizeTable( &SigmaTildeMuIndep , BlockA2->GetNObsBin() , XMaxFromFromDim , NscalenodeScale1 , NscalenodeScale2 , NSubproc );
// 	 nn3 += ReadTable  ( &SigmaTildeMuIndep , table );

// 	 ResizeTable( &SigmaTildeMuFDep , BlockA2->GetNObsBin() , XMaxFromFromDim , NscalenodeScale1 , NscalenodeScale2 , NSubproc );
// 	 nn3 += ReadTable  ( &SigmaTildeMuFDep , table );

// 	 ResizeTable( &SigmaTildeMuRDep , BlockA2->GetNObsBin() , XMaxFromFromDim , NscalenodeScale1 , NscalenodeScale2 , NSubproc );
// 	 nn3 += ReadTable  ( &SigmaTildeMuRDep , table );

// 	 ResizeTable( &SigmaRefMixed , BlockA2->GetNObsBin() , NSubproc );
// 	 nn3 += ReadTable  ( &SigmaRefMixed , table );

// 	 ResizeTable( &SigmaRef_s1 , BlockA2->GetNObsBin() , NSubproc );
// 	 nn3 += ReadTable  ( &SigmaRef_s1 , table );

// 	 ResizeTable( &SigmaRef_s2 , BlockA2->GetNObsBin() , NSubproc );
// 	 nn3 += ReadTable  ( &SigmaRef_s2 , table );
	 printf("  *  infnloBlockB::Read(). Read %d lines of NScaleDep==3 Tables.\n",nn3);

      }


   }// end of not data and not corrections

   key = 0;
   *table >> key;
   if(key != tablemagicno){
      printf("fnloBlockB::Read: At end of block found %d instead of %d.\n",key,tablemagicno);
      printf("                  You might have 'nan' in your table.\n");
      return 1;
   };
   // Put magic number back
   for(int i=0;i<(int)(log10((double)key)+1);i++){
      table->unget();
   }
   return 0;
}

int fnloBlockB::Write(ostream *table, int option){

   *table << tablemagicno << endl;
   *table << IXsectUnits << endl;
   *table << IDataFlag << endl;
   *table << IAddMultFlag << endl;
   *table << IContrFlag1 << endl;
   *table << IContrFlag2 << endl;
   //*table << IContrFlag3 << endl;	// v2.0+. for v2.1 write IContrFlag3 here, but NScaleDep only later
   *table << NScaleDep << endl;
   *table << CtrbDescript.size() << endl;
   //printf("  *  infnloBlockB::Write().  IDataFlag: %d, IAddMultFlag: %d, IContrFlag1: %d, IContrFlag2: %d, IContrFlag3: %d, NScaleDep: %d\n",
   //IDataFlag,IAddMultFlag,IContrFlag1,IContrFlag2,IContrFlag3,NScaleDep);
   for(int i=0;i<CtrbDescript.size();i++){
      *table << CtrbDescript[i] << endl;
   }
   *table << CodeDescript.size() << endl;   
   for(int i=0;i<CodeDescript.size();i++){
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
         *table << Nxtot1[i] << endl;
         for(int j=0;j<Nxtot1[i];j++){
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
         *table << ScaleDescript[i].size() << endl;
         for(int j=0;j<ScaleDescript[i].size();j++){
            *table << ScaleDescript[i][j] << endl;
         }
      }

      //! v2.1 store NScaleDep here
      //! *table << NScaleDep << endl;

      if ( NScaleDep != 3 ){
	 for(int i=0;i<NScaleDim;i++){
	    *table << Nscalevar[i] << endl;
	    *table << Nscalenode[i] << endl;
	 }
	 for(int i=0;i<NScaleDim;i++){
	    for(int j=0;j<Nscalevar[i];j++){
	       *table << ScaleFac[i][j] << endl;
	    }
	 }
	 
	int nsn = WriteTable( &ScaleNode  , table );
	//printf("  *  infnloBlockB::Write(). Wrote %d lines of ScaleNode.\n",nsn);
	int nst = WriteTable( &SigmaTilde , table , (bool)(option & DividebyNevt) , Nevt );
	//printf("  *  infnloBlockB::Write(). Wrote %d lines of SigmaTilde.\n",nst);
	printf("  *  infnloBlockB::Write(). Wrote %d lines of FNLO v2.0 tables.\n",nst+nsn);


      } // end if NScaleDep !=3.
      
      if ( NScaleDep == 3 ) {
	 int nn3 = 0;
       
	 nn3 += WriteFlexibleTable( &ScaleNode1 , table );
	 nn3 += WriteFlexibleTable( &ScaleNode2 , table );
 
 	 nn3 += WriteFlexibleTable( &SigmaTildeMuIndep, table , (bool)(option & DividebyNevt) , Nevt , true );
 	 nn3 += WriteFlexibleTable( &SigmaTildeMuFDep , table , (bool)(option & DividebyNevt) , Nevt , true );
 	 nn3 += WriteFlexibleTable( &SigmaTildeMuRDep , table , (bool)(option & DividebyNevt) , Nevt , true );
      
	 nn3 += WriteFlexibleTable( &SigmaRefMixed	, table , (bool)(option & DividebyNevt) , Nevt , true );
	 nn3 += WriteFlexibleTable( &SigmaRef_s1	, table , (bool)(option & DividebyNevt) , Nevt , true );
	 nn3 += WriteFlexibleTable( &SigmaRef_s2	, table , (bool)(option & DividebyNevt) , Nevt , true );

// 	 *table << NscalenodeScale1 << endl;
// 	 nn3 += WriteTable( &ScaleNode1 , table );
     
// 	 *table << NscalenodeScale2 << endl;
// 	 nn3 += WriteTable( &ScaleNode2 , table );
 
// 	 nn3 += WriteTable( &SigmaTildeMuIndep, table , (bool)(option & DividebyNevt) , Nevt );
// 	 nn3 += WriteTable( &SigmaTildeMuFDep , table , (bool)(option & DividebyNevt) , Nevt );
// 	 nn3 += WriteTable( &SigmaTildeMuRDep , table , (bool)(option & DividebyNevt) , Nevt );
      
// 	 nn3 += WriteTable( &SigmaRefMixed	, table , (bool)(option & DividebyNevt) , Nevt );
// 	 nn3 += WriteTable( &SigmaRef_s1	, table , (bool)(option & DividebyNevt) , Nevt );
// 	 nn3 += WriteTable( &SigmaRef_s2	, table , (bool)(option & DividebyNevt) , Nevt );

	 printf("  *  infnloBlockB::Write(). Wrote %d lines of v2.1 Tables.\n",nn3);
	  
      } // if(NScaleDep==3)
   }// end of not data and not corrections

   return 0;
}

int fnloBlockB::Copy(fnloBlockB* other){

   streambuf* streambuf = new stringbuf(ios_base::in | ios_base::out); 
   iostream* buffer = new iostream(streambuf);
   other->Write(buffer);
   *buffer << tablemagicno << endl;
   this->Read(buffer);
   delete buffer;
   delete streambuf;

   return(0);
}

void fnloBlockB::Add(fnloBlockB* other){
   double w1 = (double)Nevt / (Nevt+other->Nevt);
   double w2 = (double)other->Nevt / (Nevt+other->Nevt);
   Nevt += other->Nevt;

   if ( NScaleDep != 3 ){
      AddTableToAnotherTable( &SigmaTilde , &(other->SigmaTilde) ,w1 , w2 );
   }
   
   if ( NScaleDep == 3 ){
     AddTableToAnotherTable( &SigmaTildeMuIndep , &(other->SigmaTildeMuIndep) ,w1 , w2 );
     AddTableToAnotherTable( &SigmaTildeMuFDep , &(other->SigmaTildeMuFDep) ,w1 , w2 );
     AddTableToAnotherTable( &SigmaTildeMuRDep , &(other->SigmaTildeMuRDep) ,w1 , w2 );
     AddTableToAnotherTable( &SigmaRefMixed , &(other->SigmaRefMixed) ,w1 , w2 );
     AddTableToAnotherTable( &SigmaRef_s1 , &(other->SigmaRef_s1) ,w1 , w2 );
     AddTableToAnotherTable( &SigmaRef_s2 , &(other->SigmaRef_s2) ,w1 , w2 );
   }

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

int fnloBlockB::GetNxmax(int i){
   int nxmax = 0;
   switch (NPDFDim) {
   case 0: nxmax = Nxtot1[i];
      break;
   case 1: nxmax = ((int)pow((double)Nxtot1[i],2)+Nxtot1[i])/2;
      break;
   case 2: nxmax = Nxtot1[i]*Nxtot2[i];
      break;
   default: ;
   }
   return nxmax;
};

int fnloBlockB::GetXIndex(int Obsbin,int x1bin,int x2bin){
   int xbin = 0;
   switch (NPDFDim) {
   case 0: xbin = x1bin; // linear
      break;
   case 1: xbin = x1bin + (x2bin*(x2bin+1)/2);    // half matrix
      break;
   case 2: xbin = x1bin + x2bin * Nxtot1[Obsbin]; // full matrix
      break;
   default: ;
   }
   return xbin;
};



int fnloBlockB::GetTotalScalevars(){
   int totalscalevars=1;
   for(int scaledim=0;scaledim<NScaleDim;scaledim++){
      totalscalevars *= Nscalevar[scaledim];
   }
   return totalscalevars;
}

int fnloBlockB::GetTotalScalenodes(){
   int totalscalenodes=1;
   for(int scaledim=0;scaledim<NScaleDim;scaledim++){
      totalscalenodes *= Nscalenode[scaledim];
   }
   return totalscalenodes;
}


void fnloBlockB::ResetXsection(){
   Xsection.resize(BlockA2->GetNObsBin());
   for(int i=0;i<BlockA2->GetNObsBin();i++){
      Xsection[i] = 0.0; 
   }
}

void fnloBlockB::CalcPDFLinearComb(vector<double> pdfx1, vector<double> pdfx2, vector<double> *result){

   // indices for pdfx1 and pdfx2:
   // 0..5 = tbar, ..., ubar, dbar;
   // 6 = g;
   // 7..12 = d, u, ..., t
   
   vector<double> &pdflc = *result; // to make acces via "[]" more easy
   
   switch(IPDFdef1){
   case 2: // ep , DIS and gP
      switch(IPDFdef2){
      case 1: // DIS: determine delta,gluon,sigma
      case 2: // direct gammaP: delta,gluon,sigma

         pdflc[0] = 0.;
         for(int l=0;l<13;l++){
            double temp = (l==6 ? 0.0 : pdfx1[l]);
            if (!(l&1)) temp *= 4.;
            pdflc[0] += temp; // --- delta
         }
         pdflc[0] /= 9.;
         pdflc[1] = pdfx1[6]; // --- gluon
         if(NSubproc>2){ // only from NLO
            pdflc[2] = 0.;
            for(int l=0;l<6;l++){
               pdflc[2] += pdfx1[5-l] + pdfx1[l+7]; // --- sigma
            }
         }
         break;
      default: printf("fnloBlockB::CalcPDFLinearComb : Ipdfdef1=2, Ipdfdef2= %d not supported. Exit.\n",IPDFdef2); exit(1);
      }
      break;
   case 3:
      switch(IPDFdef2){
      case 1: // ppbar: gg   qg   gq   qr   qq   qqb   qrb 
         double B0,B,Bb;
         double A0,A,Ab;
         double D,Db;

         A0 = pdfx2[6];
         B0 = pdfx1[6];
         
         A = Ab = 0.;
         B = Bb = 0.;
         D = Db = 0.;
          for(int l=0;l<6;l++){
             A  += pdfx2[l+7];
             Ab += pdfx2[5-l];
// pp
//              B  += pdfx1[l+7];
//              Bb += pdfx1[5-l];
//              D  += pdfx1[l+7] * pdfx2[l+7] + pdfx1[5-l] * pdfx2[5-l];
//              Db += pdfx1[l+7] * pdfx2[5-l] + pdfx1[5-l] * pdfx2[l+7];
// ppbar
             Bb  += pdfx1[l+7];
             B += pdfx1[5-l];
             Db  += pdfx1[l+7] * pdfx2[l+7] + pdfx1[5-l] * pdfx2[5-l];
             D += pdfx1[l+7] * pdfx2[5-l] + pdfx1[5-l] * pdfx2[l+7];
          }         

         pdflc[0] = A0*B0; // gluon gluon
         pdflc[1] = (A + Ab)*B0; // quark gluon
         pdflc[2] = A0*(B + Bb);
         pdflc[3] = A*B + Ab*Bb - D;
         pdflc[4] = D;
         pdflc[5] = Db;
         pdflc[6] = A*Bb +Ab*B - Db;        
         break;
      default: printf("fnloBlockB::CalcPDFLinearComb :Ipdfdef1=3, Ipdfdef2= %d not supported. Exit.\n",IPDFdef2); exit(1);
      }
      break;
   case 4:
      switch(IPDFdef2){
      case 1: // resolved gammaP: gg   qg   gq   qr   qq   qqb   qrb 
         double B0,B,Bb;
         double A0,A,Ab;
         double D,Db;

         A0 = pdfx2[6];
         B0 = pdfx1[6];
         
         A = Ab = 0.;
         B = Bb = 0.;
         D = Db = 0.;
          for(int l=0;l<6;l++){
             A  += pdfx2[l+7];
             Ab += pdfx2[5-l];
             B  += pdfx1[l+7];
             Bb += pdfx1[5-l];
             D  += pdfx1[l+7] * pdfx2[l+7] + pdfx1[5-l] * pdfx2[5-l];
             Db += pdfx1[l+7] * pdfx2[5-l] + pdfx1[5-l] * pdfx2[l+7];
          }         

         pdflc[0] = A0*B0; // gluon gluon
         pdflc[1] = (A + Ab)*B0; // quark gluon
         pdflc[2] = A0*(B + Bb);
         pdflc[3] = A*B + Ab*Bb - D;
         pdflc[4] = D;
         pdflc[5] = Db;
         pdflc[6] = A*Bb +Ab*B - Db;
 
         break;
      default: printf("fnloBlockB::CalcPDFLinearComb ::Ipdfdef1=4, Ipdfdef2= %d not supported. Exit.\n",IPDFdef2); exit(1);
      }
      break;
   default: printf("fnloBlockB::CalcPDFLinearComb :Ipdfdef1= %d not supported. Exit.\n",IPDFdef1); exit(1);
   }
}

void fnloBlockB::FillPDFCache(int scalevar, void (fnloTableUser::*GetPdfs)(double x, double muf,vector<double> &xfx),fnloTableUser *tableptr){
   int scaleindex2 = 0;
   int scalevar2 = scalevar;

   if (NScaleDim>1){
      scaleindex2 = 1; // If we use multiple scales, then mu_f is by convention the second scale -> index=1 
      scalevar2 = scalevar % Nscalevar[1]; 
   }

   // linear
   if(NPDFDim == 0){
      vector<double> xfx; // PDFs of all partons
      xfx.resize(13);
       vector <double> buffer; // for resorting a pdf array
      buffer.resize(NSubproc);
     
      for(int i=0;i<BlockA2->GetNObsBin();i++){
         int nxmax = GetNxmax(i);
         for(int j=0;j<Nscalenode[scaleindex2];j++){
            for(int k=0;k<nxmax;k++){ 
               (tableptr->*GetPdfs)(XNode1[i][k],ScaleNode[i][scaleindex2][scalevar2][j],xfx);
               CalcPDFLinearComb(xfx,xfx,&buffer); //calculate linear combinations
               for(int l=0;l<NSubproc;l++){ 
                  PdfLc[i][j][k][l] = buffer[l];
               }
            }
         }
      }
   }

   // half matrix
   if(NPDFDim == 1){
      vector < vector<double> > xfx; // PDFs of all partons
      vector <double> buffer; // for resorting a pdf array
      buffer.resize(NSubproc);
      for(int i=0;i<BlockA2->GetNObsBin();i++){
         int nxmax = GetNxmax(i);  // total entries in half matrix
         int nxbins1 = Nxtot1[i]; // number of columns in half matrix
         xfx.resize(nxbins1);
         for(int j=0;j<Nscalenode[scaleindex2];j++){
            // determine all pdfs of hadron1
            for(int k=0;k<nxbins1;k++){ 
               xfx[k].resize(13);
               (tableptr->*GetPdfs)(XNode1[i][k],ScaleNode[i][scaleindex2][scalevar2][j],xfx[k]);
            }
            int x1bin = 0;
            int x2bin = 0;
            for(int k=0;k<nxmax;k++){ 
               CalcPDFLinearComb(xfx[x1bin],xfx[x2bin],&buffer); //calculate linear combinations
               for(int l=0;l<NSubproc;l++){ 
                  PdfLc[i][j][k][l] = buffer[l];
               }
               x1bin++;
               if(x1bin>x2bin){
                  x1bin = 0;
                  x2bin++;
               }
            }
         }
      }
   }
   


}

void fnloBlockB::FillPDFCache(int scalevar, void (fnloTableUser::*GetPdfs)(double x, double muf,vector<double> &xfx),
                              void (fnloTableUser::*GetPdfs2)(double x, double muf,vector<double> &xfx), fnloTableUser *tableptr){
   vector < vector<double> > xfx; // PDFs of all partons
   vector < vector<double> > xfx2;

   vector <double> buffer; // for resorting a pdf array
   buffer.resize(NSubproc);

   int scaleindex2 = 0;
   int scalevar2 = scalevar;

   if (NScaleDim>1){
      scaleindex2 = 1; // If we use multiple scales, then mu_f is by convention the second scale -> index=1 
      scalevar2 = scalevar % Nscalevar[1]; 
   }

   for(int i=0;i<BlockA2->GetNObsBin();i++){
      int nxmax = GetNxmax(i);
      for(int j=0;j<Nscalenode[scaleindex2];j++){

         // determine all pdfs of hadron1
         int nxbins1 = nxmax / Nxtot2[i];
         xfx.resize(nxbins1);
         for(int k=0;k<nxbins1;k++){ 
            xfx[k].resize(13);
            (tableptr->*GetPdfs)(XNode1[i][k],ScaleNode[i][scaleindex2][scalevar2][j],xfx[k]);
         }

         // determine all pdfs of hadron2
         int nxbins2 = nxmax / Nxtot1[i];
         xfx2.resize(nxbins2);
         for(int k=0;k<nxbins2;k++){ 
            xfx2[k].resize(13);
               (tableptr->*GetPdfs2)(XNode2[i][k],ScaleNode[i][scaleindex2][scalevar2][j],xfx2[k]);
         }         

         for(int k=0;k<nxmax;k++){ 
            int x1bin = k % Nxtot1[i];
            int x2bin = k / Nxtot1[i];
            //            printf("i=%d scalevar=%d j=%d scale=%f\n ",i,scalevar,j,ScaleNode[i][0][scalevar][j]);
            //            printf("x1=%f x2=%f gluonP=%f gluonG=%f \n",XNode1[i][x1bin],XNode2[i][x2bin],xfx[x1bin][6],xfx2[x2bin][6]);
            CalcPDFLinearComb(xfx[x1bin],xfx2[x2bin],&buffer); //calculate linear combinations
            for(int l=0;l<NSubproc;l++){ 
               PdfLc[i][j][k][l] = buffer[l];
            }
         }
      }
   }
}

void fnloBlockB::CalcXsection(double asmz, int scalevar, double rescale){

   static const double TWOPI = 6.28318530717958647692528;

   ResetXsection();

   if(IsReference()){
      for(int i=0;i<BlockA2->GetNObsBin();i++){
         for(int l=0;l<NSubproc;l++){ 
            Xsection[i] +=  SigmaTilde[i][scalevar][0][0][l];
         }
      }
      return;
   }

   if(IsLO()){
       for(int i=0;i<BlockA2->GetNObsBin();i++){
          int nxmax = GetNxmax(i);
          int scaleindex1 = 0; // Even if we use multiple scales, mu_r is by convention the first scale -> index=0
          int scalevar1 = scalevar;
          if (NScaleDim>1){          
             scalevar1 = scalevar / Nscalevar[1];
          }
          for(int j=0;j<GetTotalScalenodes();j++){
             int scalenode1 = j;
             int scalenode2 = j;
             if (NScaleDim>1){          
                scalenode1 = j / Nscalenode[1];
                scalenode2 = j % Nscalenode[1];
             }
             double alphastwopi = pow(GetAlphas(ScaleNode[i][scaleindex1][scalevar1][scalenode1],asmz)/TWOPI, Npow);
             for(int k=0;k<nxmax;k++){ 
                for(int l=0;l<NSubproc;l++){ 
                   //                   int x1bin = k % Nxtot1[i];
                   //                   int x2bin = k / Nxtot1[i];
                   //                   printf("%d %d %d : %f\n",x1bin,x2bin,k,SigmaTilde[i][scalevar][j][k][l]);
                   //                   printf(">> %d %d %d %d %g %g %g\n",i,j,k,l,SigmaTilde[i][scalevar][j][k][l],alphastwopi,PdfLc[i][j][k][l]);
                   Xsection[i] +=  SigmaTilde[i][scalevar][j][k][l] *  alphastwopi  *  PdfLc[i][scalenode2][k][l];
                }
             }
          }
       }
       return;
   }

   if(IsNLO()){
      if(fabs(1.-rescale)>1e-5){
         printf("fnloBlockB::CalcXsection: A posteriori change of scale not yet supported in usercode.\n");
         return;
      }
      double beta0=(11.*3.-2.*double(5))/3.; // = (11*CA-2*NF)/3  for NF = 5
      double mu0scale = 0.25;
//       for(int p=0;p<nscalebin;p++){ // scalebins         
//          double bweight =0.;
//          double t1 = 0.;
//          double t2 = 0.;
//          double aspow[nordmax];
//          for(int k=0;k<nxsum;k++){ // all x bins
//             for(int l=0;l<nsubproc;l++){ // subprocesses
//                double coeff = 0.;
//                switch(m){ //order in pert. series
//                case 1:
//                   coeff = weights[nbin][k][0][l][p]; // LO
//                   break;
//                case 2: 
//                   coeff = weights[nbin][k][1+fscale][l][p] + scfac*weights[nbin][k][0][l][p]  ; // NLO
//                   break;
//                default: printf("nord= %d not supported. Exit.\n",m); exit(1);
//                }
//                xsection[i][j][m-1] += coeff *  aspow[m-1]  *  pdflc[nbin][k][l][p];
      for(int i=0;i<BlockA2->GetNObsBin();i++){
         int nxmax = GetNxmax(i);
         int scaleindex1 = 0; // Even if we use multiple scales, mu_r is by convention the first scale -> index=0
         int scalevar1 = scalevar;
          if (NScaleDim>1){          
             scalevar1 = scalevar / Nscalevar[1];
          }
         for(int j=0;j<GetTotalScalenodes();j++){
             int scalenode1 = j;
             int scalenode2 = j;
             if (NScaleDim>1){          
                scalenode1 = j / Nscalenode[1];
                scalenode2 = j % Nscalenode[1];
             }
            double alphastwopi = pow(GetAlphas(ScaleNode[i][scaleindex1][scalevar1][scalenode1],asmz)/TWOPI, Npow);
            for(int k=0;k<nxmax;k++){ 
               for(int l=0;l<NSubproc;l++){                  
                  Xsection[i] +=  SigmaTilde[i][scalevar][j][k][l] *  alphastwopi  *  PdfLc[i][scalenode2][k][l];
               }
            }
         }
      }
      return;
   }
}

double fnloBlockB::GetSmallestX(int ObsBin){
   double smallestx = 1.0;

   int nxmax = GetNxmax(ObsBin);
   for(int scalevar=0;scalevar<Nscalevar[0];scalevar++ ){
      for(int j=0;j<Nscalenode[0];j++){
         for(int k=0;k<nxmax;k++){ 
            for(int l=0;l<NSubproc;l++){ 
               int x1bin = k % Nxtot1[ObsBin];
               int x2bin = k / Nxtot1[ObsBin];
               if((SigmaTilde[ObsBin][scalevar][j][k][l]!=0.)&&(XNode1[ObsBin][x1bin]<smallestx)){
                  smallestx = XNode1[ObsBin][x1bin];
               }
            }
         }
      }
   }
   return smallestx;
}

double fnloBlockB::GetSmallestX2(int ObsBin){
   double smallestx = 1.0;

   int nxmax = GetNxmax(ObsBin);
   for(int scalevar=0;scalevar<Nscalevar[0];scalevar++ ){
      for(int j=0;j<Nscalenode[0];j++){
         for(int k=0;k<nxmax;k++){ 
            for(int l=0;l<NSubproc;l++){ 
               int x1bin = k % Nxtot1[ObsBin];
               int x2bin = k / Nxtot1[ObsBin];
               if((SigmaTilde[ObsBin][scalevar][j][k][l]!=0.)&&(XNode2[ObsBin][x2bin]<smallestx)){
                  smallestx = XNode2[ObsBin][x2bin];
               }
            }
         }
      }
   }
   return smallestx;
}


//________________________________________________________________________________________________________________ //
void fnloBlockB::ResizeTable(vector<vector<vector<vector<vector<vector<vector<double > > > > > > >* v, int dim0 , int dim1, int dim2, int dim3, int dim4, int* dim5GetNxmaxFromDimI , int dim6 ){
  if ( dim0 > 0 ){
    if ( dim5GetNxmaxFromDimI[0] == 0 ) {
      v->resize(dim0);
      for ( int i= 0 ; i<dim0 ; i++){
	ResizeTable( &(v->at(i)) , dim1, dim2, dim3, dim4, GetNxmax(i), dim6 );
      }
    }
    else if ( dim5GetNxmaxFromDimI[0] != 0 ){
      cout << "Error in Resize Table. This is not yet implemented" << endl;
      exit(1);
    }
  } else{
    cout << "Error in Resize Table." << endl;
    exit(1);
  }
    
}


//________________________________________________________________________________________________________________ //
void fnloBlockB::ResizeTable(vector<vector<vector<vector<vector<vector<vector<double > > > > > > >* v, int dim0 , int dim1, int dim2, int dim3, int dim4, int dim5, int dim6 ){
  if ( dim0 > 0 ){
    v->resize(dim0);
    for ( int i= 0 ; i<dim0 ; i++){
      ResizeTable( &(v->at(i)) , dim1, dim2, dim3, dim4, dim5, dim6 );
    }
  } else{
    cout << "Error in Resize Table." << endl;
    exit(1);
  }
}


//________________________________________________________________________________________________________________ //
void fnloBlockB::ResizeTable( vector<vector<vector<vector<vector<vector<double > > > > > >* v, int dim0 , int dim1, int dim2, int dim3, int dim4, int dim5 ){
  if ( dim0 > 0 ){
    v->resize(dim0);
    for ( int i= 0 ; i<dim0 ; i++){
      ResizeTable( &(v->at(i)) , dim1, dim2, dim3, dim4, dim5 );
    }
  } else {
    cout << "Error in Resize Table." << endl;
    exit(1);
  }
}


//________________________________________________________________________________________________________________ //
void fnloBlockB::ResizeTable( vector<vector<vector<vector<vector<double > > > > >* v, int dim0 , int dim1, int dim2, int* dim3GetNxmaxFromDimI, int dim4 ){
  if ( dim0 > 0 ){
    if ( dim3GetNxmaxFromDimI[0] == 0 ) {
      v->resize(dim0);
      for ( int i= 0 ; i<dim0 ; i++){
	ResizeTable( &(v->at(i)) , dim1, dim2, GetNxmax(i), dim4 );
      }
    }
    else if ( dim3GetNxmaxFromDimI[0] != 0 ){
      cout << "Error in Resize Table. This is not yet implemented" << endl;
      exit(1);
    }
  } else{
    cout << "Error in Resize Table." << endl;
    exit(1);
  }
}


//________________________________________________________________________________________________________________ //
void fnloBlockB::ResizeTable( vector<vector<vector<vector<vector<double > > > > >* v, int dim0 , int* dim1GetNxmaxFromDimI, int dim2, int dim3, int dim4 ){
  if ( dim0 > 0 ){
    if ( dim1GetNxmaxFromDimI[0] == 0 ) {
      v->resize(dim0);
      for ( int i= 0 ; i<dim0 ; i++){
	ResizeTable( &(v->at(i)) , GetNxmax(i), dim2, dim3, dim4 );
      }
    }
    else if ( dim1GetNxmaxFromDimI[0] != 0 ){
      cout << "Error in Resize Table. This is not yet implemented" << endl;
      exit(1);
    }
  } else{
    cout << "Error in Resize Table." << endl;
    exit(1);
  }

}


//________________________________________________________________________________________________________________ //
void fnloBlockB::ResizeTable( vector<vector<vector<vector<vector<double > > > > >* v, int dim0 , int dim1, int dim2, int dim3, int dim4 ){
  if ( dim0 > 0 ){
    v->resize(dim0);
    for ( int i= 0 ; i<dim0 ; i++){
      ResizeTable( &(v->at(i)) , dim1, dim2, dim3, dim4 );
    }
  } else {
    cout << "Error in Resize Table." << endl;
    exit(1);
  }
}


//________________________________________________________________________________________________________________ //
void fnloBlockB::ResizeTable( vector<vector<vector<vector<double > > > >* v, int dim0 , int dim1, int* dim2GetNxmaxFromDimI, int dim3 ){
  if ( dim0 > 0 ){
    if ( dim2GetNxmaxFromDimI[0] == 0 ) {
      v->resize(dim0);
      for ( int i= 0 ; i<dim0 ; i++){
	ResizeTable( &(v->at(i)) , dim1, GetNxmax(i), dim3 );
      }
    }
    else if ( dim2GetNxmaxFromDimI[0] != 0 ){
      cout << "Error in Resize Table. This is not yet implemented" << endl;
      exit(1);
    }
  } else{
    cout << "Error in Resize Table." << endl;
    exit(1);
  }
  
}


//________________________________________________________________________________________________________________ //
void fnloBlockB::ResizeTable( vector<vector<vector<vector<double > > > >* v, int dim0 , int dim1, int dim2, int dim3 ){
  if ( dim0 > 0 ){
    v->resize(dim0);
    for ( int i= 0 ; i<dim0 ; i++){
      ResizeTable( &(v->at(i)) , dim1, dim2, dim3 );
    }
  } else {
    cout << "Error in Resize Table." << endl;
    exit(1);
  }
}


//________________________________________________________________________________________________________________ //
void fnloBlockB::ResizeTable( vector<vector<vector<double > > >* v, int dim0 , int* dim1GetNxmaxFromDimI, int dim2 ){
  if ( dim0 > 0 ){
    if ( dim1GetNxmaxFromDimI[0] == 0 ) {
      v->resize(dim0);
      for ( int i= 0 ; i<dim0 ; i++){
	ResizeTable( &(v->at(i)) , GetNxmax(i), dim2 );
      }
    }
    else if ( dim1GetNxmaxFromDimI[0] != 0 ){
      cout << "Error in Resize Table. This is not yet implemented" << endl;
      exit(1);
    }
  } else{
    cout << "Error in Resize Table." << endl;
    exit(1);
  }

}


//________________________________________________________________________________________________________________ //
void fnloBlockB::ResizeTable( vector<vector<vector<double > > >* v, int dim0 , int dim1, int dim2 ){
  if ( dim0 > 0 ){
    v->resize(dim0);
    for ( int i= 0 ; i<dim0 ; i++){
      ResizeTable( &(v->at(i)) , dim1, dim2 );
    }
  } else {
    cout << "Error in Resize Table." << endl;
    exit(1);
  }
}


//________________________________________________________________________________________________________________ //
void fnloBlockB::ResizeTable( vector<vector<double > >*  v, int dim0 , int* dim1GetNxmaxFromDimI ){
  if ( dim0 > 0 ){
    if ( dim1GetNxmaxFromDimI[0] == 0 ) {
      v->resize(dim0);
      for ( int i= 0 ; i<dim0 ; i++){
	ResizeTable( &(v->at(i)) , GetNxmax(i) );
      }
    }
    else if ( dim1GetNxmaxFromDimI[0] != 0 ){
      cout << "Error in Resize Table. This is not yet implemented" << endl;
      exit(1);
    }
  } else{
    cout << "Error in Resize Table." << endl;
    exit(1);
  }
}


//________________________________________________________________________________________________________________ //
void fnloBlockB::ResizeTable( vector<vector<double > >*  v, int dim0 , int dim1 ){
  if ( dim0 > 0 ){
    v->resize(dim0);
    for ( int i= 0 ; i<dim0 ; i++){
      ResizeTable( &(v->at(i)) , dim1 );
    }
  } else {
    cout << "Error in Resize Table." << endl;
    exit(1);
  }
}


//________________________________________________________________________________________________________________ //
void fnloBlockB::ResizeTable( vector<double >* v, int dim0 ){
  if ( dim0 > 0 ){
    v->resize(dim0);
  }
  else{
    cout << "Error in Resize Table." << endl;
    exit(1);
  }
}



//________________________________________________________________________________________________________________ //
int fnloBlockB::ReadTable(vector<vector<vector<vector<vector<vector<vector<double > > > > > > >* v, istream *table ){
  int nn = 0;
  for(unsigned int i0=0;i0<v->size();i0++){
    for(unsigned int i1=0;i1<v->at(i0).size();i1++){
      for(unsigned int i2=0;i2<v->at(i0)[i1].size();i2++){
 	for(unsigned int i3=0;i3<v->at(i0)[i1][i2].size();i3++){
	  for(unsigned int i4=0;i4<v->at(i0)[i1][i2][i3].size();i4++){
	    for(unsigned int i5=0;i5<v->at(i0)[i1][i2][i3][i4].size();i5++){
	      for(unsigned int i6=0;i6<v->at(i0)[i1][i2][i3][i4][i5].size();i6++){
		*table >> v->at(i0)[i1][i2][i3][i4][i5][i6];
		nn++;
	      }
	    }
	  }
	}
      }
    }
  }
  return nn;
}

//________________________________________________________________________________________________________________ //
int fnloBlockB::ReadTable(vector<vector<vector<vector<vector<vector<double > > > > > >* v, istream *table ){
  int nn = 0;
  for(unsigned int i0=0;i0<v->size();i0++){
    for(unsigned int i1=0;i1<v->at(i0).size();i1++){
      for(unsigned int i2=0;i2<v->at(i0)[i1].size();i2++){
 	for(unsigned int i3=0;i3<v->at(i0)[i1][i2].size();i3++){
	  for(unsigned int i4=0;i4<v->at(i0)[i1][i2][i3].size();i4++){
	    for(unsigned int i5=0;i5<v->at(i0)[i1][i2][i3][i4].size();i5++){
	      *table >> v->at(i0)[i1][i2][i3][i4][i5];
	      nn++;
	    }
	  }
	}
      }
    }
  }
  return nn;
}

//________________________________________________________________________________________________________________ //
int fnloBlockB::ReadTable(vector<vector<vector<vector<vector<double > > > > >* v, istream *table ){
  int nn = 0;
  for(unsigned int i0=0;i0<v->size();i0++){
    for(unsigned int i1=0;i1<v->at(i0).size();i1++){
      for(unsigned int i2=0;i2<v->at(i0)[i1].size();i2++){
 	for(unsigned int i3=0;i3<v->at(i0)[i1][i2].size();i3++){
	  for(unsigned int i4=0;i4<v->at(i0)[i1][i2][i3].size();i4++){
	    *table >> v->at(i0)[i1][i2][i3][i4];
	    nn++;
	  }
	}
      }
    }
  }
  return nn;
}

//________________________________________________________________________________________________________________ //
int fnloBlockB::ReadTable(vector<vector<vector<vector<double > > > >* v, istream *table ){
  int nn = 0;
  for(unsigned int i0=0;i0<v->size();i0++){
    for(unsigned int i1=0;i1<v->at(i0).size();i1++){
      for(unsigned int i2=0;i2<v->at(i0)[i1].size();i2++){
 	for(unsigned int i3=0;i3<v->at(i0)[i1][i2].size();i3++){
	  *table >> v->at(i0)[i1][i2][i3];
	  nn++;
	}
      }
    }
  }
  return nn;
}

//________________________________________________________________________________________________________________ //
int fnloBlockB::ReadTable(vector<vector<vector<double > > >* v, istream *table ){
  int nn = 0;
  for(unsigned int i0=0;i0<v->size();i0++){
    for(unsigned int i1=0;i1<v->at(i0).size();i1++){
      for(unsigned int i2=0;i2<v->at(i0)[i1].size();i2++){
	*table >> v->at(i0)[i1][i2];
	nn++;
      }
    }
  }
  return nn;
}

//________________________________________________________________________________________________________________ //
int fnloBlockB::ReadTable(vector<vector<double > >* v, istream *table ){
  int nn = 0;
  for(unsigned int i0=0;i0<v->size();i0++){
    for(unsigned int i1=0;i1<v->at(i0).size();i1++){
      *table >> v->at(i0)[i1];
      nn++;
    }
  }
  return nn;
}

//________________________________________________________________________________________________________________ //
int fnloBlockB::ReadTable(vector<double >* v, istream *table ){
  int nn = 0;
  for(unsigned int i0=0;i0<v->size();i0++){
    *table >> v->at(i0);
    nn++;
  }
  return nn;
}


//________________________________________________________________________________________________________________ //
int fnloBlockB::WriteTable(vector<vector<vector<vector<vector<vector<vector<double > > > > > > >* v, ostream *table , bool DivByNevt , int Nevt  ){
  int nn = 0;
  for(unsigned int i0=0;i0<v->size();i0++){
    for(unsigned int i1=0;i1<v->at(i0).size();i1++){
      for(unsigned int i2=0;i2<v->at(i0)[i1].size();i2++){
 	for(unsigned int i3=0;i3<v->at(i0)[i1][i2].size();i3++){
	  for(unsigned int i4=0;i4<v->at(i0)[i1][i2][i3].size();i4++){
	    for(unsigned int i5=0;i5<v->at(i0)[i1][i2][i3][i4].size();i5++){
	      for(unsigned int i6=0;i6<v->at(i0)[i1][i2][i3][i4][i5].size();i6++){
		if( DivByNevt && Nevt>0){
 		  *table << v->at(i0)[i1][i2][i3][i4][i5][i6] / Nevt << endl;
 		}else{
 		  *table << v->at(i0)[i1][i2][i3][i4][i5][i6] << endl;
		}
		nn++;
	      }
	    }
	  }
	}
      }
    }
  }
  return nn;
}

//________________________________________________________________________________________________________________ //
int fnloBlockB::WriteTable(vector<vector<vector<vector<vector<vector<double > > > > > >* v, ostream *table , bool DivByNevt , int Nevt ){
  int nn = 0;
  for(unsigned int i0=0;i0<v->size();i0++){
    for(unsigned int i1=0;i1<v->at(i0).size();i1++){
      for(unsigned int i2=0;i2<v->at(i0)[i1].size();i2++){
 	for(unsigned int i3=0;i3<v->at(i0)[i1][i2].size();i3++){
	  for(unsigned int i4=0;i4<v->at(i0)[i1][i2][i3].size();i4++){
	    for(unsigned int i5=0;i5<v->at(i0)[i1][i2][i3][i4].size();i5++){
	      if( DivByNevt && Nevt>0){
 		*table << v->at(i0)[i1][i2][i3][i4][i5] / Nevt << endl;
	      }else{
		*table << v->at(i0)[i1][i2][i3][i4][i5] << endl;
	      }
	      nn++;
	    }
	  }
	}
      }
    }
  }
  return nn;
}

//________________________________________________________________________________________________________________ //
int fnloBlockB::WriteTable(vector<vector<vector<vector<vector<double > > > > >* v, ostream *table , bool DivByNevt , int Nevt ){
  int nn = 0;
  for(unsigned int i0=0;i0<v->size();i0++){
    for(unsigned int i1=0;i1<v->at(i0).size();i1++){
      for(unsigned int i2=0;i2<v->at(i0)[i1].size();i2++){
 	for(unsigned int i3=0;i3<v->at(i0)[i1][i2].size();i3++){
	  for(unsigned int i4=0;i4<v->at(i0)[i1][i2][i3].size();i4++){
	    if( DivByNevt && Nevt>0){
 	      *table << v->at(i0)[i1][i2][i3][i4] / Nevt << endl;
		}else{
	      *table << v->at(i0)[i1][i2][i3][i4] << endl;
	    }
	    nn++;
	  }
	}
      }
    }
  }
  return nn;
}

//________________________________________________________________________________________________________________ //
int fnloBlockB::WriteTable(vector<vector<vector<vector<double > > > >* v, ostream *table , bool DivByNevt , int Nevt ){
  int nn = 0;
  for(unsigned int i0=0;i0<v->size();i0++){
    for(unsigned int i1=0;i1<v->at(i0).size();i1++){
      for(unsigned int i2=0;i2<v->at(i0)[i1].size();i2++){
 	for(unsigned int i3=0;i3<v->at(i0)[i1][i2].size();i3++){
	  if( DivByNevt && Nevt>0){
 	    *table << v->at(i0)[i1][i2][i3] / Nevt << endl;
	  }else{
	    *table << v->at(i0)[i1][i2][i3] << endl;
	  }
	  nn++;
	}
      }
    }
  }
  return nn;
}

//________________________________________________________________________________________________________________ //
int fnloBlockB::WriteTable(vector<vector<vector<double > > >* v, ostream *table , bool DivByNevt , int Nevt ){
  int nn = 0;
  for(unsigned int i0=0;i0<v->size();i0++){
    for(unsigned int i1=0;i1<v->at(i0).size();i1++){
      for(unsigned int i2=0;i2<v->at(i0)[i1].size();i2++){
	if( DivByNevt && Nevt>0){
	  *table << v->at(i0)[i1][i2] / Nevt << endl;
	}else{
	  *table << v->at(i0)[i1][i2] << endl;
	}
	nn++;
      }
    }
  }
  return nn;
}

//________________________________________________________________________________________________________________ //
int fnloBlockB::WriteTable(vector<vector<double > >* v, ostream *table , bool DivByNevt , int Nevt ){
  int nn = 0;
  for(unsigned int i0=0;i0<v->size();i0++){
    for(unsigned int i1=0;i1<v->at(i0).size();i1++){
      if( DivByNevt && Nevt>0){
	*table << v->at(i0)[i1] / Nevt << endl;
      }else{
	*table << v->at(i0)[i1] << endl;
      }
      nn++;
    }
  }
  return nn;
}



//________________________________________________________________________________________________________________ //
int fnloBlockB::WriteTable(vector<double >* v, ostream *table , bool DivByNevt , int Nevt ){
  int nn = 0;
  for(unsigned int i0=0;i0<v->size();i0++){
    if( DivByNevt && Nevt>0){
      *table << v->at(i0) / Nevt << endl;
    }else{
      *table << v->at(i0) << endl;
    }
    nn++;
  }
  return nn;
}


//________________________________________________________________________________________________________________ //
int fnloBlockB::WriteFlexibleTable(vector<vector<vector<vector<vector<vector<vector< double > > > > > > >* v, ostream *table , bool DivByNevt , int Nevt , bool nProcLast ){
   int nn = 1;
   *table << v->size() << endl;
   for(unsigned int i0=0;i0<v->size();i0++){
      nn += WriteFlexibleTable( &(v->at(i0)) , table , DivByNevt, Nevt , nProcLast );
   }
   return nn;
}
//________________________________________________________________________________________________________________ //
int fnloBlockB::WriteFlexibleTable(vector<vector<vector<vector<vector<vector<double > > > > > >* v, ostream *table , bool DivByNevt , int Nevt , bool nProcLast ){
   int nn = 1;
   *table << v->size() << endl;
   for(unsigned int i0=0;i0<v->size();i0++){
      nn += WriteFlexibleTable( &(v->at(i0)) , table , DivByNevt, Nevt , nProcLast );
   }
   return nn;
}
//________________________________________________________________________________________________________________ //
int fnloBlockB::WriteFlexibleTable(vector<vector<vector<vector<vector<double > > > >  >* v, ostream *table , bool DivByNevt , int Nevt , bool nProcLast ){
   int nn = 1;
   *table << v->size() << endl;
   for(unsigned int i0=0;i0<v->size();i0++){
      nn += WriteFlexibleTable( &(v->at(i0)) , table , DivByNevt, Nevt , nProcLast );
   }
   return nn;
}
//________________________________________________________________________________________________________________ //
int fnloBlockB::WriteFlexibleTable(vector<vector<vector<vector<double > >  > >* v, ostream *table , bool DivByNevt , int Nevt , bool nProcLast ){
   int nn = 1;
   *table << v->size() << endl;
   for(unsigned int i0=0;i0<v->size();i0++){
      nn += WriteFlexibleTable( &(v->at(i0)) , table , DivByNevt, Nevt , nProcLast );
   }
   return nn;
}
//________________________________________________________________________________________________________________ //
int fnloBlockB::WriteFlexibleTable(vector<vector<vector<double > > >* v, ostream *table , bool DivByNevt , int Nevt , bool nProcLast ){
   int nn = 1;
   *table << v->size() << endl;
   for(unsigned int i0=0;i0<v->size();i0++){
      nn += WriteFlexibleTable( &(v->at(i0)) , table , DivByNevt, Nevt , nProcLast );
   }
   return nn;
}
//________________________________________________________________________________________________________________ //
int fnloBlockB::WriteFlexibleTable(vector<vector<double > >* v, ostream *table , bool DivByNevt , int Nevt , bool nProcLast ){
   int nn = 1;
   *table << v->size() << endl;
   for(unsigned int i0=0;i0<v->size();i0++){
      nn += WriteFlexibleTable( &(v->at(i0)) , table , DivByNevt, Nevt , nProcLast );
   }
   return nn;
}

//________________________________________________________________________________________________________________ //
int fnloBlockB::WriteFlexibleTable(vector<double >* v, ostream *table , bool DivByNevt , int Nevt , bool nProcLast ){
   int nn = 1;
   if ( !nProcLast )*table << v->size() << endl;
   for(unsigned int i0=0;i0<v->size();i0++){
      if( DivByNevt && Nevt>0)	*table << v->at(i0) / Nevt << endl;
      else			*table << v->at(i0) << endl;
      nn++;
   }
   return nn;
}

//________________________________________________________________________________________________________________ //
int fnloBlockB::ReadFlexibleVector(vector<vector<vector<vector<vector<vector<vector<double > > > > > > >* v, istream *table , bool nProcLast ){
   int nn = 0;
   int size = 0;
   *table >> size; nn++;
   v->resize(size);
   for(unsigned int i0=0;i0<v->size();i0++){
      nn += ReadFlexibleVector(&(v->at(i0)),table,nProcLast);
   }
   return nn;
}
//________________________________________________________________________________________________________________ //
int fnloBlockB::ReadFlexibleVector(vector<vector<vector<vector<vector<vector<double > > > > > >* v, istream *table , bool nProcLast ){
   int nn = 0;
   int size = 0;
   *table >> size; nn++;
   v->resize(size);
   for(unsigned int i0=0;i0<v->size();i0++){
      nn += ReadFlexibleVector(&(v->at(i0)),table,nProcLast);
   }
   return nn;
}
//________________________________________________________________________________________________________________ //
int fnloBlockB::ReadFlexibleVector(vector<vector<vector<vector<vector<double > > > > >* v, istream *table , bool nProcLast ){
   int nn = 0;
   int size = 0;
   *table >> size; nn++;
   v->resize(size);
   for(unsigned int i0=0;i0<v->size();i0++){
      nn += ReadFlexibleVector(&(v->at(i0)),table,nProcLast);
   }
   return nn;
}
//________________________________________________________________________________________________________________ //
int fnloBlockB::ReadFlexibleVector(vector<vector<vector<vector<double > > > >* v, istream *table , bool nProcLast ){
   int nn = 0;
   int size = 0;
   *table >> size; nn++;
   v->resize(size);
   for(unsigned int i0=0;i0<v->size();i0++){
      nn += ReadFlexibleVector(&(v->at(i0)),table,nProcLast);
   }
   return nn;
}
//________________________________________________________________________________________________________________ //
int fnloBlockB::ReadFlexibleVector(vector<vector<vector<double > > >* v, istream *table , bool nProcLast ){
   int nn = 0;
   int size = 0;
   *table >> size; nn++;
   v->resize(size);
   for(unsigned int i0=0;i0<v->size();i0++){
      nn += ReadFlexibleVector(&(v->at(i0)),table,nProcLast);
   }
   return nn;
}
//________________________________________________________________________________________________________________ //
int fnloBlockB::ReadFlexibleVector(vector<vector<double > >* v, istream *table , bool nProcLast ){
   int nn = 0;
   int size = 0;
   *table >> size; nn++;
   v->resize(size);
   for(unsigned int i0=0;i0<v->size();i0++){
      nn += ReadFlexibleVector(&(v->at(i0)),table,nProcLast);
   }
   return nn;
}
//________________________________________________________________________________________________________________ //
int fnloBlockB::ReadFlexibleVector(vector<double >* v, istream *table , bool nProcLast ){
   int nn = 0;
   if ( !nProcLast ) {
      int size = 0;
      *table >> size; nn++;
      v->resize(size);
   }
   else {
      v->resize(NSubproc);
   }
   for(unsigned int i0=0;i0<v->size();i0++){
      *table >> v->at(i0);
      nn++;
   }
   return nn;
}



//________________________________________________________________________________________________________________ //
void fnloBlockB::AddTableToAnotherTable( vector<vector<vector<vector<vector<vector<vector<double > > > > > > >* vSum, vector<vector<vector<vector<vector<vector<vector<double > > > > > > >* vAdd, double w1, double w2){
  if ( vSum->size() != vAdd->size() ) {cout<<"Error in fnloBlockB::AddTableToAnotherTable. Cannot add tables with different size."<<endl; return;}
  for ( int i = 0 ; i<vSum->size() ; i++ ){
    AddTableToAnotherTable( &(vSum->at(i)), &(vAdd->at(i)), w1 , w2  );
  }
}


//________________________________________________________________________________________________________________ //
void fnloBlockB::AddTableToAnotherTable( vector<vector<vector<vector<vector<vector<double > > > > > >* vSum, vector<vector<vector<vector<vector<vector<double > > > > > >* vAdd, double w1, double w2){

  if ( vSum->size() != vAdd->size() ) {cout<<"Error in fnloBlockB::AddTableToAnotherTable. Cannot add tables with different size."<<endl; return;}
  for ( int i = 0 ; i<vSum->size() ; i++ ){
    AddTableToAnotherTable( &(vSum->at(i)), &(vAdd->at(i)), w1 , w2  );
  }
}

//________________________________________________________________________________________________________________ //


void fnloBlockB::AddTableToAnotherTable( vector<vector<vector<vector<vector<double > > > > >* vSum, vector<vector<vector<vector<vector<double > > > > >* vAdd, double w1, double w2){
  if ( vSum->size() != vAdd->size() ) {cout<<"Error in fnloBlockB::AddTableToAnotherTable. Cannot add tables with different size."<<endl; return;}
  for ( int i = 0 ; i<vSum->size() ; i++ ){
    AddTableToAnotherTable( &(vSum->at(i)), &(vAdd->at(i)), w1 , w2  );
  }
}


//________________________________________________________________________________________________________________ //
void fnloBlockB::AddTableToAnotherTable( vector<vector<vector<vector<double > > > >* vSum, vector<vector<vector<vector<double > > > >* vAdd, double w1, double w2){
  if ( vSum->size() != vAdd->size() ) {cout<<"Error in fnloBlockB::AddTableToAnotherTable. Cannot add tables with different size."<<endl; return;}
  for ( int i = 0 ; i<vSum->size() ; i++ ){
    AddTableToAnotherTable( &(vSum->at(i)), &(vAdd->at(i)), w1 , w2  );
  }
}

//________________________________________________________________________________________________________________ //

void fnloBlockB::AddTableToAnotherTable( vector<vector<vector<double > > >* vSum, vector<vector<vector<double > > >* vAdd, double w1, double w2){

  if ( vSum->size() != vAdd->size() ) {cout<<"Error in fnloBlockB::AddTableToAnotherTable. Cannot add tables with different size."<<endl; return;}
  for ( int i = 0 ; i<vSum->size() ; i++ ){
    AddTableToAnotherTable( &(vSum->at(i)), &(vAdd->at(i)), w1 , w2  );
  }
}


//________________________________________________________________________________________________________________ //
void fnloBlockB::AddTableToAnotherTable( vector<vector<double > >* vSum, vector<vector<double > >* vAdd, double w1, double w2){
  if ( vSum->size() != vAdd->size() ) {cout<<"Error in fnloBlockB::AddTableToAnotherTable. Cannot add tables with different size."<<endl; return;}
  for ( int i = 0 ; i<vSum->size() ; i++ ){
    AddTableToAnotherTable( &(vSum->at(i)), &(vAdd->at(i)), w1 , w2  );
  }
}


//________________________________________________________________________________________________________________ //
void fnloBlockB::AddTableToAnotherTable( vector<double >* vSum, vector<double >* vAdd, double w1, double w2){
  if ( vSum->size() != vAdd->size() ) {cout<<"Error in fnloBlockB::AddTableToAnotherTable. Cannot add tables with different size."<<endl; return;}
  for ( int i = 0 ; i<vSum->size() ; i++ ){
    (*vSum)[i] =  w1*(*vSum)[i] + w2*(*vAdd)[i]; 
  }
}


//________________________________________________________________________________________________________________ //


double fnloBlockB::GetAlphas(double Q, double alphasMZ){

   static const double TWOPI = 6.28318530717958647692528;
   static const double TWOPISQR = 39.47841760435743447533796;
   static const int NF = 5;
   // DB: Todo: Update this number!
   static const double MZ = 91.1882;

   double BETA0 =  (11. - 2./3.*NF); // The beta coefficients of the QCD beta function
   double BETA1 =  (51. - 19./3.*NF);
 

   // This is from NLOJET++, alpha.cc
   double res = alphasMZ;
   double b0 = BETA0/TWOPI;
   double w = 1.0 + b0*alphasMZ*log(Q/MZ);
   res /= w;
   double b1 = BETA1/TWOPISQR;
   res *= 1.0 - alphasMZ*b1/b0*log(w)/w;
   return res;
   
}


//________________________________________________________________________________________________________________ //


void fnloBlockB::Print(){
  printf("\n **************** FastNLO Table: BlockB ****************\n\n");
  printf(" B   BlockA2->GetNObsBin()         %d\n",BlockA2->GetNObsBin());
  printf(" B   IXsectUnits                   %d\n",IXsectUnits);
  printf(" B   IDataFlag                     %d\n",IDataFlag);
  printf(" B   IAddMultFlag                  %d\n",IAddMultFlag);
  printf(" B   IContrFlag1                   %d\n",IContrFlag1);
  printf(" B   IContrFlag2                   %d\n",IContrFlag2);
  printf(" B   IContrFlag3 (always 0)        %d\n",IContrFlag3);
  printf(" B   NScaleDep                     %d\n",NScaleDep);
  for(int i=0;i<CtrbDescript.size();i++){
    printf(" B   CtrbDescript[%d]               %s\n",i,CtrbDescript[i].data());
  }
  //printf(" B   NCodeDescr                    %d\n",NCodeDescr);
  for(int i=0;i<CodeDescript.size();i++){
    printf(" B   CodeDescript[%d]               %s\n",i,CodeDescript[i].data());
  }

  if(IDataFlag==1){
    printf(" B   Nuncorrel                     %d\n",Nuncorrel);
    printf(" B   Ncorrel                       %d\n",Ncorrel);
    printf(" B   NErrMatrix                    %d\n",NErrMatrix);
    printf(" B   some more output could be printed here (IDataFlag==1).\n");
  }
  if(IAddMultFlag==1){
    printf(" B   some more output could be printed here (IAddMultFlag==1).\n");
  }

  if(!(IDataFlag==1) && !(IAddMultFlag==1)){ // that's the usual case
    printf(" B   IRef                          %d\n",IRef);
    printf(" B   IScaleDep                     %d\n",IScaleDep);
    printf(" B   Nevt                          %u\n",Nevt);
    printf(" B   Nevt                          %i\n",Nevt);
    printf(" B   Nevt                          %d\n",Nevt);
    printf(" B   Nevt                          %e\n",Nevt);
    printf(" B   Nevt                          %e\n",Nevt*1.);
    printf(" B   Nevt                          %.4e\n",Nevt);
    printf(" B   Nevt                          %.e\n",Nevt*1.);
    printf(" B   Npow                          %d\n",Npow);
    printf(" B   NPDF                          %d\n",NPDF);
    if(NPDF>0){
      for(int i=0;i<NPDF;i++){
	printf(" B    - NPDFPDG[%d]                 %d\n",i,NPDFPDG[i]);
      }      
    }
    printf(" B   NPDFDim                       %d\n",NPDFDim);
    printf(" B   NFragFunc                     %d\n",NFragFunc);
    if(NFragFunc>0){
      for(int i=0;i<NFragFunc;i++){
	printf(" B    - NFFPDG[%d]               %d\n",i,NFFPDG[i]);
      }      
    } 
    printf(" B   NFFDim                        %d\n",NFFDim);
    printf(" B   NSubproc                      %d\n",NSubproc);
    printf(" B   IPDFdef1                      %d\n",IPDFdef1);
    printf(" B   IPDFdef2                      %d\n",IPDFdef2);
    printf(" B   IPDFdef3                      %d\n",IPDFdef3);
    printf(" B   Nxtot1[0-%d]             ",BlockA2->GetNObsBin());
    for(int i=0;i<BlockA2->GetNObsBin();i++){
      printf("%d ,",Nxtot1[i]);
    } 
    printf(" B   \n");

//     for(int i=0;i<BlockA2->GetNObsBin();i++){
//       printf(" B    XNode1[%d]             ",i);
//       for(int j=0;j<Nxtot1[i];j++){
// 	printf(" B   %8.4f ,",XNode1[i][j]);
//       } 
//       printf(" B   \n");
//     }
    printf(" B   if (NPDFDim==2), you could print xnodes2 here. (NPDFDim = %d)\n",NPDFDim);
    printf(" B   if (NFragFunc>0), you could print xnodes2 here. (NFragFunc = %d)\n",NFragFunc);
    printf(" B   NScales                       %d\n",NScales);
    for(int i=0;i<NScales;i++){
      printf(" B    - Iscale[%d]                  %d\n",i,Iscale[i]);
    }
    printf(" B   NScaleDim                     %d\n",NScaleDim);
    for(int i=0;i<NScaleDim;i++){
       //printf(" B    -  NscaleDescript[%d]         %d\n",i,NscaleDescript[i]);
       for(int j=0;j<ScaleDescript[i].size();j++){
	printf(" B    -  - ScaleDescript[%d][%d]     %s\n",i,j,ScaleDescript[i][j].data());
      }
       if ( NScaleDep != 3 ) {
	  printf(" B    - Nscalenode[%d]              %d\n",i,Nscalenode[i]);
	  printf(" B    - Nscalevar[%d]               %d\n",i,Nscalevar[i]);
	  for(int j=0;j<Nscalevar[i];j++){
	     printf(" B    -  - ScaleFac[%d][%d]          %6.4f\n",i,j,ScaleFac[i][j]);
	  }
       }
    }
    printf(" B   No printing of ScaleNode implemented yet.\n");
    printf(" B   No printing of SigmaTilde implemented yet.\n");
    if ( NScaleDep == 2 )  
      printf(" B   NScaleDep == 2 :              yes.\n");
    if ( NScaleDep == 2 ) {
      printf(" B   No printing of SigmaTilde2Scales, and Scale2Nodes, etc... implemented yet.\n");
    }
    if ( NScaleDep == 3 ) {
      printf(" B   NscalenodeScale1              %d\n",NscalenodeScale1);
      printf(" B   NscalenodeScale2              %d\n",NscalenodeScale2);
    }    

  }
  printf("\n *******************************************************\n\n");

  
}


//________________________________________________________________________________________________________________ //

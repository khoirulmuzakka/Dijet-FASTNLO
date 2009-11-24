{
   fnloTableUser* table;
   if(table==0){
      cout << "Make new table object." <<endl;
      //      table = new fnloTableUser("/afs/cern.ch/user/t/tkluge/h1/fastNLO/trunk/v2.0/author/hadron/output/test-hhc-born-2jet.tab");
      table = new fnloTableUser("/afs/cern.ch/user/t/tkluge/h1/fastNLO/trunk/v2.0/author/hadron/output/test-hhc-2jet.tab");
      //      table->ReadPDF("MRST98nlo.LHgrid");
      table->ReadPDF("cteq6.LHpdf");
      //      table->ReadPDF("cteq5m.LHgrid");
      table->SetPDFSet(0);
      table->SetAlphasMZ(0.118);
   }
   table->ReadTable();
   //   for(int scalevar=0;scalevar<table->GetNScaleVar();scalevar++){
      for(int scalevar=0;scalevar<1;scalevar++){
         printf("Scale variation #%d\n",scalevar+1);
      table->SetScalevar(scalevar);     
      int murdim = table->GetIscale(0);
      int murscale = table->GetScale1index(scalevar);
      printf("mu_r = %s multiplied by %4.2f  \n",table->GetScaleDescript(murdim,0).c_str(),table->GetScaleFac(murdim,murscale));
      int mufdim = table->GetIscale(1);
      int mufscale = table->GetScale2index(scalevar);
      printf("mu_f = %s multiplied by %4.2f  \n",table->GetScaleDescript(mufdim,0).c_str(),table->GetScaleFac(mufdim,mufscale));

      table->FillPDFCache();
      table->CalcXsection();
      int nbins = table->GetNObsBin();
      printf("                       +---------------------------  LO  ----------------+     +--------------------  NLO  --------------+ \n");
      printf("Bin#        pT[GeV]     reference[nB/Gev]     fastNLO[nB/GeV]     delta[%%]    reference[nB/GeV]     fastNLO[nB/GeV]     delta[%%]\n");
      for(int i =0;i<nbins;i++){
         printf("%2d      %4.1f - %4.1f      ",
                i+1,
                table->GetLoBin(i,0),
                table->GetUpBin(i,0));

         if(table->GetXsectionRef(i, 1)!=0.){
            printf("%9.3g         %9.3g            %5.2f",
                   table->GetXsectionRef(i, 1),
                   table->GetXsection(i,1),
                   100.*(table->GetXsectionRef(i, 1)-table->GetXsection(i,1))/table->GetXsectionRef(i, 1));
         }else{
            printf("                                    ");
         }
         if(table->GetXsectionRef(i, 2)!=0.){
            printf("          %9.3g         %9.3g            %5.2f\n",
                   table->GetXsectionRef(i, 2),
                   table->GetXsection(i,2),
                   100.*(table->GetXsectionRef(i, 2)-table->GetXsection(i,2))/table->GetXsectionRef(i, 2));
         }else{
            printf("\n");
         }
      }
      printf("\n\n",scalevar);

    }
 
//    int nbins = table->GetNObsBin();
//    for(int i =0;i<nbins;i++){
//       printf("Bin #%d  smallest x=%f\n",i,table->GetSmallestX(0,i));
//    }
//    for(int i =0;i<nbins;i++){
//       printf("Bin #%d  smallest x2=%f\n",i,table->GetSmallestX2(0,i));
//    }

  
}

{
   double pthigh[22];
   pthigh[0] = 5.0;
   pthigh[1] = 6.0;
   pthigh[2] = 7.0;
   pthigh[3] = 8.0;
   pthigh[4] = 9.0;
   pthigh[5] = 10.0;
   pthigh[6] = 11.0;
   pthigh[7] = 12.0;
   pthigh[8] = 13.0;
   pthigh[9] = 14.0;
   pthigh[10] = 16.0;
   pthigh[11] = 18.0;
   pthigh[12] = 20.0;
   pthigh[13] = 22.0;
   pthigh[14] = 24.0;
   pthigh[15] = 27.0;
   pthigh[16] = 30.0;
   pthigh[17] = 33.0;
   pthigh[18] = 36.0;
   pthigh[19] = 40.0;
   pthigh[20] = 44.0;
   pthigh[21] = 50.0;

   double ptbinning[5];
   ptbinning[0] = 7.0;
   ptbinning[1] = 11.0;
   ptbinning[2] = 18.0;
   ptbinning[3] = 30.0;
   ptbinning[4] = 50.0;


   TH1 *h = new TH1D("Q1","Q1",4,ptbinning);
   TH1 *h2 = new TH1D("delta","delta",4,ptbinning);

   TIncljets jets;
   jets.ReadTable("fnh2001-tst01.txt");
   jets.ReadPDF("cteq5m.LHgrid");
   jets.SetPDFSet(0);
   jets.FillPDFCache(1.);
   jets.ResetXsection();
   jets.CalcXsection(0.118,0);
   jets.Rebin();

   for(int i=0;i<4;i++){
      h->SetBinContent(i+1,jets.GetRebinned(0,i,0));
   }
         h->DrawCopy("");

//    jets.ResetReference();
//    jets.ReadReference("nlojet-reference.lo.txt",0);

//    for(int i=0;i<4;i++){
//       h2->SetBinContent(i+1,jets.GetReference(0,i,0));
//    }
//   h2->SetLineColor(2);
   //   h2->DrawCopy("same");

//    h2->Add(h,h2,1.,-1.);
//    h2->Divide(h);
//    h2->Scale(100.);
//    h2->Draw();

   //   h2->SetMinimum(-20.);
   //   h2->SetMaximum(+20.);
   
   //   gSystem->Exit(0);
}

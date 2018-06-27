int tempDraw_res(int har=1) {

   gStyle->SetErrorX(0);
   gStyle->SetOptStat(0);

   TLatex text;
   text.SetNDC(1);
   text.SetTextFont(43);
   text.SetTextSize(18);

   char name [500];
   char name1[500];

   double ryL[3][2] = {
      {-0.1, 0},
      {0   , 0},
      {-0.1, 0.1},
   };
   double ryU[3][2] = {
      {0, 0.32},
      {0, 0.1},
      {0, 0.4},
   };

   TFile* f1 = new TFile("../Output_zdcPix_flat_wtrk_doAna.root");
   TH1D* h1[2];

   sprintf(name , "h_ResCent_mix_har%d", har-1);
   sprintf(name1, "h_ResCent_mix_har%d_Clone", har-1);
   h1[0] = (TH1D*)((TH1D*)f1->Get(name))->Clone(name1);
   sprintf(name1, "h_ResCent_mix_har%d_res", har-1);
   h1[1] = (TH1D*)h1[0]->Clone(name1);
   for (int ib=0; ib<h1[1]->GetNbinsX(); ib++) {
      h1[1]->SetBinContent(ib+1, sqrt(fabs(h1[0]->GetBinContent(ib+1))));
      h1[1]->SetBinError  (ib+1, h1[0]->GetBinError(ib+1)/(sqrt(2)*h1[1]->GetBinContent(ib+1)));
   }

   std::stringstream ss_har; ss_har << har;
   for (int ih=0; ih<2; ih++) {
      formatHist1(h1[ih]);
      h1[ih]->SetMarkerStyle(21);
      h1[ih]->SetMarkerSize (1.1);
      h1[ih]->SetMarkerColor(kGreen+2);
      h1[ih]->SetLineColor  (kGreen+2);

      std::string lab1;
      if (ih==0) {
         if (har==1) lab1 = "<cos(#Psi^{A}_{" + ss_har.str() + "}" + " - #Psi^{C}_{" + ss_har.str() + "})>";
         else        lab1 = "<(" + ss_har.str() + "( #Psi^{A}_{" + ss_har.str() + "}" + " - #Psi^{C}_{" + ss_har.str() + "}))>";
      }
      h1[ih]->GetXaxis()->SetTitle("Centrality(%)");
      if (ih==0) h1[ih]->GetYaxis()->SetTitle(lab1.c_str());
      else       h1[ih]->GetYaxis()->SetTitle("Resolution");
      h1[ih]->GetXaxis()->SetTitleOffset(1.4);
      h1[ih]->GetYaxis()->SetTitleOffset(1.6);
      h1[ih]->GetYaxis()->SetRangeUser(ryL[har-1][ih], ryU[har-1][ih]);
   }

   TCanvas* c1 = new TCanvas("c_resAll", "", 1000, 550);
   Divide(c1, 2, 1, 0.1, 0.1);

   c1->cd(1);
   h1[0]->Draw();
   c1->cd(2);
   h1[1]->Draw();
   text.DrawLatex(0.7, 0.9, ("n = " + ss_har.str()).c_str());

   return 0;
}

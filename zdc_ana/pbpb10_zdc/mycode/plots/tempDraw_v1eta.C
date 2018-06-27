int tempDraw_v1eta(int centb=32) {

   gStyle->SetErrorX(0);
   gStyle->SetOptStat(0);

   TLatex text;
   text.SetNDC(1);
   text.SetTextFont(43);
   text.SetTextSize(18);

   TLegend* l1 = new TLegend(0.65, 0.76, 0.9, 0.9);
   l1->SetFillStyle(0);
   l1->SetBorderSize(0);
   l1->SetTextFont(43);
   l1->SetTextSize(18);

   char name [500];
   char name1[500];

   std::string centLabel[] = {"0-5%","5-10%","10-15%","15-20%","20-25%","25-30%","30-35%","35-40%","40-45%","45-50%","50-55%","55-60%","60-65%","65-70%","70-75%","75-80%","80-85%","85-90%","90-95%","95-100%",   "0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%",  "0-20%","20-40%","40-60%","60-80%","80-100%",  "0-40%","40-80%"};

   char* sideLabel[] = {"Using #Psi_{1}^{C}", "Using #Psi_{1}^{A}"};

   int   sty[] = {20, 25};
   float siz[] = {1.2, 1.2};
   int   col[] = {kBlue, kRed};

   TFile* f1 = new TFile("../Output_zdcPix_flat_wtrk_doAna_reb.root");
   TH1D* h1  [2];
   TF1*  fun1[2];

   for (int iside=0; iside<2; iside++) {
      sprintf(name , "hV1_eta_res_cent%d_side%d_pt10", centb, iside);
      sprintf(name1, "hV1_eta_res_cent%d_side%d_pt10_Clone", centb, iside);
      h1[iside] = (TH1D*)((TH1D*)f1->Get(name))->Clone(name1);
   }

   for (int ih=0; ih<2; ih++) {
      formatHist1(h1[ih]);
      h1[ih]->SetMarkerStyle(sty[ih]);
      h1[ih]->SetMarkerSize (siz[ih]);
      h1[ih]->SetMarkerColor(col[ih]);
      h1[ih]->SetLineColor  (col[ih]);

      h1[ih]->GetXaxis()->SetTitle("#eta");
      h1[ih]->GetYaxis()->SetTitle("v_{1}");
      h1[ih]->GetXaxis()->SetTitleOffset(1.);
      h1[ih]->GetYaxis()->SetTitleOffset(1.);
      h1[ih]->GetYaxis()->SetRangeUser(-0.005, 0.005);
  
      h1[ih]->Fit("pol1", "Q");
      fun1[ih] = (TF1*)(h1[ih]->GetListOfFunctions()->FindObject("pol1"));
      fun1[ih]->SetLineStyle(2);
      fun1[ih]->SetLineWidth(1);
      fun1[ih]->SetLineColor(col[ih]);
   }

   sprintf(name, "c_v1eta_side0_cent%d", centb);
   TCanvas* c1 = new TCanvas(name, "", 800, 500);
   sprintf(name, "c_v1eta_side1_cent%d", centb);
   TCanvas* c2 = new TCanvas(name, "", 800, 500);
   sprintf(name, "c_v1eta_side01_cent%d", centb);
   TCanvas* c3 = new TCanvas(name, "", 800, 500);

   c1->cd(0);
   h1[0]->DrawCopy();
   text.DrawLatex(0.3, 0.84, centLabel[centb].c_str());
   text.DrawLatex(0.3, 0.78, sideLabel[0]);

   c2->cd(0);
   h1[1]->DrawCopy();
   text.DrawLatex(0.3, 0.84, centLabel[centb].c_str());
   text.DrawLatex(0.3, 0.78, sideLabel[1]);
   
   c3->cd(0);
   h1[0]->Draw();
   h1[1]->Draw("same");
   text.DrawLatex(0.3, 0.84, centLabel[centb].c_str());
   l1->AddEntry(h1[0], sideLabel[0], "pl");
   l1->AddEntry(h1[1], sideLabel[1], "pl");
   l1->Draw();

   return 0;
}

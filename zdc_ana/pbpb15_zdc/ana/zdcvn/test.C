#include <iostream>
#include <fstream>
#include <iomanip>
#include <fstream>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLine.h"
#include "THStack.h"
#include "TLatex.h"
#include "string.h"

void setStyle(TH1D *h,int iter){
  h->SetMarkerStyle(20);
  h->GetXaxis()->SetRangeUser(-1.,5.);
}

void test(){
  TFile *fin=new TFile("./out/out_286665_0.root","read");
  TH1D *hNTrk,*hTrk_w,*hTrk_wei;
  TH2F *hVnW;
  hTrk_w=(TH1D*)fin->Get("hTrk_w");
  hTrk_wei=(TH1D*)fin->Get("hTrk_wei");
  hNTrk=(TH1D*)fin->Get("hNTrk");
  hVnW=(TH2F*)fin->Get("hV1_Zdc_fg_cent0_side2");

  gStyle->SetOptStat(0);
  TCanvas *c1;
  c1=new TCanvas("c1","c1",1200,800);
  c1->Divide(2,2);
  c1->cd(1);
  gPad->SetLogy();
  hTrk_w->Draw();
  setStyle(hTrk_w,0);
  c1->cd(2);
  gPad->SetLogy();
  hTrk_wei->Draw();
  setStyle(hTrk_wei,0);
  c1->cd(3);
  hNTrk->Draw();
  c1->cd(4);
  hVnW->Draw("colz");

}

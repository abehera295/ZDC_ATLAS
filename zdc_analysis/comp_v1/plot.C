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

enum {
  NDET = 16,
  NSIDE   = 2,
  NMODULE = 2,
  NGAIN   = 2,
  NCH1    = 4,
  NCH2    = 24,
  NHAR    = 4,
  //NCENT   = 20,
  NCENT   = 11,
  NCENT_1p  = 100,
  //NK      = 12,
  NSAMP   = 10,
  //STORE_DEP = 20,
  STORE_DEP = 5,
  NVZ_EP=10,
  NPT_EP=10,
  NCH_EP=2
};

enum{
  NETA=12,
  NETA2=10,
  NCH=2,
  NPT=10,
  NHAR2=2,//v2,v3
  NHAR3=5,//Q2n for v3 is v6
  NP=4,
  NK=4
};

char name[300];

void setStyle(TH1D *h,int itype1){
  h->GetYaxis()->CenterTitle(0);       h->GetXaxis()->CenterTitle(0);  //h->GetZaxis()->CenterTitle();
  h->GetYaxis()->SetTitleOffset(1.); h->GetXaxis()->SetTitleOffset(1.);
  //h->GetZaxis()->SetTitleOffset(1.5);
  h->GetXaxis()->SetTitleFont(43);    h->GetXaxis()->SetTitleSize(25);
  h->GetXaxis()->SetLabelFont(43);    h->GetXaxis()->SetLabelSize(15);
  h->GetYaxis()->SetTitleFont(43);    h->GetYaxis()->SetTitleSize(25);
  h->GetYaxis()->SetLabelFont(43);    h->GetYaxis()->SetLabelSize(15);
  //h->GetZaxis()->SetTitleFont(43);    h->GetZaxis()->SetTitleSize(17);
  //h->GetZaxis()->SetLabelFont(43);    h->GetZaxis()->SetLabelSize(15);
  h->GetYaxis()->SetNdivisions(505);
  h->GetXaxis()->SetNdivisions(505);
  int c_mat[]={4,2,4,2,6,28,36,7,8,9,46,30};
  int m_mat[]={21,22,25,25,25,21,30,27,25,28,33,20,23,29,34};
  h->SetTitle("");

  h->SetMarkerStyle(m_mat[itype1]);
  h->SetMarkerColor(c_mat[itype1]);
  h->SetLineColor(c_mat[itype1]);
  h->SetMarkerSize(1.5);

  //h->GetYaxis()->SetRangeUser(0.,5.);
}

void setStyle2(TGraphErrors *h,int itype1){
  h->GetYaxis()->CenterTitle(1);       h->GetXaxis()->CenterTitle(1);  //h->GetZaxis()->CenterTitle();
  h->GetYaxis()->SetTitleOffset(1.1); h->GetXaxis()->SetTitleOffset(1.);
  //h->GetZaxis()->SetTitleOffset(1.5);
  h->GetXaxis()->SetTitleFont(43);    h->GetXaxis()->SetTitleSize(25);
  h->GetXaxis()->SetLabelFont(43);    h->GetXaxis()->SetLabelSize(15);
  h->GetYaxis()->SetTitleFont(43);    h->GetYaxis()->SetTitleSize(25);
  h->GetYaxis()->SetLabelFont(43);    h->GetYaxis()->SetLabelSize(15);
  //h->GetZaxis()->SetTitleFont(43);    h->GetZaxis()->SetTitleSize(17);
  //h->GetZaxis()->SetLabelFont(43);    h->GetZaxis()->SetLabelSize(15);
  h->GetYaxis()->SetNdivisions(505);
  h->GetXaxis()->SetNdivisions(505);
  int c_mat[]={1,4,2,4,2,6,28,36,7,8,9,46,30};
  int m_mat[]={20,21,22,25,25,21,30,27,25,28,33,20,23,29,34};
  h->SetTitle("");

  h->SetMarkerStyle(m_mat[itype1]);
  h->SetMarkerColor(c_mat[itype1]);
  h->SetLineColor(c_mat[itype1]);
  h->SetMarkerSize(1.5);

  //h->GetYaxis()->SetRangeUser(0.,5.);
}

void plot(){
  TFile *fin[2];
  fin[0]=new TFile("./v1plots_2010.root","read");
  fin[1]=new TFile("./v1plots_2015.root","read");
  TH1D *hResCent_mix_Zdc_sub[2],*hResCent_mix_Zdc_ful[2];
  TGraphErrors *gVnEven_pt_res_allCent[2],*gVnOdd_pt_res_allCent[2];
  TGraphErrors *gVnEven_eta_res_cent[2],*gVnOdd_eta_res_cent[2];
  TGraphErrors *gVnEven_pt_res_allCent_pub[2],*gVnOdd_pt_res_allCent_pub[2];
  
  int mhar=0;
  for(int id=0;id<2;id++){
    sprintf(name,"hResCent_mix_Zdc_sub_har%d", mhar);
    hResCent_mix_Zdc_sub[id]=(TH1D*)fin[id]->Get(name);
    sprintf(name,"hResCent_mix_Zdc_ful_har%d", mhar);
    hResCent_mix_Zdc_ful[id]=(TH1D*)fin[id]->Get(name);

    sprintf(name,"gVnEven_pt_res_allCent_har%d_eta14",mhar);
    gVnEven_pt_res_allCent[id]=(TGraphErrors*)fin[id]->Get(name);
    sprintf(name,"gVnOdd_pt_res_allCent_har%d_eta15",mhar);
    gVnOdd_pt_res_allCent[id]=(TGraphErrors*)fin[id]->Get(name);

    sprintf(name,"gVnEven_eta_res_cent2_har%d",mhar);
    gVnEven_eta_res_cent[id]=(TGraphErrors*)fin[id]->Get(name);
    sprintf(name,"gVnOdd_eta_res_cent2_har%d",mhar);
    gVnOdd_eta_res_cent[id]=(TGraphErrors*)fin[id]->Get(name);

    sprintf(name,"gVnEven_pt_res_allCent_pub");
    gVnEven_pt_res_allCent_pub[id]=(TGraphErrors*)fin[id]->Get(name);
    sprintf(name,"gVnOdd_pt_res_allCent_pub");
    gVnOdd_pt_res_allCent_pub[id]=(TGraphErrors*)fin[id]->Get(name);
  }
    
  gStyle->SetOptStat(0);
  TLine *line0=new TLine(0.,0.,5.,0.);
  line0->SetLineStyle(7);
  line0->SetLineWidth(1);
  TLatex tt;
  tt.SetNDC(1);
  tt.SetTextFont(43);
  tt.SetTextSize(25);
  
  TCanvas *c1,*c2;
  TLegend *leg1,*leg2,*leg3;
  c1=new TCanvas("c1","c1",800,600);
  leg1=new TLegend(0.6,0.7,0.88,0.88);
  //leg1->SetNColumns(4);
  leg1->SetLineColor(0);

  c2=new TCanvas("c2","c2",1400,600);
  c2->Divide(2);
  leg2=new TLegend(0.2,0.65,0.6,0.88);
  //leg2->SetNColumns(4);
  leg2->SetLineColor(0);
  leg3=new TLegend(0.2,0.65,0.6,0.88);
  //leg3->SetNColumns(4);
  leg3->SetLineColor(0);
  
  
  c1->cd();
  leg1->AddEntry(hResCent_mix_Zdc_sub[0],"2010 Pb+Pb","P");
  leg1->AddEntry(hResCent_mix_Zdc_sub[1],"2015 Pb+Pb","P");
  setStyle(hResCent_mix_Zdc_sub[0],0);
  setStyle(hResCent_mix_Zdc_sub[1],1);
  hResCent_mix_Zdc_sub[0]->Draw();
  hResCent_mix_Zdc_sub[0]->GetYaxis()->SetRangeUser(0.,0.4);
  hResCent_mix_Zdc_sub[0]->GetXaxis()->SetRangeUser(0.,80);
  sprintf(name,"Res(#Psi_{%d}^{ZDC})",1);
  hResCent_mix_Zdc_sub[0]->SetYTitle(name);
  hResCent_mix_Zdc_sub[0]->SetXTitle("Centrality %");
  hResCent_mix_Zdc_sub[1]->Draw("same");
  leg1->Draw("same");

  c2->cd(1);
  leg2->AddEntry(gVnEven_pt_res_allCent_pub[0],"2010 Pb+Pb ALICE","P");
  leg2->AddEntry(gVnEven_pt_res_allCent[0],"2010 Pb+Pb ATLAS","P");
  leg2->AddEntry(gVnEven_pt_res_allCent[1],"2015 Pb+Pb ATLAS","P");
  setStyle2(gVnEven_pt_res_allCent_pub[0],0);
  setStyle2(gVnEven_pt_res_allCent[0],1);
  setStyle2(gVnEven_pt_res_allCent[1],2);
  gVnEven_pt_res_allCent[0]->Draw("AP");
  gVnEven_pt_res_allCent[0]->SetTitle("");
  gVnEven_pt_res_allCent[0]->GetYaxis()->SetRangeUser(-0.001,0.004);
  gVnEven_pt_res_allCent[0]->GetXaxis()->SetRangeUser(0.,5.);
  gVnEven_pt_res_allCent[0]->GetXaxis()->SetTitle("p_{T}");
  gVnEven_pt_res_allCent[0]->GetYaxis()->SetTitle("Even v_{1}");
  gVnEven_pt_res_allCent_pub[0]->Draw("same P");
  gVnEven_pt_res_allCent[1]->Draw("same P");
  leg2->Draw("same");
  line0->Draw("same");
  tt.DrawLatex(0.15,0.55,"|#eta|<0.8");
  tt.DrawLatex(0.15,0.45,"Cent : 5-80%");
  
  c2->cd(2);
  leg3->AddEntry(gVnOdd_pt_res_allCent_pub[0],"2010 Pb+Pb ALICE","P");
  leg3->AddEntry(gVnOdd_pt_res_allCent[0],"2010 Pb+Pb ATLAS","P");
  leg3->AddEntry(gVnOdd_pt_res_allCent[1],"2015 Pb+Pb ATLAS","P");
  setStyle2(gVnOdd_pt_res_allCent_pub[0],0);
  setStyle2(gVnOdd_pt_res_allCent[0],1);
  setStyle2(gVnOdd_pt_res_allCent[1],2);
  gVnOdd_pt_res_allCent[0]->Draw("AP");
  gVnOdd_pt_res_allCent[0]->SetTitle("");
  gVnOdd_pt_res_allCent[0]->GetYaxis()->SetRangeUser(-0.001,0.004);
  gVnOdd_pt_res_allCent[0]->GetXaxis()->SetRangeUser(0.,5.);
  gVnOdd_pt_res_allCent[0]->GetXaxis()->SetTitle("p_{T}");
  gVnOdd_pt_res_allCent[0]->GetYaxis()->SetTitle("Odd v_{1}");
  gVnOdd_pt_res_allCent_pub[0]->Draw("same P");
  gVnOdd_pt_res_allCent[1]->Draw("same P");
  leg3->Draw("same");
  line0->Draw("same");
  tt.DrawLatex(0.15,0.55,"|#eta|<0.8");
  tt.DrawLatex(0.15,0.45,"Cent : 5-80%");

    
  c1->SaveAs("./Res_v1_comp_run1run2.pdf");
  c2->SaveAs("./v1_comp_run1run2.pdf");
}

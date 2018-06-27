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
  h->GetXaxis()->SetTitleFont(43);    h->GetXaxis()->SetTitleSize(17);
  h->GetXaxis()->SetLabelFont(43);    h->GetXaxis()->SetLabelSize(15);
  h->GetYaxis()->SetTitleFont(43);    h->GetYaxis()->SetTitleSize(17);
  h->GetYaxis()->SetLabelFont(43);    h->GetYaxis()->SetLabelSize(15);
  //h->GetZaxis()->SetTitleFont(43);    h->GetZaxis()->SetTitleSize(17);
  //h->GetZaxis()->SetLabelFont(43);    h->GetZaxis()->SetLabelSize(15);
  h->GetYaxis()->SetNdivisions(505);
  h->GetXaxis()->SetNdivisions(505);
  int c_mat[]={1,4,2,2,4,6,28,36,7,8,9,46,30};
  int m_mat[]={20,21,24,25,25,21,30,27,25,28,33,20,23,29,34};
  h->SetTitle("");

  h->SetMarkerStyle(m_mat[itype1]);
  h->SetMarkerColor(c_mat[itype1]);
  h->SetLineColor(c_mat[itype1]);
  h->SetMarkerSize(1.5);

  //h->GetYaxis()->SetRangeUser(0.,5.);
}

void setStyle2(TH1D *h,int itype1){
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
  int c_mat[]={1,4,2,2,6,28,36,7,8,9,46,30};
  int m_mat[]={20,20,24,25,21,30,27,25,28,33,20,23,29,34};
  h->SetTitle("");

  h->SetMarkerStyle(m_mat[itype1]);
  h->SetMarkerColor(c_mat[itype1]);
  h->SetLineColor(c_mat[itype1]);
  h->SetMarkerSize(1.5);

  //h->GetYaxis()->SetRangeUser(0.,5.);
}

void plot_atlas_v2(){
  //TFile *fin=new TFile("../zdcvn8/Output_doAna_vn.root","read");
  //TFile *fin=new TFile("../zdcvn8/Output_doAna_vn_old.root","read");
  //TFile *fin=new TFile("../zdcvn10/Output_doAna_vn_mix.root","read");
  //TFile *fin=new TFile("../zdcvn11/Output_doAna_vn_mix2.root","read");
  TFile *fin=new TFile("../zdcvn/Output_doAna_vn_nomix.root","read");
  TH1D *hVn_cum2_cent[NHAR2],*hVn_cum4_cent[NHAR2];
  TH1D *hVn_sub_cum2_cent[NHAR2],*hVn_sub_cum4_cent[NHAR2];
  TH1D *hVn_FCal_cent_res[NHAR],*hVn_Zdc_cent_res[NHAR];
  TH1D *hVn_FCal_cent_res1,*hVn_Zdc_cent_res1;
  TH1D *hVn_FCal_cent_res2,*hVn_Zdc_cent_res2;
  TH1D *hVn_FCal_pt_res,*hVn_Zdc_pt_res;
  
  for (int ihar=0; ihar<NHAR2; ihar++) {
    sprintf(name,"hVn_cum2_cent_har%d",ihar);
    hVn_cum2_cent[ihar]=(TH1D*)fin->Get(name);
    sprintf(name,"hVn_cum4_cent_har%d",ihar);
    hVn_cum4_cent[ihar]=(TH1D*)fin->Get(name);

    sprintf(name,"hVn_sub_cum2_cent_har%d",ihar);
    hVn_sub_cum2_cent[ihar]=(TH1D*)fin->Get(name);
    sprintf(name,"hVn_sub_cum4_cent_har%d",ihar);
    hVn_sub_cum4_cent[ihar]=(TH1D*)fin->Get(name);
  }
  
  for(int ihar=0; ihar<NHAR; ihar++) {
    sprintf(name,"hV%d_FCal_cent_res_side%d_pt%d_eta%d",ihar+1,NSIDE,NPT,NETA);
    hVn_FCal_cent_res[ihar]=(TH1D*)fin->Get(name);
    sprintf(name,"hV%d_Zdc_cent_raw_side%d_pt%d_eta%d",ihar+1,NSIDE,NPT,NETA);
    hVn_Zdc_cent_res[ihar]=(TH1D*)fin->Get(name);
  }
  
  int iside=2,ihar=0,ipt=NPT,ieta1=NETA,ieta2=NETA-5;
  int icent=2;
  sprintf(name,"hV%d_FCal_cent_res_side%d_pt%d_eta%d",ihar+1,iside,ipt,ieta1);
  hVn_FCal_cent_res1=(TH1D*)fin->Get(name);
  sprintf(name,"hV%d_Zdc_cent_res_side%d_pt%d_eta%d",ihar+2,iside,ipt,ieta1);
  hVn_Zdc_cent_res1=(TH1D*)fin->Get(name);

  sprintf(name,"hV%d_FCal_cent_res_side%d_pt%d_eta%d",ihar+1,iside,ipt,ieta2);
  hVn_FCal_cent_res2=(TH1D*)fin->Get(name);
  sprintf(name,"hV%d_Zdc_cent_res_side%d_pt%d_eta%d",ihar+2,iside,ipt,ieta2);
  hVn_Zdc_cent_res2=(TH1D*)fin->Get(name);

  sprintf(name,"hV%d_FCal_pt_res_cent%d_side%d_eta%d",ihar+1,icent,iside,NETA);
  hVn_FCal_pt_res=(TH1D*)fin->Get(name);
  sprintf(name,"hV%d_Zdc_pt_res_cent%d_side%d_eta%d",ihar+2,icent,iside,NETA);
  hVn_Zdc_pt_res=(TH1D*)fin->Get(name);
  
  
  
  gStyle->SetOptStat(0);
  
  TCanvas *c;
  c=new TCanvas("c","c",900,600);
  TLegend *leg;
  leg=new TLegend(0.3,0.25,0.55,0.55);
  //leg->SetNColumns(3);
  leg->SetLineColor(0);

  setStyle(hVn_FCal_cent_res1,0);
  setStyle(hVn_Zdc_cent_res1,1);
  setStyle(hVn_sub_cum2_cent[ihar],2);
  setStyle(hVn_sub_cum4_cent[ihar],3);
  //setStyle(hVn_cum2_cent[ihar],1);
  //setStyle(hVn_cum4_cent[ihar],3);
  //setStyle(hVn_FCal_cent_res2,4);
  //setStyle(hVn_Zdc_cent_res2,5);

  leg->AddEntry(hVn_FCal_cent_res1,"v_{2} EP FCal","P");
  leg->AddEntry(hVn_Zdc_cent_res1,"v_{2} EP Zdc","P");
  leg->AddEntry(hVn_sub_cum2_cent[ihar],"v_{2}#left{2#right} 2SubEv","P");
  leg->AddEntry(hVn_sub_cum4_cent[ihar],"v_{2}#left{4#right} 2SubEv","P");
  //leg->AddEntry(hVn_cum2_cent[ihar],"v_{2}#left{2#right} Standard","P");
  //leg->AddEntry(hVn_cum4_cent[ihar],"v_{2}#left{4#right} Standard","P");


  c->cd();
  hVn_FCal_cent_res1->Draw();
  hVn_FCal_cent_res1->GetYaxis()->SetRangeUser(0.,0.16);
  hVn_FCal_cent_res1->SetYTitle("v_{2}");
  hVn_FCal_cent_res1->SetXTitle("Centrality %");
  //hVn_cum2_cent[ihar]->Draw("same");
  hVn_sub_cum2_cent[ihar]->Draw("same");
  //hVn_cum4_cent[ihar]->Draw("same");
  hVn_sub_cum4_cent[ihar]->Draw("same");
  hVn_Zdc_cent_res1->Draw("same");
  //hVn_FCal_cent_res2->Draw("same");
  //hVn_Zdc_cent_res2->Draw("same");
  leg->Draw("same");

  c->SaveAs("./v2_comp.pdf");

  /*
  TCanvas *c3;
  c3=new TCanvas("c3","c3",800,600);
  TLegend *leg3;
  leg3=new TLegend(0.15,0.65,0.45,0.88);
  leg3->SetLineColor(0);

  for(int ihar=0; ihar<NHAR-1; ihar++) {
    sprintf(name,"v_{%d}",ihar+2);
    leg3->AddEntry(hVn_FCal_cent_res[ihar],name,"P");
    setStyle2(hVn_FCal_cent_res[ihar],ihar);
  }

  c3->cd();
  hVn_FCal_cent_res[0]->Draw();
  hVn_FCal_cent_res[0]->GetXaxis()->SetTitle("Centrality");
  hVn_FCal_cent_res[0]->GetYaxis()->SetTitle("v_{n}");
  //hVn_FCal_cent_res[0]->GetXaxis()->SetRangeUser(0,80);
  hVn_FCal_cent_res[0]->GetYaxis()->SetRangeUser(-0.01,0.2);
  for(int ihar=0; ihar<NHAR-1; ihar++) {
    hVn_FCal_cent_res[ihar]->Draw("same");
  }

  leg3->Draw("same");
  */
}

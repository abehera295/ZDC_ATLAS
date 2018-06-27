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
  int c_mat[]={1,2,4,6,28,36,7,8,9,46,30};
  int m_mat[]={24,24,24,25,25,21,30,27,25,28,33,20,23,29,34};
  h->SetTitle("");

  h->SetMarkerStyle(m_mat[itype1]);
  h->SetMarkerColor(c_mat[itype1]);
  
  h->SetMarkerSize(1.2);

  //h->GetYaxis()->SetRangeUser(0.,5.);
}

void test_mixedEv(){
  TFile *fin=new TFile("./zdcvn10/Output_doAna_vn_mix.root","read");
  TH1D *hDphi_FCal_fg,*hDphi_FCal_bg,*hDphi_FCal_rat;
  TH1D *hDphi_Zdc_fg,*hDphi_Zdc_bg,*hDphi_Zdc_rat;
  int icent=3,ihar=0;
  sprintf(name,"hDphi_Zdc_fg_cent%d_har%d",icent,ihar);
  hDphi_Zdc_fg=(TH1D*)fin->Get(name);
  sprintf(name,"hDphi_Zdc_bg_cent%d_har%d",icent,ihar);
  hDphi_Zdc_bg=(TH1D*)fin->Get(name);
  sprintf(name,"hDphi_Zdc_rat_cent%d_har%d",icent,ihar);
  hDphi_Zdc_rat=(TH1D*)fin->Get(name);

  cout<<"NBins = "<<hDphi_Zdc_rat->GetNbinsX()<<endl;
  int irebin=10;
  hDphi_Zdc_fg->Rebin(irebin);
  hDphi_Zdc_bg->Rebin(irebin);
  hDphi_Zdc_rat->Rebin(irebin);
  cout<<"NBins after rebin = "<<hDphi_Zdc_rat->GetNbinsX()<<endl;

  //SetStyle
  setStyle(hDphi_Zdc_fg,0);
  setStyle(hDphi_Zdc_bg,1);
  setStyle(hDphi_Zdc_rat,2);
  gStyle->SetOptStat(0);
  
  TCanvas *c;
  c=new TCanvas("c","c",1200,600);
  c->Divide(2);
  TLegend *leg;
  leg=new TLegend(0.4,0.75,0.7,0.88);
  //leg->SetNColumns(3);
  leg->SetLineColor(0);
  leg->AddEntry(hDphi_Zdc_fg,"Same Ev","P");
  leg->AddEntry(hDphi_Zdc_rat,"Same Ev - Mixed Ev","P");
  
  c->cd(1);
  hDphi_Zdc_fg->Draw();
  hDphi_Zdc_fg->GetYaxis()->SetRangeUser(30e3,70e3);
  hDphi_Zdc_fg->SetXTitle("#Psi_{1}^{A}-#Psi_{1}^{C}");
  hDphi_Zdc_rat->Draw("same");
  leg->Draw("same");

  c->cd(2);
  hDphi_Zdc_bg->Draw();
  hDphi_Zdc_bg->GetYaxis()->SetRangeUser(500e3,700e3);
  hDphi_Zdc_bg->SetXTitle("#Psi_{1}^{A}-#Psi_{1}^{C}");
}

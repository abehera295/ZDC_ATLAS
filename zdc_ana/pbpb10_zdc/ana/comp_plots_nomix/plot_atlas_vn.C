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
  int m_mat[]={20,20,25,25,25,21,30,27,25,28,33,20,23,29,34};
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

void plot_atlas_vn(){
  TFile *fin=new TFile("../zdcvn11/Output_doAna_vn_mix3.root","read");
  TH1D *hVn_FCal_cent_res[NHAR],*hVn_Zdc_cent_res[NHAR];
  
  for(int ihar=0; ihar<NHAR; ihar++) {
    sprintf(name,"hV%d_FCal_cent_res_side%d_pt%d_eta%d",ihar+1,2,NPT,NETA);
    hVn_FCal_cent_res[ihar]=(TH1D*)fin->Get(name);
    sprintf(name,"hV%d_Zdc_cent_res_side%d_pt%d_eta%d",ihar+1,2,NPT,NETA);
    hVn_Zdc_cent_res[ihar]=(TH1D*)fin->Get(name);
  }
  
    
  gStyle->SetOptStat(0);
  TLine *line0=new TLine(0.,0.,80.,0.);
  line0->SetLineStyle(7);
  line0->SetLineWidth(1);
  
  TCanvas *c;
  c=new TCanvas("c","c",900,600);
  TLegend *leg;
  leg=new TLegend(0.15,0.85,0.88,0.8);
  leg->SetNColumns(4);
  leg->SetLineColor(0);

  setStyle(hVn_FCal_cent_res[0],0);
  setStyle(hVn_Zdc_cent_res[1],1);
  setStyle(hVn_FCal_cent_res[1],2);
  setStyle(hVn_Zdc_cent_res[2],3);

  leg->AddEntry(hVn_FCal_cent_res[0],"v_{2}(#Psi_{2}^{FCal})","P");
  leg->AddEntry(hVn_Zdc_cent_res[1],"v_{2}(#Psi_{1}^{Zdc})","P");
  leg->AddEntry(hVn_FCal_cent_res[1],"v_{3}(#Psi_{3}^{FCal})","P");
  leg->AddEntry(hVn_Zdc_cent_res[2],"v_{3}(#Psi_{1}^{Zdc})","P");


  c->cd();
  hVn_FCal_cent_res[0]->Draw();
  hVn_FCal_cent_res[0]->GetYaxis()->SetRangeUser(-0.01,0.2);
  hVn_FCal_cent_res[0]->GetXaxis()->SetRangeUser(0.,80);
  hVn_FCal_cent_res[0]->SetYTitle("v_{n}");
  hVn_FCal_cent_res[0]->SetXTitle("Centrality %");
  hVn_Zdc_cent_res[1]->Draw("same");
  hVn_FCal_cent_res[1]->Draw("same");
  hVn_Zdc_cent_res[2]->Draw("same");
  leg->Draw("same");
  line0->Draw("same");
  c->SaveAs("./vn_comp.pdf");

}

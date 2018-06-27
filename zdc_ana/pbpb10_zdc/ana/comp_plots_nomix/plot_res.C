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

void plot_res(){
  TFile *fin=new TFile("../zdcvn11/Output_doAna_vn_mix3.root","read");
  TH1D *hResCent_mix_FCal_sub[NHAR],*hResCent_mix_FCal_ful[NHAR];
  TH1D *hResCent_mix_Zdc_sub[NHAR],*hResCent_mix_Zdc_ful[NHAR];
  
  for(int ihar=0; ihar<NHAR; ihar++) {
    sprintf(name,"hResCent_mix_FCal_sub_har%d", ihar);
    hResCent_mix_FCal_sub[ihar]=(TH1D*)fin->Get(name);
    sprintf(name,"hResCent_mix_FCal_ful_har%d", ihar);
    hResCent_mix_FCal_ful[ihar]=(TH1D*)fin->Get(name);

    sprintf(name,"hResCent_mix_Zdc_sub_har%d", ihar);
    hResCent_mix_Zdc_sub[ihar]=(TH1D*)fin->Get(name);
    sprintf(name,"hResCent_mix_Zdc_ful_har%d", ihar);
    hResCent_mix_Zdc_ful[ihar]=(TH1D*)fin->Get(name);
  }
  
    
  gStyle->SetOptStat(0);
  TLine *line0=new TLine(0.,0.,80.,0.);
  line0->SetLineStyle(7);
  line0->SetLineWidth(1);
  
  TCanvas *c[NHAR][2];
  TLegend *leg[NHAR][2];
  for(int ihar=0;ihar<NHAR;ihar++){
    for(int i=0;i<2;i++){
      sprintf(name,"c_har%d_%d",ihar,i);
      c[ihar][i]=new TCanvas(name,"",900,600);
      leg[ihar][i]=new TLegend(0.7,0.7,0.88,0.88);
      //leg[ihar][i]->SetNColumns(4);
      leg[ihar][i]->SetLineColor(0);
      if(i==0){//FCal
	leg[ihar][i]->AddEntry(hResCent_mix_FCal_sub[ihar],"Half Det","P");
	leg[ihar][i]->AddEntry(hResCent_mix_FCal_ful[ihar],"Full Det","P");
      }
      if(i==1){//Zdc
	leg[ihar][i]->AddEntry(hResCent_mix_Zdc_sub[ihar],"Half Det","P");
	leg[ihar][i]->AddEntry(hResCent_mix_Zdc_ful[ihar],"Full Det","P");
      }
    }
  }


  for(int ihar=0;ihar<NHAR;ihar++){
    setStyle(hResCent_mix_FCal_sub[ihar],0);
    setStyle(hResCent_mix_FCal_ful[ihar],1);
    setStyle(hResCent_mix_Zdc_sub[ihar],0);
    setStyle(hResCent_mix_Zdc_ful[ihar],1);

    c[ihar][0]->cd();
    hResCent_mix_FCal_ful[ihar]->Draw();
    hResCent_mix_FCal_ful[ihar]->GetYaxis()->SetRangeUser(0.,1.);
    hResCent_mix_FCal_ful[ihar]->GetXaxis()->SetRangeUser(0.,80);
    sprintf(name,"Res(%d#Psi_{%d}^{FCAL})",ihar+2,ihar+2);
    hResCent_mix_FCal_ful[ihar]->SetYTitle(name);
    hResCent_mix_FCal_ful[ihar]->SetXTitle("Centrality %");
    hResCent_mix_FCal_sub[ihar]->Draw("same");
    leg[ihar][0]->Draw("same");

    c[ihar][1]->cd();
    hResCent_mix_Zdc_ful[ihar]->Draw();
    hResCent_mix_Zdc_ful[ihar]->GetXaxis()->SetRangeUser(0.,80);
    hResCent_mix_Zdc_ful[ihar]->GetYaxis()->SetRangeUser(0.,0.5);
    if(ihar>0) hResCent_mix_Zdc_ful[ihar]->GetYaxis()->SetRangeUser(0.,0.1);
    sprintf(name,"Res(#Psi_{%d}^{ZDC})",1);
    if(ihar>0) sprintf(name,"Res(%d#Psi_{%d}^{ZDC})",ihar+1,1);
    hResCent_mix_Zdc_ful[ihar]->SetYTitle(name);
    hResCent_mix_Zdc_ful[ihar]->SetXTitle("Centrality %");
    hResCent_mix_Zdc_sub[ihar]->Draw("same");
    leg[ihar][1]->Draw("same");
    sprintf(name,"Res_FCal_v%d.pdf",ihar+2);
    c[ihar][0]->SaveAs(name);
    sprintf(name,"Res_Zdc_v%d.pdf",ihar+1);
    c[ihar][1]->SaveAs(name);
  }
}

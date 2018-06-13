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

//Variables and Hists
char name[300];
static double PI = acos(-1.0);
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
  STORE_DEP = 3,
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

int cent_mat[]={0,5,10,20,30,40,50,60,70,80,90,100};
double cent_mat_double[]={0,5,10,20,30,40,50,60,70,80,90,100};
double pt_mat[]={0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 7.0, 9.0, 12.0};


TProfile *hMean_pt[NCENT][3][NSIDE+1];
TH1D *hVn_Zdc_eta_raw[NCENT][NSIDE+1][NHAR][NPT+1],*hVn_Zdc_eta_raw_fg[NCENT][NSIDE+1][NHAR][NPT+1],*hVn_Zdc_eta_raw_bg[NCENT][NSIDE+1][NHAR][NPT+1];
  TH1D *hVn_Zdc_eta_res[NCENT][NSIDE+1][NHAR][NPT+1],*hVn_Zdc_eta_res_fg[NCENT][NSIDE+1][NHAR][NPT+1],*hVn_Zdc_eta_res_bg[NCENT][NSIDE+1][NHAR
][NPT+1];
  TH1D *hVn_Zdc_pt_raw[NCENT][NSIDE+1][NHAR][NETA+4],*hVn_Zdc_pt_raw_fg[NCENT][NSIDE+1][NHAR][NETA+4],*hVn_Zdc_pt_raw_bg[NCENT][NSIDE+1][NHAR]
[NETA+4];
  TH1D *hVn_Zdc_pt_res[NCENT][NSIDE+1][NHAR][NETA+4],*hVn_Zdc_pt_res_fg[NCENT][NSIDE+1][NHAR][NETA+4],*hVn_Zdc_pt_res_bg[NCENT][NSIDE+1][NHAR]
[NETA+4];
  TH1D *hVn_Zdc_cent_raw[NSIDE+1][NHAR][NPT+1][NETA+4],*hVn_Zdc_cent_raw_fg[NSIDE+1][NHAR][NPT+1][NETA+4],*hVn_Zdc_cent_raw_bg[NSIDE+1][NHAR][
NPT+1][NETA+4];
  TH1D *hVn_Zdc_cent_res[NSIDE+1][NHAR][NPT+1][NETA+4],*hVn_Zdc_cent_res_fg[NSIDE+1][NHAR][NPT+1][NETA+4],*hVn_Zdc_cent_res_bg[NSIDE+1][NHAR][
NPT+1][NETA+4];
TH1D *hVn_Zdc_pt_raw_allCent[NSIDE+1][NHAR][NETA+4],*hVn_Zdc_pt_raw_allCent_fg[NSIDE+1][NHAR][NETA+4],*hVn_Zdc_pt_raw_allCent_bg[NSIDE+1][NHAR][NETA+4];
TH1D *hVn_Zdc_pt_res_allCent[NSIDE+1][NHAR][NETA+4],*hVn_Zdc_pt_res_allCent_fg[NSIDE+1][NHAR][NETA+4],*hVn_Zdc_pt_res_allCent_bg[NSIDE+1][NHAR][NETA+4];

TH1D *hVnEven_cent_raw[NHAR][NPT+1][NETA+4],*hVnOdd_cent_raw[NHAR][NPT+1][NETA+4];
TH1D *hVnEven_cent_res[NHAR][NPT+1][NETA+4],*hVnOdd_cent_res[NHAR][NPT+1][NETA+4];
TH1D *hVnEven_eta_raw[NCENT][NHAR][NPT+1],*hVnOdd_eta_raw[NCENT][NHAR][NPT+1];
TH1D *hVnEven_eta_res[NCENT][NHAR][NPT+1],*hVnOdd_eta_res[NCENT][NHAR][NPT+1];


TH1D *hVnEven_pt_raw_allCent[NHAR][NETA+4],*hVnOdd_pt_raw_allCent[NHAR][NETA+4];
TH1D *hVnEven_pt_res_allCent[NHAR][NETA+4],*hVnOdd_pt_res_allCent[NHAR][NETA+4];

TGraphErrors *gVnEven_pt_raw_allCent[NHAR][NETA+4];
TGraphErrors *gVnEven_pt_res_allCent[NHAR][NETA+4];
TGraphErrors *gVnOdd_pt_raw_allCent[NHAR][NETA+4];
TGraphErrors *gVnOdd_pt_res_allCent[NHAR][NETA+4];

void setStyle(TH1D *h,int itype1){
  h->GetYaxis()->CenterTitle(0);       h->GetXaxis()->CenterTitle(0);  //h->GetZaxis()->Cen\
terTitle();
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
  int c_mat[]={1,2,4,6,28,36,7,8,9,46,30};
  int m_mat[]={24,30,27,25,28,33,20,23,29,34};
  h->SetTitle("");

  h->SetMarkerStyle(24);
  h->SetMarkerSize(1.5);
  h->SetMarkerColor(c_mat[itype1]);
  h->SetLineColor(c_mat[itype1]);

  //h->GetYaxis()->SetRangeUser(0.,5.);
}

void setStyle2(TGraphErrors *h,int itype1){
  h->GetYaxis()->CenterTitle(0);       h->GetXaxis()->CenterTitle(0);  //h->GetZaxis()->Cen \
terTitle();
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
  int c_mat[]={1,2,1,2,4,6,28,36,7,8,9,46,30};
  int m_mat[]={20,24,21,21,30,27,25,28,33,20,23,29,34};
  h->SetTitle("");

  h->SetMarkerStyle(m_mat[itype1]);
  h->SetMarkerSize(1.5);
  h->SetMarkerColor(c_mat[itype1]);
  h->SetLineColor(c_mat[itype1]);

  //h->GetYaxis()->SetRangeUser(0.,5.);
}


void plot_alice2(){
  TFile *fin=new TFile("./zdcvn8/Output_doAna_vn2.root","read");
  TProfile *hRes_Zdc_cos_fg[NCENT],*hRes_Zdc_sin_fg[NCENT];
  for(int icent=0;icent<NCENT;icent++){
    sprintf(name, "hRes_Zdc_cos_fg_cent%d", icent);
    hRes_Zdc_cos_fg[icent]=(TProfile*)fin->Get(name);
    sprintf(name, "hRes_Zdc_sin_fg_cent%d", icent);
    hRes_Zdc_sin_fg[icent]=(TProfile*)fin->Get(name);
  }

  TH1D *hRes_Zdc_cos_psi1_atlas,*hRes_Zdc_sin_psi1_atlas,*hRes_Zdc_psi1_atlas;
  sprintf(name,"hRes_Zdc_cos_psi1_atlas");
  hRes_Zdc_cos_psi1_atlas=new TH1D(name,"",NCENT,cent_mat_double);
  sprintf(name,"hRes_Zdc_sin_psi1_atlas");
  hRes_Zdc_sin_psi1_atlas=new TH1D(name,"",NCENT,cent_mat_double);
  sprintf(name,"hRes_Zdc_psi1_atlas");
  hRes_Zdc_psi1_atlas=new TH1D(name,"",NCENT,cent_mat_double);

  double y,y1,y2,dy,dy1,dy2;
  for(int icent=0;icent<NCENT;icent++){
    y1=hRes_Zdc_cos_fg[icent]->GetBinContent(1);
    dy1=hRes_Zdc_cos_fg[icent]->GetBinError(1);
    hRes_Zdc_cos_psi1_atlas->SetBinContent(icent+1,y1);
    hRes_Zdc_cos_psi1_atlas->SetBinError(icent+1,dy1);

    y2=hRes_Zdc_sin_fg[icent]->GetBinContent(1);
    dy2=hRes_Zdc_sin_fg[icent]->GetBinError(1);
    hRes_Zdc_sin_psi1_atlas->SetBinContent(icent+1,y2);
    hRes_Zdc_sin_psi1_atlas->SetBinError(icent+1,dy2);

    y=sqrt(fabs((y1+y2)/2.));
    dy=0.5*pow(fabs((y1+y2)/2.),-0.5)*0.5*sqrt(pow(dy1,2)+pow(dy2,2));
    hRes_Zdc_psi1_atlas->SetBinContent(icent+1,y);
    hRes_Zdc_psi1_atlas->SetBinError(icent+1,dy);
  }
  
  TH1D *hRes_Zdc_cos_psi1_alice,*hRes_Zdc_sin_psi1_alice,*hRes_Zdc_psi1_alice;
  sprintf(name,"hRes_Zdc_cos_psi1_alice");
  hRes_Zdc_cos_psi1_alice=new TH1D(name,"",NCENT,cent_mat_double);
  sprintf(name,"hRes_Zdc_sin_psi1_alice");
  hRes_Zdc_sin_psi1_alice=new TH1D(name,"",NCENT,cent_mat_double);
  sprintf(name,"hRes_Zdc_psi1_alice");
  hRes_Zdc_psi1_alice=new TH1D(name,"",NCENT,cent_mat_double);
  
  //Get data for Alice published results
  ifstream inf;
  inf.open("./pubdata/data_res_alice.txt");
  for(int icent=1;icent<9;icent++){
    double Z[5];
    inf>>Z[0]>>Z[1]>>Z[2]>>Z[3]>>Z[4];
    y1=100*Z[1],dy1=100*Z[2];
    y2=100*Z[3],dy2=100*Z[4];

    hRes_Zdc_cos_psi1_alice->SetBinContent(icent+1,y1);
    hRes_Zdc_cos_psi1_alice->SetBinError(icent+1,dy1);

    hRes_Zdc_sin_psi1_alice->SetBinContent(icent+1,y2);
    hRes_Zdc_sin_psi1_alice->SetBinError(icent+1,dy2);

    y=sqrt(fabs((y1+y2)/2.));
    dy=0.5*pow(fabs((y1+y2)/2.),-0.5)*0.5*sqrt(pow(dy1,2)+pow(dy2,2));
    hRes_Zdc_psi1_alice->SetBinContent(icent+1,y);
    hRes_Zdc_psi1_alice->SetBinError(icent+1,dy);
  }
  inf.close();


  //Plot
  gStyle->SetOptStat(0);
  TLine *line0=new TLine(0.,0.,100.,0.);
  line0->SetLineStyle(7);
  line0->SetLineWidth(1);

  TCanvas *c;
  c=new TCanvas("c","c",1200,600);
  c->Divide(2);
  
  TLegend *leg1;
  leg1=new TLegend(0.25,0.7,0.5,0.85);
  //leg1->SetNColumns(2);
  leg1->SetLineColor(0);
  leg1->AddEntry(hRes_Zdc_cos_psi1_alice,"<Q_{x}^{t}Q_{x}^{p}> ALICE","P");
  leg1->AddEntry(hRes_Zdc_sin_psi1_alice,"<Q_{y}^{t}Q_{y}^{p}>ALICE","P");
  setStyle(hRes_Zdc_cos_psi1_alice,0);
  setStyle(hRes_Zdc_sin_psi1_alice,1);

  c->cd(1);
  hRes_Zdc_cos_psi1_alice->Draw();
  hRes_Zdc_cos_psi1_alice->GetYaxis()->SetRangeUser(-1,0.);
  //hRes_Zdc_cos_psi1_alice->GetXaxis()->SetRangeUser(0.,5.);
  hRes_Zdc_cos_psi1_alice->GetXaxis()->SetTitle("Centrality");
  hRes_Zdc_cos_psi1_alice->GetYaxis()->SetTitle("<Q_{x}^{t}Q_{x}^{p}>");
  hRes_Zdc_sin_psi1_alice->Draw("same");
  leg1->Draw("same");

  TLegend *leg2;
  leg2=new TLegend(0.25,0.7,0.5,0.85);
  //leg2->SetNColumns(2);
  leg2->SetLineColor(0);
  leg2->AddEntry(hRes_Zdc_cos_psi1_atlas,"<Q_{x}^{t}Q_{x}^{p}> ATLAS","P");
  leg2->AddEntry(hRes_Zdc_sin_psi1_atlas,"<Q_{y}^{t}Q_{y}^{p}>ATLAS","P");
  setStyle(hRes_Zdc_cos_psi1_atlas,0);
  setStyle(hRes_Zdc_sin_psi1_atlas,1);
  c->cd(2);
  hRes_Zdc_cos_psi1_atlas->Draw();
  hRes_Zdc_cos_psi1_atlas->GetYaxis()->SetRangeUser(-0.055,0.);
  //hRes_Zdc_cos_psi1_atlas->GetXaxis()->SetRangeUser(0.,5.);
  hRes_Zdc_cos_psi1_atlas->GetXaxis()->SetTitle("Centrality");
  hRes_Zdc_cos_psi1_atlas->GetYaxis()->SetTitle("<Q_{x}^{t}Q_{x}^{p}>");
  hRes_Zdc_sin_psi1_atlas->Draw("same");
  leg2->Draw("same");


  TCanvas *c2;
  c2=new TCanvas("c2","c2",1200,600);
  c2->Divide(2);
  
  TLegend *leg3;
  leg3=new TLegend(0.15,0.7,0.85,0.85);
  //leg3->SetNColumns(2);
  leg3->SetLineColor(0);
  leg3->AddEntry(hRes_Zdc_psi1_alice,"sqrt((<Q_{x}^{t}Q_{x}^{p}>+<Q_{y}^{t}Q_{y}^{p}>)/2) ALICE","P");
  leg3->AddEntry(hRes_Zdc_psi1_atlas,"sqrt((<Q_{x}^{t}Q_{x}^{p}>+<Q_{y}^{t}Q_{y}^{p}>)/2) ATLAS","P");
  setStyle(hRes_Zdc_psi1_alice,0);
  setStyle(hRes_Zdc_psi1_atlas,1);

  c2->cd(1);
  hRes_Zdc_psi1_alice->Draw();
  hRes_Zdc_psi1_alice->GetYaxis()->SetRangeUser(0.,1.5);
  //hRes_Zdc_psi1_alice->GetXaxis()->SetRangeUser(0.,5.);
  hRes_Zdc_psi1_alice->GetXaxis()->SetTitle("Centrality");
  hRes_Zdc_psi1_alice->GetYaxis()->SetTitle("<Q_{x}^{t}Q_{x}^{p}>");
  hRes_Zdc_psi1_atlas->Draw("same");
  leg3->Draw("same");
  //line0->Draw("same");

  TH1D *hR;
  hR=(TH1D*)hRes_Zdc_psi1_alice->Clone();
  hR->Divide(hRes_Zdc_psi1_atlas);

  c2->cd(2);
  hR->Draw();
  hR->GetYaxis()->SetRangeUser(2.,8.);
  hR->GetXaxis()->SetTitle("Centrality");
  hR->GetYaxis()->SetTitle("Ratio");
  
  //hR->GetYaxis()->SetRangeUser(0.,0.3);

  //Save Canvas:
  //c->SaveAs("v1_pt_old.pdf");
  
}



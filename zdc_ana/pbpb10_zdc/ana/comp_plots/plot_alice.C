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

#define withBG

//Variables and Hists
char name[300];
static double PI = acos(-1.0);
const int NCENT=10;
const int NHAR=4;
const int NPT=10;
const int NSIDE=2;
const int NETA=12;
const double pt_mat[]={0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 7.0, 9.0, 12.0};


TProfile *hMean_pt[NCENT][3][NSIDE+1],*hMean_eta[NCENT][3][NSIDE+1];
TH1D *hResCent_mix_FCal_sub[NHAR],*hResCent_mix_FCal_ful[NHAR];
TH1D *hResCent_mix_Zdc_sub[NHAR],*hResCent_mix_Zdc_ful[NHAR];

TH1D *hVn_Zdc_eta_raw[NCENT][NSIDE+1][NHAR][NPT+2],*hVn_Zdc_eta_raw_fg[NCENT][NSIDE+1][NHAR][NPT+2],*hVn_Zdc_eta_raw_bg[NCENT][NSIDE+1][NHAR][NPT+2];
  TH1D *hVn_Zdc_eta_res[NCENT][NSIDE+1][NHAR][NPT+2],*hVn_Zdc_eta_res_fg[NCENT][NSIDE+1][NHAR][NPT+2],*hVn_Zdc_eta_res_bg[NCENT][NSIDE+1][NHAR
][NPT+2];
  TH1D *hVn_Zdc_pt_raw[NCENT][NSIDE+1][NHAR][NETA+4],*hVn_Zdc_pt_raw_fg[NCENT][NSIDE+1][NHAR][NETA+4],*hVn_Zdc_pt_raw_bg[NCENT][NSIDE+1][NHAR]
[NETA+4];
  TH1D *hVn_Zdc_pt_res[NCENT][NSIDE+1][NHAR][NETA+4],*hVn_Zdc_pt_res_fg[NCENT][NSIDE+1][NHAR][NETA+4],*hVn_Zdc_pt_res_bg[NCENT][NSIDE+1][NHAR]
[NETA+4];
  TH1D *hVn_Zdc_cent_raw[NSIDE+1][NHAR][NPT+2][NETA+4],*hVn_Zdc_cent_raw_fg[NSIDE+1][NHAR][NPT+2][NETA+4],*hVn_Zdc_cent_raw_bg[NSIDE+1][NHAR][
NPT+2][NETA+4];
  TH1D *hVn_Zdc_cent_res[NSIDE+1][NHAR][NPT+2][NETA+4],*hVn_Zdc_cent_res_fg[NSIDE+1][NHAR][NPT+2][NETA+4],*hVn_Zdc_cent_res_bg[NSIDE+1][NHAR][
NPT+2][NETA+4];

TH1D *hVn_Zdc_pt_raw_allCent[NSIDE+1][NHAR][NETA+4],*hVn_Zdc_pt_raw_allCent_fg[NSIDE+1][NHAR][NETA+4],*hVn_Zdc_pt_raw_allCent_bg[NSIDE+1][NHAR][NETA+4];
TH1D *hVn_Zdc_pt_res_allCent[NSIDE+1][NHAR][NETA+4],*hVn_Zdc_pt_res_allCent_fg[NSIDE+1][NHAR][NETA+4],*hVn_Zdc_pt_res_allCent_bg[NSIDE+1][NHAR][NETA+4];

TH1D *hVnEven_cent_raw[NHAR][NPT+2][NETA+4],*hVnOdd_cent_raw[NHAR][NPT+2][NETA+4];
TH1D *hVnEven_cent_res[NHAR][NPT+2][NETA+4],*hVnOdd_cent_res[NHAR][NPT+2][NETA+4];
TH1D *hVnEven_eta_raw[NCENT][NHAR][NPT+2],*hVnOdd_eta_raw[NCENT][NHAR][NPT+2];
TH1D *hVnEven_eta_res[NCENT][NHAR][NPT+2],*hVnOdd_eta_res[NCENT][NHAR][NPT+2];


TProfile *hMean_pt_allCent[NSIDE+1],*hMean_eta_allCent[NSIDE+1];
TH1D *hVnEven_pt_raw_allCent[NHAR][NETA+4],*hVnOdd_pt_raw_allCent[NHAR][NETA+4];
TH1D *hVnEven_pt_res_allCent[NHAR][NETA+4],*hVnOdd_pt_res_allCent[NHAR][NETA+4];

TGraphErrors *gVnEven_pt_raw_allCent[NHAR][NETA+4];
TGraphErrors *gVnEven_pt_res_allCent[NHAR][NETA+4];
TGraphErrors *gVnOdd_pt_raw_allCent[NHAR][NETA+4];
TGraphErrors *gVnOdd_pt_res_allCent[NHAR][NETA+4];

TGraphErrors *gVnEven_eta_raw[NCENT][NHAR][NPT+2];
TGraphErrors *gVnEven_eta_res[NCENT][NHAR][NPT+2];
TGraphErrors *gVnOdd_eta_raw[NCENT][NHAR][NPT+2];
TGraphErrors *gVnOdd_eta_res[NCENT][NHAR][NPT+2];


/*
  enum {
  NDET = 16,
  NSIDE   = 2,
  NMODULE = 2,
  NGAIN   = 2,
  NCH1    = 4,
  NCH2    = 24,
  NHAR    = 4,
  //NCENT   = 20,
  NCENT   = 10,
  NCENT_1p  = 100,
  NK      = 12,
  NSAMP   = 10,
  //STORE_DEP = 20,
  STORE_DEP = 5,
  NVZ_EP=10,
  NPT_EP=10,
  NCH_EP=2
  };
*/
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
  int m_mat[]={20,24,21,21,24,24,25,28,33,20,23,29,34};
  h->SetTitle("");

  h->SetMarkerStyle(m_mat[itype1]);
  h->SetMarkerSize(1.5);
  h->SetMarkerColor(c_mat[itype1]);
  h->SetLineColor(c_mat[itype1]);

  //h->GetYaxis()->SetRangeUser(0.,5.);
}


void FindEvenOdd(TH1D *hEven,TH1D *hOdd,TH1D *h1,TH1D *h2){
  //hEven=(TH1D*)h1->Clone();
  //hOdd=(TH1D*)h2->Clone();
  hEven->Reset();
  hOdd->Reset();
  int N=h1->GetNbinsX();
  double y1,y2,dy1,dy2,yE,yO,dyE,dyO;
  for(int i=0;i<N;i++){
    y1=h1->GetBinContent(i+1);
    y2=h2->GetBinContent(i+1);
    dy1=h1->GetBinError(i+1);
    dy2=h2->GetBinError(i+1);
    yE=-0.5 * (y1+y2);
    yO=-0.5 * (y2-y1);
    dyE=0.5 * sqrt(dy1*dy1+dy2*dy2);
    dyO=0.5 * sqrt(dy1*dy1+dy2*dy2);
    hEven->SetBinContent(i+1,yE);
    hEven->SetBinError(i+1,dyE);
    hOdd->SetBinContent(i+1,yO);
    hOdd->SetBinError(i+1,dyO);
  }
}


void plot_alice(){
  
#ifdef withBG
  TFile *fin=new TFile("../zdcvn11/Output_doAna_vn_mix3.root","read");
#endif
#ifndef withBG
  TFile *fin=new TFile("../zdcvn11/Output_doAna_vn_nomix.root","read");
#endif
  for (int icent=0; icent<NCENT; icent++) {
    for (int iside=0; iside<NSIDE+1; iside++) {
      sprintf(name, "hMean_pt_cent%d_ch%d_side%d",icent,2,iside);
      hMean_pt[icent][2][iside]=(TProfile*)fin->Get(name);
      sprintf(name, "hMean_eta_cent%d_ch%d_side%d",icent,2,iside);
      hMean_eta[icent][2][iside]=(TProfile*)fin->Get(name);
    }
  }

  for (int iside=0; iside<NSIDE+1; iside++) {
    sprintf(name, "hMean_pt_allCent_side%d",iside);
    hMean_pt_allCent[iside]=(TProfile*)hMean_pt[1][2][iside]->Clone(name);
    sprintf(name, "hMean_eta_allCent_side%d",iside);
    hMean_eta_allCent[iside]=(TProfile*)hMean_eta[1][2][iside]->Clone(name);
    for (int icent=2; icent<9; icent++) {
      hMean_pt_allCent[iside]->Add(hMean_pt[icent][2][iside],1);
      hMean_eta_allCent[iside]->Add(hMean_eta[icent][2][iside],1);
    }
  }
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
  
  for (int ihar=0; ihar<NHAR; ihar++) {
    for (int iside=0; iside<NSIDE+1; iside++) {
      for (int ieta=0; ieta<NETA+4; ieta++) {
	sprintf(name, "hV%d_Zdc_pt_raw_allCent_side%d_eta%d", ihar+1,iside,ieta);
	hVn_Zdc_pt_raw_allCent[iside][ihar][ieta]=(TH1D*)fin->Get(name);
	sprintf(name, "hV%d_Zdc_pt_res_allCent_side%d_eta%d", ihar+1,iside,ieta);
	hVn_Zdc_pt_res_allCent[iside][ihar][ieta]=(TH1D*)fin->Get(name);

	sprintf(name, "hV%d_Zdc_cent_raw_side%d_pt%d_eta%d", ihar+1,iside,NPT,ieta);
	hVn_Zdc_cent_raw[iside][ihar][NPT][ieta]=(TH1D*)fin->Get(name);
	sprintf(name, "hV%d_Zdc_cent_res_side%d_pt%d_eta%d", ihar+1,iside,NPT,ieta);
	hVn_Zdc_cent_res[iside][ihar][NPT][ieta]=(TH1D*)fin->Get(name);
      }
      for(int icent=0;icent<NCENT;icent++){
	sprintf(name, "hV%d_Zdc_eta_raw_cent%d_side%d_pt%d", ihar+1,icent,iside,NPT);
	hVn_Zdc_eta_raw[icent][iside][ihar][NPT]=(TH1D*)fin->Get(name);
	sprintf(name, "hV%d_Zdc_eta_res_cent%d_side%d_pt%d", ihar+1,icent,iside,NPT);
	hVn_Zdc_eta_res[icent][iside][ihar][NPT]=(TH1D*)fin->Get(name);
      }
    }
  }
  
  for (int ihar=0; ihar<NHAR; ihar++) {
    for (int ieta=0; ieta<NETA+4; ieta++) {

      //Raw
      sprintf(name, "hV%dEven_pt_raw_allCent_eta%d", ihar+1,  ieta);
      hVnEven_pt_raw_allCent[ihar][ieta]=(TH1D*)hVn_Zdc_pt_raw_allCent[0][ihar][ieta]->Clone(name);
      sprintf(name, "hV%dOdd_pt_raw_allCent_eta%d", ihar+1,  ieta);
      hVnOdd_pt_raw_allCent[ihar][ieta]=(TH1D*)hVn_Zdc_pt_raw_allCent[0][ihar][ieta]->Clone(name);
      FindEvenOdd(hVnEven_pt_raw_allCent[ihar][ieta],hVnOdd_pt_raw_allCent[ihar][ieta],hVn_Zdc_pt_raw_allCent[0][ihar][ieta],hVn_Zdc_pt_raw_allCent[1][ihar][ieta]);
      
      //Resolution corrected
      sprintf(name, "hV%dEven_pt_res_allCent_eta%d", ihar+1,  ieta);
      hVnEven_pt_res_allCent[ihar][ieta]=(TH1D*)hVn_Zdc_pt_res_allCent[0][ihar][ieta]->Clone(name);
      sprintf(name, "hV%dOdd_pt_res_allCent_eta%d", ihar+1,  ieta);
      hVnOdd_pt_res_allCent[ihar][ieta]=(TH1D*)hVn_Zdc_pt_res_allCent[0][ihar][ieta]->Clone(name);
      FindEvenOdd(hVnEven_pt_res_allCent[ihar][ieta],hVnOdd_pt_res_allCent[ihar][ieta],hVn_Zdc_pt_res_allCent[0][ihar][ieta],hVn_Zdc_pt_res_allCent[1][ihar][ieta]);

      //As a function of centrality
      sprintf(name, "hV%dEven_cent_raw_pt%d_eta%d", ihar+1, NPT,ieta);
      hVnEven_cent_raw[ihar][NPT][ieta]=(TH1D*)hVn_Zdc_cent_raw[0][ihar][NPT][ieta]->Clone(name);
      sprintf(name, "hV%dOdd_cent_raw_pt%d_eta%d", ihar+1, NPT,ieta);
      hVnOdd_cent_raw[ihar][NPT][ieta]=(TH1D*)hVn_Zdc_cent_raw[0][ihar][NPT][ieta]->Clone(name);
      FindEvenOdd(hVnEven_cent_raw[ihar][NPT][ieta],hVnOdd_cent_raw[ihar][NPT][ieta],hVn_Zdc_cent_raw[0][ihar][NPT][ieta],hVn_Zdc_cent_raw[1][ihar][NPT][ieta]);
	
      sprintf(name, "hV%dEven_cent_res_pt%d_eta%d", ihar+1, NPT,ieta);
      hVnEven_cent_res[ihar][NPT][ieta]=(TH1D*)hVn_Zdc_cent_res[0][ihar][NPT][ieta]->Clone(name);
      sprintf(name, "hV%dOdd_cent_res_pt%d_eta%d", ihar+1, NPT,ieta);
      hVnOdd_cent_res[ihar][NPT][ieta]=(TH1D*)hVn_Zdc_cent_res[0][ihar][NPT][ieta]->Clone(name);
      FindEvenOdd(hVnEven_cent_res[ihar][NPT][ieta],hVnOdd_cent_res[ihar][NPT][ieta],hVn_Zdc_cent_res[0][ihar][NPT][ieta],hVn_Zdc_cent_res[1][ihar][NPT][ieta]);
    }
  }
  
  //As a function of eta
  for (int ihar=0; ihar<NHAR; ihar++) {
    for(int icent=0;icent<NCENT;icent++){
      sprintf(name, "hV%dEven_cent_raw_cent%d_pt%d", ihar+1,icent,NPT);
      hVnEven_eta_raw[icent][ihar][NPT]=(TH1D*)hVn_Zdc_eta_raw[icent][0][ihar][NPT]->Clone(name);
      sprintf(name, "hV%dOdd_cent_raw_cent%d_pt%d", ihar+1,icent,NPT);
      hVnOdd_eta_raw[icent][ihar][NPT]=(TH1D*)hVn_Zdc_eta_raw[icent][0][ihar][NPT]->Clone(name);
      FindEvenOdd(hVnEven_eta_raw[icent][ihar][NPT],hVnOdd_eta_raw[icent][ihar][NPT],hVn_Zdc_eta_raw[icent][0][ihar][NPT],hVn_Zdc_eta_raw[icent][1][ihar][NPT]);

      sprintf(name, "hV%dEven_cent_res_cent%d_pt%d", ihar+1,icent,NPT);
      hVnEven_eta_res[icent][ihar][NPT]=(TH1D*)hVn_Zdc_eta_res[icent][0][ihar][NPT]->Clone(name);
      sprintf(name, "hV%dOdd_cent_res_cent%d_pt%d", ihar+1,icent,NPT);
      hVnOdd_eta_res[icent][ihar][NPT]=(TH1D*)hVn_Zdc_eta_res[icent][0][ihar][NPT]->Clone(name);
      FindEvenOdd(hVnEven_eta_res[icent][ihar][NPT],hVnOdd_eta_res[icent][ihar][NPT],hVn_Zdc_eta_res[icent][0][ihar][NPT],hVn_Zdc_eta_res[icent][1][ihar][NPT]);
    }
  }

  //Get TGraphs from TH1D
  double y,y1,y2,yc,dy,dy1,dy2,dyc,yc2;
  double X[100],dX[100],Y[100],dY[100];
  int N;
  int ihar=0,iside=2;
  //As a function of pt
  for(int ieta=0;ieta<NETA+4;ieta++){
    N=hVnEven_pt_res_allCent[ihar][ieta]->GetNbinsX();
    for (int ipt=0; ipt<N; ipt++) {
      y1=hVnEven_pt_res_allCent[ihar][ieta]->GetBinContent(ipt+1);
      dy1=hVnEven_pt_res_allCent[ihar][ieta]->GetBinError(ipt+1);
      yc=y1;
      dyc=dy1;
      X[ipt]=hMean_pt_allCent[iside]->GetBinContent(ipt+1);
      Y[ipt]=yc;
      dY[ipt]=dyc;
    }
    gVnEven_pt_res_allCent[ihar][ieta]=new TGraphErrors(N,X,Y,0,dY);
    sprintf(name,"gVnEven_pt_res_allCent_har%d_eta%d",ihar,ieta);
    gVnEven_pt_res_allCent[ihar][ieta]->SetName(name);
  }

  for(int ieta=0;ieta<NETA+4;ieta++){
    N=hVnOdd_pt_res_allCent[ihar][ieta]->GetNbinsX();
    for (int ipt=0; ipt<N; ipt++) {
      y1=hVnOdd_pt_res_allCent[ihar][ieta]->GetBinContent(ipt+1);
      dy1=hVnOdd_pt_res_allCent[ihar][ieta]->GetBinError(ipt+1);
      yc=y1;
      dyc=dy1;
      X[ipt]=hMean_pt_allCent[iside]->GetBinContent(ipt+1);
      Y[ipt]=yc;
      dY[ipt]=dyc;
    }
    gVnOdd_pt_res_allCent[ihar][ieta]=new TGraphErrors(N,X,Y,0,dY);
    sprintf(name,"gVnOdd_pt_res_allCent_har%d_eta%d",ihar,ieta);
    gVnOdd_pt_res_allCent[ihar][ieta]->SetName(name);
  }

  //As a function of eta
  for(int icent=0;icent<NCENT;icent++){
    N=hVnEven_eta_res[icent][ihar][NPT]->GetNbinsX();
    for (int ieta=0; ieta<N; ieta++) {
      y1=hVnEven_eta_res[icent][ihar][NPT]->GetBinContent(ieta+1);
      dy1=hVnEven_eta_res[icent][ihar][NPT]->GetBinError(ieta+1);
      yc=y1;
      dyc=dy1;
      X[ieta]=hMean_eta[icent][2][iside]->GetBinContent(ieta+1);
      Y[ieta]=yc;
      dY[ieta]=dyc;
    }
    gVnEven_eta_res[icent][ihar][NPT]=new TGraphErrors(N,X,Y,0,dY);
    sprintf(name,"gVnEven_eta_res_cent%d_har%d_pt%d",icent,ihar,NPT);
    gVnEven_eta_res[icent][ihar][NPT]->SetName(name);
  }

  for(int icent=0;icent<NCENT;icent++){
    N=hVnOdd_eta_res[icent][ihar][NPT]->GetNbinsX();
    for (int ieta=0; ieta<N; ieta++) {
      y1=hVnOdd_eta_res[icent][ihar][NPT]->GetBinContent(ieta+1);
      dy1=hVnOdd_eta_res[icent][ihar][NPT]->GetBinError(ieta+1);
      yc=y1;
      dyc=dy1;
      X[ieta]=hMean_eta[icent][2][iside]->GetBinContent(ieta+1);
      Y[ieta]=yc;
      dY[ieta]=dyc;
    }
    gVnOdd_eta_res[icent][ihar][NPT]=new TGraphErrors(N,X,Y,0,dY);
    sprintf(name,"gVnOdd_eta_res_cent%d_har%d_pt%d",icent,ihar,NPT);
    gVnOdd_eta_res[icent][ihar][NPT]->SetName(name);
  }

  //Get graphs for Alice published results
  ifstream inf;
  inf.open("../pubdata/data_v1_even_alice.txt");
  for(int i=0;i<10;i++){
    double Z[6];
    inf>>Z[0]>>Z[1]>>Z[2]>>Z[3]>>Z[4]>>Z[5];
    X[i]=Z[0];
    Y[i]=Z[3];
    dY[i]=Z[4];
  }
  inf.close();
  TGraphErrors *gVnEven_pt_res_allCent_pub;
  gVnEven_pt_res_allCent_pub=new TGraphErrors(10,X,Y,0,dY);
  gVnEven_pt_res_allCent_pub->SetName("gVnEven_pt_res_allCent_pub");
  
  inf.open("../pubdata/data_v1_odd_alice.txt");
  for(int i=0;i<10;i++){
    double Z[6];
    inf>>Z[0]>>Z[1]>>Z[2]>>Z[3]>>Z[4]>>Z[5];
    X[i]=Z[0];
    Y[i]=Z[3];
    dY[i]=Z[4];
  }
  inf.close();
  TGraphErrors *gVnOdd_pt_res_allCent_pub;
  gVnOdd_pt_res_allCent_pub=new TGraphErrors(10,X,Y,0,dY);
  gVnOdd_pt_res_allCent_pub->SetName("gVnOdd_pt_res_allCent_pub");

  //As a function of eta
  //10-20% Centrality
  inf.open("../pubdata/data_v1_odd_eta_alice.txt");
  for(int i=0;i<10;i++){
    double Z[5];
    inf>>Z[0]>>Z[1]>>Z[2]>>Z[3]>>Z[4];
    X[i]=Z[0];
    Y[i]=Z[3];
    dY[i]=Z[4];
  }
  inf.close();
  TGraphErrors *gVnOdd_eta_res_cent2_pub;
  gVnOdd_eta_res_cent2_pub=new TGraphErrors(10,X,Y,0,dY);
  gVnOdd_eta_res_cent2_pub->SetName("gVnOdd_eta_res_cent2_pub");
  
  inf.open("../pubdata/data_v1_even_eta_alice.txt");
  for(int i=0;i<10;i++){
    double Z[5];
    inf>>Z[0]>>Z[1]>>Z[2]>>Z[3]>>Z[4];
    X[i]=Z[0];
    Y[i]=Z[3];
    dY[i]=Z[4];
  }
  inf.close();
  TGraphErrors *gVnEven_eta_res_cent2_pub;
  gVnEven_eta_res_cent2_pub=new TGraphErrors(10,X,Y,0,dY);
  gVnEven_eta_res_cent2_pub->SetName("gVnEven_eta_res_cent2_pub");

  //30-40% Centrality
  inf.open("../pubdata/data_v1_odd_eta2_alice.txt");
  for(int i=0;i<10;i++){
    double Z[5];
    inf>>Z[0]>>Z[1]>>Z[2]>>Z[3]>>Z[4];
    X[i]=Z[0];
    Y[i]=Z[3];
    dY[i]=Z[4];
  }
  inf.close();
  TGraphErrors *gVnOdd_eta_res_cent4_pub;
  gVnOdd_eta_res_cent4_pub=new TGraphErrors(10,X,Y,0,dY);
  gVnOdd_eta_res_cent4_pub->SetName("gVnOdd_eta_res_cent4_pub");

  inf.open("../pubdata/data_v1_even_eta2_alice.txt");
  for(int i=0;i<10;i++){
    double Z[5];
    inf>>Z[0]>>Z[1]>>Z[2]>>Z[3]>>Z[4];
    X[i]=Z[0];
    Y[i]=Z[3];
    dY[i]=Z[4];
  }
  inf.close();
  TGraphErrors *gVnEven_eta_res_cent4_pub;
  gVnEven_eta_res_cent4_pub=new TGraphErrors(10,X,Y,0,dY);
  gVnEven_eta_res_cent4_pub->SetName("gVnEven_eta_res_cent4_pub");
  
  //As a function of cent
  inf.open("../pubdata/data_v1_odd_cent_alice.txt");
  for(int i=0;i<10;i++){
    double Z[5];
    inf>>Z[0]>>Z[1]>>Z[2]>>Z[3]>>Z[4];
    X[i]=Z[0];
    Y[i]=Z[3];
    dY[i]=Z[4];
  }
  inf.close();
  TGraphErrors *gVnOdd_cent_res_pub;
  gVnOdd_cent_res_pub=new TGraphErrors(10,X,Y,0,dY);
  gVnOdd_cent_res_pub->SetName("gVnOdd_cent_res_pub");

  inf.open("../pubdata/data_v1_even_cent_alice.txt");
  for(int i=0;i<10;i++){
    double Z[5];
    inf>>Z[0]>>Z[1]>>Z[2]>>Z[3]>>Z[4];
    X[i]=Z[0];
    Y[i]=Z[3];
    dY[i]=Z[4];
  }
  inf.close();
  TGraphErrors *gVnEven_cent_res_pub;
  gVnEven_cent_res_pub=new TGraphErrors(10,X,Y,0,dY);
  gVnEven_cent_res_pub->SetName("gVnEven_cent_res_pub");



  //Plot
  //NETA+2 is eta normal used for even vn and NETA+3 is eta flipped (vn=-vn for eta<0) used for odd vn
  int ieta1=NETA+2,ieta2=NETA+3;
  gStyle->SetOptStat(0);
  TCanvas *c;
  c=new TCanvas("c","c",1200,600);
  c->Divide(2);
  
  TLine *line0=new TLine(0.,0.,5.,0.);
  line0->SetLineStyle(7);
  line0->SetLineWidth(1);

  TLegend *leg1;
  leg1=new TLegend(0.15,0.7,0.5,0.88);
  //leg1->SetNColumns(2);
  leg1->SetLineColor(0);
  leg1->AddEntry(gVnEven_pt_res_allCent_pub,"v_{1} Even (ALICE)","P");
  leg1->AddEntry(gVnEven_pt_res_allCent[ihar][ieta1],"v_{1} Even (ATLAS)","P");
  setStyle2(gVnEven_pt_res_allCent_pub,0);
  setStyle2(gVnEven_pt_res_allCent[ihar][ieta1],1);

  c->cd(1);
  gVnEven_pt_res_allCent[ihar][ieta1]->Draw("AP");
  gVnEven_pt_res_allCent[ihar][ieta1]->SetTitle("");
  gVnEven_pt_res_allCent[ihar][ieta1]->GetYaxis()->SetRangeUser(-0.002,0.004);
  gVnEven_pt_res_allCent[ihar][ieta1]->GetXaxis()->SetRangeUser(0.,5.);
  gVnEven_pt_res_allCent[ihar][ieta1]->GetXaxis()->SetTitle("p_{T}");
  gVnEven_pt_res_allCent[ihar][ieta1]->GetYaxis()->SetTitle("v_{1}");
  gVnEven_pt_res_allCent_pub->Draw("P");
  //hVn_Zdc_pt_res_allCent[2][ihar][ieta1]->Draw("same");
  leg1->Draw("same");
  line0->Draw("same");

  TLegend *leg2;
  leg2=new TLegend(0.15,0.7,0.5,0.88);
  //leg2->SetNColumns(2);
  leg2->SetLineColor(0);
  leg2->AddEntry(gVnOdd_pt_res_allCent_pub,"v_{1} Odd (ALICE)","P");
  leg2->AddEntry(gVnOdd_pt_res_allCent[ihar][ieta2],"v_{1} Odd (ATLAS)","P");
  setStyle2(gVnOdd_pt_res_allCent_pub,0);
  setStyle2(gVnOdd_pt_res_allCent[ihar][ieta2],1);
  c->cd(2);
  gVnOdd_pt_res_allCent[ihar][ieta2]->Draw("AP");
  gVnOdd_pt_res_allCent[ihar][ieta2]->SetTitle("");
  gVnOdd_pt_res_allCent[ihar][ieta2]->GetYaxis()->SetRangeUser(-0.002,0.004);
  gVnOdd_pt_res_allCent[ihar][ieta2]->GetXaxis()->SetRangeUser(0.,5.);
  gVnOdd_pt_res_allCent[ihar][ieta2]->GetXaxis()->SetTitle("p_{T}");
  gVnOdd_pt_res_allCent[ihar][ieta2]->GetYaxis()->SetTitle("v_{1}");
  gVnOdd_pt_res_allCent_pub->Draw("P");
  //hVn_Zdc_pt_res_allCent[2][ihar][ieta2]->Draw("same");
  leg2->Draw("same");
  line0->Draw("same");
  
  
  int ipt=NPT;
  TCanvas *c2;
  c2=new TCanvas("c2","c2",800,600);
  //c2->Divide(2);
  
  TLegend *leg3;
  leg3=new TLegend(0.15,0.65,0.55,0.88);
  //leg3->SetNColumns(2);
  leg3->SetLineColor(0);
  leg3->AddEntry( hVnEven_cent_res[ihar][ipt][ieta1],"v_{1} Even ATLAS (p_{T}>0.5GeV)","P");
  leg3->AddEntry(hVnOdd_cent_res[ihar][ipt][ieta2],"v_{1} Odd ATLAS (p_{T}>0.5GeV)","P");
  leg3->AddEntry(gVnEven_cent_res_pub,"v_{1} Even ALICE (p_{T}>0.15GeV)","P");
  leg3->AddEntry(gVnOdd_cent_res_pub,"v_{1} Odd ALICE (p_{T}>0.15GeV)","P");

  TLine *line0_2=new TLine(0.,0.,100.,0.);
  line0_2->SetLineStyle(7);
  line0_2->SetLineWidth(1);
  setStyle(hVnEven_cent_res[ihar][ipt][ieta1],0);
  setStyle(hVnOdd_cent_res[ihar][ipt][ieta2],1);
  setStyle2(gVnEven_cent_res_pub,2);
  setStyle2(gVnOdd_cent_res_pub,3);
	    
  c2->cd();
  hVnEven_cent_res[ihar][ipt][ieta1]->Draw();
  hVnEven_cent_res[ihar][ipt][ieta1]->GetXaxis()->SetTitle("Centrality");
  hVnEven_cent_res[ihar][ipt][ieta1]->GetYaxis()->SetTitle("v_{1}");
  hVnEven_cent_res[ihar][ipt][ieta1]->GetYaxis()->SetRangeUser(-1e-3,1e-3);
  hVnOdd_cent_res[ihar][ipt][ieta2]->Draw("same");
  gVnEven_cent_res_pub->Draw("P");
  gVnOdd_cent_res_pub->Draw("P");
  leg3->Draw("same");
  line0_2->Draw("same");


  int icent1=2,icent2=4;
  TCanvas *c3;
  c3=new TCanvas("c3","c3",1200,600);
  c3->Divide(2);
  
  TLine *line0_3=new TLine(-2.5,0.,2.5,0.);
  line0_3->SetLineStyle(7);
  line0_3->SetLineWidth(1);


  TLegend *leg4;
  leg4=new TLegend(0.45,0.65,0.88,0.88);
  //leg4->SetNColumns(2);
  leg4->SetLineColor(0);
  leg4->AddEntry( gVnEven_eta_res[icent1][ihar][ipt],"v_{1} Even ATLAS (p_{T}>0.5GeV)","P");
  leg4->AddEntry(gVnOdd_eta_res[icent1][ihar][ipt],"v_{1} Odd ATLAS (p_{T}>0.5GeV)","P");
  leg4->AddEntry(gVnEven_eta_res_cent2_pub,"v_{1} Even ALICE (p_{T}>0.15GeV)","P");
  leg4->AddEntry(gVnOdd_eta_res_cent2_pub,"v_{1} Odd ALICE (p_{T}>0.15GeV)","P");
  
  setStyle(hVnEven_eta_res[icent1][ihar][ipt],2);
  setStyle(hVnOdd_eta_res[icent1][ihar][ipt],3);
  setStyle2(gVnEven_eta_res[icent1][ihar][ipt],4);
  setStyle2(gVnOdd_eta_res[icent1][ihar][ipt],5);
  setStyle2(gVnEven_eta_res_cent2_pub,2);
  setStyle2(gVnOdd_eta_res_cent2_pub,3);

  c3->cd(1);
  gVnOdd_eta_res[icent1][ihar][ipt]->Draw("AP");
  //hVnOdd_eta_res[icent1][ihar][ipt]->GetXaxis()->SetRangeUser(-0.8,0.8);
  gVnOdd_eta_res[icent1][ihar][ipt]->GetXaxis()->SetTitle("#eta");
  gVnOdd_eta_res[icent1][ihar][ipt]->GetYaxis()->SetTitle("v_{1}");
  gVnEven_eta_res[icent1][ihar][ipt]->Draw("same P");
  gVnOdd_eta_res_cent2_pub->Draw("P");
  gVnEven_eta_res_cent2_pub->Draw("P");
  leg4->Draw("same");
  line0_3->Draw("same");

  TLegend *leg5;
  leg5=new TLegend(0.45,0.65,0.88,0.88);
  //leg5->SetNColumns(2);
  leg5->SetLineColor(0);
  leg5->AddEntry( gVnEven_eta_res[icent2][ihar][ipt],"v_{1} Even ATLAS (p_{T}>0.5GeV)","P");
  leg5->AddEntry(gVnOdd_eta_res[icent2][ihar][ipt],"v_{1} Odd ATLAS (p_{T}>0.5GeV)","P");
  leg5->AddEntry(gVnEven_eta_res_cent4_pub,"v_{1} Even ALICE (p_{T}>0.15GeV)","P");
  leg5->AddEntry(gVnOdd_eta_res_cent4_pub,"v_{1} Odd ALICE (p_{T}>0.15GeV)","P");

  setStyle(hVnEven_eta_res[icent2][ihar][ipt],2);
  setStyle(hVnOdd_eta_res[icent2][ihar][ipt],3);
  setStyle2(gVnEven_eta_res[icent2][ihar][ipt],4);
  setStyle2(gVnOdd_eta_res[icent2][ihar][ipt],5);
  setStyle2(gVnEven_eta_res_cent4_pub,2);
  setStyle2(gVnOdd_eta_res_cent4_pub,3);

  c3->cd(2);
  gVnOdd_eta_res[icent2][ihar][ipt]->Draw("AP");
  //gVnOdd_eta_res[icent2][ihar][ipt]->GetXaxis()->SetRangeUser(-0.8,0.8);
  gVnOdd_eta_res[icent2][ihar][ipt]->GetXaxis()->SetTitle("#eta");
  gVnOdd_eta_res[icent2][ihar][ipt]->GetYaxis()->SetTitle("v_{1}");
  gVnEven_eta_res[icent2][ihar][ipt]->Draw("same P");
  gVnOdd_eta_res_cent4_pub->Draw("P");
  gVnEven_eta_res_cent4_pub->Draw("P");
  leg5->Draw("same");
  line0_3->Draw("same");
  

  //Save Canvas:
  /*
  c->SaveAs("v1_pt_old.pdf");
  c2->SaveAs("v1_cent_old.pdf");
  c3->SaveAs("v1_eta_old.pdf");
  */
  c->SaveAs("v1_pt.pdf");
  c2->SaveAs("v1_cent.pdf");
  c3->SaveAs("v1_eta.pdf");

  TFile *fout=new TFile("./v1plots_2010.root","recreate");
  fout->cd();
  hResCent_mix_Zdc_sub[ihar]->Write();
  hResCent_mix_Zdc_ful[ihar]->Write();
  gVnEven_pt_res_allCent_pub->Write();
  gVnOdd_pt_res_allCent_pub->Write();
  gVnEven_pt_res_allCent[ihar][ieta1]->Write();
  gVnOdd_pt_res_allCent[ihar][ieta2]->Write();
  gVnOdd_eta_res[icent1][ihar][ipt]->Write();
  gVnEven_eta_res[icent1][ihar][ipt]->Write();
  hVnEven_cent_res[ihar][ipt][ieta1]->Write();
  hVnOdd_cent_res[ihar][ipt][ieta2]->Write();
  fout->Close();
}



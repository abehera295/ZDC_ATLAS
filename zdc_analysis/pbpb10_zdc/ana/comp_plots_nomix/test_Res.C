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

int cent_mat[]={0,5,10,20,30,40,50,60,70,80,90,100};
double cent_mat_double[]={0,5,10,20,30,40,50,60,70,80,90,100};

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
  h->SetLineColor(c_mat[itype1]);
  
  h->SetMarkerSize(1.2);

  //h->GetYaxis()->SetRangeUser(0.,5.);
}

void test_Res(){
  TFile *fin=new TFile("../zdcvn11/Output_res.root","read");
  TProfile *hResCent_FCal_fg,*hResCent_FCal_bg,*hResCent_FCal_rat;
  TProfile *hResCent_Zdc_fg,*hResCent_Zdc_bg,*hResCent_Zdc_rat;
  TH1D *hResCent_FCal,*hResCent_Zdc;
  int icent=3,ihar=0;
  sprintf(name,"hResCent_Zdc_fg_har%d",ihar);
  hResCent_Zdc_fg=(TProfile*)fin->Get(name);
  sprintf(name,"hResCent_Zdc_bg_har%d",ihar);
  hResCent_Zdc_bg=(TProfile*)fin->Get(name);

  int N=hResCent_Zdc_fg->GetNbinsX();
  hResCent_Zdc=new TH1D("hResCent_Zdc","",NCENT,0,NCENT);
  for(int i=0;i<N;i++){
    double y1=hResCent_Zdc_fg->GetBinContent(i+1);
    double y2=hResCent_Zdc_bg->GetBinContent(i+1);
    double dy1=hResCent_Zdc_fg->GetBinError(i+1);
    double dy2=hResCent_Zdc_bg->GetBinError(i+1);
    double y=y1-y2;
    double dy=sqrt(pow(dy1,2)+pow(dy2,2));
    hResCent_Zdc->SetBinContent(i+1,y);
    hResCent_Zdc->SetBinError(i+1,dy);
    cout<<i+1<<"\t"<<hResCent_Zdc_fg->GetBinContent(i+1)<<"\t"<<hResCent_Zdc_fg->GetBinEntries(i+1)<<endl;
  }
  //SetStyle
  setStyle(hResCent_Zdc_fg,0);
  setStyle(hResCent_Zdc_bg,1);
  setStyle(hResCent_Zdc,2);
  gStyle->SetOptStat(0);
  
  TCanvas *c;
  c=new TCanvas("c","c",1200,600);
  c->Divide(2);
  TLegend *leg;
  leg=new TLegend(0.4,0.75,0.7,0.88);
  //leg->SetNColumns(3);
  leg->SetLineColor(0);
  leg->AddEntry(hResCent_Zdc_fg,"Same Ev","P");
  leg->AddEntry(hResCent_Zdc,"Same Ev - Mixed Ev","P");
  
  c->cd(1);
  hResCent_Zdc_fg->Draw();
  //hResCent_Zdc_fg->GetYaxis()->SetRangeUser(30e3,70e3);
  hResCent_Zdc_fg->SetXTitle("#Psi_{1}^{A}-#Psi_{1}^{C}");
  hResCent_Zdc->Draw("same");
  leg->Draw("same");

  c->cd(2);
  hResCent_Zdc_bg->Draw();
  //hResCent_Zdc_bg->GetYaxis()->SetRangeUser(500e3,700e3);
  hResCent_Zdc_bg->SetXTitle("#Psi_{1}^{A}-#Psi_{1}^{C}");
}

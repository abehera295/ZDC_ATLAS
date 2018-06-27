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
  int m_mat[]={20,21,33,24,25,27,25,28,33,20,23,29,34};
  h->SetTitle("");

  h->SetMarkerStyle(m_mat[itype1]);
  h->SetMarkerSize(1.5);
  h->SetMarkerColor(c_mat[itype1]);
  h->SetLineColor(c_mat[itype1]);
  //h->SetMarkerColor(4);
  //h->SetLineColor(4);

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
  int c_mat[]={1,2,4,1,2,4,6,28,36,7,8,9,46,30};
  int m_mat[]={24,25,27,20,21,33,30,27,25,28,33,20,23,29,34};
  h->SetTitle("");

  h->SetMarkerStyle(m_mat[itype1]);
  h->SetMarkerSize(1.5);
  h->SetMarkerColor(c_mat[itype1]);
  h->SetLineColor(c_mat[itype1]);
  //h->SetMarkerColor(1);
  //h->SetLineColor(1);

  //h->GetYaxis()->SetRangeUser(0.,5.);
}


void plot_atlas_v2_pt(){
  TProfile *hMean_pt[NCENT][NCH+1][NSIDE+1];
  TH1D *hVn_FCal_cent_res[NHAR],*hVn_Zdc_cent_res[NHAR];
  TH1D *hVn_FCal_pt_res[NHAR][NCENT],*hVn_Zdc_pt_res[NHAR][NCENT];
  TH1D *hVn_FCal_cent_res1,*hVn_Zdc_cent_res1;
  TH1D *hVn_FCal_cent_res2,*hVn_Zdc_cent_res2;
  TGraphErrors *gVn_FCal_pt_res[NHAR][NCENT],*gVn_Zdc_pt_res[NHAR][NCENT];

  /*
  TFile *fin1=new TFile("./zdcvn11/Output_SP.root","read");
  for (int icent=0; icent<NCENT; icent++) {
    for (int iside=0; iside<NSIDE+1; iside++) {
      sprintf(name, "hMean_pt_cent%d_ch%d_side%d",icent,2,iside);
      hMean_pt[icent][2][iside]=(TProfile*)fin1->Get(name);
    }
  }
  */
  //TFile *fin=new TFile("../zdcvn8/Output_doAna_vn.root","read");
  //TFile *fin=new TFile("../zdcvn8/Output_doAna_vn_old.root","read");
  //TFile *fin=new TFile("../zdcvn10/Output_doAna_vn_mix.root","read");
  //TFile *fin=new TFile("../zdcvn11/Output_doAna_vn_mix2.root","read");
  TFile *fin=new TFile("../zdcvn2/Output_doAna_vn_nomix.root","read");
  for(int ihar=0; ihar<NHAR; ihar++) {
    sprintf(name,"hV%d_FCal_cent_res_side%d_pt%d_eta%d",ihar+1,0,NPT,NETA);
    hVn_FCal_cent_res[ihar]=(TH1D*)fin->Get(name);
    sprintf(name,"hV%d_Zdc_cent_res_side%d_pt%d_eta%d",ihar+1,0,NPT,NETA);
    hVn_Zdc_cent_res[ihar]=(TH1D*)fin->Get(name);
    for(int icent=0;icent<NCENT;icent++){
      sprintf(name,"hV%d_FCal_pt_res_cent%d_side%d_eta%d",ihar+1,icent,1,NETA);
      hVn_FCal_pt_res[ihar][icent]=(TH1D*)fin->Get(name);
      sprintf(name,"hV%d_Zdc_pt_res_cent%d_side%d_eta%d",ihar+1,icent,1,NETA);
      hVn_Zdc_pt_res[ihar][icent]=(TH1D*)fin->Get(name);
    }
  }

  for (int icent=0; icent<NCENT; icent++) {
    for (int iside=0; iside<NSIDE+1; iside++) {
      sprintf(name, "hMean_pt_cent%d_ch%d_side%d",icent,2,iside);
      hMean_pt[icent][2][iside]=(TProfile*)fin->Get(name);
    }
  }
  
  ifstream inf;
  double X[100],dX[100],Y[100],dY[100];
  int N;
  //Make graphs from hists
  for(int ihar=0; ihar<NHAR; ihar++) {
    for(int icent=0;icent<NCENT;icent++){
      N=hVn_FCal_pt_res[ihar][icent]->GetNbinsX();
      for(int i=0;i<N;i++){
	X[i]=hMean_pt[icent][2][2]->GetBinContent(i+1);
	Y[i]=hVn_FCal_pt_res[ihar][icent]->GetBinContent(i+1);
	dY[i]=hVn_FCal_pt_res[ihar][icent]->GetBinError(i+1);
      }
      gVn_FCal_pt_res[ihar][icent]=new TGraphErrors(N,X,Y,0,dY);
    }
  }
  

  //Cent 0-5%
  inf.open("../pubdata/data_v2_cent0.txt");
  for(int i=0;i<25;i++){
    char data[100];
    double Z[5];
    inf>>Z[0]>>Z[1]>>Z[2]>>Z[3]>>Z[4];
    X[i]=Z[0];
    Y[i]=Z[3];
    dY[i]=Z[4];
  }
  inf.close();
  TGraphErrors *gVn_FCal_pt_pub_cent0;
  gVn_FCal_pt_pub_cent0=new TGraphErrors(25,X,Y,0,dY);

  //Cent 10-20%
  inf.open("../pubdata/data_v2_cent2.txt");
  for(int i=0;i<25;i++){
    char data[100];
    double Z[5];
    inf>>Z[0]>>Z[1]>>Z[2]>>Z[3]>>Z[4];
    X[i]=Z[0];
    Y[i]=Z[3];
    dY[i]=Z[4];
  }
  inf.close();
  TGraphErrors *gVn_FCal_pt_pub_cent2;
  gVn_FCal_pt_pub_cent2=new TGraphErrors(25,X,Y,0,dY);

  //Cent 30-40%
  inf.open("../pubdata/data_v2_cent4.txt");
  for(int i=0;i<25;i++){
    char data[100];
    double Z[5];
    inf>>Z[0]>>Z[1]>>Z[2]>>Z[3]>>Z[4];
    X[i]=Z[0];
    Y[i]=Z[3];
    dY[i]=Z[4];
  }
  inf.close();
  TGraphErrors *gVn_FCal_pt_pub_cent4;
  gVn_FCal_pt_pub_cent4=new TGraphErrors(25,X,Y,0,dY);

  
  //Plot
  gStyle->SetOptStat(0);

  int ihar=1,icent=2;
  TCanvas *c;
  c=new TCanvas("c","c",800,600);
  TLegend *leg;
  leg=new TLegend(0.15,0.75,0.85,0.88);
  leg->SetNColumns(2);
  leg->SetLineColor(0);
  leg->AddEntry(gVn_FCal_pt_pub_cent0,"v_{2} (0-5%%) FCal Pub","P");
  leg->AddEntry(gVn_FCal_pt_res[ihar][0],"v_{2} (0-5%%)  FCal","P");
  leg->AddEntry(gVn_FCal_pt_pub_cent2,"v_{2} (10-20%%) FCal Pub","P");
  leg->AddEntry(gVn_FCal_pt_res[ihar][2],"v_{2} (10-20%%)  FCal","P");
  leg->AddEntry(gVn_FCal_pt_pub_cent4,"v_{2} (30-40%%) FCal Pub","P");
  leg->AddEntry(gVn_FCal_pt_res[ihar][4],"v_{2} (30-40%%)  FCal","P");

  
  //leg->AddEntry(hVn_Zdc_pt_res[ihar][icent],"v_{2} EP Zdc","P");
  setStyle(hVn_FCal_pt_res[ihar][0],0);
  setStyle(hVn_FCal_pt_res[ihar][2],1);
  setStyle(hVn_FCal_pt_res[ihar][4],2);

  setStyle2(gVn_FCal_pt_pub_cent0,0);
  setStyle2(gVn_FCal_pt_pub_cent2,1);
  setStyle2(gVn_FCal_pt_pub_cent4,2);
  setStyle2(gVn_FCal_pt_res[ihar][0],3);
  setStyle2(gVn_FCal_pt_res[ihar][2],4);
  setStyle2(gVn_FCal_pt_res[ihar][4],5);
  
  
  c->cd();
  gVn_FCal_pt_pub_cent0->Draw("AP");
  gVn_FCal_pt_pub_cent0->GetYaxis()->SetRangeUser(0.,0.3);
  gVn_FCal_pt_pub_cent0->GetXaxis()->SetTitle("p_{T}");
  gVn_FCal_pt_pub_cent0->GetYaxis()->SetTitle("v_{2}");
  gVn_FCal_pt_res[ihar][0]->Draw("P");
  gVn_FCal_pt_pub_cent2->Draw("P");
  gVn_FCal_pt_res[ihar][2]->Draw("P");
  gVn_FCal_pt_pub_cent4->Draw("P");
  gVn_FCal_pt_res[ihar][4]->Draw("P");
  
  //hVn_Zdc_pt_res[ihar][icent]->Draw("same");
  leg->Draw("same");
  
  c->SaveAs("v2_pt_comp.pdf");  
  
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

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
const int NCENT=11;
const int NHAR=4;
const int NPT=10;
const int NSIDE=2;
const int NETA=12;
const double pt_mat[]={0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 7.0, 9.0, 12.0};

TProfile *hRes_FCal_fg[NHAR];
TProfile *hRes_Zdc_fg[NHAR];
TProfile *hMean_pt[NCENT][3][NSIDE+1];
TProfile *hVn_FCal_pt_fg[NCENT][NSIDE+1][NHAR];
TProfile *hVn_FCal_eta_fg[NCENT][NSIDE+1][NHAR];
TH1D *hVnc_FCal_pt_fg[NCENT][NSIDE+1][NHAR];
TH1D *hVnc_FCal_eta_fg[NCENT][NSIDE+1][NHAR];

TProfile *hVn_Zdc_pt_fg[NCENT][NSIDE+1][NHAR];
TProfile *hVn_Zdc_eta_fg[NCENT][NSIDE+1][NHAR];
TH1D *hVnc_Zdc_pt_fg[NCENT][NSIDE+1][NHAR];
TH1D *hVnc_Zdc_eta_fg[NCENT][NSIDE+1][NHAR];

TH1D *hVn_FCal_pt_res1,*hVn_Zdc_pt_res1;
TH1D *hVn_FCal_pt_res2,*hVn_Zdc_pt_res2;

TGraphErrors *gVnc_FCal_pt_fg[NCENT][NSIDE+1][NHAR];
TGraphErrors *gVn_FCal_pt_res1,*gVn_Zdc_pt_res2;
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
  h->GetXaxis()->SetTitleFont(43);    h->GetXaxis()->SetTitleSize(17);
  h->GetXaxis()->SetLabelFont(43);    h->GetXaxis()->SetLabelSize(15);
  h->GetYaxis()->SetTitleFont(43);    h->GetYaxis()->SetTitleSize(17);
  h->GetYaxis()->SetLabelFont(43);    h->GetYaxis()->SetLabelSize(15);
  //h->GetZaxis()->SetTitleFont(43);    h->GetZaxis()->SetTitleSize(17);
  //h->GetZaxis()->SetLabelFont(43);    h->GetZaxis()->SetLabelSize(15);
  h->GetYaxis()->SetNdivisions(505);
  h->GetXaxis()->SetNdivisions(505);
  int c_mat[]={1,2,4,6,28,36,7,8,9,46,30};
  int m_mat[]={24,30,27,25,28,33,20,23,29,34};
  h->SetTitle("");

  if(itype1==0 || itype1==1) h->SetMarkerStyle(20);
  if(itype1==2 || itype1==3) h->SetMarkerStyle(21);
  if(itype1==4 || itype1==5) h->SetMarkerStyle(24);
  h->SetMarkerSize(1.);
  h->SetMarkerColor(c_mat[itype1+1]);
  h->SetLineColor(c_mat[itype1+1]);

  //h->GetYaxis()->SetRangeUser(0.,5.);
}

void test2(){
  //TFile *fin=new TFile("./out/out_0_1_0.root","read");
  //TFile *fin=new TFile("./Output_3pv.root","read");
  TFile *fin=new TFile("./Output_SP.root","read");
  for (int icent=0; icent<NCENT; icent++) {
    for (int iside=0; iside<NSIDE+1; iside++) {
      sprintf(name, "hMean_pt_cent%d_ch%d_side%d",icent,2,iside);
       hMean_pt[icent][2][iside]=(TProfile*)fin->Get(name);
      for (int ihar=0; ihar<NHAR; ihar++) {
	sprintf(name, "hVn_FCal_pt_fg_cent%d_side%d_har%d",icent,iside,ihar);
	hVn_FCal_pt_fg[icent][iside][ihar] = (TProfile*)fin->Get(name);
	sprintf(name, "hVn_FCal_eta_fg_cent%d_side%d_har%d",icent,iside,ihar);
	hVn_FCal_eta_fg[icent][iside][ihar] = (TProfile*)fin->Get(name);
	/*
	sprintf(name, "hVnc_FCal_pt_fg_cent%d_side%d_har%d",icent,iside,ihar);
	hVnc_FCal_pt_fg[icent][iside][ihar]=(TProfile*)hVn_FCal_pt_fg[icent][iside][ihar]->Clone(name);
	hVnc_FCal_pt_fg[icent][iside][ihar]->Reset();
	*/

	sprintf(name, "hVn_Zdc_pt_fg_cent%d_side%d_har%d",icent,iside,ihar);
	hVn_Zdc_pt_fg[icent][iside][ihar] = (TProfile*)fin->Get(name);
	sprintf(name, "hVn_Zdc_eta_fg_cent%d_side%d_har%d",icent,iside,ihar);
	hVn_Zdc_eta_fg[icent][iside][ihar] = (TProfile*)fin->Get(name);
	/*
	sprintf(name, "hVnc_Zdc_pt_fg_cent%d_side%d_har%d",icent,iside,ihar);
	hVnc_Zdc_pt_fg[icent][iside][ihar]=(TProfile*)hVn_Zdc_pt_fg[icent][iside][ihar]->Clone(name);
	hVnc_Zdc_pt_fg[icent][iside][ihar]->Reset();
	*/
      }
    }
  }
  //Init TH1D
  for (int icent=0; icent<NCENT; icent++) {
    for (int iside=0; iside<NSIDE+1; iside++) {
      for (int ihar=0; ihar<NHAR; ihar++) {
	sprintf(name, "hVnc_FCal_pt_fg_cent%d_side%d_har%d",icent,iside,ihar);
        hVnc_FCal_pt_fg[icent][iside][ihar]=new TH1D(name,"",NPT,pt_mat);
	hVnc_FCal_pt_fg[icent][iside][ihar]->Sumw2();
	sprintf(name, "hVnc_Zdc_pt_fg_cent%d_side%d_har%d",icent,iside,ihar);
        hVnc_Zdc_pt_fg[icent][iside][ihar]=new	TH1D(name,"",NPT,pt_mat);
	hVnc_Zdc_pt_fg[icent][iside][ihar]->Sumw2();
      }
    }
  }
  
  for (int ihar=0; ihar<NHAR; ihar++) {
    sprintf(name, "hRes_FCal_fg_har%d",ihar);
    hRes_FCal_fg[ihar]= (TProfile*)fin->Get(name);
    sprintf(name, "hRes_Zdc_fg_har%d",ihar);
    hRes_Zdc_fg[ihar]= (TProfile*)fin->Get(name);
  }
  
  double y,y1,y2,yc,dy,dy1,dy2,dyc,yc2;
  double X[100],dX[100],Y[100],dY[100];
  for (int icent=0; icent<NCENT; icent++) {
    for (int iside=0; iside<NSIDE+1; iside++) {
      for (int ihar=0; ihar<NHAR; ihar++) {
	int N=hVn_FCal_pt_fg[icent][iside][ihar]->GetNbinsX();
	for (int ipt=0; ipt<N; ipt++) {
	  y1=hVn_FCal_pt_fg[icent][iside][ihar]->GetBinContent(ipt+1);
	  dy1=hVn_FCal_pt_fg[icent][iside][ihar]->GetBinError(ipt+1);
	  //yc=sqrt(fabs(y))*y/fabs(y);
	  //yc=sqrt(fabs(y));
	  //dyc=0.5*dy/sqrt(fabs(y));
	  y2=hRes_FCal_fg[ihar]->GetBinContent(icent+1);
	  yc=y1/sqrt(fabs(y2));
	  dyc=pow(fabs(y2),-1./2.)*dy1-0.5*y1*dy2*pow(fabs(y2),-3./2.);
	  hVnc_FCal_pt_fg[icent][iside][ihar]->SetBinContent(ipt+1,yc);
	  hVnc_FCal_pt_fg[icent][iside][ihar]->SetBinError(ipt+1,dyc);
	  //yc2=hVnc_FCal_pt_fg[icent][iside][ihar]->GetBinContent(ipt+1);
	  //cout<<"\t"<<y<<"\t"<<yc<<"\t"<<yc2<<endl;
	  X[ipt]=hMean_pt[icent][2][iside]->GetBinContent(ipt+1);
	  Y[ipt]=yc;
	  dY[ipt]=dyc;
	  if(icent==1 && iside==0 && ihar==0) cout<<X[ipt]<<"\t-\t"<<Y[ipt]<<endl;
	}
	gVnc_FCal_pt_fg[icent][iside][ihar]=new TGraphErrors(N,X,Y,0,dY);
      }
    }
  }

  for (int icent=0; icent<NCENT; icent++) {
    for (int iside=0; iside<NSIDE+1; iside++) {
      for (int ihar=0; ihar<NHAR; ihar++) {
	int N=hVn_Zdc_pt_fg[icent][iside][ihar]->GetNbinsX();
	for (int ipt=0; ipt<N; ipt++) {
	  y1=hVn_Zdc_pt_fg[icent][iside][ihar]->GetBinContent(ipt+1);
	  dy1=hVn_Zdc_pt_fg[icent][iside][ihar]->GetBinError(ipt+1);
	  //yc=sqrt(fabs(y))*y/fabs(y);
	  //yc=sqrt(fabs(y));
	  //dyc=0.5*dy/sqrt(fabs(y));
	  y2=hRes_Zdc_fg[ihar]->GetBinContent(icent+1);
	  yc=y1/sqrt(fabs(y2));
	  dyc=pow(fabs(y2),-1./2.)*dy1-0.5*y1*dy2*pow(fabs(y2),-3./2.);
	  hVnc_Zdc_pt_fg[icent][iside][ihar]->SetBinContent(ipt+1,yc);
	  hVnc_Zdc_pt_fg[icent][iside][ihar]->SetBinError(ipt+1,dyc);
	  //yc2=hVnc_Zdc_pt_fg[icent][iside][ihar]->GetBinContent(ipt+1);
	  //cout<<"\t"<<y<<"\t"<<yc<<"\t"<<yc2<<endl;
	}
      }
    }
  }

  //Make ratio histograms
  TH1D *hR[NCENT][NSIDE+1][NHAR];
  for (int icent=0; icent<NCENT; icent++) {
    for (int iside=0; iside<NSIDE+1; iside++) {
      for (int ihar=0; ihar<NHAR-1; ihar++) {
	hR[icent][iside][ihar]=(TH1D*)hVnc_Zdc_pt_fg[icent][iside][ihar+1]->Clone();
	hR[icent][iside][ihar]->Divide(hVnc_FCal_pt_fg[icent][iside][ihar]);
      }
    }
  }

  ifstream inf;
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
  TGraphErrors *gVn_FCal_pt_pub;
  gVn_FCal_pt_pub=new TGraphErrors(25,X,Y,0,dY);
  gVn_FCal_pt_pub->SetMarkerStyle(20);
  gVn_FCal_pt_pub->SetMarkerColor(1);

  int icent=2;
  int iside=0;
  int ihar=0;
  int ieta=NETA;
  TFile *fin2=new TFile("./Output_doAna.root","read");
  sprintf(name,"hV%d_FCal_pt_res_cent%d_side%d_eta%d",ihar+1,icent,iside,ieta);
  hVn_FCal_pt_res1=(TH1D*)fin2->Get(name);
  sprintf(name,"hV%d_Zdc_pt_res_cent%d_side%d_eta%d",ihar+2,icent,iside,ieta);
  hVn_Zdc_pt_res1=(TH1D*)fin2->Get(name);

  sprintf(name,"hV%d_FCal_pt_res_cent%d_side%d_eta%d",ihar+1,icent,NSIDE,ieta);
  hVn_FCal_pt_res2=(TH1D*)fin2->Get(name);
  sprintf(name,"hV%d_Zdc_pt_res_cent%d_side%d_eta%d",ihar+2,icent,NSIDE,ieta);
  hVn_Zdc_pt_res2=(TH1D*)fin2->Get(name);

  int N=hVn_FCal_pt_res2->GetNbinsX();
  for (int ipt=0; ipt<N; ipt++) {
    y=hVn_FCal_pt_res2->GetBinContent(ipt+1);
    dy=hVn_FCal_pt_res2->GetBinError(ipt+1);
    
    X[ipt]=hMean_pt[icent][2][NSIDE]->GetBinContent(ipt+1);
    Y[ipt]=y;
    dY[ipt]=dy;
  }
  gVn_FCal_pt_res2=new TGraphErrors(N,X,Y,0,dY);

  /*
  cout<<"Enter icent :"<<endl;
  cin>>icent;
  cout<<"Enter ihar :"<<endl;
  cin>>ihar;
  */
  setStyle(hVnc_FCal_pt_fg[icent][iside][ihar],0);
  setStyle(hVnc_Zdc_pt_fg[icent][iside][ihar+1],1);
  setStyle(hVn_FCal_pt_res1,2);
  setStyle(hVn_Zdc_pt_res1,3);
  setStyle(hVn_FCal_pt_res2,2);
  
  cout<<"Res_FCal="<<sqrt(fabs(hRes_FCal_fg[ihar]->GetBinContent(icent+1)))<<endl;
  cout<<"Res_Zdc="<<sqrt(fabs(hRes_Zdc_fg[ihar+1]->GetBinContent(icent+1)))<<endl;

  gStyle->SetOptStat(0);
  TCanvas *c;
  c=new TCanvas("c","c",1000,400);
  c->Divide(2);
  TLegend *leg;
  leg=new TLegend(0.15,0.6,0.45,0.88);
  leg->SetLineColor(0);
  
  leg->AddEntry(hVnc_FCal_pt_fg[icent][iside][ihar],"v_{2} SP FCal","P");
  //leg->AddEntry(hVnc_Zdc_pt_fg[icent][iside][ihar+1],"v_{2} SP Zdc","P");
  leg->AddEntry(hVn_FCal_pt_res1,"v_{2} EP FCal","P");
  //leg->AddEntry(hVn_Zdc_pt_res1,"v_{2} EP Zdc","P");
  
  c->cd(1);
  hVnc_FCal_pt_fg[icent][iside][ihar]->Draw();
  //hVnc_FCal_pt_fg[icent][iside][ihar]->GetXaxis()->SetRangeUser(0.5,5.);
  hVnc_FCal_pt_fg[icent][iside][ihar]->GetYaxis()->SetRangeUser(0.,0.2);
  //hVnc_Zdc_pt_fg[icent][iside][ihar+1]->Draw("same");
  hVn_FCal_pt_res1->Draw("same");
  //hVn_Zdc_pt_res1->Draw("same");
  //hVn_Zdc_pt_fg[icent][iside][ihar]->GetYaxis()->SetRangeUser(0.,0.5);
  leg->Draw("same");
  c->cd(2);
  hR[icent][iside][ihar]->Draw();
  hR[icent][iside][ihar]->GetXaxis()->SetRangeUser(0.5,5.);
  hR[icent][iside][ihar]->GetYaxis()->SetRangeUser(0.,1.5);

  TCanvas *c2;
  c2=new TCanvas("c2","c2",800,600);
  TLegend *leg2;
  leg2=new TLegend(0.2,0.2,0.45,0.5);
  leg2->SetLineColor(0);
  leg2->AddEntry(gVn_FCal_pt_pub,"v_{2} EP FCal Pub","P");
  leg2->AddEntry(gVn_FCal_pt_res2,"v_{2} EP FCal","P");
  //leg2->AddEntry(gVnc_FCal_pt_fg[icent][iside][ihar],"v_{2} SP FCal","P");
  gVn_FCal_pt_res2->Draw("AP");
  gVn_FCal_pt_res2->SetTitle("");
  gVn_FCal_pt_res2->GetYaxis()->SetRangeUser(0.,0.2);
  gVn_FCal_pt_res2->GetXaxis()->SetTitle("p_{T}");
  gVn_FCal_pt_res2->GetYaxis()->SetTitle("v_{2}");
  gVn_FCal_pt_res2->SetMarkerStyle(24);
  gVn_FCal_pt_res2->SetMarkerColor(2);
  //gVnc_FCal_pt_fg[icent][iside][ihar]->Draw("P");
  gVnc_FCal_pt_fg[icent][iside][ihar]->SetMarkerStyle(24);
  gVnc_FCal_pt_fg[icent][iside][ihar]->SetMarkerColor(4);
  gVn_FCal_pt_pub->Draw("P");
  leg2->Draw("same");

}

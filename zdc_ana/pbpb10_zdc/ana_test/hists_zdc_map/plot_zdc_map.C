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

enum{
  MAXTRK=100000,
  NDET=16,
  NHAR=4,
  NSIDE=2,
  NMODULE=2,
  NGAIN=2,
};

char name[200];

void Divide_Pad(TCanvas *c,TPad *p[],int nrow,int ncol){
  int N=nrow*ncol;
  double xmin,ymin,xmax,ymax;
  xmin=0.0,ymin=0.0;
  xmax=1.-0.0,ymax=1.-0.0;
  double dx,dy,xw,yw,xw1,yw1,xwn,ywn;
  dx=0.1,dy=0.1;
  xw=(1.-xmin)/ncol;
  yw=(1.-ymin)/nrow;
  double MX_min[10][10],MX_max[10][10];
  double MY_min[10][10],MY_max[10][10];
  for(int icol=0;icol<ncol;icol++){
    for(int irow=nrow-1;irow>=0;irow--){
      xw1=xmax/ncol+0.025;
      yw1=ymax/nrow+0.034;
      xw=(xmax-xw1)/(ncol-1);
      yw=(ymax-yw1)/(nrow-1);
      xwn=xw;//xmax/ncol+0.01;
      ywn=yw;//ymax/nrow+0.01;
      xwn=xw+0.005;
      ywn=yw+0.015;
      xw=(xmax-xw1-xwn)/(ncol-2);
      yw=(ymax-yw1-ywn)/(nrow-2);

      if(icol==0){
        MX_min[irow][icol]=xmin;
        MX_max[irow][icol]=xmin+xw1;
      }
      if(icol>0 && icol!=ncol-1){
        MX_min[irow][icol]=xmin+xw1+(icol-1)*xw;
        MX_max[irow][icol]=xmin+xw1+(icol)*xw;
      }
      if(icol==ncol-1){
        MX_min[irow][icol]=MX_max[irow][icol-1];//xmin+xw1+(icol-1)*xw;
        MX_max[irow][icol]=xmax;
      }

      if(irow==nrow-1){
        MY_min[irow][icol]=ymin;
        MY_max[irow][icol]=ymin+yw1;
      }
      if(irow!=nrow-1 && irow!=0){
        MY_min[irow][icol]=ymin+yw1+(nrow-2-irow)*yw;
        MY_max[irow][icol]=ymin+yw1+(nrow-2-irow+1)*yw;
      }
      if(irow==0){
        MY_min[irow][icol]=MY_max[irow+1][icol];//ymin+yw1+(nrow-2-irow)*yw;
        MY_max[irow][icol]=ymax;
      }
    }

    int cnt=0;
    for(int irow=0;irow<nrow;irow++){
      for(int icol=0;icol<ncol;icol++){
        sprintf(name,"p_%d",cnt);
        //cout<<MX_min[irow][icol]<<",\t"<<MX_max[irow][icol]<<",\t"<<MY_min[irow][icol]<<",\t"<<MY_max[irow][icol]<<endl;

        p[cnt]=new TPad(name,name,MX_min[irow][icol],MY_min[irow][icol],MX_max[irow][icol],MY_max[irow][icol]);
        c->cd();
        p[cnt]->SetRightMargin(0);
        p[cnt]->SetLeftMargin(0);
        p[cnt]->SetTopMargin(0);
        p[cnt]->SetBottomMargin(0);
        if(icol==0) p[cnt]->SetLeftMargin(0.17);
        if(icol==ncol-1) p[cnt]->SetRightMargin(0.05);
        if(irow==0)       p[cnt]->SetTopMargin(0.1);
        if(irow==nrow-1) p[cnt]->SetBottomMargin(0.2);
        p[cnt]->SetTicks(1,1);
        p[cnt]->Draw();
        cnt++;
      }
    }
  }

}

void setStyle1(TH1 *h,int itype1){
  h->GetYaxis()->CenterTitle();       h->GetXaxis()->CenterTitle();  //h->GetZaxis()->CenterTitle();
  h->GetYaxis()->SetTitleOffset(2.5); h->GetXaxis()->SetTitleOffset(2.5);
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

  h->SetMarkerStyle(24);
  h->SetMarkerSize(1.);
  h->SetMarkerColor(c_mat[itype1]);
  h->SetLineColor(c_mat[itype1]);

  //h->GetYaxis()->SetRangeUser(0.,5.);
}

void plot_zdc_map(){
  TH2F *hEta_Layer;
  TH2F *hPixel_Energy[NSIDE+1][20];

  TFile *fin;
  fin=new TFile("./hists_zdc.root","read");

  hEta_Layer=(TH2F*)fin->Get("hEta_Layer");
  for(int iside=0;iside<NSIDE+1;iside++){
    for(int icent=0;icent<20;icent++){
      sprintf(name,"hPixel_Energy_side%d_cent%d",iside,icent);
      hPixel_Energy[iside][icent]=(TH2F*)fin->Get(name);
      setStyle1(hPixel_Energy[iside][icent],2);
    }
  }


  //Plot
  gStyle->SetOptStat(0);
  TLatex text;
  text.SetNDC(1);
  text.SetTextFont(43);
  text.SetTextSize(18);
  const int NCan1=3;
  TCanvas *c[NCan1];
  TPad *pad[NCan1][6];
  TLegend *leg[NCan1];
  for(int ican=0;ican<NCan1;ican++){
    sprintf(name,"c_%d",ican);
      c[ican]=new TCanvas(name,name,1200,800);
      leg[ican]=new TLegend(0.45,0.55,0.88,0.88);
      leg[ican]->SetLineColor(0);
  }

  int icent_mat[]={0,4,6,8,10,15};
  int cent_mat[]={0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100};
  for(int ican=0;ican<NCan1;ican++){
    int iside=ican;
    Divide_Pad(c[ican],pad[ican],2,3);
    for(int i=0;i<6;i++){
      int icent=icent_mat[i];
      c[ican]->cd();
      pad[ican][i]->cd();
      hPixel_Energy[iside][icent]->Draw("colz");
      if(i==5) hPixel_Energy[iside][icent]->SetXTitle("Grid Point : X-axis");
      if(i==0) hPixel_Energy[iside][icent]->SetYTitle("Grid Point : Y-axis");
      hPixel_Energy[iside][icent]->GetXaxis()->SetRangeUser(-1.,9.);
      hPixel_Energy[iside][icent]->GetYaxis()->SetRangeUser(-1.,12.);
      
      sprintf(name,"%d-%d%%",cent_mat[icent],cent_mat[icent+1]);
      text.DrawLatex(0.2, 0.2,name);
    }
    
  }
  
}

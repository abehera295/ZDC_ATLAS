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
const int NCENT=10;
const int NHAR=4;
TH1D *hRes_flat_FCal[NCENT][NHAR];
TH1D *hResCent_flat_FCal[NHAR];
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
int find_Min(TH1D *h){
  int N=h->GetNbinsX();
  int imin=0;
  double y,ymin;
  for(int i=0;i<N;i++){
    y=h->GetBinContent(i+1);
    if(i==0) ymin=y;
    if(y<=ymin) ymin=y;
  }
  for(int i=0;i<N;i++){
    y=h->GetBinContent(i+1);
    if(y==ymin) imin=i;
  }
  return imin;
}

double cal_Res(double x){
  double res=x*sqrt(PI)/2.;
  res*=exp(-x*x/2.);
  //res*=TMath::Bessell0(x*x/2.)+TMath::Bessell1(x*x/2.);
  res*=TMath::BesselI0(x*x/2.)+TMath::BesselI1(x*x/2.);
  return res;
}

void test(){
  TFile *fin=new TFile("./Output_vn.root","read");
  for (int icent=0; icent<NCENT; icent++) {
    for (int ihar=0; ihar<NHAR; ihar++) {
      sprintf(name, "hRes_flat_FCal_cent%d_har%d", icent, ihar);
      hRes_flat_FCal[icent][ihar] = (TH1D*)fin->Get(name);
    }
  }

  for (int ihar=0; ihar<NHAR; ihar++) {
    sprintf(name,"hResCent_flat_FCal_har%d",ihar);
    hResCent_flat_FCal[ihar]=new TH1D (name,"",NCENT,0.-0.5,NCENT-0.5);
    for (int icent=0; icent<NCENT; icent++) {
      //For FCal
      //FCal has v2,v3,v4,v5
      hResCent_flat_FCal[ihar]->SetBinContent(icent+1, hRes_flat_FCal[icent][ihar]->GetMean());
      hResCent_flat_FCal[ihar]->SetBinError(icent+1, hRes_flat_FCal[icent][ihar]->GetMeanError());  
    }
  }

  for (int ihar=0; ihar<NHAR; ihar++) {
    for (int icent=0; icent<NCENT; icent++) {
      double y1=hResCent_flat_FCal[ihar]->GetBinContent(icent+1);
      double dy1=hResCent_flat_FCal[ihar]->GetBinError(icent+1);
      double y2=sqrt(fabs(y1));
      double dy2=0.5*dy1/sqrt(fabs(y1));
      hResCent_flat_FCal[ihar]->SetBinContent(icent+1,y2);
      hResCent_flat_FCal[ihar]->SetBinError(icent+1,dy2);
    }
  }
  
  TH1D *hX[NHAR],*hDRes[NCENT];
  for (int ihar=0; ihar<NHAR; ihar++) {
    sprintf(name,"hX_har%d",ihar);
    hX[ihar]=new TH1D(name,"",NCENT,0.-0.5,NCENT-0.5);
    hX[ihar]->Sumw2();
  }
  for (int icent=0; icent<NCENT; icent++) {
    sprintf(name,"hDRes_cent%d",icent);
    hDRes[icent]=new TH1D(name,"",10000,0.-0.5,10000.-0.5);
  }
  double X[10000],DRes[10000];
  double dx=10./10000.;
  int ihar=0;
  for (int icent=0; icent<NCENT; icent++) {
    double y1=hResCent_flat_FCal[ihar]->GetBinContent(icent+1);
    double dy1=hResCent_flat_FCal[ihar]->GetBinError(icent+1);
    for(int i=0;i<10000;i++){
      X[i]=0.+i*dx;
      //double Res=cal_Res(X[i]);
      double Res=(X[i]*sqrt(PI)/2.) * (exp(-X[i]*X[i]/2.)) * (TMath::BesselI0(X[i]*X[i]/2.)+TMath::BesselI1(X[i]*X[i]/2.)) ;
      hDRes[icent]->SetBinContent(i+1,fabs(y1-Res));
    }
    int idres_min=hDRes[icent]->GetMinimumBin();
    hX[ihar]->SetBinContent(icent+1,X[idres_min]);
    hX[ihar]->SetBinError(icent+1,dy1);
  }

  TCanvas *c;
  c=new TCanvas("c","c",1200,500);
  c->Divide(2);
  c->cd(1);
  hResCent_flat_FCal[ihar]->Draw();
  c->cd(2);
  hX[ihar]->Draw();

  TCanvas *c2;
  c2=new TCanvas("c2","c2",1200,900);
  c2->Divide(3,3);
  for(int i=0;i<9;i++){
    c2->cd(i+1);
    hDRes[i]->Draw();
  }

}

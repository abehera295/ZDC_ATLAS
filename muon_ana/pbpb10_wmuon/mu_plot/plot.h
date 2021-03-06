#ifndef plot_h
#define plot_h

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

static double PI = acos(-1.0);
static double PI2 = 2.*acos(-1.0);
enum {
  NDET=16,
  NSIDE=2,
  NMODULE=2,
  NGAIN=2,
  NCH1=4,
  NCH2=24,
  NHAR=4,
  //NCENT=20,
  NCENT=7,
  NCENT_1p=100,
  //NK=12,
  NSAMP=10,
  //STORE_DEP=20,
  STORE_DEP=5,
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
  //NK=4,
  NK=12,
  NZ=8,//25mm bin
};

enum{
  NPT_TrkWei=9,
  NCENT_TrkWei=14,
  NETA_TrkWei=50,
  NPHI_TrkWei=64,
  MPT=5,
};

char name[200];
int cent_mat[8]={0,10,20,30,40,60,80,100};
double cent_mat_double[8]={0,10,20,30,40,60,80,100};
const int NPT_Mu=9;
const double pt_mat_Mu[10]={4.,4.5,5.,5.5,6.,7.,8.,10.,12.,14.};
//const int NETA_Mu=8;
//const double eta_mat_Mu[9]={-2.,-1.5,-1.,-0.5,0.,0.5,1.,1.5,2.};
const int NETA_Mu=4;
const double eta_mat_Mu[5]={-2.,-1.,0.,1.,2.};

//Muon hists
TH1D *hMu_n_raw,*hMu_pt_raw,*hMu_phi_raw,*hMu_eta_raw,*hMu_eloss_all_raw;
TH1D *hMu_ms_phi_raw,*hMu_ms_theta_raw,*hMu_ms_qoverp_raw;
TH1D *hMu_n,*hMu_pt,*hMu_phi,*hMu_eta,*hMu_eloss_all;
TH1D *hMu_ms_phi,*hMu_ms_theta,*hMu_ms_qoverp;
TH1D *hMu_dphi_yield_FCal[NCENT][NSIDE+1][NPT_Mu][NHAR];
TH2F *hMu_dphi_eloss_FCal[NCENT][NSIDE+1][NPT_Mu][NHAR];
TH1D *hMu_dphi_yield_Zdc[NCENT][NSIDE+1][NETA_Mu][NHAR];
TH2F *hMu_dphi_eloss_Zdc[NCENT][NSIDE+1][NETA_Mu][NHAR];

class plot{
 public:
  //Functions
  plot();
  void init();
  void Divide_Pad(TCanvas *c,TPad *p[],int nrow,int ncol);
  void plot1_mu();
  void plot1_recombine();
  void plot2_fit();
  void plot2_readHists();
  void plot2_doFit();
  void plot2_drawFit();
  
  //Variables
  TFile *fin,*fin2;
  TLatex text;
};
#endif


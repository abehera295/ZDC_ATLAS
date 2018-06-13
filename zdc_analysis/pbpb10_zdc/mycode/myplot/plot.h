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
  NSIDE   = 2,
  NMODULE = 2,
  NGAIN   = 2,
  NCH1    = 4,
  NCH2    = 24,
  NHAR    = 3,
  NCENT   = 20,
  NCENT_1p  = 100,
  NK      = 12,
  NSAMP   = 10,
  STORE_DEP = 20
};

enum{
  NPT=10,
  NETA=12,
  NCENT_NEW=37
};

const double pt_mat[]={0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 7.0, 9.0, 12.0};
const double eta_mat[]={-2.5, -2.0, -1.5, -1.0, -0.5, -0.25, 0.0, 0.25, 0.5, 1.0, 1.5, 2.0, 2.5};
string Cent_mat[]={"0-5%","5-10%","10-15%","15-20%","20-25%","25-30%","30-35%","35-40%","40-45%","45-50%","50-55%","55-60%","60-65%","65-70%","70-75%","75-80%","80-85%","85-90%","90-95%","95-100%",};
string Cent_mat_new[]={"0-5%","5-10%","10-15%","15-20%","20-25%","25-30%","30-35%","35-40%","40-45%","45-50%","50-55%","55-60%","60-65%","65-70%","70-75%","75-80%","80-85%","85-90%","90-95%","95-100%","0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%","0-20%","20-40%","40-60%","60-80%","80-100%","0-40%","40-80%",};

char name[20];

TH1D *hPsi[NCENT][NSIDE+1][NHAR];
TH1D *hPsiC[NCENT][NSIDE+1][NHAR];
TH1D *hPsiF[NCENT][NSIDE+1][NHAR];

TH1D *hRes      [NCENT][NHAR];
TH1D *hRes_recn [NCENT][NHAR];
TH1D *hRes_flat [NCENT][NHAR];

TH1D *hResCent      [NHAR];
TH1D *hResCent_recn [NHAR];
TH1D *hResCent_flat [NHAR];
TH1D *hResCent_mix  [NHAR];

TH1D *hDphi_fg [NCENT][NHAR];
TH1D *hDphi_bg [NCENT][NHAR];
TH1D *hDphi_rat[NCENT][NHAR];
/*
std::vector<std::vector<TH1D*> > hVn         [NCENT][NSIDE+1][NHAR];
std::vector<TH1D*>               hVn_eta_raw [NCENT][NSIDE+1][NHAR];
std::vector<TH1D*>               hVn_eta_res [NCENT][NSIDE+1][NHAR];
std::vector<TH1D*>               hVn_pt_raw  [NCENT][NSIDE+1][NHAR];
std::vector<TH1D*>               hVn_pt_res  [NCENT][NSIDE+1][NHAR];
std::vector<std::vector<TH1D*> > hVn_cent_raw[NSIDE+1][NHAR];
std::vector<std::vector<TH1D*> > hVn_cent_res[NSIDE+1][NHAR];

std::vector<std::vector<TH1D*> > hVnS         [NCENT][NSIDE+1][NHAR];
std::vector<TH1D*>               hVnS_eta_raw [NCENT][NSIDE+1][NHAR];
std::vector<TH1D*>               hVnS_eta_res [NCENT][NSIDE+1][NHAR];
std::vector<TH1D*>               hVnS_pt_raw  [NCENT][NSIDE+1][NHAR];
std::vector<TH1D*>               hVnS_pt_res  [NCENT][NSIDE+1][NHAR];
std::vector<std::vector<TH1D*> > hVnS_cent_raw[NSIDE+1][NHAR];
std::vector<std::vector<TH1D*> > hVnS_cent_res[NSIDE+1][NHAR];
*/
TH1D *hVn_eta_raw[NCENT][NSIDE+1][NHAR][NPT],*hVnS_eta_raw[NCENT][NSIDE+1][NHAR][NPT];
TH1D *hVn_eta_res[NCENT][NSIDE+1][NHAR][NPT],*hVnS_eta_res[NCENT][NSIDE+1][NHAR][NPT];
TH1D *hVn_pt_res[NCENT][NSIDE+1][NHAR][NETA],*hVnS_pt_res[NCENT][NSIDE+1][NHAR][NETA];

TH1D *hV1_eta_res[NCENT][NSIDE+1];
TProfile *h_prof_sumx[NSIDE+1][NHAR];
TProfile *h_prof_sumy[NSIDE+1][NHAR];
TH2D *h_Q2D [NCENT][NSIDE+1][NHAR];

TProfile *h_prof_flatCos[NCENT][NSIDE+1][NHAR];
TProfile *h_prof_flatSin[NCENT][NSIDE+1][NHAR];

class plot{
 public:
  //Functions
  plot();
  void init();
  void setStyle1(TH1D *h,int i);
  void setStyle2(TH1D *h,int i);
  void setStyle3(TH1D *h,int i);
  double FindAvg(TH1D *h);
  void plot1_Psi();
  void plot2_Res();
  void plot3_v1_eta();

  //Variables
  TFile *fin,*fin2;
  TLatex text;
};
#endif

plot::plot(){

}

double plot::FindAvg(TH1D *h){
  int N=h->GetNbinsX();
  double sum=0.;
  for(int i=0;i<N;i++){
    sum+=h->GetBinContent(i+1);
  }
  return sum/double(N);
}

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
  NCENT=11,
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

int cent_mat[]={0,5,10,20,30,40,50,60,70,80,90,100};
double cent_mat_double[]={0,5,10,20,30,40,50,60,70,80,90,100};
const double pt_mat[]={0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 7.0, 9.0, 12.0};
const double eta_mat[]={-2.5, -2.0, -1.5, -1.0, -0.5, -0.25, 0.0, 0.25, 0.5, 1.0, 1.5, 2.0, 2.5};
string Cent_mat[]={"0-5%","5-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%",};
//string Cent_mat[]={"0-5%","5-10%","10-15%","15-20%","20-25%","25-30%","30-35%","35-40%","40-45%","45-50%","50-55%","55-60%","60-65%","65-70%","70-75%","75-80%","80-85%","85-90%","90-95%","95-100%",};
string Cent_mat_new[]={"0-5%","5-10%","10-15%","15-20%","20-25%","25-30%","30-35%","35-40%","40-45%","45-50%","50-55%","55-60%","60-65%","65-70%","70-75%","75-80%","80-85%","85-90%","90-95%","95-100%","0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%","0-20%","20-40%","40-60%","60-80%","80-100%","0-40%","40-80%",};

char name[20];

//EP Correction
TH1D* hPsiZ_raw[NCENT][NSIDE+1][NHAR];
TH1D* hPsiZ_rec[NCENT][NSIDE+1][NHAR];
TH1D* hPsiZ_flat[NCENT][NSIDE+1][NHAR];


//New
//For FCal
TH1D *hPsi_FCal[NCENT][NSIDE+1][NHAR];
TH1D *hRes_flat_FCal[NCENT][NHAR];
TH1D *hResCent_flat_FCal[NHAR];
TH1D *hResCent_bg_FCal[NHAR];
TH1D *hResCent_flat_FCal_sub[NHAR];
TH1D *hResCent_flat_FCal_ful[NHAR];
TH1D *hResXCent_flat_FCal_sub[NHAR];
TH1D *hResXCent_flat_FCal_ful[NHAR];
TH1D *hResCent_mix_FCal[NHAR];
TH1D *hResCent_mix_FCal_sub[NHAR];
TH1D *hResCent_mix_FCal_ful[NHAR];
TH1D *hResXCent_mix_FCal_sub[NHAR];
TH1D *hResXCent_mix_FCal_ful[NHAR];
TH1D *hDphi_FCal_fg[NCENT][NHAR];
TH1D *hDphi_FCal_bg[NCENT][NHAR];
TH1D *hDphi_FCal_rat[NCENT][NHAR];
TProfile *hRes_FCal_cos_fg[NCENT],*hRes_FCal_cos_bg[NCENT];
TProfile *hRes_FCal_sin_fg[NCENT],*hRes_FCal_sin_bg[NCENT];
TProfile *hResCent_FCal_fg[NHAR],*hResCent_FCal_bg[NHAR];

TProfile2D *hVn_FCal[NCENT][NSIDE+1][NHAR],*hVnS_FCal[NCENT][NSIDE+1][NHAR];
TProfile2D *hVn_FCal_fg[NCENT][NSIDE+1][NHAR],*hVnS_FCal_fg[NCENT][NSIDE+1][NHAR];
TProfile2D *hVn_FCal_bg[NCENT][NSIDE+1][NHAR],*hVnS_FCal_bg[NCENT][NSIDE+1][NHAR];

TH1D *hVn_FCal_eta_raw[NCENT][NSIDE+1][NHAR][NPT+1],*hVn_FCal_eta_raw_fg[NCENT][NSIDE+1][NHAR][NPT+1],*hVn_FCal_eta_raw_bg[NCENT][NSIDE+1][NHAR][NPT+1];
TH1D *hVn_FCal_eta_res[NCENT][NSIDE+1][NHAR][NPT+1],*hVn_FCal_eta_res_fg[NCENT][NSIDE+1][NHAR][NPT+1],*hVn_FCal_eta_res_bg[NCENT][NSIDE+1][NHAR][NPT+1];
TH1D *hVn_FCal_pt_raw[NCENT][NSIDE+1][NHAR][NETA+4],*hVn_FCal_pt_raw_fg[NCENT][NSIDE+1][NHAR][NETA+4],*hVn_FCal_pt_raw_bg[NCENT][NSIDE+1][NHAR][NETA+4];
TH1D *hVn_FCal_pt_res[NCENT][NSIDE+1][NHAR][NETA+4],*hVn_FCal_pt_res_fg[NCENT][NSIDE+1][NHAR][NETA+4],*hVn_FCal_pt_res_bg[NCENT][NSIDE+1][NHAR][NETA+4];
TH1D *hVn_FCal_cent_raw[NSIDE+1][NHAR][NPT+1][NETA+4],*hVn_FCal_cent_raw_fg[NSIDE+1][NHAR][NPT+1][NETA+4],*hVn_FCal_cent_raw_bg[NSIDE+1][NHAR][NPT+1][NETA+4];
TH1D *hVn_FCal_cent_res[NSIDE+1][NHAR][NPT+1][NETA+4],*hVn_FCal_cent_res_fg[NSIDE+1][NHAR][NPT+1][NETA+4],*hVn_FCal_cent_res_bg[NSIDE+1][NHAR][NPT+1][NETA+4];

TH1D *hVnS_FCal_eta_raw[NCENT][NSIDE+1][NHAR][NPT+1],*hVnS_FCal_eta_raw_fg[NCENT][NSIDE+1][NHAR][NPT+1],*hVnS_FCal_eta_raw_bg[NCENT][NSIDE+1][NHAR][NPT+1];
TH1D *hVnS_FCal_eta_res[NCENT][NSIDE+1][NHAR][NPT+1],*hVnS_FCal_eta_res_fg[NCENT][NSIDE+1][NHAR][NPT+1],*hVnS_FCal_eta_res_bg[NCENT][NSIDE+1][NHAR][NPT+1];
TH1D *hVnS_FCal_pt_raw[NCENT][NSIDE+1][NHAR][NETA+4],*hVnS_FCal_pt_raw_fg[NCENT][NSIDE+1][NHAR][NETA+4],*hVnS_FCal_pt_raw_bg[NCENT][NSIDE+1][NHAR][NETA+4];
TH1D *hVnS_FCal_pt_res[NCENT][NSIDE+1][NHAR][NETA+4],*hVnS_FCal_pt_res_fg[NCENT][NSIDE+1][NHAR][NETA+4],*hVnS_FCal_pt_res_bg[NCENT][NSIDE+1][NHAR][NETA+4];
TH1D *hVnS_FCal_cent_raw[NSIDE+1][NHAR][NPT+1][NETA+4],*hVnS_FCal_cent_raw_fg[NSIDE+1][NHAR][NPT+1][NETA+4],*hVnS_FCal_cent_raw_bg[NSIDE+1][NHAR][NPT+1][NETA+4];
TH1D *hVnS_FCal_cent_res[NSIDE+1][NHAR][NPT+1][NETA+4],*hVnS_FCal_cent_res_fg[NSIDE+1][NHAR][NPT+1][NETA+4],*hVnS_FCal_cent_res_bg[NSIDE+1][NHAR][NPT+1][NETA+4];

//In SP  
TProfile *hRes_FCal_fg[NHAR],*hRes_FCal_bg[NHAR];
TProfile *hVn_FCal_cent_fg[NSIDE+1][NHAR],*hVn_FCal_cent_bg[NSIDE+1][NHAR];
TProfile *hVn_FCal_pt_fg[NCENT][NSIDE+1][NHAR],*hVn_FCal_pt_bg[NCENT][NSIDE+1][NHAR];
TProfile *hVn_FCal_eta_fg[NCENT][NSIDE+1][NHAR],*hVn_FCal_eta_bg[NCENT][NSIDE+1][NHAR];
TProfile2D *hC112_FCal_pt_same_fg[NCENT],*hC112_FCal_pt_same_bg[NCENT];
TProfile2D *hC112_FCal_eta_same_fg[NCENT],*hC112_FCal_eta_same_bg[NCENT];
TProfile2D *hC112_FCal_pt_opp_fg[NCENT],*hC112_FCal_pt_opp_bg[NCENT];
TProfile2D *hC112_FCal_eta_opp_fg[NCENT],*hC112_FCal_eta_opp_bg[NCENT];

//For Zdc
TH1D *hPsi_Zdc[NCENT][NSIDE+1][NHAR];
TH1D *hRes_flat_Zdc[NCENT][NHAR];
TH1D *hResCent_flat_Zdc[NHAR];
TH1D *hResCent_bg_Zdc[NHAR];
TH1D *hResCent_flat_Zdc_sub[NHAR];
TH1D *hResCent_flat_Zdc_ful[NHAR];
TH1D *hResXCent_flat_Zdc_sub[NHAR];
TH1D *hResXCent_flat_Zdc_ful[NHAR];
TH1D *hResCent_mix_Zdc[NHAR];
TH1D *hResCent_mix_Zdc_sub[NHAR];
TH1D *hResCent_mix_Zdc_ful[NHAR];
TH1D *hResXCent_mix_Zdc_sub[NHAR];
TH1D *hResXCent_mix_Zdc_ful[NHAR];
TH1D *hDphi_Zdc_fg[NCENT][NHAR];
TH1D *hDphi_Zdc_bg[NCENT][NHAR];
TH1D *hDphi_Zdc_rat[NCENT][NHAR];
TProfile *hRes_Zdc_cos_fg[NCENT],*hRes_Zdc_cos_bg[NCENT];
TProfile *hRes_Zdc_sin_fg[NCENT],*hRes_Zdc_sin_bg[NCENT];
TProfile *hResCent_Zdc_fg[NHAR],*hResCent_Zdc_bg[NHAR];
  
TProfile2D *hVn_Zdc[NCENT][NSIDE+1][NHAR],*hVnS_Zdc[NCENT][NSIDE+1][NHAR];
TProfile2D *hVn_Zdc_fg[NCENT][NSIDE+1][NHAR],*hVnS_Zdc_fg[NCENT][NSIDE+1][NHAR];
TProfile2D *hVn_Zdc_bg[NCENT][NSIDE+1][NHAR],*hVnS_Zdc_bg[NCENT][NSIDE+1][NHAR];

TH1D *hVn_Zdc_eta_raw[NCENT][NSIDE+1][NHAR][NPT+1],*hVn_Zdc_eta_raw_fg[NCENT][NSIDE+1][NHAR][NPT+1],*hVn_Zdc_eta_raw_bg[NCENT][NSIDE+1][NHAR][NPT+1];
TH1D *hVn_Zdc_eta_res[NCENT][NSIDE+1][NHAR][NPT+1],*hVn_Zdc_eta_res_fg[NCENT][NSIDE+1][NHAR][NPT+1],*hVn_Zdc_eta_res_bg[NCENT][NSIDE+1][NHAR][NPT+1];
TH1D *hVn_Zdc_pt_raw[NCENT][NSIDE+1][NHAR][NETA+4],*hVn_Zdc_pt_raw_fg[NCENT][NSIDE+1][NHAR][NETA+4],*hVn_Zdc_pt_raw_bg[NCENT][NSIDE+1][NHAR][NETA+4];
TH1D *hVn_Zdc_pt_res[NCENT][NSIDE+1][NHAR][NETA+4],*hVn_Zdc_pt_res_fg[NCENT][NSIDE+1][NHAR][NETA+4],*hVn_Zdc_pt_res_bg[NCENT][NSIDE+1][NHAR][NETA+4];
TH1D *hVn_Zdc_cent_raw[NSIDE+1][NHAR][NPT+1][NETA+4],*hVn_Zdc_cent_raw_fg[NSIDE+1][NHAR][NPT+1][NETA+4],*hVn_Zdc_cent_raw_bg[NSIDE+1][NHAR][NPT+1][NETA+4];
TH1D *hVn_Zdc_cent_res[NSIDE+1][NHAR][NPT+1][NETA+4],*hVn_Zdc_cent_res_fg[NSIDE+1][NHAR][NPT+1][NETA+4],*hVn_Zdc_cent_res_bg[NSIDE+1][NHAR][NPT+1][NETA+4];

TH1D *hVnS_Zdc_eta_raw[NCENT][NSIDE+1][NHAR][NPT+1],*hVnS_Zdc_eta_raw_fg[NCENT][NSIDE+1][NHAR][NPT+1],*hVnS_Zdc_eta_raw_bg[NCENT][NSIDE+1][NHAR][NPT+1];
TH1D *hVnS_Zdc_eta_res[NCENT][NSIDE+1][NHAR][NPT+1],*hVnS_Zdc_eta_res_fg[NCENT][NSIDE+1][NHAR][NPT+1],*hVnS_Zdc_eta_res_bg[NCENT][NSIDE+1][NHAR][NPT+1];
TH1D *hVnS_Zdc_pt_raw[NCENT][NSIDE+1][NHAR][NETA+4],*hVnS_Zdc_pt_raw_fg[NCENT][NSIDE+1][NHAR][NETA+4],*hVnS_Zdc_pt_raw_bg[NCENT][NSIDE+1][NHAR][NETA+4];
TH1D *hVnS_Zdc_pt_res[NCENT][NSIDE+1][NHAR][NETA+4],*hVnS_Zdc_pt_res_fg[NCENT][NSIDE+1][NHAR][NETA+4],*hVnS_Zdc_pt_res_bg[NCENT][NSIDE+1][NHAR][NETA+4];
TH1D *hVnS_Zdc_cent_raw[NSIDE+1][NHAR][NPT+1][NETA+4],*hVnS_Zdc_cent_raw_fg[NSIDE+1][NHAR][NPT+1][NETA+4],*hVnS_Zdc_cent_raw_bg[NSIDE+1][NHAR][NPT+1][NETA+4];
TH1D *hVnS_Zdc_cent_res[NSIDE+1][NHAR][NPT+1][NETA+4],*hVnS_Zdc_cent_res_fg[NSIDE+1][NHAR][NPT+1][NETA+4],*hVnS_Zdc_cent_res_bg[NSIDE+1][NHAR][NPT+1][NETA+4];

//For 10-80%
TH1D *hDphi_Zdc_allCent_fg[NHAR];
TH1D *hDphi_Zdc_allCent_bg[NHAR];
TH1D *hDphi_Zdc_allCent_rat[NHAR];
TH1D *hRes_mix_Zdc_allCent,*hRes_flat_Zdc_allCent,*hRes_bg_Zdc_allCent;
TProfile *hRes_Zdc_cos_allCent_fg,*hRes_Zdc_cos_allCent_bg;
TProfile *hRes_Zdc_sin_allCent_fg,*hRes_Zdc_sin_allCent_bg;
  
TProfile2D *hVn_Zdc_allCent[NSIDE+1][NHAR],*hVnS_Zdc_allCent[NSIDE+1][NHAR];
TProfile2D *hVn_Zdc_allCent_fg[NSIDE+1][NHAR],*hVnS_Zdc_allCent_fg[NSIDE+1][NHAR];
TProfile2D *hVn_Zdc_allCent_bg[NSIDE+1][NHAR],*hVnS_Zdc_allCent_bg[NSIDE+1][NHAR];
  
TH1D *hVn_Zdc_pt_raw_allCent[NSIDE+1][NHAR][NETA+4],*hVn_Zdc_pt_raw_allCent_fg[NSIDE+1][NHAR][NETA+4],*hVn_Zdc_pt_raw_allCent_bg[NSIDE+1][NHAR][NETA+4];
TH1D *hVn_Zdc_pt_res_allCent[NSIDE+1][NHAR][NETA+4],*hVn_Zdc_pt_res_allCent_fg[NSIDE+1][NHAR][NETA+4],*hVn_Zdc_pt_res_allCent_bg[NSIDE+1][NHAR][NETA+4];
 
//In SC
TProfile *hRes_Zdc_fg[NHAR],*hRes_Zdc_bg[NHAR];
TProfile *hVn_Zdc_cent_fg[NSIDE+1][NHAR],*hVn_Zdc_cent_bg[NSIDE+1][NHAR];
TProfile *hVn_Zdc_pt_fg[NCENT][NSIDE+1][NHAR],*hVn_Zdc_pt_bg[NCENT][NSIDE+1][NHAR];
TProfile *hVn_Zdc_eta_fg[NCENT][NSIDE+1][NHAR],*hVn_Zdc_eta_bg[NCENT][NSIDE+1][NHAR];
TProfile2D *hC112_Zdc_pt_same_fg[NCENT],*hC112_Zdc_pt_same_bg[NCENT];
TProfile2D *hC112_Zdc_eta_same_fg[NCENT],*hC112_Zdc_eta_same_bg[NCENT];
TProfile2D *hC112_Zdc_pt_opp_fg[NCENT],*hC112_Zdc_pt_opp_bg[NCENT];
TProfile2D *hC112_Zdc_eta_opp_fg[NCENT],*hC112_Zdc_eta_opp_bg[NCENT];

//initiallise temp hists
TH2F *hVn_FCal_tmp_fg[NCENT][NSIDE+1][NHAR],*hVn_FCal_tmp_bg[NCENT][NSIDE+1][NHAR][5];
TH2F *hVnS_FCal_tmp_fg[NCENT][NSIDE+1][NHAR],*hVnS_FCal_tmp_bg[NCENT][NSIDE+1][NHAR][5];
TH2F *hVn_Zdc_tmp_fg[NCENT][NSIDE+1][NHAR],*hVn_Zdc_tmp_bg[NCENT][NSIDE+1][NHAR][5];
TH2F *hVnS_Zdc_tmp_fg[NCENT][NSIDE+1][NHAR],*hVnS_Zdc_tmp_bg[NCENT][NSIDE+1][NHAR][5];
TH2F *hVnW_tmp_fg[NCENT][NSIDE+1][NHAR],*hVnW_tmp_bg[NCENT][NSIDE+1][NHAR][5];

//Mean pt and eta distributions
TProfile *hMean_pt[NCENT][NCH+1][NSIDE+1];
TProfile *hMean_eta[NCENT][NCH+1][NSIDE+1];
//Cumulant Histograms
TProfile *hCrn2_cent_fg[NHAR],*hCrn2_cent_bg[NHAR];
TProfile *hCrn4_cent_fg[NHAR],*hCrn4_cent_bg[NHAR];
TProfile *hCrn2_pt_fg[NCENT][NHAR],*hCrn2_pt_bg[NCENT][NHAR];
TProfile *hCrn4_pt_fg[NCENT][NHAR],*hCrn4_pt_bg[NCENT][NHAR];

TProfile *hCrn2_sub_cent_fg[NHAR],*hCrn2_sub_cent_bg[NHAR];
TProfile *hCrn4_sub_cent_fg[NHAR],*hCrn4_sub_cent_bg[NHAR];
TProfile *hCrn2_sub_pt_fg[NCENT][NHAR],*hCrn2_sub_pt_bg[NCENT][NHAR];
TProfile *hCrn4_sub_pt_fg[NCENT][NHAR],*hCrn4_sub_pt_bg[NCENT][NHAR];

TH1D *hCn_cum2_cent[NHAR],*hCn_cum4_cent[NHAR];
TH1D *hVn_cum2_cent[NHAR],*hVn_cum4_cent[NHAR];
TH1D *hCn_cum2_pt[NCENT][NHAR],*hCn_cum4_pt[NCENT][NHAR];
TH1D *hVn_cum2_pt[NCENT][NHAR],*hVn_cum4_pt[NCENT][NHAR];

TH1D *hCn_sub_cum2_cent[NHAR],*hCn_sub_cum4_cent[NHAR];
TH1D *hVn_sub_cum2_cent[NHAR],*hVn_sub_cum4_cent[NHAR];
TH1D *hCn_sub_cum2_pt[NCENT][NHAR],*hCn_sub_cum4_pt[NCENT][NHAR];
TH1D *hVn_sub_cum2_pt[NCENT][NHAR],*hVn_sub_cum4_pt[NCENT][NHAR];


TH1D* h_Qx [NCENT][NSIDE+1][NHAR];
TH1D* h_QxC[NCENT][NSIDE+1][NHAR];
TH1D* h_Qy [NCENT][NSIDE+1][NHAR];
TH1D* h_QyC[NCENT][NSIDE+1][NHAR];

TProfile* h_prof_sumx[NSIDE+1][NHAR];
TProfile* h_prof_sumy[NSIDE+1][NHAR];
TH2F* h_Q2D [NCENT][NSIDE+1][NHAR];

TProfile* h_prof_flatCos[NCENT][NSIDE+1][NHAR];
TProfile* h_prof_flatSin[NCENT][NSIDE+1][NHAR];


//For plotting
TH1D *hPsi[NCENT][NSIDE+1][NHAR];
TH1D *hPsiC[NCENT][NSIDE+1][NHAR];
TH1D *hPsiF[NCENT][NSIDE+1][NHAR];

TH1D *hVnEven_eta_raw[NCENT][NHAR][NPT+1],*hVnOdd_eta_raw[NCENT][NHAR][NPT+1];
TH1D *hVnEven_eta_res[NCENT][NHAR][NPT+1],*hVnOdd_eta_res[NCENT][NHAR][NPT+1];
TH1D *hVnEven_pt_raw[NCENT][NHAR][NETA+2],*hVnOdd_pt_raw[NCENT][NHAR][NETA+2];
TH1D *hVnEven_pt_res[NCENT][NHAR][NETA+2],*hVnOdd_pt_res[NCENT][NHAR][NETA+2];
TH1D *hVnEven_cent_raw[NHAR][NPT+1][NETA+2],*hVnOdd_cent_raw[NHAR][NPT+1][NETA+2];
TH1D *hVnEven_cent_res[NHAR][NPT+1][NETA+2],*hVnOdd_cent_res[NHAR][NPT+1][NETA+2];

TH1D *hVnEven_pt_raw_allCent[NHAR][NETA+2],*hVnOdd_pt_raw_allCent[NHAR][NETA+2];
TH1D *hVnEven_pt_res_allCent[NHAR][NETA+2],*hVnOdd_pt_res_allCent[NHAR][NETA+2];

TH1D *hV1_eta_res[NCENT][NSIDE+1];

class plot{
 public:
  //Functions
  plot();
  void init();
  void setStyle1(TH1D *h,int i);
  void setStyle2(TH1D *h,int i);
  void setStyle3(TH1D *h,int i);
  void setStyle4(TH1D *h,int i1,int i2);
  void setStyle5(TH1D *h,int i);
  void setStyle6(TH1D *h,int i,int i2);
  void setStyle7(TH1D *h,int i,int i2);
  void setStyle8(TProfile *h,int i,int i2);
  void setStyle9(TH1D *h,int i,int i2);
  double FindAvg(TH1D *h);
  void Divide_Pad(TCanvas *c,TPad *p[],int nrow,int ncol);
  void FindEvenOdd(TH1D *hE,TH1D *hO,TH1D *h1,TH1D *h2);
  void plot1_Psi();
  void plot2_Res();
  void plot3_v1_eta();
  void plot4_v1_pt();
  void plot5_v2_eta();
  void plot6_v2_pt();
  void plot7_v2_pt_FCal();
  void plot8_v2_pt();
  void plot9_v2_pt_comp();
  void write_allCent();
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
  int cnt=0;
  for(int i=0;i<N;i++){
    sum+=h->GetBinContent(i+1);
    if(h->GetBinContent(i+1)!=0.) cnt++;
  }
  return sum/double(cnt);
}

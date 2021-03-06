//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Oct  8 14:24:28 2014 by ROOT version 5.34/21
// from TTree HeavyIonD3PD/HeavyIonD3PD
// found on file: ../../rawOut/Output_zdcPico_sp24_170482.root
//////////////////////////////////////////////////////////
#ifndef extractor_flat_h
#define extractor_flat_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TMath.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <string>
#include <stdexcept>
#include "TProfile.h"
#include "TProfile2D.h"
#include "TComplex.h"

#include "Tool.cxx"

#define MAXTRK 100000

using namespace std;

//FCAL_ET cuts for 0 to 98% centrality in 2010 and 2011 PbPbruns
/*
static double centcuts[]={
  3.4432,3.2977,3.1648,3.0400,2.9222,2.8101,2.7033,2.6014,2.5038,2.4104,
  2.3208,2.2343,2.1516,2.0715,1.9941,1.9194,1.8473,1.7776,1.7100,1.6448,
  1.5815,1.5203,1.4610,1.4034,1.3479,1.2939,1.2417,1.1910,1.1418,1.0942,
  1.0482,1.0036,0.9604,0.9185,0.8780,0.8388,0.8009,0.7642,0.7288,0.6946,
  0.6615,0.6295,0.5987,0.5689,0.5401,0.5126,0.4859,0.4604,0.4357,0.4120,
  0.3892,0.3674,0.3465,0.3265,0.3073,0.2889,0.2714,0.2546,0.2387,0.2235,
  0.2090,0.1953,0.1822,0.1698,0.1582,0.1471,0.1366,0.1268,0.1175,0.1088,
  0.1005,0.0928,0.0856,0.0789,0.0725,0.0666,0.0611,0.0559,0.0511,0.0467,
  0.0425,0.0386,0.0350,0.0316,0.0285,0.0256,0.0229,0.0203,0.0180,0.0158,
  0.0137,0.0118,0.0100,0.0083,0.0067,0.0051,0.0033,-0.0085,   -1,    -2
};
*/
//FCAL_ET cuts for 0 to 80% centrality in 2015 PbPb runs
static double centcuts[]={4.26258,4.08137,3.91763,3.7635,3.61844,3.48077,3.34945,3.22397,3.10407,2.98931,2.87864,2.77237,2.66999,2.57162,2.47658,2.38468,2.29572,2.21002,2.12711,2.04651,1.96859,1.89316,1.81997,1.74932,1.68058,1.61434,1.55005,1.48744,1.42719,1.36875,1.31197,1.25693,1.20373,1.15214,1.10211,1.05367,1.0068,0.961609,0.917795,0.87541,0.834538,0.795018,0.756791,0.719896,0.684377,0.65018,0.617108,0.585275,0.554569,0.525092,0.49675,0.46959,0.443549,0.418573,0.394518,0.371561,0.349697,0.328744,0.308686,0.289595,0.27137,0.25407,0.237615,0.22199,0.207148,0.193096,0.179776,0.167193,0.155307,0.14414,0.133573,0.123657,0.114352,0.105619,0.097388,0.089723,0.082548,0.075838,0.06956,0.063719};

static double PI = acos(-1.0);
//mycode for track variable conversion
static double PI2 = 2.*acos(-1.0);
enum{
  PHIM=0x1ff, PHIS=9,
  ETAM=0x1ff, ETAS=9,
  CHM =0x1,   CHS=1,
  PTM =0x3f,  PTS=6,
  QUALM=0x7,  QUALS=3
};
static double PIslice = PI2/(PHIM+1);
static double ETAslice = 5./500;
static double ETAoff = 2.5;
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

class Event_EP{
 public:
  Event_EP(int );
  ~Event_EP(){;}
  //void clear();
  float Psi_FCal[NSIDE+1][NHAR];
  float Psi_Zdc[NSIDE+1][NHAR];

  int id;
  int icent;
  int ivz;
  //save tracks for pair cuts.
};

class Event{
 public:
  Event(int );
  ~Event(){;}
  //void clear();
  float Qx_ID[NCH+1][NPT+1][NSIDE+1][NETA2+1][NHAR];
  float Qy_ID[NCH+1][NPT+1][NSIDE+1][NETA2+1][NHAR];
  float Qw_ID[NCH+1][NPT+1][NSIDE+1][NETA2+1];

  float Qx_ID_cum[NPT+1][NSIDE+1][NHAR3][NK];
  float Qy_ID_cum[NPT+1][NSIDE+1][NHAR3][NK];
  float Qw_ID_cum[NPT+1][NSIDE+1][NP][NK];
  
  float Qx_FCal[NSIDE+1][NHAR];
  float Qy_FCal[NSIDE+1][NHAR];
  float Qw_FCal[NSIDE+1];
  float Qx_Zdc[NSIDE+1][NHAR];
  float Qy_Zdc[NSIDE+1][NHAR];
  float Qw_Zdc[NSIDE+1];

  int id;
  int icent;
  int ivz;
  //save tracks for pair cuts.
};

typedef Event_EP* EVENT_PTR_EP;
typedef Event* EVENT_PTR;

class extractor_flat {
 public :
  enum{
    //MAXPOOL=100*100,
    MAXPOOL=100,
  };
  int Pool_size;
  int count1[NCENT],count2[NCENT];
  vector< EVENT_PTR_EP > Pool_EP[NCENT][NVZ_EP];
  vector< EVENT_PTR > Pool[NCENT];
  //vector< EVENT_PTR > poolback[MAXPOOL];
  
  //Variables
  int from,to,NEv_id,NEv,RNum;
  string outname;
  //int NETA;
  //int NPT;
  int nevents;
  int centb;
  int vzb;
  int centb_1p;
  int ntrkQ;
  std::vector<int> sampSel;
  //std::vector<double> Psi_FCal_store[NCENT_1p][NSIDE][NHAR];
  //std::vector<double> Psi_Zdc_store[NCENT_1p][NSIDE][NHAR];
  std::vector<float> diffPtBins;
  std::vector<float> diffEtaBins;

  //Q-vectors
  float Psi_FCal[NSIDE+1][NHAR];
  float Psi_Zdc[NSIDE+1][NHAR];

  float Qx_ID[NCH+1][NPT+1][NSIDE+1][NETA2+1][NHAR];
  float Qy_ID[NCH+1][NPT+1][NSIDE+1][NETA2+1][NHAR];
  float Qw_ID[NCH+1][NPT+1][NSIDE+1][NETA2+1];

  float Qx_ID_cum[NPT+1][NSIDE+1][NHAR3][NK];
  float Qy_ID_cum[NPT+1][NSIDE+1][NHAR3][NK];
  float Qw_ID_cum[NPT+1][NSIDE+1][NP][NK];
  
  //EP Correction
  int    evt[NZ][NCENT_1p];
  double qmean[NZ][NCENT_1p];//fcal mean et
  float Qwmean[NDET][NZ][NCENT_1p];
  float Qxmean[NDET][NZ][NCENT_1p][NHAR];
  float Qymean[NDET][NZ][NCENT_1p][NHAR];
  float Qx_FCal[NSIDE+1][NHAR];
  float Qy_FCal[NSIDE+1][NHAR];
  float Qw_FCal[NSIDE+1];
  float Qx_Zdc[NSIDE+1][NHAR];
  float Qy_Zdc[NSIDE+1][NHAR];
  float Qw_Zdc[NSIDE+1];

  
  TFile *fOut_zdc_psi,*fOut_vn,*fOut_vn2,*fOut_doAna,*fOut_wei,*fOut_SP,*fOut_cum,*fOut_res;

  TH1D* hEt;
  TH1D* hNtrk;
  TH1D* hCent;
  TH2F* hEt_ntrk;
  TH1D* h_eTot    [NSIDE];
  TH1D* h_ePix    [NSIDE];
  TH2F* hEt_eTot  [NSIDE];
  TH2F* hEt_ePix  [NSIDE];
  TH2F* hNtrk_eTot[NSIDE];
  TH2F* hNtrk_ePix[NSIDE];

  //FCal EP Correction
  TH2F* hevts;
  TH2F* hetaphi[NZ][NCENT][NCH][NPT];
  TProfile *hQw[NDET][NZ],*hEt_vz[NZ],*hQw_m[NDET][NZ];
  TProfile *hQx[NDET][NZ][NHAR];
  TProfile *hQy[NDET][NZ][NHAR];
  TProfile *hQxS[NDET][NZ][NHAR];
  TProfile *hQyS[NDET][NZ][NHAR];
  TH1D *hqx[NDET][NZ][NHAR];
  TH1D *hqy[NDET][NZ][NHAR];
  TH1D *hqxS[NDET][NZ][NHAR];
  TH1D *hqyS[NDET][NZ][NHAR];
  //ZDC EP Correction
  TProfile* h_prof_sumx[NSIDE+1][NHAR];
  TProfile* h_prof_sumy[NSIDE+1][NHAR];
  TH2F* h_Q2D [NCENT][NSIDE+1][NHAR];
  TProfile* h_prof_flatCos[NCENT_1p][NSIDE+1][NHAR];
  TProfile* h_prof_flatSin[NCENT_1p][NSIDE+1][NHAR];


  TH1D* h_Qx [NCENT][NSIDE+1][NHAR];
  TH1D* h_QxC[NCENT][NSIDE+1][NHAR];
  TH1D* h_Qy [NCENT][NSIDE+1][NHAR];
  TH1D* h_QyC[NCENT][NSIDE+1][NHAR];



  TH1D* hPsiZ_raw[NCENT][NSIDE+1][NHAR];
  TH1D* hPsiZ_rec[NCENT][NSIDE+1][NHAR];
  TH1D* hPsiZ_flat[NCENT][NSIDE+1][NHAR];
  
  TH1D* hRes      [NCENT][NHAR];
  TH1D* hRes_recn [NCENT][NHAR];
  TH1D* hRes_flat [NCENT][NHAR];

  TH1D* hResCent      [NHAR];
  TH1D* hResCent_recn [NHAR];
  TH1D* hResCent_flat [NHAR];
  TH1D* hResCent_mix  [NHAR];
  
  TH1D* hDphi_fg [NCENT][NHAR];
  TH1D* hDphi_bg [NCENT][NHAR];
  TH1D* hDphi_rat[NCENT][NHAR];


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



  TH1D *hNTrk;
  //Track histograms
  TH1D *hPt,*hEta,*hPhi,*hPhi_wei,*hCharge,*hTrk_w,*hTrk_wei,*hTrk_eff,*hWeight;
  TProfile *hTrk_eff2[NCENT][NPT];
  TH2F *hEtaPhi[NCENT][NVZ_EP][NPT_EP][NCH_EP],*hEtaPhi2[NCENT][NVZ_EP][NPT_EP][NCH_EP],*hEtaPhi_wei[NCENT][NVZ_EP][NPT_EP][NCH_EP];
  Tool *effTool;
      
  extractor_flat(){};
  extractor_flat(string filelist,int rno,int iev,int nev);
  virtual ~extractor_flat();

  /*
    void Exec_zdc_psi(int entSel = 1);
    void Init_zdc_psi();
    void Run_zdc_psi();
    void SaveHistos_zdc_psi();

    void Exec_wei(int entSel = 1);
    void Init_wei();
    void Run_wei();
    void SaveHistos_wei();

    void FillCalib();
    void FillFlattening();

  */
  void Exec_res(int entSel = 1);
  void ReadHistos_res();
  void Init_res();
  void Run_res();
  void FillGlobal_res();
  void SaveHistos_res();

  void Exec_vn(int entSel = 1);
  void ReadHistos_vn();
  void Init_vn();
  void Run_vn();
  void FillGlobal_vn();
  void FillTrack_vn(float Psi_FCal[3][4],float Psi_Zdc[3][4]);
  void SaveHistos_vn();

  void Exec_doAna(int entSel = 1);
  void ReadHistos_doAna();
  void Init_doAna();
  void Run_doAna();
  void doAna_ResFactor(TH1D *hRes_ful,TH1D *hResX_ful,TH1D *hResX_sub,TH1D *hRes);
  void doAna_FCal();
  void doAna_Zdc();
  void doAna_Cum();
  void SaveHistos_doAna();

  /*
  void Exec_vn2(int entSel = 1);
  void ReadHistos_vn2();
  void Init_vn2();
  void Run_vn2();
  void FillGlobal_vn2();
  void FillTrack_vn2(float Psi_FCal[3][4],float Psi_Zdc[3][4]);
  void SaveHistos_vn2();

  void Exec_SP(int entSel = 1);
  void ReadHistos_SP();
  void Init_SP();
  void Run_SP();
  void GetTrackQn_SP();
  void FillHisto_SP(int id0);
  void Fill_2PCfg_SP(Event* ev1, Event *ev2);
  void Fill_3PCfg_SP(Event* ev1, Event *ev2, Event *ev3);
  void SaveHistos_SP();

  void Exec_cum(int entSel = 1);
  void ReadHistos_cum();
  void Init_cum();
  void Run_cum();
  void GetTrackQn_cum();
  void FillHisto_cum(int id0);
  void Fill_2PCfg_cum(Event* ev1);
  void Fill_4PCfg_cum(Event* ev1);
  void SaveHistos_cum();
  */
  //void FillTrackVn(std::vector<std::vector<double> > PsiF);
  //void FillTrackVn(double PsiF[3][4]);
  void get_qVec_zdc(std::vector<std::vector<double> >& Qx, std::vector<std::vector<double> >& Qy, std::vector<double>& Qw);//new
  void get_qVec(std::vector<std::vector<double> >& Qx, std::vector<std::vector<double> >& Qy, std::vector<double>& Qw);//old
  std::vector<double> discFourrier(const TH1D* h1, int ihar, TFile* f1=NULL);
  float get_pt(int ptbin);
  int  get_centb (double et );
  int  get_vzb (double vz );
  int  get_ptBin (double pt );
  int  get_etaBin(double eta);
  int get_ptBin_TrkWei(int ptb);

  void set_ptBins ();
  void set_etaBins();
  void DivideHists(TH1D *h,TH1D *hnum,TH1D *hden);
  
  //virtual Int_t    Cut(Long64_t entry);
  //virtual Int_t    GetEntry(Long64_t entry);
  //virtual Long64_t LoadTree(Long64_t entry);
  //virtual void     Init(TTree *tree);
  //virtual Bool_t   Notify();
  //virtual void     Show(Long64_t entry = -1);

  TChain          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain
  char name[500];
      
  // Declaration of leaf types
  UInt_t          RunNumber;
  UInt_t          EventNumber;
  UInt_t          lbn;
  Int_t           vx_n;
  Float_t         vx_x;
  Float_t         vx_y;
  Float_t         vx_z;
  Float_t         Fcal_Et;
  Float_t         Fcal_Et_p;
  Float_t         Fcal_Et_n;
  Int_t           Centrality;
  Int_t           nloose;
  Int_t           ntight;
  Int_t           trk_n;
  UInt_t          trk_data[MAXTRK];   //[Mytrk_n]
  Float_t         Qw[16];
  Float_t         Qx[16][7];
  Float_t         Qy[16][7];
  Float_t         Qw_zdc[3];
  Float_t         Qx_zdc[3][7];
  Float_t         Qy_zdc[3][7];
  Int_t           Zdc_n;
  vector<float>   *Zdc_Energy_LG;
  vector<float>   *Zdc_Energy_HG;
  vector<unsigned int> *Zdc_Id;
  vector<int>     *Zdc_Side;
  vector<int>     *Zdc_Type;
  vector<int>     *Zdc_Module;
  vector<int>     *Zdc_Channel;

  //Extra
  Int_t           mbtime_countA;
  Int_t           mbtime_countC;
  Float_t         trk_pt[MAXTRK];   //[Mytrk_n]
  Float_t         trk_eta[MAXTRK];   //[Mytrk_n]
  Float_t         trk_phi0_wrt_PV[MAXTRK];   //[Mytrk_n]
  Int_t           trk_Quality2[MAXTRK];   //[Mytrk_n]
  Int_t           trk_charge[MAXTRK];   //[Mytrk_n]



  //extra mycode
  Int_t trk_phibin[MAXTRK];
  Int_t trk_etabin[MAXTRK];
  Int_t trk_ptbin[MAXTRK];

  Float_t trk_w[MAXTRK];
  Float_t trk_wei[MAXTRK];
  Float_t trk_eff[MAXTRK];
  Float_t trk_phi[MAXTRK];
      
};

#endif


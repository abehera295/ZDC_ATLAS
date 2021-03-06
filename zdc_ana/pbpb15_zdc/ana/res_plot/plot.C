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
  NSIDE   = 2,
  NMODULE = 2,
  NGAIN   = 2,
  NCH1    = 4,
  NCH2    = 24,
  NHAR    = 3,
  //NCENT   = 20,
  NCENT   = 11,
  NCENT_1p  = 100,
  NK      = 12,
  NSAMP   = 10,
  STORE_DEP = 20
};

enum{
  NETA=12,
  NETA2=10,
  NCH=2,
  NPT=10
};

const double pt_mat[]={0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 7.0, 9.0, 12.0};
const double eta_mat[]={-2.5, -2.0, -1.5, -1.0, -0.5, -0.25, 0.0, 0.25, 0.5, 1.0, 1.5, 2.0, 2.5};
string Cent_mat[]={"0-5%","5-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%",};
//string Cent_mat[]={"0-5%","5-10%","10-15%","15-20%","20-25%","25-30%","30-35%","35-40%","40-45%","45-50%","50-55%","55-60%","60-65%","65-70%","70-75%","75-80%","80-85%","85-90%","90-95%","95-100%",};
string Cent_mat_new[]={"0-5%","5-10%","10-15%","15-20%","20-25%","25-30%","30-35%","35-40%","40-45%","45-50%","50-55%","55-60%","60-65%","65-70%","70-75%","75-80%","80-85%","85-90%","90-95%","95-100%","0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%","0-20%","20-40%","40-60%","60-80%","80-100%","0-40%","40-80%",};
double cent_mat_double[]={0,5,10,20,30,40,50,60,70,80,90,100};

char name[20];
const double PI=acos(-1.);

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

std::vector<double> discFourrier(const TH1D* h1, int ihar) {

  std::vector<double> vn;
  if(h1->GetEntries() == 0) {
    //std::cout << "discFourrier ERROR : Empty histogram " << h1->GetName() << std::endl;
    vn.push_back(0.);
    vn.push_back(0.);
    return vn;
  }

  double p = 0, ep = 0;
  double bwid    = h1->GetBinWidth(1)/2;
  int    Nbins   = h1->GetNbinsX();

  double Sum=0,Sum_err=0;//Find vn
  int b1,b2;
  b1=h1->FindFirstBinAbove(0.);
  b2=h1->FindLastBinAbove(0.);
  for(int phibin=1; phibin<=Nbins;phibin++){
    //for(int phibin=b1+1; phibin<=b2-1;phibin++){
    float binval = h1->GetBinContent(phibin);
    float phi    = h1->GetBinCenter(phibin);
    Sum         += binval;
    Sum_err     += pow(h1->GetBinError(phibin),2);
    //p += cos( (ihar)*phi) * binval;
    p += cos( (ihar+1)*phi) * binval;
  }
  double SumAvg = Sum/Nbins;
  //double ErrAvg = sqrt(Sum_err)/Nbins;
  p /=Sum;

  for(int phibin=1; phibin<=Nbins;phibin++){//Find Errors
  //for(int phibin=b1+1; phibin<=b2-1;phibin++){//Find Errors
    float binerr = h1->GetBinError(phibin);
    float phi    = h1->GetBinCenter(phibin);
    //ep  += pow((cos((ihar)*phi)-p)*binerr,2);
    ep  += pow((cos((ihar+1)*phi)-p)*binerr,2);
  }

  //Correct for Bin Width
  double corr = (ihar+1)*bwid/sin((ihar+1)*bwid);
  corr=1.;
  p*=corr;
  ep=sqrt(ep)/Sum*corr;

  vn.push_back(p);
  vn.push_back(ep);
  /*
  if (f1 != NULL) {//Write Histogram with Harmonics
    std::string hName = h1->GetName();
    hName += "_fourhar";
    TH1D* hFour = (TH1D*)h1->Clone(hName.c_str());
    hFour->SetMarkerStyle(33);
    hFour->SetMarkerSize (1.3);
    hFour->SetMarkerColor(1);
    hFour->SetLineColor  (1);

    double hLowEdge = h1->GetBinLowEdge(1);
    double hUpEdge  = h1->GetBinLowEdge(h1->GetNbinsX()+1);

    TF1* fHar = (TF1*)NULL;
    std::stringstream      ss_har; ss_har << ihar;
    std::string fName    = "FourrierHar" + ss_har.str();
    std::string fFormula = "[0]*(1+2*[1]*cos(" + ss_har.str() + "*x))";

    fHar = new TF1(fName.c_str(), fFormula.c_str(), hLowEdge, hUpEdge);
    fHar ->SetLineWidth(1);
    fHar ->SetLineStyle(5);
    fHar ->SetLineColor(ihar%8+2);
    hFour->Fit(fHar,"QR");
    fHar ->SetParameter(1, vn.at(0));
    fHar ->SetParError (1, vn.at(1));
    hFour->GetListOfFunctions()->Add(fHar);

    f1->cd(0); if (!(f1->GetListOfKeys()->Contains(hFour->GetName()))) hFour->Write();
    hFour->GetListOfFunctions()->Remove(fHar);

    delete fHar;
  }
  */
  return vn;
}
  
void setStyle1(TH1D *h,int itype1){
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

double FindAvg(TH1D *h){
  int N=h->GetNbinsX();
  double sum=0.;
  int cnt=0;
  for(int i=0;i<N;i++){
    sum+=h->GetBinContent(i+1);
    cnt++;
  }
  return sum/double(N);
}

double FindAvg2(TH1D *h){
  int N=h->GetNbinsX();
  double sum=0.;
  int b1,b2;
  b1=h->FindBin(-3.);
  b2=h->FindBin(3.);
  int cnt=0;
  for(int i=b1;i<b2;i++){
    sum+=h->GetBinContent(i+1);
    cnt++;
  }
  return sum/double(cnt);
}

void GetTH1D(TProfile *h1,TH1D *h2){
  int N=h1->GetNbinsX();
  for(int i=0;i<N;i++){
    double y=h1->GetBinContent(i+1);
    double dy=h1->GetBinError(i+1);
    h2->SetBinContent(i+1,y);
    h2->SetBinError(i+1,dy);
  }
}

void plot(){
  //Read hists
  static const int NCENT2=100;
  TH1D *hPsi2[NCENT2][NSIDE+1][NHAR];
  TH1D *hPsiC2[NCENT2][NSIDE+1][NHAR];
  TH1D *hPsiF2[NCENT2][NSIDE+1][NHAR];


  TH1D *hPsi[NCENT][NSIDE+1][NHAR];
  TH1D *hPsiC[NCENT][NSIDE+1][NHAR];
  TH1D *hPsiF[NCENT][NSIDE+1][NHAR];

  TH1D *hDphi_FCal_fg[NCENT][NHAR];
  TH1D *hDphi_FCal_bg[NCENT][NHAR];
  TH1D *hDphi_FCal_rat[NCENT][NHAR];

  TH1D *hCorrCent_fg_FCal[NHAR],*hCorrCent_bg_FCal[NHAR],*hCorrCent_FCal[NHAR];
  TH1D *hCorrCent_fg_FCal2[NHAR],*hCorrCent_bg_FCal2[NHAR],*hCorrCent_FCal2[NHAR];
  TH1D *hCorrCent_fg_Zdc[NHAR],*hCorrCent_bg_Zdc[NHAR],*hCorrCent_Zdc[NHAR];
  TH1D *hCorrCent_fg_Zdc2[NHAR],*hCorrCent_bg_Zdc2[NHAR],*hCorrCent_Zdc2[NHAR];

  TH1D *hResCent_fg_FCal[NHAR],*hResCent_bg_FCal[NHAR],*hResCent_FCal[NHAR];
  TH1D *hResCent_fg_FCal2[NHAR],*hResCent_bg_FCal2[NHAR],*hResCent_FCal2[NHAR];
  TH1D *hResCent_fg_Zdc[NHAR],*hResCent_bg_Zdc[NHAR],*hResCent_Zdc[NHAR];
  TH1D *hResCent_fg_Zdc2[NHAR],*hResCent_bg_Zdc2[NHAR],*hResCent_Zdc2[NHAR];

  TH1D *hDphi_Zdc_fg[NCENT][NHAR];
  TH1D *hDphi_Zdc_bg[NCENT][NHAR];
  TH1D *hDphi_Zdc_rat[NCENT][NHAR];
  
  TFile *fin1,*fin2;
  TLatex text;
  text.SetNDC(1);
  text.SetTextFont(43);
  text.SetTextSize(18);

  //1. Get hPsi
  /*
  fin1=new TFile("../ana/hists_step3/EP287843_flattened.root","read");
  //fin1=new TFile("../zdcvn/Output_res.root","read");
  for (int icent=0; icent<NCENT2; icent++) {
    for (int iside=0; iside<NSIDE+1; iside++) {
      for (int ihar=0; ihar<NHAR; ihar++) {
	sprintf(name, "hPsi_zdc_cent%d_side%d_har%d", icent, iside, ihar);
	//sprintf(name, "hPsiZ_raw_cent%d_side%d_har%d", icent, iside, ihar);
	hPsi2[icent][iside][ihar]=(TH1D*)fin1->Get(name);
	sprintf(name, "hPsiC_zdc_cent%d_side%d_har%d", icent, iside, ihar);
	//sprintf(name, "hPsiZ_rec_cent%d_side%d_har%d", icent, iside, ihar);
	hPsiC2[icent][iside][ihar]=(TH1D*)fin1->Get(name);
	sprintf(name, "hPsiF_zdc_cent%d_side%d_har%d", icent, iside, ihar);
	//sprintf(name, "hPsiZ_flat_cent%d_side%d_har%d", icent, iside, ihar);
	hPsiF2[icent][iside][ihar]=(TH1D*)fin1->Get(name);
      }
    }
  }

  //Add Centralities to make 100 cents to 10 cents
  for (int iside=0; iside<NSIDE+1; iside++) {
    for (int ihar=0; ihar<NHAR; ihar++) {
      int cnt1=0,cnt2=0;
      for(int icent=0;icent<NCENT2;icent++){
	if(cnt1==0){
	  hPsi[cnt2][iside][ihar]=(TH1D*)hPsi2[icent][iside][ihar]->Clone();
	  hPsiC[cnt2][iside][ihar]=(TH1D*)hPsiC2[icent][iside][ihar]->Clone();
	  hPsiF[cnt2][iside][ihar]=(TH1D*)hPsiF2[icent][iside][ihar]->Clone();
	  sprintf(name,"hPsi_cent%d_side%d_har%d",cnt2,iside,ihar);
	  hPsi[cnt2][iside][ihar]->SetName(name);
	  sprintf(name,"hPsiC_cent%d_side%d_har%d",cnt2,iside,ihar);
	  hPsiC[cnt2][iside][ihar]->SetName(name);
	  sprintf(name,"hPsiF_cent%d_side%d_har%d",cnt2,iside,ihar);
	  hPsiF[cnt2][iside][ihar]->SetName(name);
	}
	else{
	  hPsi[cnt2][iside][ihar]->Add(hPsi2[icent][iside][ihar]);
	  hPsiC[cnt2][iside][ihar]->Add(hPsiC2[icent][iside][ihar]);
	  hPsiF[cnt2][iside][ihar]->Add(hPsiF2[icent][iside][ihar]);
	}
	cnt1++;
	if(icent<=10){
	  if(cnt1==5) {cnt1=0; cnt2++;}
	}
	else{
	  if(cnt1==10) {cnt1=0; cnt2++;}
	}
      }
    }
  }
  */
  //1. Get hPsi
  fin1=new TFile("../zdcvn/Output_res.root","read");
  for (int icent=0; icent<NCENT; icent++) {
    for (int iside=0; iside<NSIDE+1; iside++) {
      for (int ihar=0; ihar<NHAR; ihar++) {
	sprintf(name, "hPsiZ_raw_cent%d_side%d_har%d", icent, iside, ihar);
	hPsi[icent][iside][ihar]=(TH1D*)fin1->Get(name);
	sprintf(name, "hPsiZ_rec_cent%d_side%d_har%d", icent, iside, ihar);
	hPsiC[icent][iside][ihar]=(TH1D*)fin1->Get(name);
	sprintf(name, "hPsiZ_flat_cent%d_side%d_har%d", icent, iside, ihar);
	hPsiF[icent][iside][ihar]=(TH1D*)fin1->Get(name);
      }
    }
  }
  
  for (int icent=0; icent<NCENT; icent++) {
    for (int iside=0; iside<NSIDE+1; iside++) {
      for (int ihar=0; ihar<NHAR; ihar++) {
	//Scale
	double int1=hPsi[icent][iside][ihar]->Integral();
	if(int1!=0.) hPsi[icent][iside][ihar]->Scale(1./int1);
	int1=hPsiC[icent][iside][ihar]->Integral();
	if(int1!=0.) hPsiC[icent][iside][ihar]->Scale(1./int1);
	int1=hPsiF[icent][iside][ihar]->Integral();
	if(int1!=0.) hPsiF[icent][iside][ihar]->Scale(1./int1);
	double mean0=FindAvg(hPsiF[icent][iside][ihar]);
	if(mean0!=0.) hPsi[icent][iside][ihar]->Scale(1./mean0);
	if(mean0!=0.) hPsiC[icent][iside][ihar]->Scale(1./mean0);
	if(mean0!=0.) hPsiF[icent][iside][ihar]->Scale(1./mean0);
	//Style
	setStyle1(hPsi[icent][iside][ihar],0);
	setStyle1(hPsiC[icent][iside][ihar],1);
	setStyle1(hPsiF[icent][iside][ihar],2);
      }
    }
  }

  //2. Get hDPsi
  fin2=new TFile("../zdcvn/Output_res.root","read");
  for (int icent=0; icent<NCENT; icent++) {
    for (int ihar=0; ihar<NHAR; ihar++) {
      sprintf(name, "hDphi_FCal_fg_cent%d_har%d", icent, ihar);
      hDphi_FCal_fg[icent][ihar] = (TH1D*)fin2->Get(name);
      sprintf(name, "hDphi_FCal_bg_cent%d_har%d", icent, ihar);
      hDphi_FCal_bg[icent][ihar] = (TH1D*)fin2->Get(name);

      /*
      sprintf(name, "hDphi_FCal_rat_cent%d_har%d", icent, ihar);
      hDphi_FCal_rat[icent][ihar] = (TH1D*)hDphi_FCal_fg[icent][ihar]->Clone(name);
      hDphi_FCal_rat[icent][ihar]->Divide(hDphi_FCal_bg[icent][ihar]);
      
      //Rebin
      hDphi_FCal_fg[icent][ihar]->Rebin(8);
      hDphi_FCal_bg[icent][ihar]->Rebin(8);
      hDphi_FCal_rat[icent][ihar]->Rebin(8);
      */
      //Scale
      double int1=hDphi_FCal_fg[icent][ihar]->Integral();
      if(int1!=0.) hDphi_FCal_fg[icent][ihar]->Scale(1./int1);
      int1=hDphi_FCal_bg[icent][ihar]->Integral();
      if(int1!=0.) hDphi_FCal_bg[icent][ihar]->Scale(1./int1);
      //hDphi_FCal_rat[icent][ihar]->Scale(1./hDphi_FCal_rat[icent][ihar]->Integral());
      double mean0=FindAvg2(hDphi_FCal_bg[icent][ihar]);
      if(mean0!=0.) hDphi_FCal_fg[icent][ihar]->Scale(1./mean0);
      if(mean0!=0.) hDphi_FCal_bg[icent][ihar]->Scale(1./mean0);
      //hDphi_FCal_rat[icent][ihar]->Scale(1./mean0);
      //Style
      setStyle1(hDphi_FCal_fg[icent][ihar],0);
      setStyle1(hDphi_FCal_bg[icent][ihar],1);
      //setStyle1(hDphi_FCal_rat[icent][ihar],2);

      sprintf(name, "hDphi_Zdc_fg_cent%d_har%d", icent, ihar);
      hDphi_Zdc_fg[icent][ihar] = (TH1D*)fin2->Get(name);
      sprintf(name, "hDphi_Zdc_bg_cent%d_har%d", icent, ihar);
      hDphi_Zdc_bg[icent][ihar] = (TH1D*)fin2->Get(name);

      /*
      sprintf(name, "hDphi_Zdc_rat_cent%d_har%d", icent, ihar);
      hDphi_Zdc_rat[icent][ihar] = (TH1D*)hDphi_Zdc_fg[icent][ihar]->Clone(name);
      hDphi_Zdc_rat[icent][ihar]->Divide(hDphi_Zdc_bg[icent][ihar]);
      */
      /*
      //Rebin
      hDphi_Zdc_fg[icent][ihar]->Rebin(8);
      hDphi_Zdc_bg[icent][ihar]->Rebin(8);
      hDphi_Zdc_rat[icent][ihar]->Rebin(8);
      */
      //Scale
      int1=hDphi_Zdc_fg[icent][ihar]->Integral();
      if(int1!=0.) hDphi_Zdc_fg[icent][ihar]->Scale(1./int1);
      int1=hDphi_Zdc_bg[icent][ihar]->Integral();
      if(int1!=0.) hDphi_Zdc_bg[icent][ihar]->Scale(1./int1);
      //hDphi_Zdc_rat[icent][ihar]->Scale(1./hDphi_Zdc_rat[icent][ihar]->Integral());
      mean0=FindAvg2(hDphi_Zdc_bg[icent][ihar]);
      if(mean0!=0.) hDphi_Zdc_fg[icent][ihar]->Scale(1./mean0);
      if(mean0!=0.) hDphi_Zdc_bg[icent][ihar]->Scale(1./mean0);
      //hDphi_Zdc_rat[icent][ihar]->Scale(1./mean0);
      //Style
      setStyle1(hDphi_Zdc_fg[icent][ihar],0);
      setStyle1(hDphi_Zdc_bg[icent][ihar],1);
      //setStyle1(hDphi_Zdc_rat[icent][ihar],2);
    }
  }

  //3. Get hResCent
  //Init
  for(int ihar=0;ihar<NHAR;ihar++){
    //Correlation
    sprintf(name,"hCorrCent_fg_FCal_har%d",ihar);
    hCorrCent_fg_FCal[ihar]=new TH1D(name,"",NCENT,cent_mat_double);
    sprintf(name,"hCorrCent_bg_FCal_har%d",ihar);
    hCorrCent_bg_FCal[ihar]=new TH1D(name,"",NCENT,cent_mat_double);
    
    sprintf(name,"hCorrCent_fg_FCal2_har%d",ihar);
    hCorrCent_fg_FCal2[ihar]=new TH1D(name,"",NCENT,cent_mat_double);
    sprintf(name,"hCorrCent_bg_FCal2_har%d",ihar);
    hCorrCent_bg_FCal2[ihar]=new TH1D(name,"",NCENT,cent_mat_double);

    sprintf(name,"hCorrCent_fg_Zdc_har%d",ihar);
    hCorrCent_fg_Zdc[ihar]=new TH1D(name,"",NCENT,cent_mat_double);
    sprintf(name,"hCorrCent_bg_Zdc_har%d",ihar);
    hCorrCent_bg_Zdc[ihar]=new TH1D(name,"",NCENT,cent_mat_double);

    sprintf(name,"hCorrCent_fg_Zdc2_har%d",ihar);
    hCorrCent_fg_Zdc2[ihar]=new TH1D(name,"",NCENT,cent_mat_double);
    sprintf(name,"hCorrCent_bg_Zdc2_har%d",ihar);
    hCorrCent_bg_Zdc2[ihar]=new TH1D(name,"",NCENT,cent_mat_double);

    //Resolution
    sprintf(name,"hResCent_fg_FCal_har%d",ihar);
    hResCent_fg_FCal[ihar]=new TH1D(name,"",NCENT,cent_mat_double);
    sprintf(name,"hResCent_bg_FCal_har%d",ihar);
    hResCent_bg_FCal[ihar]=new TH1D(name,"",NCENT,cent_mat_double);
    
    sprintf(name,"hResCent_fg_FCal2_har%d",ihar);
    hResCent_fg_FCal2[ihar]=new TH1D(name,"",NCENT,cent_mat_double);
    sprintf(name,"hResCent_bg_FCal2_har%d",ihar);
    hResCent_bg_FCal2[ihar]=new TH1D(name,"",NCENT,cent_mat_double);

    sprintf(name,"hResCent_fg_Zdc_har%d",ihar);
    hResCent_fg_Zdc[ihar]=new TH1D(name,"",NCENT,cent_mat_double);
    sprintf(name,"hResCent_bg_Zdc_har%d",ihar);
    hResCent_bg_Zdc[ihar]=new TH1D(name,"",NCENT,cent_mat_double);

    sprintf(name,"hResCent_fg_Zdc2_har%d",ihar);
    hResCent_fg_Zdc2[ihar]=new TH1D(name,"",NCENT,cent_mat_double);
    sprintf(name,"hResCent_bg_Zdc2_har%d",ihar);
    hResCent_bg_Zdc2[ihar]=new TH1D(name,"",NCENT,cent_mat_double);
  }
  for(int ihar=0;ihar<NHAR;ihar++){
    for (int icent=0; icent<NCENT; icent++) {
      double y,dy;
      std::vector<double> res_vec_FCal_fg = discFourrier(hDphi_FCal_fg[icent][ihar],0);
      hCorrCent_fg_FCal[ihar]->SetBinContent(icent+1, res_vec_FCal_fg.at(0) );
      hCorrCent_fg_FCal[ihar]->SetBinError(icent+1, res_vec_FCal_fg.at(1) );
      y=0.,dy=0.;
      hResCent_fg_FCal[ihar]->SetBinContent(icent+1, sqrt(fabs(res_vec_FCal_fg.at(0))) );
      if(res_vec_FCal_fg.at(0)!=0.) dy=0.5*pow(fabs(res_vec_FCal_fg.at(0)),-0.5)*res_vec_FCal_fg.at(1);
      hResCent_fg_FCal[ihar]->SetBinError(icent+1, fabs(dy) );

      std::vector<double> res_vec_FCal_bg = discFourrier(hDphi_FCal_bg[icent][ihar],0);
      hCorrCent_bg_FCal[ihar]->SetBinContent(icent+1, res_vec_FCal_bg.at(0) );
      hCorrCent_bg_FCal[ihar]->SetBinError(icent+1, res_vec_FCal_bg.at(1) );
      y=0.,dy=0.;
      hResCent_bg_FCal[ihar]->SetBinContent(icent+1, sqrt(fabs(res_vec_FCal_bg.at(0))) );
      if(res_vec_FCal_bg.at(0)!=0.) dy=0.5*pow(fabs(res_vec_FCal_bg.at(0)),-0.5)*res_vec_FCal_bg.at(1);
      hResCent_bg_FCal[ihar]->SetBinError(icent+1, fabs(dy) );
      
      std::vector<double> res_vec_Zdc_fg = discFourrier(hDphi_Zdc_fg[icent][ihar],0);
      hCorrCent_fg_Zdc[ihar]->SetBinContent(icent+1, res_vec_Zdc_fg.at(0) );
      hCorrCent_fg_Zdc[ihar]->SetBinError(icent+1, res_vec_Zdc_fg.at(1) );
      y=0.,dy=0.;
      hResCent_fg_Zdc[ihar]->SetBinContent(icent+1, sqrt(fabs(res_vec_Zdc_fg.at(0))) );
      if(res_vec_Zdc_fg.at(0)!=0.) dy=0.5*pow(fabs(res_vec_Zdc_fg.at(0)),-0.5)*res_vec_Zdc_fg.at(1);
      hResCent_fg_Zdc[ihar]->SetBinError(icent+1, fabs(dy) );
      
      std::vector<double> res_vec_Zdc_bg = discFourrier(hDphi_Zdc_bg[icent][ihar],0);
      hCorrCent_bg_Zdc[ihar]->SetBinContent(icent+1, res_vec_Zdc_bg.at(0) );
      hCorrCent_bg_Zdc[ihar]->SetBinError(icent+1, res_vec_Zdc_bg.at(1) );
      y=0.,dy=0.;
      hResCent_bg_Zdc[ihar]->SetBinContent(icent+1, sqrt(fabs(res_vec_Zdc_bg.at(0))) );
      if(res_vec_Zdc_bg.at(0)!=0.) dy=0.5*pow(fabs(res_vec_Zdc_bg.at(0)),-0.5)*res_vec_Zdc_bg.at(1);
      hResCent_bg_Zdc[ihar]->SetBinError(icent+1, fabs(dy) );
    }
    sprintf(name,"hCorrCent_FCal_har%d",ihar);
    hCorrCent_FCal[ihar]=(TH1D*)hCorrCent_fg_FCal[ihar]->Clone(name);
    hCorrCent_FCal[ihar]->Add(hCorrCent_bg_FCal[ihar],-1.);
    
    sprintf(name,"hCorrCent_Zdc_har%d",ihar);
    hCorrCent_Zdc[ihar]=(TH1D*)hCorrCent_fg_Zdc[ihar]->Clone(name);
    hCorrCent_Zdc[ihar]->Add(hCorrCent_bg_Zdc[ihar],-1.);

    sprintf(name,"hResCent_FCal_har%d",ihar);
    hResCent_FCal[ihar]=(TH1D*)hResCent_fg_FCal[ihar]->Clone(name);
    //hResCent_FCal[ihar]->Add(hResCent_bg_FCal[ihar],-1.);
    hResCent_FCal[ihar]->Reset();
    
    sprintf(name,"hResCent_Zdc_har%d",ihar);
    hResCent_Zdc[ihar]=(TH1D*)hResCent_fg_Zdc[ihar]->Clone(name);
    //hResCent_Zdc[ihar]->Add(hResCent_bg_Zdc[ihar],-1.);
    hResCent_Zdc[ihar]->Reset();
    
    for (int icent=0; icent<NCENT; icent++) {
      double y,y1,y2,dy,dy1,dy2;
      y1=hCorrCent_FCal[ihar]->GetBinContent(icent+1);
      dy1=hCorrCent_FCal[ihar]->GetBinError(icent+1);
      y=0.,dy=0.;
      y=sqrt(fabs(y1));
      if(y1!=0.) dy=0.5*pow(fabs(y1),-0.5)*dy1;
      hResCent_FCal[ihar]->SetBinContent(icent+1, y );
      hResCent_FCal[ihar]->SetBinError(icent+1, fabs(dy) );

      y1=hCorrCent_Zdc[ihar]->GetBinContent(icent+1);
      dy1=hCorrCent_Zdc[ihar]->GetBinError(icent+1);
      y=0.,dy=0.;
      y=sqrt(fabs(y1));
      if(y1!=0.) dy=0.5*pow(fabs(y1),-0.5)*dy1;
      hResCent_Zdc[ihar]->SetBinContent(icent+1, y );
      hResCent_Zdc[ihar]->SetBinError(icent+1, fabs(dy) );
    }
  }
  
  //4. Get hResCent using TProfile
  //Read TProfiles
  TProfile *htmp;
  for (int ihar=0; ihar<NHAR; ihar++) {
    //Correlation
    sprintf(name, "hResCent_FCal_fg_har%d",ihar);
    htmp= (TProfile*)fin2->Get(name);
    GetTH1D(htmp,hCorrCent_fg_FCal2[ihar]);
    sprintf(name, "hResCent_FCal_bg_har%d",ihar);
    htmp= (TProfile*)fin2->Get(name);
    GetTH1D(htmp,hCorrCent_bg_FCal2[ihar]);
    sprintf(name,"hCorrCent_fg_FCal2_ihar%d",ihar);
    hCorrCent_FCal2[ihar]=(TH1D*)hCorrCent_fg_FCal2[ihar]->Clone(name);
    hCorrCent_FCal2[ihar]->Add(hCorrCent_bg_FCal2[ihar],-1.);
    
    
    sprintf(name, "hResCent_Zdc_fg_har%d",ihar);
    htmp= (TProfile*)fin2->Get(name);
    GetTH1D(htmp,hCorrCent_fg_Zdc2[ihar]);
    sprintf(name, "hResCent_Zdc_bg_har%d",ihar);
    htmp= (TProfile*)fin2->Get(name);
    GetTH1D(htmp,hCorrCent_bg_Zdc2[ihar]);
    sprintf(name,"hCorrCent_fg_Zdc2_ihar%d",ihar);
    hCorrCent_Zdc2[ihar]=(TH1D*)hCorrCent_fg_Zdc2[ihar]->Clone(name);
    hCorrCent_Zdc2[ihar]->Add(hCorrCent_bg_Zdc2[ihar],-1.);

    //Resolution
     for (int icent=0; icent<NCENT; icent++) {
       double y,y1,y2,dy,dy1,dy2;
       y1=hCorrCent_fg_FCal2[ihar]->GetBinContent(icent+1);
       dy1=hCorrCent_fg_FCal2[ihar]->GetBinError(icent+1);
       y=0.,dy=0.;
       y=sqrt(fabs(y1));
       if(y1!=0.) dy=0.5*pow(fabs(y1),-0.5)*dy1;
       hResCent_fg_FCal2[ihar]->SetBinContent(icent+1, y );
       hResCent_fg_FCal2[ihar]->SetBinError(icent+1, fabs(dy) );

       y1=hCorrCent_bg_FCal2[ihar]->GetBinContent(icent+1);
       dy1=hCorrCent_bg_FCal2[ihar]->GetBinError(icent+1);
       y=0.,dy=0.;
       y=sqrt(fabs(y1));
       if(y1!=0.) dy=0.5*pow(fabs(y1),-0.5)*dy1;
       hResCent_bg_FCal2[ihar]->SetBinContent(icent+1, y );
       hResCent_bg_FCal2[ihar]->SetBinError(icent+1, fabs(dy) );

       y1=hCorrCent_fg_Zdc2[ihar]->GetBinContent(icent+1);
       dy1=hCorrCent_fg_Zdc2[ihar]->GetBinError(icent+1);
       y=0.,dy=0.;
       y=sqrt(fabs(y1));
       if(y1!=0.) dy=0.5*pow(fabs(y1),-0.5)*dy1;
       hResCent_fg_Zdc2[ihar]->SetBinContent(icent+1, y );
       hResCent_fg_Zdc2[ihar]->SetBinError(icent+1, fabs(dy) );

       y1=hCorrCent_bg_Zdc2[ihar]->GetBinContent(icent+1);
       dy1=hCorrCent_bg_Zdc2[ihar]->GetBinError(icent+1);
       y=0.,dy=0.;
       y=sqrt(fabs(y1));
       if(y1!=0.) dy=0.5*pow(fabs(y1),-0.5)*dy1;
       hResCent_bg_Zdc2[ihar]->SetBinContent(icent+1, y );
       hResCent_bg_Zdc2[ihar]->SetBinError(icent+1, fabs(dy) );
     }
     sprintf(name,"hResCent_FCal2_har%d",ihar);
     hResCent_FCal2[ihar]=(TH1D*)hResCent_fg_FCal2[ihar]->Clone(name);
     //hResCent_FCal2[ihar]->Add(hResCent_bg_FCal2[ihar],-1.);
     hResCent_FCal2[ihar]->Reset();

     sprintf(name,"hResCent_Zdc2_har%d",ihar);
     hResCent_Zdc2[ihar]=(TH1D*)hResCent_fg_Zdc2[ihar]->Clone(name);
     //hResCent_Zdc2[ihar]->Add(hResCent_bg_Zdc2[ihar],-1.);
     hResCent_Zdc2[ihar]->Reset();
     
     for (int icent=0; icent<NCENT; icent++) {
       double y,y1,y2,dy,dy1,dy2;
       y1=hCorrCent_FCal2[ihar]->GetBinContent(icent+1);
       dy1=hCorrCent_FCal2[ihar]->GetBinError(icent+1);
       y=0.,dy=0.;
       y=sqrt(fabs(y1));
       if(y1!=0.) dy=0.5*pow(fabs(y1),-0.5)*dy1;
       hResCent_FCal2[ihar]->SetBinContent(icent+1, y );
       hResCent_FCal2[ihar]->SetBinError(icent+1, fabs(dy) );

       y1=hCorrCent_Zdc2[ihar]->GetBinContent(icent+1);
       dy1=hCorrCent_Zdc2[ihar]->GetBinError(icent+1);
       y=0.,dy=0.;
       y=sqrt(fabs(y1));
       if(y1!=0.)dy=0.5*pow(fabs(y1),-0.5)*dy1;
       hResCent_Zdc2[ihar]->SetBinContent(icent+1, y );
       hResCent_Zdc2[ihar]->SetBinError(icent+1, fabs(dy) );
     }

  }
  
  
  //Plot hPsi
  gStyle->SetOptStat(0);
  const int NCan1=5;
  const int NHar1=1;
  TCanvas *c1[NHar1][NCan1];
  TLegend *leg1[NHar1][NCan1];
  TPad *pad1[NHar1][NCan1][6];
  TLegend *leg;
  for(int ihar=0;ihar<NHar1;ihar++){
    for(int i=0;i<NCan1;i++){
      sprintf(name,"c1_%d_%d",ihar,i);
      c1[ihar][i]=new TCanvas(name,name,1100,700);
      leg1[ihar][i]=new TLegend(0.45,0.55,0.88,0.88);
      leg1[ihar][i]->SetLineColor(0);
    }
  }
  //const int icent_mat[]={0,4,6,8,10,15};
  int icent_mat[]={0,2,3,4,5,7};
  for(int ihar=0;ihar<NHar1;ihar++){
    //Plot for ZDC-sideC
    int i=0;
    Divide_Pad(c1[ihar][i],pad1[ihar][i],2,3);
    for(int j=0;j<6;j++){
      c1[ihar][i]->cd();
      pad1[ihar][i][j]->cd();
      hPsi[icent_mat[j]][0][ihar]->Draw();
      sprintf(name,"#Psi_{%d}^{C}",ihar+1);
      hPsi[icent_mat[j]][0][ihar]->SetXTitle(name);
      hPsi[icent_mat[j]][0][ihar]->GetYaxis()->SetRangeUser(0.,3.);
      hPsiC[icent_mat[j]][0][ihar]->Draw("same");
      hPsiF[icent_mat[j]][0][ihar]->Draw("same");
      if(j==0){
	leg1[ihar][i]->AddEntry(hPsi[icent_mat[j]][0][ihar],"Raw distribution","P");
	leg1[ihar][i]->AddEntry(hPsiC[icent_mat[j]][0][ihar],"after Recentering","P");
	leg1[ihar][i]->AddEntry(hPsiF[icent_mat[j]][0][ihar],"after Flattening","P");
	leg1[ihar][i]->Draw("same");
	text.DrawLatex(0.2, 0.7,"Side C");
      }
      sprintf(name,"%s",Cent_mat[icent_mat[j]].c_str());
      text.DrawLatex(0.2, 0.8,name);
    }
  
    //Plot for ZDC-sideA
    i=1;
    Divide_Pad(c1[ihar][i],pad1[ihar][i],2,3);
    for(int j=0;j<6;j++){
      c1[ihar][i]->cd();
      pad1[ihar][i][j]->cd();
      hPsi[icent_mat[j]][1][ihar]->Draw();
      sprintf(name,"#Psi_{%d}^{A}",ihar+1);
      hPsi[icent_mat[j]][1][ihar]->SetXTitle(name);
      hPsi[icent_mat[j]][1][ihar]->GetYaxis()->SetRangeUser(0.,3.);
      hPsiC[icent_mat[j]][1][ihar]->Draw("same");
      hPsiF[icent_mat[j]][1][ihar]->Draw("same");
      if(j==0){
	leg1[ihar][i]->AddEntry(hPsi[icent_mat[j]][1][ihar],"Raw distribution","P");
	leg1[ihar][i]->AddEntry(hPsiC[icent_mat[j]][1][ihar],"after Recentering","P");
	leg1[ihar][i]->AddEntry(hPsiF[icent_mat[j]][1][ihar],"after Flattening","P");
	leg1[ihar][i]->Draw("same");
	text.DrawLatex(0.2, 0.7,"Side A");
      }
      sprintf(name,"%s",Cent_mat[icent_mat[j]].c_str());
      text.DrawLatex(0.2, 0.8,name);
    }
    //Plot for ZDC-Side A+C
    i=2;
    Divide_Pad(c1[ihar][i],pad1[ihar][i],2,3);
    for(int j=0;j<6;j++){
      c1[ihar][i]->cd();
      pad1[ihar][i][j]->cd();
      hPsi[icent_mat[j]][2][ihar]->Draw();
      sprintf(name,"#Psi_{%d}^{A+C}",ihar+1);
      hPsi[icent_mat[j]][2][ihar]->SetXTitle(name);
      hPsi[icent_mat[j]][2][ihar]->GetYaxis()->SetRangeUser(0.,3.);
      hPsiC[icent_mat[j]][2][ihar]->Draw("same");
      hPsiF[icent_mat[j]][2][ihar]->Draw("same");
      if(j==0){
	leg1[ihar][i]->AddEntry(hPsi[icent_mat[j]][2][ihar],"Raw distribution","P");
	leg1[ihar][i]->AddEntry(hPsiC[icent_mat[j]][2][ihar],"after Recentering","P");
	leg1[ihar][i]->AddEntry(hPsiF[icent_mat[j]][2][ihar],"after Flattening","P");
	leg1[ihar][i]->Draw("same");
	text.DrawLatex(0.2, 0.7,"Side A+C");
      }
      sprintf(name,"%s",Cent_mat[icent_mat[j]].c_str());
      text.DrawLatex(0.2, 0.8,name);
    }

    //Plot Correlation ZDC
    i=3;
    Divide_Pad(c1[ihar][i],pad1[ihar][i],2,3);

    for(int j=0;j<6;j++){
      c1[ihar][i]->cd();
      pad1[ihar][i][j]->cd();
      hDphi_Zdc_fg[icent_mat[j]][ihar]->Draw();
      hDphi_Zdc_fg[icent_mat[j]][ihar]->GetXaxis()->SetRangeUser(-PI,PI);
      hDphi_Zdc_fg[icent_mat[j]][ihar]->GetYaxis()->SetRangeUser(0.5,1.5);
      if(ihar==0) sprintf(name,"#Psi_{%d}^{A}-#Psi_{%d}^{C}",ihar+1,ihar+1);
      if(ihar==1){
	sprintf(name,"%d(#Psi_{%d}^{A}-#Psi_{%d}^{C})",ihar+1,ihar,ihar);
	hDphi_Zdc_fg[icent_mat[j]][ihar]->GetYaxis()->SetRangeUser(0.8,1.2);

	leg=new TLegend(0.2,0.15,0.88,0.3);
	leg->SetLineColor(0);
	leg->SetNColumns(2);
	if(j==0){
	  leg->AddEntry(hDphi_Zdc_fg[icent_mat[j]][ihar],"same event","P");
	  leg->AddEntry(hDphi_Zdc_bg[icent_mat[j]][ihar],"mixed event","P");
	  leg->AddEntry(hDphi_Zdc_rat[icent_mat[j]][ihar],"correlation","P");
	}
      }
      hDphi_Zdc_fg[icent_mat[j]][ihar]->SetXTitle(name);
      hDphi_Zdc_bg[icent_mat[j]][ihar]->Draw("same");
      //hDphi_Zdc_rat[icent_mat[j]][ihar]->Draw("same");
      if(j==0){
	leg1[ihar][i]->AddEntry(hDphi_Zdc_fg[icent_mat[j]][ihar],"same event","P");
	leg1[ihar][i]->AddEntry(hDphi_Zdc_bg[icent_mat[j]][ihar],"mixed event","P");
	leg1[ihar][i]->AddEntry(hDphi_Zdc_rat[icent_mat[j]][ihar],"correlation","P");
	if(ihar==0) leg1[ihar][i]->Draw("same");
	if(ihar==1) leg->Draw("same");
	//text.DrawLatex(0.2, 0.8,"Side A");
      }
      sprintf(name,"%s",Cent_mat[icent_mat[j]].c_str());
      text.DrawLatex(0.2, 0.8,name);
    }
    
    i=4;
    Divide_Pad(c1[ihar][i],pad1[ihar][i],2,3);
    for(int j=0;j<6;j++){
      c1[ihar][i]->cd();
      pad1[ihar][i][j]->cd();
      hDphi_FCal_fg[icent_mat[j]][ihar]->Draw();
      hDphi_FCal_fg[icent_mat[j]][ihar]->GetXaxis()->SetRangeUser(-PI,PI);
      hDphi_FCal_fg[icent_mat[j]][ihar]->GetYaxis()->SetRangeUser(0.5,1.5);
      sprintf(name,"%d(#Psi_{%d}^{A}-#Psi_{%d}^{C})",ihar+2,ihar+2,ihar+2);
      hDphi_FCal_fg[icent_mat[j]][ihar]->SetXTitle(name);
      hDphi_FCal_bg[icent_mat[j]][ihar]->Draw("same");
      //hDphi_FCal_rat[icent_mat[j]][ihar]->Draw("same");
      if(j==0){
	leg1[ihar][i]->AddEntry(hDphi_FCal_fg[icent_mat[j]][ihar],"same event","P");
	leg1[ihar][i]->AddEntry(hDphi_FCal_bg[icent_mat[j]][ihar],"mixed event","P");
	leg1[ihar][i]->AddEntry(hDphi_FCal_rat[icent_mat[j]][ihar],"correlation","P");
	leg1[ihar][i]->Draw("same");
	//text.DrawLatex(0.2, 0.8,"Side A");
      }
      sprintf(name,"%s",Cent_mat[icent_mat[j]].c_str());
      text.DrawLatex(0.2, 0.8,name);
    }
    
    //Save Canvases
    sprintf(name,"./plots/Psi%d_sideC.pdf",ihar+1);
    c1[ihar][0]->SaveAs(name);
    sprintf(name,"./plots/Psi%d_sideA.pdf",ihar+1);
    c1[ihar][1]->SaveAs(name);
    sprintf(name,"./plots/Psi%d_sideA+C.pdf",ihar+1);
    c1[ihar][2]->SaveAs(name);
    sprintf(name,"./plots/Zdc_dPsi%d.pdf",ihar+1);
    c1[ihar][3]->SaveAs(name);
    sprintf(name,"./plots/FCal_dPsi%d.pdf",ihar+2);
    c1[ihar][4]->SaveAs(name);
    
  }

  
  //Plot hRes
  int ihar=0;
  TCanvas *c2[4];
  for(int i=0;i<4;i++){
    sprintf(name,"c2_%d",i);
    c2[i]=new TCanvas(name,"",1200,600);
    c2[i]->Divide(2);
  }
  TLegend *leg2[4];
  for(int i=0;i<4;i++){
    leg2[i]=new TLegend(0.2,0.3,0.5,0.5);
    if(i==1) leg2[i]=new TLegend(0.55,0.3,0.88,0.5);
    leg2[i]->SetLineColor(0);
  }
  
  setStyle1(hCorrCent_FCal[0],2);
  setStyle1(hCorrCent_FCal2[0],2);
  setStyle1(hCorrCent_fg_FCal[0],1);
  setStyle1(hCorrCent_fg_FCal2[0],2);
  setStyle1(hCorrCent_bg_FCal[0],1);
  setStyle1(hCorrCent_bg_FCal2[0],2);

  setStyle1(hCorrCent_Zdc[0],2);
  setStyle1(hCorrCent_Zdc2[0],2);
  setStyle1(hCorrCent_fg_Zdc[0],1);
  setStyle1(hCorrCent_fg_Zdc2[0],2);
  setStyle1(hCorrCent_bg_Zdc[0],1);
  setStyle1(hCorrCent_bg_Zdc2[0],2);

  
  setStyle1(hResCent_FCal[0],2);
  setStyle1(hResCent_FCal2[0],2);
  setStyle1(hResCent_fg_FCal[0],1);
  setStyle1(hResCent_fg_FCal2[0],2);
  setStyle1(hResCent_bg_FCal[0],1);
  setStyle1(hResCent_bg_FCal2[0],2);

  setStyle1(hResCent_Zdc[0],2);
  setStyle1(hResCent_Zdc2[0],2);
  setStyle1(hResCent_fg_Zdc[0],1);
  setStyle1(hResCent_fg_Zdc2[0],2);
  setStyle1(hResCent_bg_Zdc[0],1);
  setStyle1(hResCent_bg_Zdc2[0],2);
  

  //Correlation
  c2[0]->cd(1);
  hCorrCent_fg_FCal[0]->Draw();
  hCorrCent_FCal[0]->Draw("same");
  hCorrCent_FCal[0]->SetMarkerSize(1.8);
  hCorrCent_fg_FCal[0]->SetXTitle("Centrality");
  hCorrCent_fg_FCal[0]->SetYTitle("<cos2(#Psi_{2}^{A}-#Psi_{2}^{C})>");
  hCorrCent_fg_FCal[0]->GetXaxis()->SetTitleOffset(1.5);
  hCorrCent_fg_FCal[0]->GetYaxis()->SetTitleOffset(1.8);
  leg2[0]->AddEntry(hCorrCent_fg_FCal[0],"Same Ev","P");
  leg2[0]->AddEntry(hCorrCent_FCal[0],"Same-Mixed","P");
  leg2[0]->Draw("same");
    
  c2[0]->cd(2);
  hCorrCent_bg_FCal[0]->Draw();
  hCorrCent_bg_FCal2[0]->Draw("same");
  hCorrCent_bg_FCal2[0]->SetMarkerSize(1.8);
  //hCorrCent_bg_FCal[0]->GetYaxis()->SetRangeUser(0.,0.06);
  hCorrCent_bg_FCal[0]->SetXTitle("Centrality");
  hCorrCent_bg_FCal[0]->SetYTitle("Mixed Ev BG");
  hCorrCent_bg_FCal[0]->GetXaxis()->SetTitleOffset(1.5);
  hCorrCent_bg_FCal[0]->GetYaxis()->SetTitleOffset(1.8);

  c2[1]->cd(1);
  hCorrCent_fg_Zdc[0]->Draw();
  hCorrCent_Zdc[0]->Draw("same");
  hCorrCent_Zdc[0]->SetMarkerSize(1.8);
  hCorrCent_fg_Zdc[0]->GetYaxis()->SetRangeUser(-0.035,0.);
  hCorrCent_fg_Zdc[0]->SetXTitle("Centrality");
  hCorrCent_fg_Zdc[0]->SetYTitle("<cos(#Psi_{1}^{A}-#Psi_{1}^{C})>");
  hCorrCent_fg_Zdc[0]->GetXaxis()->SetTitleOffset(1.5);
  hCorrCent_fg_Zdc[0]->GetYaxis()->SetTitleOffset(1.8);
  leg2[1]->AddEntry(hCorrCent_fg_Zdc[0],"Same Ev","P");
  leg2[1]->AddEntry(hCorrCent_Zdc[0],"Same-Mixed","P");
  leg2[1]->Draw("same");
  
  c2[1]->cd(2);
  hCorrCent_bg_Zdc[0]->Draw();
  hCorrCent_bg_Zdc2[0]->Draw("same");
  hCorrCent_bg_Zdc2[0]->SetMarkerSize(1.8);
  //hCorrCent_bg_Zdc[0]->GetYaxis()->SetRangeUser(0.,0.015);
  hCorrCent_bg_Zdc[0]->SetXTitle("Centrality");
  hCorrCent_bg_Zdc[0]->SetYTitle("Mixed Ev BG");
  hCorrCent_bg_Zdc[0]->GetXaxis()->SetTitleOffset(1.5);
  hCorrCent_bg_Zdc[0]->GetYaxis()->SetTitleOffset(1.8);

  //Resolution
  c2[2]->cd(1);
  hResCent_fg_FCal[0]->Draw();
  hResCent_FCal[0]->Draw("same");
  hResCent_FCal[0]->SetMarkerSize(1.8);
  hResCent_fg_FCal[0]->SetXTitle("Centrality");
  hResCent_fg_FCal[0]->SetYTitle("Resolution");
  hResCent_fg_FCal[0]->GetXaxis()->SetTitleOffset(1.5);
  hResCent_fg_FCal[0]->GetYaxis()->SetTitleOffset(1.8);
  leg2[2]->AddEntry(hResCent_fg_FCal[0],"Same Ev","P");
  leg2[2]->AddEntry(hResCent_FCal[0],"Same-Mixed","P");
  leg2[2]->Draw("same");

  c2[2]->cd(2);
  hResCent_bg_FCal[0]->Draw();
  hResCent_bg_FCal2[0]->Draw("same");
  hResCent_bg_FCal2[0]->SetMarkerSize(1.8);
  hResCent_bg_FCal[0]->SetXTitle("Centrality");
  hResCent_bg_FCal[0]->SetYTitle("Mixed Ev BG");
  hResCent_bg_FCal[0]->GetXaxis()->SetTitleOffset(1.5);
  hResCent_bg_FCal[0]->GetYaxis()->SetTitleOffset(1.8);

  c2[3]->cd(1);
  hResCent_fg_Zdc[0]->Draw();
  hResCent_Zdc[0]->Draw("same");
  hResCent_Zdc[0]->SetMarkerSize(1.8);
  hResCent_fg_Zdc[0]->SetXTitle("Centrality");
  hResCent_fg_Zdc[0]->SetYTitle("Resolution");
  hResCent_fg_Zdc[0]->GetXaxis()->SetTitleOffset(1.5);
  hResCent_fg_Zdc[0]->GetYaxis()->SetTitleOffset(1.8);
  leg2[3]->AddEntry(hResCent_fg_Zdc[0],"Same Ev","P");
  leg2[3]->AddEntry(hResCent_Zdc[0],"Same-Mixed","P");
  leg2[3]->Draw("same");

  c2[3]->cd(2);
  hResCent_bg_Zdc[0]->Draw();
  hResCent_bg_Zdc2[0]->Draw("same");
  hResCent_bg_Zdc2[0]->SetMarkerSize(1.8);
  hResCent_bg_Zdc[0]->SetXTitle("Centrality");
  hResCent_bg_Zdc[0]->SetYTitle("Mixed Ev BG");
  hResCent_bg_Zdc[0]->GetXaxis()->SetTitleOffset(1.5);
  hResCent_bg_Zdc[0]->GetYaxis()->SetTitleOffset(1.8);

  sprintf(name,"./plots/Corr_FCal_v%d.pdf",ihar+2);
  c2[0]->SaveAs(name);
  sprintf(name,"./plots/Corr_Zdc_v%d.pdf",ihar+1);
  c2[1]->SaveAs(name);
  sprintf(name,"./plots/Res_FCal_v%d.pdf",ihar+2);
  c2[2]->SaveAs(name);
  sprintf(name,"./plots/Res_Zdc_v%d.pdf",ihar+1);
  c2[3]->SaveAs(name);
  
}

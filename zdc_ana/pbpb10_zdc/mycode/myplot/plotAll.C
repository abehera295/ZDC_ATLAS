#include "./plot.h"
#include "./plot1_Psi.C"
#include "./plot2_Res.C"
#include "./plot3_v1_eta.C"

void plotAll(){
  plot *a=new plot();
  a->init();
  a->plot1_Psi();
  a->plot2_Res();
  a->plot3_v1_eta();
}

void plot::init(){
  fin=new TFile("../Output_zdcPix_flat_wtrk_local_wv2_sin.root");
  
  //Read Psi Histograms
  for (int icent=0; icent<NCENT; icent++) {
    for (int iside=0; iside<NSIDE+1; iside++) {
      for (int ihar=0; ihar<NHAR; ihar++) {
	sprintf(name, "h_Psi_cent%d_side%d_har%d", icent, iside, ihar);
	hPsi[icent][iside][ihar]=(TH1D*)fin->Get(name);
	sprintf(name, "h_PsiC_cent%d_side%d_har%d", icent, iside, ihar);
	hPsiC[icent][iside][ihar]=(TH1D*)fin->Get(name);
	sprintf(name, "h_PsiF_cent%d_side%d_har%d", icent, iside, ihar);
	hPsiF[icent][iside][ihar]=(TH1D*)fin->Get(name);
	//Scale
	hPsi[icent][iside][ihar]->Scale(1./hPsi[icent][iside][ihar]->Integral());
        hPsiC[icent][iside][ihar]->Scale(1./hPsiC[icent][iside][ihar]->Integral());
        hPsiF[icent][iside][ihar]->Scale(1./hPsiF[icent][iside][ihar]->Integral());
	double mean0=FindAvg(hPsiF[icent][iside][ihar]);
	hPsi[icent][iside][ihar]->Scale(1./mean0);
	hPsiC[icent][iside][ihar]->Scale(1./mean0);
	hPsiF[icent][iside][ihar]->Scale(1./mean0);
	//Style
	setStyle1(hPsi[icent][iside][ihar],0);
	setStyle1(hPsiC[icent][iside][ihar],1);
	setStyle1(hPsiF[icent][iside][ihar],2);

      }
    }
  }

  for (int icent=0; icent<NCENT; icent++) {
    for (int ihar=0; ihar<NHAR; ihar++) {
      sprintf(name, "h_Dphi_fg_cent%d_har%d", icent, ihar);
      hDphi_fg[icent][ihar]=(TH1D*)fin->Get(name);
      sprintf(name, "h_Dphi_bg_cent%d_har%d", icent, ihar);
      hDphi_bg[icent][ihar]=(TH1D*)fin->Get(name);
      sprintf(name, "h_Dphi_rat_cent%d_har%d", icent, ihar);
      hDphi_rat[icent][ihar]=(TH1D*)fin->Get(name);
      //Rebin
      hDphi_fg[icent][ihar]->Rebin(8);
      hDphi_bg[icent][ihar]->Rebin(8);
      hDphi_rat[icent][ihar]->Rebin(8);
      //Scale
      hDphi_fg[icent][ihar]->Scale(1./hDphi_fg[icent][ihar]->Integral());
      hDphi_bg[icent][ihar]->Scale(1./hDphi_bg[icent][ihar]->Integral());
      hDphi_rat[icent][ihar]->Scale(1./hDphi_rat[icent][ihar]->Integral());
      double mean0=FindAvg(hDphi_bg[icent][ihar]);
      hDphi_fg[icent][ihar]->Scale(1./mean0);
      hDphi_bg[icent][ihar]->Scale(1./mean0);
      hDphi_rat[icent][ihar]->Scale(1./mean0);
      //Style
      setStyle1(hDphi_fg[icent][ihar],0);
      setStyle1(hDphi_bg[icent][ihar],1);
      setStyle1(hDphi_rat[icent][ihar],2);
    }
  }
  
  //Read Res histograms
  for (int ihar=0; ihar<NHAR; ihar++) {
    sprintf(name,"h_ResCent_har%d",ihar);
    hResCent[ihar]=(TH1D*)fin->Get(name);
    sprintf(name,"h_ResCent_recn_har%d",ihar);
    hResCent_recn[ihar]=(TH1D*)fin->Get(name);
    sprintf(name,"h_ResCent_flat_har%d",ihar);
    hResCent_flat[ihar]=(TH1D*)fin->Get(name);
    sprintf(name,"h_ResCent_mix_har%d",ihar);
    hResCent_mix[ihar]=(TH1D*)fin->Get(name);

 }

  fin2=new TFile("../Output_zdcPix_flat_wtrk_local_wv2_doAna_reb.root","read");
  //Read vn Histograms
  for (int icent=0; icent<NCENT_NEW; icent++) {
    for (int iside=0; iside<NSIDE+1; iside++) {
      sprintf(name , "hV1_eta_res_cent%d_side%d_pt10", icent, iside);
      hV1_eta_res[icent][iside]=(TH1D*)fin2->Get(name);
      for (int ihar=0; ihar<NHAR; ihar++) {
	for (int ipt=0; ipt<NPT+1; ipt++) {
	  sprintf(name, "hV%d_eta_raw_cent%d_side%d_pt%d", ihar+1, icent, iside, ipt);
	  hVn_eta_raw[icent][iside][ihar][ipt]=(TH1D*)fin2->Get(name);
	  sprintf(name, "hV%dS_eta_raw_cent%d_side%d_pt%d", ihar+1, icent, iside, ipt);
	  hVnS_eta_raw[icent][iside][ihar][ipt]=(TH1D*)fin2->Get(name);

	  sprintf(name, "hV%d_eta_res_cent%d_side%d_pt%d", ihar+1, icent, iside, ipt);
	  hVn_eta_res[icent][iside][ihar][ipt]=(TH1D*)fin2->Get(name);
	  sprintf(name, "hV%dS_eta_res_cent%d_side%d_pt%d", ihar+1, icent, iside, ipt);
	  hVnS_eta_res[icent][iside][ihar][ipt]=(TH1D*)fin2->Get(name);
	}
      }
    }
  }
  
  gStyle->SetOptStat(0);

  text.SetNDC(1);
  text.SetTextFont(43);
  text.SetTextSize(18);
}


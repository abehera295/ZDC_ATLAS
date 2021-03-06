const int NCENT2=5;
const double cent_mat2[NCENT2+1]={0.,10.,20.,30.,40.,60.};
TF1 *func_FCal[NCENT2][NPT_Mu2][NSIDE+1];
TF1 *func_Zdc[NCENT2][NETA_Mu][NSIDE+1];
double vn_mat[50][50];
double vn_e_mat[50][50];
TH1D *hV2_FCal_pt[NCENT2][NSIDE+1],*hV2_FCal_cent[NPT_Mu2][NSIDE+1];
TH1D *hV1_Zdc_eta[NCENT2][NSIDE+1],*hV1_Zdc_cent[NETA_Mu][NSIDE+1];

TProfile *hResCent_FCal_fg[NHAR],*hResCent_FCal_bg[NHAR];
TProfile *hResCent_Zdc_fg[NHAR],*hResCent_Zdc_bg[NHAR];
TH1D *hResCent_FCal[NHAR],*hResCent_Zdc[NHAR];


void setStyle2(TH1D *h,int itype1){
  h->GetYaxis()->CenterTitle(0);       h->GetXaxis()->CenterTitle(0);  //h->GetZaxis()->CenterTitle();
  h->GetYaxis()->SetTitleOffset(2.); h->GetXaxis()->SetTitleOffset(2.);
  //h->GetZaxis()->SetTitleOffset(1.5);
  h->GetXaxis()->SetTitleFont(43);    h->GetXaxis()->SetTitleSize(18);
  h->GetXaxis()->SetLabelFont(43);    h->GetXaxis()->SetLabelSize(15);
  h->GetYaxis()->SetTitleFont(43);    h->GetYaxis()->SetTitleSize(18);
  h->GetYaxis()->SetLabelFont(43);    h->GetYaxis()->SetLabelSize(15);
  //h->GetZaxis()->SetTitleFont(43);    h->GetZaxis()->SetTitleSize(20);
  //h->GetZaxis()->SetLabelFont(43);    h->GetZaxis()->SetLabelSize(15);
  h->GetYaxis()->SetNdivisions(505);
  h->GetXaxis()->SetNdivisions(505);
  int c_mat[]={4,1,2,4,6,28,36,7,8,9,46,30};
  int m_mat[]={20,30,27,25,28,33,20,23,29,34};
  h->SetTitle("");
  h->SetMarkerStyle(m_mat[itype1]);
  h->SetMarkerSize(1.);
  h->SetMarkerColor(c_mat[itype1]);
  h->SetLineColor(c_mat[itype1]);
  //h->GetYaxis()->SetRangeUser(0.,5.);
}

void fixTitle2(TH1D *h){
  h->GetYaxis()->SetTitleOffset(1.2); h->GetXaxis()->SetTitleOffset(1.2);
  h->GetYaxis()->SetTitleSize(20.); h->GetXaxis()->SetTitleSize(20.);
}


void plot::plot2_fit(){
  plot2_readHists();
  plot2_doFit();
  plot2_drawFit();
}

void Subt_bg(TH1D *h,TProfile *hfg,TProfile *hbg){
  int N=h->GetNbinsX();
  double x1,dx1,y1,dy1,x2,dx2,y2,dy2,x,dx,y,dy;
  for(int ibin=0;ibin<N;ibin++){
    y1=hfg->GetBinContent(ibin+1);
    dy1=hfg->GetBinError(ibin+1);
    y2=hbg->GetBinContent(ibin+1);
    dy2=hbg->GetBinError(ibin+1);
    y=y1-y2;
    dy=sqrt(pow(dy1,2)+pow(dy2,2));
    h->SetBinContent(ibin+1,y);
    h->SetBinError(ibin+1,dy);
  }
  
}

void plot::plot2_readHists(){
  TFile *fin=new TFile("./step1.root","read");
  for(int icent=0;icent<NCENT2;icent++){
    for (int iside=0; iside<NSIDE+1; iside++){
      for (int ihar=0; ihar<NHAR; ihar++){
	for(int ipt=0;ipt<NPT_Mu2;ipt++){
	  sprintf(name,"hMu_dphi_yield_FCal2_cent%d_side%d_pt%d_har%d",icent,iside,ipt,ihar);
	  hMu_dphi_yield_FCal2[icent][iside][ipt][ihar]=(TH1D*)fin->Get(name);
	}
	for(int ieta=0;ieta<NETA_Mu;ieta++){
	  sprintf(name,"hMu_dphi_yield_Zdc_cent%d_side%d_eta%d_har%d",icent,iside,ieta,ihar);
	  hMu_dphi_yield_Zdc[icent][iside][ieta][ihar]=(TH1D*)fin->Get(name);
	}
      }
    }
  }
  TFile *fin2=new TFile("../zdcvn/Output_res.root","read");
  for (int ihar=0; ihar<NHAR; ihar++) {
    sprintf(name, "hResCent_FCal_fg_har%d", ihar);
    hResCent_FCal_fg[ihar]=(TProfile*)fin2->Get(name);
    sprintf(name, "hResCent_FCal_bg_har%d", ihar);
    hResCent_FCal_bg[ihar]=(TProfile*)fin2->Get(name);

    sprintf(name, "hResCent_Zdc_fg_har%d", ihar);
    hResCent_Zdc_fg[ihar]=(TProfile*)fin2->Get(name);
    sprintf(name, "hResCent_Zdc_bg_har%d", ihar);
    hResCent_Zdc_bg[ihar]=(TProfile*)fin2->Get(name);

    sprintf(name, "hResCent_FCal_har%d", ihar);
    hResCent_FCal[ihar]=new TH1D(name,"",NCENT,cent_mat_double);
    sprintf(name, "hResCent_Zdc_har%d", ihar);
    hResCent_Zdc[ihar]=new TH1D(name,"",NCENT,cent_mat_double);
    Subt_bg(hResCent_FCal[ihar],hResCent_FCal_fg[ihar],hResCent_FCal_bg[ihar]);
    Subt_bg(hResCent_Zdc[ihar],hResCent_Zdc_fg[ihar],hResCent_Zdc_bg[ihar]);
  }
  //Make Hists for v2
  for (int iside=0; iside<NSIDE+1; iside++){
    for(int icent=0;icent<NCENT2;icent++){
      sprintf(name,"hV2_FCal_pt_cent%d_side%d",icent,iside);
      hV2_FCal_pt[icent][iside]=new TH1D(name,"",NPT_Mu2,pt_mat_Mu2);
      sprintf(name,"hV1_Zdc_eta_cent%d_side%d",icent,iside);
      hV1_Zdc_eta[icent][iside]=new TH1D(name,"",NETA_Mu,eta_mat_Mu);
    }
    for(int ipt=0;ipt<NPT_Mu2;ipt++){
      sprintf(name,"hV2_FCal_cent_pt%d_side%d",ipt,iside);
      hV2_FCal_cent[ipt][iside]=new TH1D(name,"",NCENT2,cent_mat2);
    }
    for(int ieta=0;ieta<NETA_Mu;ieta++){
      sprintf(name,"hV1_Zdc_cent_eta%d_side%d",ieta,iside);
      hV1_Zdc_cent[ieta][iside]=new TH1D(name,"",NCENT2,cent_mat2);
    }
  }
  
}

void plot::plot2_doFit(){
  int jside,jhar;
  double res,res_e;
  for(int icent=0;icent<NCENT2;icent++){
    for (int iside=0; iside<NSIDE+1; iside++){
      for(int ipt=0;ipt<NPT_Mu2;ipt++){
	jhar=1;
	sprintf(name,"func_FCal_cent%d_pt%d_side%d",icent,ipt,iside);
	func_FCal[icent][ipt][iside]=new TF1(name,"[0]*(1.+2.*[1]*cos(x))",0.,PI);
	hMu_dphi_yield_FCal2[icent][iside][ipt][jhar]->Fit(name,"R");
	hMu_dphi_yield_FCal2[icent][iside][ipt][jhar]->Fit(name,"R");
	double y=hResCent_FCal[jhar]->GetBinContent(icent+1);
	double dy=hResCent_FCal[jhar]->GetBinError(icent+1);
	res=sqrt(fabs(y));
        res_e=0.5*pow(fabs(y),-0.5)*dy;
	if(iside==2){
	  res=sqrt(fabs(2.*y));
	  res_e=0.5*pow(fabs(2.*y),-0.5)*2.*dy;
	}
	//vn_mat[icent][ipt]=func_FCal[icent][ipt][iside]->GetParameter(1);
	//vn_e_mat[icent][ipt]=func_FCal[icent][ipt][iside]->GetParError(1);
	vn_mat[icent][ipt]=func_FCal[icent][ipt][iside]->GetParameter(1)/res;
	vn_e_mat[icent][ipt]=sqrt(pow(func_FCal[icent][ipt][iside]->GetParError(1),2)+pow(res_e,2));
	hV2_FCal_pt[icent][iside]->SetBinContent(ipt+1,vn_mat[icent][ipt]);
	hV2_FCal_pt[icent][iside]->SetBinError(ipt+1,vn_e_mat[icent][ipt]);
	hV2_FCal_cent[ipt][iside]->SetBinContent(icent+1,vn_mat[icent][ipt]);
	hV2_FCal_cent[ipt][iside]->SetBinError(icent+1,vn_e_mat[icent][ipt]);
      }
      for(int ieta=0;ieta<NETA_Mu;ieta++){
	jhar=0;
	sprintf(name,"func_Zdc_cent%d_eta%d_side%d",icent,ieta,iside);
	func_Zdc[icent][ieta][iside]=new TF1(name,"[0]*(1.+2.*[1]*cos(x))",0.,PI);
	func_Zdc[icent][ieta][iside]->SetParLimits(0,0.,1e5);
	func_Zdc[icent][ieta][iside]->SetParLimits(1,-1.,1.);
	hMu_dphi_yield_Zdc[icent][iside][ieta][jhar]->Fit(name,"R");
	hMu_dphi_yield_Zdc[icent][iside][ieta][jhar]->Fit(name,"R");
	double y=hResCent_Zdc[jhar]->GetBinContent(icent+1);
	double dy=hResCent_Zdc[jhar]->GetBinError(icent+1);
	res=sqrt(fabs(y));
        res_e=0.5*pow(fabs(y),-0.5)*dy;
	if(iside==2){
          res=sqrt(fabs(2.*y));
          res_e=0.5*pow(fabs(2.*y),-0.5)*2.*dy;
	}
	//vn_mat[icent][ieta]=func_Zdc[icent][ieta][iside]->GetParameter(1)/res;
	//vn_e_mat[icent][ieta]=sqrt(pow(func_Zdc[icent][ieta][iside]->GetParError(1),2)+pow(res_e,2));
	vn_mat[icent][ieta]=func_Zdc[icent][ieta][iside]->GetParameter(1);
	vn_e_mat[icent][ieta]=func_Zdc[icent][ieta][iside]->GetParError(1);
	hV1_Zdc_eta[icent][iside]->SetBinContent(ieta+1,vn_mat[icent][ieta]);
	hV1_Zdc_eta[icent][iside]->SetBinError(ieta+1,vn_e_mat[icent][ieta]);
	hV1_Zdc_cent[ieta][iside]->SetBinContent(icent+1,vn_mat[icent][ieta]);
	hV1_Zdc_cent[ieta][iside]->SetBinError(icent+1,vn_e_mat[icent][ieta]);
      }
    }
  }

}

void plot::plot2_drawFit(){
  int icent_mat[3]={0,2,4};//0-10%,20-30%,40-60%
  string centP_mat[3]={"0-10%","20-30%","40-60%"};
  int ipt_mat[2]={0,3};
  string ptR_mat[2]={"4< p_{T} < 5 GeV","8 < p_{T} < 10 GeV"};
  int ieta_mat[3]={0,1,2};
  string etaR_mat[3]={"-2 < #eta < -1.5","-1.5 < #eta < -1.","-1 < #eta < 0.5"};

  //Plot Yield
  TCanvas *c2;
  TLegend *leg2[2];
  TPad *pad2[6];
  sprintf(name,"c2");
  c2=new TCanvas(name,name,1200,800);
  c2->Divide(3,2);
  int jdet,jhar;
  //For FCal
  jdet=0;jhar=1;//v2
  //Divide_Pad(c2[jdet],pad2[jdet],2,3);
  for(int i=0;i<3;i++){
    int icent=icent_mat[i];
    int ipt;
    c2->cd(i+1);
    //pad2[jdet][i]->cd();
    ipt=ipt_mat[0];
    setStyle1(hMu_dphi_yield_FCal2[icent][2][ipt][jhar],0);
    //hMu_dphi_yield_FCal2[icent][2][ipt][jhar]->Scale(1./1000.);
    hMu_dphi_yield_FCal2[icent][2][ipt][jhar]->Draw();
    hMu_dphi_yield_FCal2[icent][2][ipt][jhar]->SetXTitle("2(#phi-#Psi_{2}FCal)");
    hMu_dphi_yield_FCal2[icent][2][ipt][jhar]->SetYTitle("Signal Counts");
    //hMu_dphi_yield_FCal2[icent][2][ipt][jhar]->GetYaxis()->SetRangeUser(2000.,15000);
    sprintf(name,"%s",centP_mat[i].c_str());
    text.DrawLatex(0.2,0.4,name);
    sprintf(name,"%s",ptR_mat[0].c_str());
    text.DrawLatex(0.2,0.2,name);
    text.DrawLatex(0.2,0.3,"|#eta|<2");
    
    c2->cd(i+3+1);
    //pad2[jdet][i+3]->cd();
    ipt=ipt_mat[1];
    setStyle1(hMu_dphi_yield_FCal2[icent][2][ipt][jhar],0);
    //hMu_dphi_yield_FCal2[icent][2][ipt][jhar]->Scale(1./1000.);
    hMu_dphi_yield_FCal2[icent][2][ipt][jhar]->Draw();
    hMu_dphi_yield_FCal2[icent][2][ipt][jhar]->SetXTitle("2(#phi-#Psi_{2}FCal)");
    hMu_dphi_yield_FCal2[icent][2][ipt][jhar]->SetYTitle("Signal Counts");
    //hMu_dphi_yield_FCal2[icent][2][ipt][jhar]->GetYaxis()->SetRangeUser(0.,2000);
    sprintf(name,"%s",centP_mat[i].c_str());
    text.DrawLatex(0.2,0.4,name);
    sprintf(name,"%s",ptR_mat[1].c_str());
    text.DrawLatex(0.2,0.2,name);
    text.DrawLatex(0.2,0.3,"|#eta|<2");
  }
  
  TCanvas *c3[2];
  for(int ican=0;ican<2;ican++){
    sprintf(name,"c3_%d",ican);
    c3[ican]=new TCanvas(name,name,1200,800);
    c3[ican]->Divide(2,2);
  }
  
  int jside;
  //For Zdc
  jside=0;jhar=0;//v1
  //Divide_Pad(c2[jdet],pad2[jdet],2,3);
  for(int i=0;i<4;i++){
    int ieta=i;
    int icent;
    c3[jside]->cd(i+1);
    //pad2[jdet][i]->cd();
    icent=icent_mat[0];
    setStyle1(hMu_dphi_yield_Zdc[icent][jside][ieta][jhar],0);
    //hMu_dphi_yield_Zdc[icent][jside][ieta][jhar]->Scale(1./1000.);
    hMu_dphi_yield_Zdc[icent][jside][ieta][jhar]->Draw();
    hMu_dphi_yield_Zdc[icent][jside][ieta][jhar]->SetXTitle("#phi-#Psi_{1}Zdc");
    hMu_dphi_yield_Zdc[icent][jside][ieta][jhar]->SetYTitle("Signal Counts");
    //hMu_dphi_yield_Zdc[icent][2][ieta][jhar]->GetYaxis()->SetRangeUser(2000.,15000);
    sprintf(name,"%1.1f < #eta < %1.1f",eta_mat_Mu[i],eta_mat_Mu[i+1]);
    text.DrawLatex(0.2,0.4,name);
    sprintf(name,"%s",centP_mat[0].c_str());
    text.DrawLatex(0.2,0.2,name);
    text.DrawLatex(0.2,0.3,"p_{T}>4");
  }
  
  jside=1;jhar=0;//v1
  for(int i=0;i<4;i++){
    int ieta=i;
    int icent;
    c3[jside]->cd(i+1);
    //pad2[jdet][i]->cd();
    icent=icent_mat[0];
    setStyle1(hMu_dphi_yield_Zdc[icent][jside][ieta][jhar],0);
    //hMu_dphi_yield_Zdc[icent][jside][ieta][jhar]->Scale(1./1000.);
    hMu_dphi_yield_Zdc[icent][jside][ieta][jhar]->Draw();
    hMu_dphi_yield_Zdc[icent][jside][ieta][jhar]->SetXTitle("#phi-#Psi_{1}Zdc");
    hMu_dphi_yield_Zdc[icent][jside][ieta][jhar]->SetYTitle("Signal Counts");
    //hMu_dphi_yield_Zdc[icent][2][ieta][jhar]->GetYaxis()->SetRangeUser(2000.,15000);
    sprintf(name,"%1.1f < #eta < %1.1f",eta_mat_Mu[i],eta_mat_Mu[i+1]);
    text.DrawLatex(0.2,0.4,name);
    sprintf(name,"%s",centP_mat[0].c_str());
    text.DrawLatex(0.2,0.2,name);
    text.DrawLatex(0.2,0.3,"p_{T}>4");
  }
  
  c2->SaveAs("./plots/Yield_FCal_fit.pdf");
  c3[0]->SaveAs("./plots/Yield_Zdc_sideC_fit.pdf");
  c3[1]->SaveAs("./plots/Yield_Zdc_sideA_fit.pdf");

  
  //Plot v2
  TCanvas *c4[2];
  TLegend *leg3[2];
  TPad *pad3[2][6];
  for(int ican=0;ican<2;ican++){
    sprintf(name,"c4_%d",ican);
    c4[ican]=new TCanvas(name,name,1000,800);
    c4[ican]->Divide(2,2);
  }
  //For FCal
  jside=2;
  setStyle2(hV2_FCal_pt[0][jside],0);
  setStyle2(hV2_FCal_pt[2][jside],0);
  setStyle2(hV2_FCal_cent[0][jside],0);
  setStyle2(hV2_FCal_cent[2][jside],0);
  c4[0]->cd(1);
  hV2_FCal_pt[0][jside]->Draw();
  hV2_FCal_pt[0][jside]->GetYaxis()->SetRangeUser(0.,0.1);
  hV2_FCal_pt[0][jside]->SetXTitle("p_{T}");
  hV2_FCal_pt[0][jside]->SetYTitle("v_{2}");
  text.DrawLatex(0.2,0.3,"0-10%");
  c4[0]->cd(2);
  hV2_FCal_pt[2][jside]->Draw();
  hV2_FCal_pt[2][jside]->GetYaxis()->SetRangeUser(0.,0.1);
  hV2_FCal_pt[2][jside]->SetXTitle("p_{T}");
  hV2_FCal_pt[2][jside]->SetYTitle("v_{2}");
  text.DrawLatex(0.2,0.3,"20-30%");
  c4[0]->cd(3);
  hV2_FCal_cent[0][jside]->Draw();
  hV2_FCal_cent[0][jside]->GetYaxis()->SetRangeUser(0.,0.1);
  hV2_FCal_cent[0][jside]->SetXTitle("Centrality");
  hV2_FCal_cent[0][jside]->SetYTitle("v_{2}");
  text.DrawLatex(0.6,0.3,"4 < p_{T} < 5");
  c4[0]->cd(4);
  hV2_FCal_cent[2][jside]->Draw();
  hV2_FCal_cent[2][jside]->GetYaxis()->SetRangeUser(0.,0.1);
  hV2_FCal_cent[2][jside]->SetXTitle("Centrality");
  hV2_FCal_cent[2][jside]->SetYTitle("v_{2}");
  text.DrawLatex(0.6,0.3,"6 < p_{T} < 8");
  
  //For Zdc
  jside=0;//Side-C
  setStyle2(hV1_Zdc_eta[0][jside],0);
  setStyle2(hV1_Zdc_eta[2][jside],0);
  c4[1]->cd(1);
  hV1_Zdc_eta[0][jside]->Draw();
  hV1_Zdc_eta[0][jside]->GetYaxis()->SetRangeUser(-0.01,0.01);
  hV1_Zdc_eta[0][jside]->SetXTitle("#eta");
  hV1_Zdc_eta[0][jside]->SetYTitle("v_{1}");
  text.DrawLatex(0.2,0.8,"SideC");
  text.DrawLatex(0.2,0.75,"0-10%");
  c4[1]->cd(2);
  hV1_Zdc_eta[2][jside]->Draw();
  hV1_Zdc_eta[2][jside]->GetYaxis()->SetRangeUser(-0.01,0.01);
  hV1_Zdc_eta[2][jside]->SetXTitle("#eta");
  hV1_Zdc_eta[2][jside]->SetYTitle("v_{1}");
  text.DrawLatex(0.2,0.8,"SideC");
  text.DrawLatex(0.2,0.75,"20-30%");

  jside=1;//Side-A
  setStyle2(hV1_Zdc_eta[0][jside],0);
  setStyle2(hV1_Zdc_eta[2][jside],0);
  c4[1]->cd(3);
  hV1_Zdc_eta[0][jside]->Draw();
  hV1_Zdc_eta[0][jside]->GetYaxis()->SetRangeUser(-0.01,0.01);
  hV1_Zdc_eta[0][jside]->SetXTitle("Centrality");
  hV1_Zdc_eta[0][jside]->SetYTitle("v_{1}");
  text.DrawLatex(0.2,0.8,"SideA");
  text.DrawLatex(0.2,0.75,"0-10%");
  c4[1]->cd(4);
  hV1_Zdc_eta[2][jside]->Draw();
  hV1_Zdc_eta[2][jside]->GetYaxis()->SetRangeUser(-0.01,0.01);
  hV1_Zdc_eta[2][jside]->SetXTitle("Centrality");
  hV1_Zdc_eta[2][jside]->SetYTitle("v_{1}");
  text.DrawLatex(0.2,0.8,"SideA");
  text.DrawLatex(0.2,0.75,"20-30%");
  
  c4[0]->SaveAs("./plots/v2_FCal.pdf");
  c4[1]->SaveAs("./plots/v1_Zdc.pdf");
  
  //Write hists to root file
  TFile *fout1=new TFile("fit.root","recreate");
  fout1->cd();
  int jside=2;
  for(int icent=0;icent<NCENT2;icent++){
    for(int ipt=0;ipt<NPT_Mu2;ipt++){
      hMu_dphi_yield_FCal2[icent][jside][ipt][1]->Write();
    }
    for(int ieta=0;ieta<NETA_Mu;ieta++){
      hMu_dphi_yield_Zdc[icent][jside][ieta][0]->Write();
    }
  }
  fout1->Close();

  TFile *fout2=new TFile("vn.root","recreate");
  fout2->cd();
  for (int iside=0; iside<NSIDE+1; iside++){
    for(int icent=0;icent<NCENT2;icent++){
      hV2_FCal_pt[icent][iside]->Write();
      hV1_Zdc_eta[icent][iside]->Write();
    }
    for(int ipt=0;ipt<NPT_Mu2;ipt++){
      hV2_FCal_cent[ipt][iside]->Write();
    }
    for(int ieta=0;ieta<NETA_Mu;ieta++){
      hV1_Zdc_cent[ieta][iside]->Write();
    }
  }
  fout2->Close();

}

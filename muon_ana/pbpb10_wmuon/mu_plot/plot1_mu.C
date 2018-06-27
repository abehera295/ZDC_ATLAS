//Recombine Cent bins and pT bins
const int NPT_Mu2=5;//4-5,5-6,6-8,8-10,10-14
const double pt_mat_Mu2[NPT_Mu2+1]={4.,5.,6.,8.,10.,14.};

TH1D *hMu_dphi_yield_FCal2[NCENT][NSIDE+1][NPT_Mu2][NHAR];
TH2F *hMu_dphi_eloss_FCal2[NCENT][NSIDE+1][NPT_Mu2][NHAR];
TH1D *hMu_eloss_FCal2[NCENT][NSIDE+1][NPT_Mu2][NHAR][4];
//TH1D *hMu_dphi_yield_Zdc2[NCENT][NSIDE+1][NPT_Mu2][NHAR];
//TH2F *hMu_dphi_eloss_Zdc2[NCENT][NSIDE+1][NPT_Mu2][NHAR];
//TH1D *hMu_eloss_Zdc2[NCENT][NSIDE+1][NETA_Mu2][NHAR][4];

void setStyle1(TH1D *h,int itype1){
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

void fixTitle1(TH1D *h){
  h->GetYaxis()->SetTitleOffset(1.); h->GetXaxis()->SetTitleOffset(1.);
  h->GetYaxis()->SetTitleSize(25.); h->GetXaxis()->SetTitleSize(25.);
}

void get_1Dhist(TH1D *h1d,TH2F *h2d){
  int Nx=h2d->GetNbinsX();

}


void plot::plot1_mu(){
  cout<<"Total Events = "<<hMu_n->GetEntries()<<endl;
  cout<<"Total muon tracks = "<<hMu_pt->GetEntries()<<endl;

  //Plot Track plots
  TCanvas *c1[6];
  for(int ican=0;ican<6;ican++){
    sprintf(name,"c1_%d",ican);
    c1[ican]=new TCanvas(name,name,800,600);
  }

  setStyle1(hMu_n_raw,0);
  setStyle1(hMu_pt_raw,0);
  setStyle1(hMu_eta_raw,0);
  fixTitle1(hMu_n_raw);
  fixTitle1(hMu_pt_raw);
  fixTitle1(hMu_eta_raw);
  setStyle1(hMu_n,0);
  setStyle1(hMu_pt,0);
  setStyle1(hMu_eta,0);
  fixTitle1(hMu_n);
  fixTitle1(hMu_pt);
  fixTitle1(hMu_eta);

  c1[0]->cd();
  gPad->SetLogy();
  hMu_n_raw->Draw();
  hMu_n_raw->SetXTitle("N_{Ch}^{muons}");
  hMu_n_raw->GetXaxis()->SetRangeUser(0.,10.);
  sprintf(name,"Entries = %1.2e",hMu_n_raw->GetEntries());
  text.DrawLatex(0.7,0.8,name);
  
  c1[1]->cd();
  gPad->SetLogy();
  hMu_pt_raw->Draw();
  hMu_pt_raw->Rebin(40);//Make 25 bins for 0-25 pT
  hMu_pt_raw->SetXTitle("p_{T}");
  sprintf(name,"Entries = %1.2e",hMu_pt_raw->GetEntries());
  text.DrawLatex(0.7,0.8,name);

  c1[2]->cd();
  hMu_eta_raw->Draw();
  hMu_eta_raw->SetXTitle("#eta");
  hMu_eta_raw->GetXaxis()->SetRangeUser(-3.,3.);
  sprintf(name,"Entries = %1.2e",hMu_eta_raw->GetEntries());
  text.DrawLatex(0.7,0.8,name);

  c1[3]->cd();
  gPad->SetLogy();
  hMu_n->Draw();
  hMu_n->SetXTitle("N_{Ch}^{muons}");
  hMu_n->GetXaxis()->SetRangeUser(0.,10.);
  sprintf(name,"Entries = %1.2e",hMu_n->GetEntries());
  text.DrawLatex(0.7,0.8,name);

  c1[4]->cd();
  gPad->SetLogy();
  hMu_pt->Draw();
  hMu_pt->Rebin(40);//Make 25 bins for 0-25 pT
  hMu_pt->SetXTitle("p_{T}");
  sprintf(name,"Entries = %1.2e",hMu_pt->GetEntries());
  text.DrawLatex(0.7,0.8,name);

  c1[5]->cd();
  hMu_eta->Draw();
  hMu_eta->SetXTitle("#eta");
  hMu_eta->GetXaxis()->SetRangeUser(-3.,3.);
  sprintf(name,"Entries = %1.2e",hMu_eta->GetEntries());
  text.DrawLatex(0.7,0.8,name);

  plot1_recombine();

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
  
  TFile *fout=new TFile("step1.root","recreate");
  fout->cd();
  for(int icent=0;icent<NCENT;icent++){
    for (int iside=0; iside<NSIDE+1; iside++){
      for (int ihar=0; ihar<NHAR; ihar++){
	for(int ipt=0;ipt<NPT_Mu2;ipt++){
	  hMu_dphi_yield_FCal2[icent][iside][ipt][ihar]->Write();
	}
	for(int ieta=0;ieta<NETA_Mu;ieta++){
	  hMu_dphi_yield_Zdc[icent][iside][ieta][ihar]->Write();
	}
      }
    }
  }
  fout->Close();

  //Save Plots as pdf
  c1[0]->SaveAs("./plots/Mu_n_raw.pdf");
  c1[1]->SaveAs("./plots/Mu_pt_raw.pdf");
  c1[2]->SaveAs("./plots/Mu_eta_raw.pdf");
  c1[3]->SaveAs("./plots/Mu_n.pdf");
  c1[4]->SaveAs("./plots/Mu_pt.pdf");
  c1[5]->SaveAs("./plots/Mu_eta.pdf");
  c2->SaveAs("./plots/Yield_FCal.pdf");
  c3[0]->SaveAs("./plots/Yield_Zdc_sideC.pdf");
  c3[1]->SaveAs("./plots/Yield_Zdc_sideA.pdf");
}


void plot::plot1_recombine(){
  //Add hists
  TH1D *htmp[NPT_Mu];
  TH2F *htmp2[NPT_Mu];
  
  //Recombine Yield plots
  //For FCal
  for(int icent=0;icent<NCENT;icent++){
    for (int iside=0; iside<NSIDE+1; iside++){
      for (int ihar=0; ihar<NHAR; ihar++){
	for(int ipt=0;ipt<NPT_Mu;ipt++){
	  htmp[ipt]=(TH1D*)hMu_dphi_yield_FCal[icent][iside][ipt][ihar]->Clone();
	}
	htmp[0]->Add(htmp[1],1);
	htmp[2]->Add(htmp[3],1);
	htmp[4]->Add(htmp[5],1);
	htmp[7]->Add(htmp[8],1);
	sprintf(name,"hMu_dphi_yield_FCal2_cent%d_side%d_pt0_har%d",icent,iside,ihar);
	hMu_dphi_yield_FCal2[icent][iside][0][ihar]=(TH1D*)htmp[0]->Clone(name);
	sprintf(name,"hMu_dphi_yield_FCal2_cent%d_side%d_pt1_har%d",icent,iside,ihar);
	hMu_dphi_yield_FCal2[icent][iside][1][ihar]=(TH1D*)htmp[2]->Clone(name);
	sprintf(name,"hMu_dphi_yield_FCal2_cent%d_side%d_pt2_har%d",icent,iside,ihar);
	hMu_dphi_yield_FCal2[icent][iside][2][ihar]=(TH1D*)htmp[4]->Clone(name);
	sprintf(name,"hMu_dphi_yield_FCal2_cent%d_side%d_pt3_har%d",icent,iside,ihar);
	hMu_dphi_yield_FCal2[icent][iside][3][ihar]=(TH1D*)htmp[6]->Clone(name);
	sprintf(name,"hMu_dphi_yield_FCal2_cent%d_side%d_pt4_har%d",icent,iside,ihar);
	hMu_dphi_yield_FCal2[icent][iside][4][ihar]=(TH1D*)htmp[7]->Clone(name);
      }
    }
  }

  /*
  //For Zdc
  for(int icent=0;icent<NCENT;icent++){
    for (int iside=0; iside<NSIDE+1; iside++){
      for (int ihar=0; ihar<NHAR; ihar++){
	for(int ipt=0;ipt<NPT_Mu;ipt++){
	  htmp[ipt]=(TH1D*)hMu_dphi_yield_Zdc[icent][iside][ipt][ihar]->Clone();
	}
	htmp[0]->Add(htmp[1],1);
	htmp[2]->Add(htmp[3],1);
	htmp[4]->Add(htmp[5],1);
	htmp[7]->Add(htmp[8],1);
	sprintf(name,"hMu_dphi_yield_Zdc2_cent%d_side%d_pt0_har%d",icent,iside,ihar);
	hMu_dphi_yield_Zdc2[icent][iside][0][ihar]=(TH1D*)htmp[0]->Clone(name);
	sprintf(name,"hMu_dphi_yield_Zdc2_cent%d_side%d_pt1_har%d",icent,iside,ihar);
	hMu_dphi_yield_Zdc2[icent][iside][1][ihar]=(TH1D*)htmp[2]->Clone(name);
	sprintf(name,"hMu_dphi_yield_Zdc2_cent%d_side%d_pt2_har%d",icent,iside,ihar);
	hMu_dphi_yield_Zdc2[icent][iside][2][ihar]=(TH1D*)htmp[4]->Clone(name);
	sprintf(name,"hMu_dphi_yield_Zdc2_cent%d_side%d_pt3_har%d",icent,iside,ihar);
	hMu_dphi_yield_Zdc2[icent][iside][3][ihar]=(TH1D*)htmp[6]->Clone(name);
	sprintf(name,"hMu_dphi_yield_Zdc2_cent%d_side%d_pt4_har%d",icent,iside,ihar);
	hMu_dphi_yield_Zdc2[icent][iside][4][ihar]=(TH1D*)htmp[7]->Clone(name);
      }
    }
  }
  */
  //Recombine Momentum imbalance plots
  //For FCal
  for(int icent=0;icent<NCENT;icent++){
    for (int iside=0; iside<NSIDE+1; iside++){
      for (int ihar=0; ihar<NHAR; ihar++){
	for(int ipt=0;ipt<NPT_Mu;ipt++){
	  htmp2[ipt]=(TH2F*)hMu_dphi_eloss_FCal[icent][iside][ipt][ihar]->Clone();
	}
	htmp2[0]->Add(htmp2[1],1);
	htmp2[2]->Add(htmp2[3],1);
	htmp2[4]->Add(htmp2[5],1);
	htmp2[7]->Add(htmp2[8],1);
	sprintf(name,"hMu_dphi_eloss_FCal2_cent%d_side%d_pt0_har%d",icent,iside,ihar);
	hMu_dphi_eloss_FCal2[icent][iside][0][ihar]=(TH2F*)htmp2[0]->Clone(name);
	sprintf(name,"hMu_dphi_eloss_FCal2_cent%d_side%d_pt1_har%d",icent,iside,ihar);
	hMu_dphi_eloss_FCal2[icent][iside][1][ihar]=(TH2F*)htmp2[2]->Clone(name);
	sprintf(name,"hMu_dphi_eloss_FCal2_cent%d_side%d_pt2_har%d",icent,iside,ihar);
	hMu_dphi_eloss_FCal2[icent][iside][2][ihar]=(TH2F*)htmp2[4]->Clone(name);
	sprintf(name,"hMu_dphi_eloss_FCal2_cent%d_side%d_pt3_har%d",icent,iside,ihar);
	hMu_dphi_eloss_FCal2[icent][iside][3][ihar]=(TH2F*)htmp2[6]->Clone(name);
	sprintf(name,"hMu_dphi_eloss_FCal2_cent%d_side%d_pt4_har%d",icent,iside,ihar);
	hMu_dphi_eloss_FCal2[icent][iside][4][ihar]=(TH2F*)htmp2[7]->Clone(name);
      }
    }
  }
  /*
  //For Zdc
  for(int icent=0;icent<NCENT;icent++){
    for (int iside=0; iside<NSIDE+1; iside++){
      for (int ihar=0; ihar<NHAR; ihar++){
	for(int ipt=0;ipt<NPT_Mu;ipt++){
	  htmp2[ipt]=(TH2F*)hMu_dphi_eloss_Zdc[icent][iside][ipt][ihar]->Clone();
	}
	htmp2[0]->Add(htmp2[1],1);
	htmp2[2]->Add(htmp2[3],1);
	htmp2[4]->Add(htmp2[5],1);
	htmp2[7]->Add(htmp2[8],1);
	sprintf(name,"hMu_dphi_eloss_Zdc2_cent%d_side%d_pt0_har%d",icent,iside,ihar);
	hMu_dphi_eloss_Zdc2[icent][iside][0][ihar]=(TH2F*)htmp2[0]->Clone(name);
	sprintf(name,"hMu_dphi_eloss_Zdc2_cent%d_side%d_pt1_har%d",icent,iside,ihar);
	hMu_dphi_eloss_Zdc2[icent][iside][1][ihar]=(TH2F*)htmp2[2]->Clone(name);
	sprintf(name,"hMu_dphi_eloss_Zdc2_cent%d_side%d_pt2_har%d",icent,iside,ihar);
	hMu_dphi_eloss_Zdc2[icent][iside][2][ihar]=(TH2F*)htmp2[4]->Clone(name);
	sprintf(name,"hMu_dphi_eloss_Zdc2_cent%d_side%d_pt3_har%d",icent,iside,ihar);
	hMu_dphi_eloss_Zdc2[icent][iside][3][ihar]=(TH2F*)htmp2[6]->Clone(name);
	sprintf(name,"hMu_dphi_eloss_Zdc2_cent%d_side%d_pt4_har%d",icent,iside,ihar);
	hMu_dphi_eloss_Zdc2[icent][iside][4][ihar]=(TH2F*)htmp2[7]->Clone(name);
      }
    }
  }
  */
}

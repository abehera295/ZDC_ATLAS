
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


//Recombine Cent bins and pT bins
const int NCENT2=4;//0-10,10-20,20-40,40-60
const int NPT_Mu2=5;//4-5,5-6,6-8,8-10,10-14
const double cent_mat2[NCENT2+1]={0.,10.,20.,40.,60.};
const double pt_mat_Mu2[NPT_Mu2+1]={4.,5.,6.,8.,10.,14.};

TH1D *hMu_dphi_yield_FCal2[NCENT2][NSIDE+1][NPT_Mu2][NHAR];
TH1D *hMu_dphi_yield_Zdc2[NCENT2][NSIDE+1][NPT_Mu2][NHAR];
TH2F *hMu_dphi_eloss_FCal2[NCENT2][NSIDE+1][NPT_Mu2][NHAR];
TH2F *hMu_dphi_eloss_Zdc2[NCENT2][NSIDE+1][NPT_Mu2][NHAR];

void plot::plot1_mu(){
  cout<<"Total Events = "<<hMu_n->GetEntries()<<endl;
  cout<<"Total muon tracks = "<<hMu_pt->GetEntries()<<endl;

  //Plot Track plots
  TCanvas *c1[3];
  for(int ican=0;ican<3;ican++){
    sprintf(name,"c1_%d",ican);
    c1[ican]=new TCanvas(name,name,800,600);
  }

  setStyle1(hMu_n,0);
  setStyle1(hMu_pt,0);
  setStyle1(hMu_eta,0);
  fixTitle1(hMu_n);
  fixTitle1(hMu_pt);
  fixTitle1(hMu_eta);
  c1[0]->cd();
  gPad->SetLogy();
  hMu_n->Draw();
  hMu_n->SetXTitle("N_{Ch}^{muons}");
  hMu_n->GetXaxis()->SetRangeUser(0.,10.);

  c1[1]->cd();
  gPad->SetLogy();
  hMu_pt->Draw();
  hMu_pt->SetXTitle("p_{T}");
  hMu_pt->GetXaxis()->SetRangeUser(0.,25.);

  c1[2]->cd();
  hMu_eta->Draw();
  hMu_eta->SetXTitle("#eta");
  hMu_eta->GetXaxis()->SetRangeUser(-3.,3.);

  plot1_recombine();
  int icent_mat[3]={0,2,3};//0-10%,20-40%,40-60%
  string centP_mat[3]={"0-10%","20-40%","40-60%"};
  int ipt_mat[2]={0,3};
  string ptR_mat[2]={"4< p_{T} < 5 GeV","8 < p_{T} < 10 GeV"};
  //Plot Yield
  
  TCanvas *c2[2];
  TLegend *leg2[2];
  TPad *pad2[2][6];
  for(int ican=0;ican<2;ican++){
    sprintf(name,"c2_%d",ican);
    c2[ican]=new TCanvas(name,name,1200,800);
    c2[ican]->Divide(3,2);
  }
  int jdet,jhar;
  //For FCal
  jdet=0;jhar=1;//v2
  //Divide_Pad(c2[jdet],pad2[jdet],2,3);
  for(int i=0;i<3;i++){
    int icent=icent_mat[i];
    int ipt;
    c2[jdet]->cd(i+1);
    //pad2[jdet][i]->cd();
    ipt=ipt_mat[0];
    setStyle1(hMu_dphi_yield_FCal2[icent][2][ipt][jhar],0);
    hMu_dphi_yield_FCal2[icent][2][ipt][jhar]->Scale(1./1000.);
    hMu_dphi_yield_FCal2[icent][2][ipt][jhar]->Draw();
    hMu_dphi_yield_FCal2[icent][2][ipt][jhar]->SetXTitle("2(#phi-#Psi_{2}FCal)");
    hMu_dphi_yield_FCal2[icent][2][ipt][jhar]->SetYTitle("Signal Counts/1000");
    //hMu_dphi_yield_FCal2[icent][2][ipt][jhar]->GetYaxis()->SetRangeUser(2000.,15000);
    sprintf(name,"%s",centP_mat[i].c_str());
    text.DrawLatex(0.25,0.3,name);
    sprintf(name,"%s",ptR_mat[0].c_str());
    text.DrawLatex(0.2,0.2,name);

    c2[jdet]->cd(i+3+1);
    //pad2[jdet][i+3]->cd();
    ipt=ipt_mat[1];
    setStyle1(hMu_dphi_yield_FCal2[icent][2][ipt][jhar],0);
    hMu_dphi_yield_FCal2[icent][2][ipt][jhar]->Scale(1./1000.);
    hMu_dphi_yield_FCal2[icent][2][ipt][jhar]->Draw();
    hMu_dphi_yield_FCal2[icent][2][ipt][jhar]->SetXTitle("2(#phi-#Psi_{2}FCal)");
    hMu_dphi_yield_FCal2[icent][2][ipt][jhar]->SetYTitle("Signal Counts/1000");
    //hMu_dphi_yield_FCal2[icent][2][ipt][jhar]->GetYaxis()->SetRangeUser(0.,2000);
    sprintf(name,"%s",centP_mat[i].c_str());
    text.DrawLatex(0.25,0.3,name);
    sprintf(name,"%s",ptR_mat[1].c_str());
    text.DrawLatex(0.2,0.2,name);
  }

  //For Zdc
  jdet=1;jhar=0;//v1
  //Divide_Pad(c2[jdet],pad2[jdet],2,3);
  for(int i=0;i<3;i++){
    int icent=icent_mat[i];
    int ipt;
    c2[jdet]->cd(i+1);
    //pad2[jdet][i]->cd();
    ipt=ipt_mat[0];
    setStyle1(hMu_dphi_yield_Zdc2[icent][2][ipt][jhar],0);
    hMu_dphi_yield_Zdc2[icent][2][ipt][jhar]->Scale(1./1000.);
    hMu_dphi_yield_Zdc2[icent][2][ipt][jhar]->Draw();
    hMu_dphi_yield_Zdc2[icent][2][ipt][jhar]->SetXTitle("#phi-#Psi_{1}Zdc");
    hMu_dphi_yield_Zdc2[icent][2][ipt][jhar]->SetYTitle("Signal Counts/1000");
    //hMu_dphi_yield_Zdc2[icent][2][ipt][jhar]->GetYaxis()->SetRangeUser(2000.,15000);
    sprintf(name,"%s",centP_mat[i].c_str());
    text.DrawLatex(0.25,0.3,name);
    sprintf(name,"%s",ptR_mat[0].c_str());
    text.DrawLatex(0.2,0.2,name);

    c2[jdet]->cd(i+3+1);
    //pad2[jdet][i+3]->cd();
    ipt=ipt_mat[1];
    setStyle1(hMu_dphi_yield_Zdc2[icent][2][ipt][jhar],0);
    hMu_dphi_yield_Zdc2[icent][2][ipt][jhar]->Scale(1./1000.);
    hMu_dphi_yield_Zdc2[icent][2][ipt][jhar]->Draw();
    hMu_dphi_yield_Zdc2[icent][2][ipt][jhar]->SetXTitle("#phi-#Psi_{1}Zdc");
    hMu_dphi_yield_Zdc2[icent][2][ipt][jhar]->SetYTitle("Signal Counts/1000");
    //hMu_dphi_yield_Zdc2[icent][2][ipt][jhar]->GetYaxis()->SetRangeUser(0.,2000);
    sprintf(name,"%s",centP_mat[i].c_str());
    text.DrawLatex(0.25,0.3,name);
    sprintf(name,"%s",ptR_mat[1].c_str());
    text.DrawLatex(0.2,0.2,name);
  }

  TFile *fout=new TFile("step1.root","recreate");
  fout->cd();
  for(int icent=0;icent<NCENT2;icent++){
    for (int iside=0; iside<NSIDE+1; iside++){
      for(int ipt=0;ipt<NPT_Mu2;ipt++){
        for (int ihar=0; ihar<NHAR; ihar++){
	  hMu_dphi_yield_FCal2[icent][iside][ipt][ihar]->Write();
	  hMu_dphi_yield_Zdc2[icent][iside][ipt][ihar]->Write();
	}
      }
    }
  }
  cout<<"NCENT2="<<NCENT2<<endl;
  fout->Close();

}


void plot::plot1_recombine(){
  //Add hists
  TH1D *htmp[NPT_Mu];
  TH2F *htmp2[NPT_Mu];
  int jcent;
  
  //Recombine Yield plots
  //For FCal
  for (int iside=0; iside<NSIDE+1; iside++){
    for (int ihar=0; ihar<NHAR; ihar++){
      //1. Cent:0-10
      jcent=0;
      for(int ipt=0;ipt<NPT_Mu;ipt++){
	htmp[ipt]=(TH1D*)hMu_dphi_yield_FCal[0][iside][ipt][ihar]->Clone();
	htmp[ipt]->Add(hMu_dphi_yield_FCal[1][iside][ipt][ihar],1);
      }
      htmp[0]->Add(htmp[1],1);
      htmp[2]->Add(htmp[3],1);
      htmp[4]->Add(htmp[5],1);
      htmp[7]->Add(htmp[8],1);
      sprintf(name,"hMu_dphi_yield_FCal2_cent%d_side%d_pt0_har%d",jcent,iside,ihar);
      hMu_dphi_yield_FCal2[jcent][iside][0][ihar]=(TH1D*)htmp[0]->Clone(name);
      sprintf(name,"hMu_dphi_yield_FCal2_cent%d_side%d_pt1_har%d",jcent,iside,ihar);
      hMu_dphi_yield_FCal2[jcent][iside][1][ihar]=(TH1D*)htmp[2]->Clone(name);
      sprintf(name,"hMu_dphi_yield_FCal2_cent%d_side%d_pt2_har%d",jcent,iside,ihar);
      hMu_dphi_yield_FCal2[jcent][iside][2][ihar]=(TH1D*)htmp[4]->Clone(name);
      sprintf(name,"hMu_dphi_yield_FCal2_cent%d_side%d_pt3_har%d",jcent,iside,ihar);
      hMu_dphi_yield_FCal2[jcent][iside][3][ihar]=(TH1D*)htmp[6]->Clone(name);
      sprintf(name,"hMu_dphi_yield_FCal2_cent%d_side%d_pt4_har%d",jcent,iside,ihar);
      hMu_dphi_yield_FCal2[jcent][iside][4][ihar]=(TH1D*)htmp[7]->Clone(name);

      //2. Cent:10-20
      jcent=1;
      for(int ipt=0;ipt<NPT_Mu;ipt++){
	htmp[ipt]=(TH1D*)hMu_dphi_yield_FCal[2][iside][ipt][ihar]->Clone();
      }
      htmp[0]->Add(htmp[1],1);
      htmp[2]->Add(htmp[3],1);
      htmp[4]->Add(htmp[5],1);
      htmp[7]->Add(htmp[8],1);
      sprintf(name,"hMu_dphi_yield_FCal2_cent%d_side%d_pt0_har%d",jcent,iside,ihar);
      hMu_dphi_yield_FCal2[jcent][iside][0][ihar]=(TH1D*)htmp[0]->Clone(name);
      sprintf(name,"hMu_dphi_yield_FCal2_cent%d_side%d_pt1_har%d",jcent,iside,ihar);
      hMu_dphi_yield_FCal2[jcent][iside][1][ihar]=(TH1D*)htmp[2]->Clone(name);
      sprintf(name,"hMu_dphi_yield_FCal2_cent%d_side%d_pt2_har%d",jcent,iside,ihar);
      hMu_dphi_yield_FCal2[jcent][iside][2][ihar]=(TH1D*)htmp[4]->Clone(name);
      sprintf(name,"hMu_dphi_yield_FCal2_cent%d_side%d_pt3_har%d",jcent,iside,ihar);
      hMu_dphi_yield_FCal2[jcent][iside][3][ihar]=(TH1D*)htmp[6]->Clone(name);
      sprintf(name,"hMu_dphi_yield_FCal2_cent%d_side%d_pt4_har%d",jcent,iside,ihar);
      hMu_dphi_yield_FCal2[jcent][iside][4][ihar]=(TH1D*)htmp[7]->Clone(name);

      //3. Cent:20-40
      jcent=2;
      for(int ipt=0;ipt<NPT_Mu;ipt++){
	htmp[ipt]=(TH1D*)hMu_dphi_yield_FCal[3][iside][ipt][ihar]->Clone();
	htmp[ipt]->Add(hMu_dphi_yield_FCal[4][iside][ipt][ihar],1);
      }
      htmp[0]->Add(htmp[1],1);
      htmp[2]->Add(htmp[3],1);
      htmp[4]->Add(htmp[5],1);
      htmp[7]->Add(htmp[8],1);
      sprintf(name,"hMu_dphi_yield_FCal2_cent%d_side%d_pt0_har%d",jcent,iside,ihar);
      hMu_dphi_yield_FCal2[jcent][iside][0][ihar]=(TH1D*)htmp[0]->Clone(name);
      sprintf(name,"hMu_dphi_yield_FCal2_cent%d_side%d_pt1_har%d",jcent,iside,ihar);
      hMu_dphi_yield_FCal2[jcent][iside][1][ihar]=(TH1D*)htmp[2]->Clone(name);
      sprintf(name,"hMu_dphi_yield_FCal2_cent%d_side%d_pt2_har%d",jcent,iside,ihar);
      hMu_dphi_yield_FCal2[jcent][iside][2][ihar]=(TH1D*)htmp[4]->Clone(name);
      sprintf(name,"hMu_dphi_yield_FCal2_cent%d_side%d_pt3_har%d",jcent,iside,ihar);
      hMu_dphi_yield_FCal2[jcent][iside][3][ihar]=(TH1D*)htmp[6]->Clone(name);
      sprintf(name,"hMu_dphi_yield_FCal2_cent%d_side%d_pt4_har%d",jcent,iside,ihar);
      hMu_dphi_yield_FCal2[jcent][iside][4][ihar]=(TH1D*)htmp[7]->Clone(name);

      //3. Cent:40-60
      jcent=3;
      for(int ipt=0;ipt<NPT_Mu;ipt++){
	htmp[ipt]=(TH1D*)hMu_dphi_yield_FCal[5][iside][ipt][ihar]->Clone();
	htmp[ipt]->Add(hMu_dphi_yield_FCal[6][iside][ipt][ihar],1);
      }
      htmp[0]->Add(htmp[1],1);
      htmp[2]->Add(htmp[3],1);
      htmp[4]->Add(htmp[5],1);
      htmp[7]->Add(htmp[8],1);
      sprintf(name,"hMu_dphi_yield_FCal2_cent%d_side%d_pt0_har%d",jcent,iside,ihar);
      hMu_dphi_yield_FCal2[jcent][iside][0][ihar]=(TH1D*)htmp[0]->Clone(name);
      sprintf(name,"hMu_dphi_yield_FCal2_cent%d_side%d_pt1_har%d",jcent,iside,ihar);
      hMu_dphi_yield_FCal2[jcent][iside][1][ihar]=(TH1D*)htmp[2]->Clone(name);
      sprintf(name,"hMu_dphi_yield_FCal2_cent%d_side%d_pt2_har%d",jcent,iside,ihar);
      hMu_dphi_yield_FCal2[jcent][iside][2][ihar]=(TH1D*)htmp[4]->Clone(name);
      sprintf(name,"hMu_dphi_yield_FCal2_cent%d_side%d_pt3_har%d",jcent,iside,ihar);
      hMu_dphi_yield_FCal2[jcent][iside][3][ihar]=(TH1D*)htmp[6]->Clone(name);
      sprintf(name,"hMu_dphi_yield_FCal2_cent%d_side%d_pt4_har%d",jcent,iside,ihar);
      hMu_dphi_yield_FCal2[jcent][iside][4][ihar]=(TH1D*)htmp[7]->Clone(name);

    }
  }

  //For Zdc
  for (int iside=0; iside<NSIDE+1; iside++){
    for (int ihar=0; ihar<NHAR; ihar++){
      //1. Cent:0-10
      jcent=0;
      for(int ipt=0;ipt<NPT_Mu;ipt++){
	htmp[ipt]=(TH1D*)hMu_dphi_yield_Zdc[0][iside][ipt][ihar]->Clone();
	htmp[ipt]->Add(hMu_dphi_yield_Zdc[1][iside][ipt][ihar],1);
      }
      htmp[0]->Add(htmp[1],1);
      htmp[2]->Add(htmp[3],1);
      htmp[4]->Add(htmp[5],1);
      htmp[7]->Add(htmp[8],1);
      sprintf(name,"hMu_dphi_yield_Zdc2_cent%d_side%d_pt0_har%d",jcent,iside,ihar);
      hMu_dphi_yield_Zdc2[jcent][iside][0][ihar]=(TH1D*)htmp[0]->Clone(name);
      sprintf(name,"hMu_dphi_yield_Zdc2_cent%d_side%d_pt1_har%d",jcent,iside,ihar);
      hMu_dphi_yield_Zdc2[jcent][iside][1][ihar]=(TH1D*)htmp[2]->Clone(name);
      sprintf(name,"hMu_dphi_yield_Zdc2_cent%d_side%d_pt2_har%d",jcent,iside,ihar);
      hMu_dphi_yield_Zdc2[jcent][iside][2][ihar]=(TH1D*)htmp[4]->Clone(name);
      sprintf(name,"hMu_dphi_yield_Zdc2_cent%d_side%d_pt3_har%d",jcent,iside,ihar);
      hMu_dphi_yield_Zdc2[jcent][iside][3][ihar]=(TH1D*)htmp[6]->Clone(name);
      sprintf(name,"hMu_dphi_yield_Zdc2_cent%d_side%d_pt4_har%d",jcent,iside,ihar);
      hMu_dphi_yield_Zdc2[jcent][iside][4][ihar]=(TH1D*)htmp[7]->Clone(name);

      //2. Cent:10-20
      jcent=1;
      for(int ipt=0;ipt<NPT_Mu;ipt++){
	htmp[ipt]=(TH1D*)hMu_dphi_yield_Zdc[2][iside][ipt][ihar]->Clone();
      }
      htmp[0]->Add(htmp[1],1);
      htmp[2]->Add(htmp[3],1);
      htmp[4]->Add(htmp[5],1);
      htmp[7]->Add(htmp[8],1);
      sprintf(name,"hMu_dphi_yield_Zdc2_cent%d_side%d_pt0_har%d",jcent,iside,ihar);
      hMu_dphi_yield_Zdc2[jcent][iside][0][ihar]=(TH1D*)htmp[0]->Clone(name);
      sprintf(name,"hMu_dphi_yield_Zdc2_cent%d_side%d_pt1_har%d",jcent,iside,ihar);
      hMu_dphi_yield_Zdc2[jcent][iside][1][ihar]=(TH1D*)htmp[2]->Clone(name);
      sprintf(name,"hMu_dphi_yield_Zdc2_cent%d_side%d_pt2_har%d",jcent,iside,ihar);
      hMu_dphi_yield_Zdc2[jcent][iside][2][ihar]=(TH1D*)htmp[4]->Clone(name);
      sprintf(name,"hMu_dphi_yield_Zdc2_cent%d_side%d_pt3_har%d",jcent,iside,ihar);
      hMu_dphi_yield_Zdc2[jcent][iside][3][ihar]=(TH1D*)htmp[6]->Clone(name);
      sprintf(name,"hMu_dphi_yield_Zdc2_cent%d_side%d_pt4_har%d",jcent,iside,ihar);
      hMu_dphi_yield_Zdc2[jcent][iside][4][ihar]=(TH1D*)htmp[7]->Clone(name);

      //3. Cent:20-40
      jcent=2;
      for(int ipt=0;ipt<NPT_Mu;ipt++){
	htmp[ipt]=(TH1D*)hMu_dphi_yield_Zdc[3][iside][ipt][ihar]->Clone();
	htmp[ipt]->Add(hMu_dphi_yield_Zdc[4][iside][ipt][ihar],1);
      }
      htmp[0]->Add(htmp[1],1);
      htmp[2]->Add(htmp[3],1);
      htmp[4]->Add(htmp[5],1);
      htmp[7]->Add(htmp[8],1);
      sprintf(name,"hMu_dphi_yield_Zdc2_cent%d_side%d_pt0_har%d",jcent,iside,ihar);
      hMu_dphi_yield_Zdc2[jcent][iside][0][ihar]=(TH1D*)htmp[0]->Clone(name);
      sprintf(name,"hMu_dphi_yield_Zdc2_cent%d_side%d_pt1_har%d",jcent,iside,ihar);
      hMu_dphi_yield_Zdc2[jcent][iside][1][ihar]=(TH1D*)htmp[2]->Clone(name);
      sprintf(name,"hMu_dphi_yield_Zdc2_cent%d_side%d_pt2_har%d",jcent,iside,ihar);
      hMu_dphi_yield_Zdc2[jcent][iside][2][ihar]=(TH1D*)htmp[4]->Clone(name);
      sprintf(name,"hMu_dphi_yield_Zdc2_cent%d_side%d_pt3_har%d",jcent,iside,ihar);
      hMu_dphi_yield_Zdc2[jcent][iside][3][ihar]=(TH1D*)htmp[6]->Clone(name);
      sprintf(name,"hMu_dphi_yield_Zdc2_cent%d_side%d_pt4_har%d",jcent,iside,ihar);
      hMu_dphi_yield_Zdc2[jcent][iside][4][ihar]=(TH1D*)htmp[7]->Clone(name);

      //3. Cent:40-60
      jcent=3;
      for(int ipt=0;ipt<NPT_Mu;ipt++){
	htmp[ipt]=(TH1D*)hMu_dphi_yield_Zdc[5][iside][ipt][ihar]->Clone();
	htmp[ipt]->Add(hMu_dphi_yield_Zdc[6][iside][ipt][ihar],1);
      }
      htmp[0]->Add(htmp[1],1);
      htmp[2]->Add(htmp[3],1);
      htmp[4]->Add(htmp[5],1);
      htmp[7]->Add(htmp[8],1);
      sprintf(name,"hMu_dphi_yield_Zdc2_cent%d_side%d_pt0_har%d",jcent,iside,ihar);
      hMu_dphi_yield_Zdc2[jcent][iside][0][ihar]=(TH1D*)htmp[0]->Clone(name);
      sprintf(name,"hMu_dphi_yield_Zdc2_cent%d_side%d_pt1_har%d",jcent,iside,ihar);
      hMu_dphi_yield_Zdc2[jcent][iside][1][ihar]=(TH1D*)htmp[2]->Clone(name);
      sprintf(name,"hMu_dphi_yield_Zdc2_cent%d_side%d_pt2_har%d",jcent,iside,ihar);
      hMu_dphi_yield_Zdc2[jcent][iside][2][ihar]=(TH1D*)htmp[4]->Clone(name);
      sprintf(name,"hMu_dphi_yield_Zdc2_cent%d_side%d_pt3_har%d",jcent,iside,ihar);
      hMu_dphi_yield_Zdc2[jcent][iside][3][ihar]=(TH1D*)htmp[6]->Clone(name);
      sprintf(name,"hMu_dphi_yield_Zdc2_cent%d_side%d_pt4_har%d",jcent,iside,ihar);
      hMu_dphi_yield_Zdc2[jcent][iside][4][ihar]=(TH1D*)htmp[7]->Clone(name);

    }
  }


  //Recombine Momentum imbalance plots
  //For FCal
  for (int iside=0; iside<NSIDE+1; iside++){
    for (int ihar=0; ihar<NHAR; ihar++){
      //1. Cent:0-10
      jcent=0;
      for(int ipt=0;ipt<NPT_Mu;ipt++){
	htmp2[ipt]=(TH2F*)hMu_dphi_eloss_FCal[0][iside][ipt][ihar]->Clone();
	htmp2[ipt]->Add(hMu_dphi_eloss_FCal[1][iside][ipt][ihar],1);
      }
      htmp2[0]->Add(htmp2[1],1);
      htmp2[2]->Add(htmp2[3],1);
      htmp2[4]->Add(htmp2[5],1);
      htmp2[7]->Add(htmp2[8],1);
      sprintf(name,"hMu_dphi_eloss_FCal2_cent%d_side%d_pt0_har%d",jcent,iside,ihar);
      hMu_dphi_eloss_FCal2[jcent][iside][0][ihar]=(TH2F*)htmp2[0]->Clone(name);
      sprintf(name,"hMu_dphi_eloss_FCal2_cent%d_side%d_pt1_har%d",jcent,iside,ihar);
      hMu_dphi_eloss_FCal2[jcent][iside][1][ihar]=(TH2F*)htmp2[2]->Clone(name);
      sprintf(name,"hMu_dphi_eloss_FCal2_cent%d_side%d_pt2_har%d",jcent,iside,ihar);
      hMu_dphi_eloss_FCal2[jcent][iside][2][ihar]=(TH2F*)htmp2[4]->Clone(name);
      sprintf(name,"hMu_dphi_eloss_FCal2_cent%d_side%d_pt3_har%d",jcent,iside,ihar);
      hMu_dphi_eloss_FCal2[jcent][iside][3][ihar]=(TH2F*)htmp2[6]->Clone(name);
      sprintf(name,"hMu_dphi_eloss_FCal2_cent%d_side%d_pt4_har%d",jcent,iside,ihar);
      hMu_dphi_eloss_FCal2[jcent][iside][4][ihar]=(TH2F*)htmp2[7]->Clone(name);

      //2. Cent:10-20
      jcent=1;
      for(int ipt=0;ipt<NPT_Mu;ipt++){
	htmp2[ipt]=(TH2F*)hMu_dphi_eloss_FCal[2][iside][ipt][ihar]->Clone();
      }
      htmp2[0]->Add(htmp2[1],1);
      htmp2[2]->Add(htmp2[3],1);
      htmp2[4]->Add(htmp2[5],1);
      htmp2[7]->Add(htmp2[8],1);
      sprintf(name,"hMu_dphi_eloss_FCal2_cent%d_side%d_pt0_har%d",jcent,iside,ihar);
      hMu_dphi_eloss_FCal2[jcent][iside][0][ihar]=(TH2F*)htmp2[0]->Clone(name);
      sprintf(name,"hMu_dphi_eloss_FCal2_cent%d_side%d_pt1_har%d",jcent,iside,ihar);
      hMu_dphi_eloss_FCal2[jcent][iside][1][ihar]=(TH2F*)htmp2[2]->Clone(name);
      sprintf(name,"hMu_dphi_eloss_FCal2_cent%d_side%d_pt2_har%d",jcent,iside,ihar);
      hMu_dphi_eloss_FCal2[jcent][iside][2][ihar]=(TH2F*)htmp2[4]->Clone(name);
      sprintf(name,"hMu_dphi_eloss_FCal2_cent%d_side%d_pt3_har%d",jcent,iside,ihar);
      hMu_dphi_eloss_FCal2[jcent][iside][3][ihar]=(TH2F*)htmp2[6]->Clone(name);
      sprintf(name,"hMu_dphi_eloss_FCal2_cent%d_side%d_pt4_har%d",jcent,iside,ihar);
      hMu_dphi_eloss_FCal2[jcent][iside][4][ihar]=(TH2F*)htmp2[7]->Clone(name);

      //3. Cent:20-40
      jcent=2;
      for(int ipt=0;ipt<NPT_Mu;ipt++){
	htmp2[ipt]=(TH2F*)hMu_dphi_eloss_FCal[3][iside][ipt][ihar]->Clone();
	htmp2[ipt]->Add(hMu_dphi_eloss_FCal[4][iside][ipt][ihar],1);
      }
      htmp2[0]->Add(htmp2[1],1);
      htmp2[2]->Add(htmp2[3],1);
      htmp2[4]->Add(htmp2[5],1);
      htmp2[7]->Add(htmp2[8],1);
      sprintf(name,"hMu_dphi_eloss_FCal2_cent%d_side%d_pt0_har%d",jcent,iside,ihar);
      hMu_dphi_eloss_FCal2[jcent][iside][0][ihar]=(TH2F*)htmp2[0]->Clone(name);
      sprintf(name,"hMu_dphi_eloss_FCal2_cent%d_side%d_pt1_har%d",jcent,iside,ihar);
      hMu_dphi_eloss_FCal2[jcent][iside][1][ihar]=(TH2F*)htmp2[2]->Clone(name);
      sprintf(name,"hMu_dphi_eloss_FCal2_cent%d_side%d_pt2_har%d",jcent,iside,ihar);
      hMu_dphi_eloss_FCal2[jcent][iside][2][ihar]=(TH2F*)htmp2[4]->Clone(name);
      sprintf(name,"hMu_dphi_eloss_FCal2_cent%d_side%d_pt3_har%d",jcent,iside,ihar);
      hMu_dphi_eloss_FCal2[jcent][iside][3][ihar]=(TH2F*)htmp2[6]->Clone(name);
      sprintf(name,"hMu_dphi_eloss_FCal2_cent%d_side%d_pt4_har%d",jcent,iside,ihar);
      hMu_dphi_eloss_FCal2[jcent][iside][4][ihar]=(TH2F*)htmp2[7]->Clone(name);

      //3. Cent:40-60
      jcent=3;
      for(int ipt=0;ipt<NPT_Mu;ipt++){
	htmp2[ipt]=(TH2F*)hMu_dphi_eloss_FCal[5][iside][ipt][ihar]->Clone();
	htmp2[ipt]->Add(hMu_dphi_eloss_FCal[6][iside][ipt][ihar],1);
      }
      htmp2[0]->Add(htmp2[1],1);
      htmp2[2]->Add(htmp2[3],1);
      htmp2[4]->Add(htmp2[5],1);
      htmp2[7]->Add(htmp2[8],1);
      sprintf(name,"hMu_dphi_eloss_FCal2_cent%d_side%d_pt0_har%d",jcent,iside,ihar);
      hMu_dphi_eloss_FCal2[jcent][iside][0][ihar]=(TH2F*)htmp2[0]->Clone(name);
      sprintf(name,"hMu_dphi_eloss_FCal2_cent%d_side%d_pt1_har%d",jcent,iside,ihar);
      hMu_dphi_eloss_FCal2[jcent][iside][1][ihar]=(TH2F*)htmp2[2]->Clone(name);
      sprintf(name,"hMu_dphi_eloss_FCal2_cent%d_side%d_pt2_har%d",jcent,iside,ihar);
      hMu_dphi_eloss_FCal2[jcent][iside][2][ihar]=(TH2F*)htmp2[4]->Clone(name);
      sprintf(name,"hMu_dphi_eloss_FCal2_cent%d_side%d_pt3_har%d",jcent,iside,ihar);
      hMu_dphi_eloss_FCal2[jcent][iside][3][ihar]=(TH2F*)htmp2[6]->Clone(name);
      sprintf(name,"hMu_dphi_eloss_FCal2_cent%d_side%d_pt4_har%d",jcent,iside,ihar);
      hMu_dphi_eloss_FCal2[jcent][iside][4][ihar]=(TH2F*)htmp2[7]->Clone(name);

    }
  }

  //For Zdc
  for (int iside=0; iside<NSIDE+1; iside++){
    for (int ihar=0; ihar<NHAR; ihar++){
      //1. Cent:0-10
      jcent=0;
      for(int ipt=0;ipt<NPT_Mu;ipt++){
	htmp2[ipt]=(TH2F*)hMu_dphi_eloss_Zdc[0][iside][ipt][ihar]->Clone();
	htmp2[ipt]->Add(hMu_dphi_eloss_Zdc[1][iside][ipt][ihar],1);
      }
      htmp2[0]->Add(htmp2[1],1);
      htmp2[2]->Add(htmp2[3],1);
      htmp2[4]->Add(htmp2[5],1);
      htmp2[7]->Add(htmp2[8],1);
      sprintf(name,"hMu_dphi_eloss_Zdc2_cent%d_side%d_pt0_har%d",jcent,iside,ihar);
      hMu_dphi_eloss_Zdc2[jcent][iside][0][ihar]=(TH2F*)htmp2[0]->Clone(name);
      sprintf(name,"hMu_dphi_eloss_Zdc2_cent%d_side%d_pt1_har%d",jcent,iside,ihar);
      hMu_dphi_eloss_Zdc2[jcent][iside][1][ihar]=(TH2F*)htmp2[2]->Clone(name);
      sprintf(name,"hMu_dphi_eloss_Zdc2_cent%d_side%d_pt2_har%d",jcent,iside,ihar);
      hMu_dphi_eloss_Zdc2[jcent][iside][2][ihar]=(TH2F*)htmp2[4]->Clone(name);
      sprintf(name,"hMu_dphi_eloss_Zdc2_cent%d_side%d_pt3_har%d",jcent,iside,ihar);
      hMu_dphi_eloss_Zdc2[jcent][iside][3][ihar]=(TH2F*)htmp2[6]->Clone(name);
      sprintf(name,"hMu_dphi_eloss_Zdc2_cent%d_side%d_pt4_har%d",jcent,iside,ihar);
      hMu_dphi_eloss_Zdc2[jcent][iside][4][ihar]=(TH2F*)htmp2[7]->Clone(name);

      //2. Cent:10-20
      jcent=1;
      for(int ipt=0;ipt<NPT_Mu;ipt++){
	htmp2[ipt]=(TH2F*)hMu_dphi_eloss_Zdc[2][iside][ipt][ihar]->Clone();
      }
      htmp2[0]->Add(htmp2[1],1);
      htmp2[2]->Add(htmp2[3],1);
      htmp2[4]->Add(htmp2[5],1);
      htmp2[7]->Add(htmp2[8],1);
      sprintf(name,"hMu_dphi_eloss_Zdc2_cent%d_side%d_pt0_har%d",jcent,iside,ihar);
      hMu_dphi_eloss_Zdc2[jcent][iside][0][ihar]=(TH2F*)htmp2[0]->Clone(name);
      sprintf(name,"hMu_dphi_eloss_Zdc2_cent%d_side%d_pt1_har%d",jcent,iside,ihar);
      hMu_dphi_eloss_Zdc2[jcent][iside][1][ihar]=(TH2F*)htmp2[2]->Clone(name);
      sprintf(name,"hMu_dphi_eloss_Zdc2_cent%d_side%d_pt2_har%d",jcent,iside,ihar);
      hMu_dphi_eloss_Zdc2[jcent][iside][2][ihar]=(TH2F*)htmp2[4]->Clone(name);
      sprintf(name,"hMu_dphi_eloss_Zdc2_cent%d_side%d_pt3_har%d",jcent,iside,ihar);
      hMu_dphi_eloss_Zdc2[jcent][iside][3][ihar]=(TH2F*)htmp2[6]->Clone(name);
      sprintf(name,"hMu_dphi_eloss_Zdc2_cent%d_side%d_pt4_har%d",jcent,iside,ihar);
      hMu_dphi_eloss_Zdc2[jcent][iside][4][ihar]=(TH2F*)htmp2[7]->Clone(name);

      //3. Cent:20-40
      jcent=2;
      for(int ipt=0;ipt<NPT_Mu;ipt++){
	htmp2[ipt]=(TH2F*)hMu_dphi_eloss_Zdc[3][iside][ipt][ihar]->Clone();
	htmp2[ipt]->Add(hMu_dphi_eloss_Zdc[4][iside][ipt][ihar],1);
      }
      htmp2[0]->Add(htmp2[1],1);
      htmp2[2]->Add(htmp2[3],1);
      htmp2[4]->Add(htmp2[5],1);
      htmp2[7]->Add(htmp2[8],1);
      sprintf(name,"hMu_dphi_eloss_Zdc2_cent%d_side%d_pt0_har%d",jcent,iside,ihar);
      hMu_dphi_eloss_Zdc2[jcent][iside][0][ihar]=(TH2F*)htmp2[0]->Clone(name);
      sprintf(name,"hMu_dphi_eloss_Zdc2_cent%d_side%d_pt1_har%d",jcent,iside,ihar);
      hMu_dphi_eloss_Zdc2[jcent][iside][1][ihar]=(TH2F*)htmp2[2]->Clone(name);
      sprintf(name,"hMu_dphi_eloss_Zdc2_cent%d_side%d_pt2_har%d",jcent,iside,ihar);
      hMu_dphi_eloss_Zdc2[jcent][iside][2][ihar]=(TH2F*)htmp2[4]->Clone(name);
      sprintf(name,"hMu_dphi_eloss_Zdc2_cent%d_side%d_pt3_har%d",jcent,iside,ihar);
      hMu_dphi_eloss_Zdc2[jcent][iside][3][ihar]=(TH2F*)htmp2[6]->Clone(name);
      sprintf(name,"hMu_dphi_eloss_Zdc2_cent%d_side%d_pt4_har%d",jcent,iside,ihar);
      hMu_dphi_eloss_Zdc2[jcent][iside][4][ihar]=(TH2F*)htmp2[7]->Clone(name);

      //3. Cent:40-60
      jcent=3;
      for(int ipt=0;ipt<NPT_Mu;ipt++){
	htmp2[ipt]=(TH2F*)hMu_dphi_eloss_Zdc[5][iside][ipt][ihar]->Clone();
	htmp2[ipt]->Add(hMu_dphi_eloss_Zdc[6][iside][ipt][ihar],1);
      }
      htmp2[0]->Add(htmp2[1],1);
      htmp2[2]->Add(htmp2[3],1);
      htmp2[4]->Add(htmp2[5],1);
      htmp2[7]->Add(htmp2[8],1);
      sprintf(name,"hMu_dphi_eloss_Zdc2_cent%d_side%d_pt0_har%d",jcent,iside,ihar);
      hMu_dphi_eloss_Zdc2[jcent][iside][0][ihar]=(TH2F*)htmp2[0]->Clone(name);
      sprintf(name,"hMu_dphi_eloss_Zdc2_cent%d_side%d_pt1_har%d",jcent,iside,ihar);
      hMu_dphi_eloss_Zdc2[jcent][iside][1][ihar]=(TH2F*)htmp2[2]->Clone(name);
      sprintf(name,"hMu_dphi_eloss_Zdc2_cent%d_side%d_pt2_har%d",jcent,iside,ihar);
      hMu_dphi_eloss_Zdc2[jcent][iside][2][ihar]=(TH2F*)htmp2[4]->Clone(name);
      sprintf(name,"hMu_dphi_eloss_Zdc2_cent%d_side%d_pt3_har%d",jcent,iside,ihar);
      hMu_dphi_eloss_Zdc2[jcent][iside][3][ihar]=(TH2F*)htmp2[6]->Clone(name);
      sprintf(name,"hMu_dphi_eloss_Zdc2_cent%d_side%d_pt4_har%d",jcent,iside,ihar);
      hMu_dphi_eloss_Zdc2[jcent][iside][4][ihar]=(TH2F*)htmp2[7]->Clone(name);

    }
  }

}

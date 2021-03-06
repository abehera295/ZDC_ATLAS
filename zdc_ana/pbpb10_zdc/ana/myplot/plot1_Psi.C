
void plot::setStyle1(TH1D *h,int itype1){
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

void plot::plot1_Psi(){
  //Read hists
  static const int NCENT2=100;
  TH1D *hPsi2[NCENT2][NSIDE+1][NHAR];
  TH1D *hPsiC2[NCENT2][NSIDE+1][NHAR];
  TH1D *hPsiF2[NCENT2][NSIDE+1][NHAR];

  TFile *fin3;
  fin3=new TFile("../hists_step3/EP169927_flattened.root","read");
  for (int icent=0; icent<NCENT2; icent++) {
    for (int iside=0; iside<NSIDE+1; iside++) {
      for (int ihar=0; ihar<NHAR; ihar++) {
	sprintf(name, "hPsi_zdc_cent%d_side%d_har%d", icent, iside, ihar);
	hPsi2[icent][iside][ihar]=(TH1D*)fin3->Get(name);
	sprintf(name, "hPsiC_zdc_cent%d_side%d_har%d", icent, iside, ihar);
	hPsiC2[icent][iside][ihar]=(TH1D*)fin3->Get(name);
	sprintf(name, "hPsiF_zdc_cent%d_side%d_har%d", icent, iside, ihar);
	hPsiF2[icent][iside][ihar]=(TH1D*)fin3->Get(name);
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
  
  for (int icent=0; icent<NCENT; icent++) {
    for (int iside=0; iside<NSIDE+1; iside++) {
      for (int ihar=0; ihar<NHAR; ihar++) {
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
  
  const int NCan1=5;
  const int NHar1=2;
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
      hPsi[icent_mat[j]][0][ihar]->GetYaxis()->SetRangeUser(0.,5.);
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
      hPsi[icent_mat[j]][1][ihar]->GetYaxis()->SetRangeUser(0.,5.);
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
      hPsi[icent_mat[j]][2][ihar]->GetYaxis()->SetRangeUser(0.,5.);
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
	hDphi_Zdc_fg[icent_mat[j]][ihar]->GetYaxis()->SetRangeUser(0.95,1.05);
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
      hDphi_Zdc_rat[icent_mat[j]][ihar]->Draw("same");
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
      hDphi_FCal_rat[icent_mat[j]][ihar]->Draw("same");
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

}

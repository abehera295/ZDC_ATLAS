
void plot::setStyle1(TH1D *h,int itype1){
  h->GetYaxis()->CenterTitle();       h->GetXaxis()->CenterTitle();  //h->GetZaxis()->CenterTitle();
  h->GetYaxis()->SetTitleOffset(1.2); h->GetXaxis()->SetTitleOffset(1.5);
  //h->GetZaxis()->SetTitleOffset(1.5);
  h->GetXaxis()->SetTitleFont(43);    h->GetXaxis()->SetTitleSize(17);
  h->GetXaxis()->SetLabelFont(43);    h->GetXaxis()->SetLabelSize(15);
  h->GetYaxis()->SetTitleFont(43);    h->GetYaxis()->SetTitleSize(17);
  h->GetYaxis()->SetLabelFont(43);    h->GetYaxis()->SetLabelSize(15);
  //h->GetZaxis()->SetTitleFont(43);    h->GetZaxis()->SetTitleSize(17);
  //h->GetZaxis()->SetLabelFont(43);    h->GetZaxis()->SetLabelSize(15);
  h->GetYaxis()->SetNdivisions(505);
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
  const int NCan=3;
  TCanvas *c[NCan];
  TLegend *leg[NCan];
  for(int i=0;i<NCan;i++){
    sprintf(name,"c%d",i);
    c[i]=new TCanvas(name,name,1100,700);
    c[i]->Divide(3,2);
    leg[i]=new TLegend(0.45,0.55,0.88,0.88);
    leg[i]->SetLineColor(0);
  }
  const int icent_mat[]={0,4,6,8,10,15};
  //Plot for ZDC-sideC
  int i=0;
  for(int j=0;j<6;j++){
    c[i]->cd(j+1);
    if(j<3) gPad->SetBottomMargin(0);
    if(j>2) gPad->SetTopMargin(0);
    if(j!=0 && j!=3) gPad->SetLeftMargin(0);
    gPad->SetRightMargin(0);
    hPsi[icent_mat[j]][0][0]->Draw();
    hPsi[icent_mat[j]][0][0]->SetXTitle("#Psi_{1}^{C}");
    hPsi[icent_mat[j]][0][0]->GetYaxis()->SetRangeUser(0.,5.);
    hPsiC[icent_mat[j]][0][0]->Draw("same");
    hPsiF[icent_mat[j]][0][0]->Draw("same");
    if(j==0){
      leg[i]->AddEntry(hPsi[icent_mat[j]][0][0],"Raw distribution","PL");
      leg[i]->AddEntry(hPsiC[icent_mat[j]][0][0],"after Recentering","PL");
      leg[i]->AddEntry(hPsiF[icent_mat[j]][0][0],"after Flattening","PL");
      leg[i]->Draw("same");
      text.DrawLatex(0.2, 0.7,"Side C");
    }
    sprintf(name,"%s",Cent_mat[icent_mat[j]].c_str());
    text.DrawLatex(0.2, 0.8,name);
  }
  //Plot for ZDC-sideA
  i=1;
  for(int j=0;j<6;j++){
    c[i]->cd(j+1);
    if(j<3) gPad->SetBottomMargin(0);
    if(j>2) gPad->SetTopMargin(0);
    if(j!=0 && j!=3) gPad->SetLeftMargin(0);
    gPad->SetRightMargin(0);
    hPsi[icent_mat[j]][1][0]->Draw();
    hPsi[icent_mat[j]][1][0]->SetXTitle("#Psi_{1}^{A}");
    hPsi[icent_mat[j]][1][0]->GetYaxis()->SetRangeUser(0.,5.);
    hPsiC[icent_mat[j]][1][0]->Draw("same");
    hPsiF[icent_mat[j]][1][0]->Draw("same");
    if(j==0){
      leg[i]->AddEntry(hPsi[icent_mat[j]][1][0],"Raw distribution","PL");
      leg[i]->AddEntry(hPsiC[icent_mat[j]][1][0],"after Recentering","PL");
      leg[i]->AddEntry(hPsiF[icent_mat[j]][1][0],"after Flattening","PL");
      leg[i]->Draw("same");
      text.DrawLatex(0.2, 0.7,"Side A");
    }
    sprintf(name,"%s",Cent_mat[icent_mat[j]].c_str());
    text.DrawLatex(0.2, 0.8,name);
  }

  //Plot Correlation
  i=2;
  for(int j=0;j<6;j++){
    c[i]->cd(j+1);
    if(j<3) gPad->SetBottomMargin(0);
    if(j>2) gPad->SetTopMargin(0);
    if(j!=0 && j!=3) gPad->SetLeftMargin(0);
    gPad->SetRightMargin(0);
    hDphi_fg[icent_mat[j]][0]->Draw();
    hDphi_fg[icent_mat[j]][0]->GetYaxis()->SetRangeUser(0.5,1.5);
    hDphi_fg[icent_mat[j]][0]->SetXTitle("#Psi_{1}^{A}-#Psi_{1}^{C}");
    hDphi_bg[icent_mat[j]][0]->Draw("same");
    hDphi_rat[icent_mat[j]][0]->Draw("same");
    if(j==0){
      leg[i]->AddEntry(hDphi_fg[icent_mat[j]][0],"same event","PL");
      leg[i]->AddEntry(hDphi_bg[icent_mat[j]][0],"mixed event","PL");
      leg[i]->AddEntry(hDphi_rat[icent_mat[j]][0],"correlation","PL");
      leg[i]->Draw("same");
      //text.DrawLatex(0.2, 0.8,"Side A");
    }
    sprintf(name,"%s",Cent_mat[icent_mat[j]].c_str());
    text.DrawLatex(0.2, 0.8,name);
  }

  //Save Canvases
  c[0]->SaveAs("Psi_sideC.pdf");
  c[1]->SaveAs("Psi_sideA.pdf");
  c[2]->SaveAs("dPsi.pdf");
}

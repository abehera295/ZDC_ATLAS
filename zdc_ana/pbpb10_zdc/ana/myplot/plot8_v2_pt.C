void plot::setStyle8(TProfile *h,int itype1,int itype2){
h->GetYaxis()->CenterTitle(0);       h->GetXaxis()->CenterTitle(0);  //h->GetZaxis()->CenterTitle\
();
  h->GetYaxis()->SetTitleOffset(2.5); h->GetXaxis()->SetTitleOffset(2.5);
  //h->GetZaxis()->SetTitleOffset(1.5);
  h->GetXaxis()->SetTitleFont(43);    h->GetXaxis()->SetTitleSize(20);
  h->GetXaxis()->SetLabelFont(43);    h->GetXaxis()->SetLabelSize(17);
  h->GetYaxis()->SetTitleFont(43);    h->GetYaxis()->SetTitleSize(20);
  h->GetYaxis()->SetLabelFont(43);    h->GetYaxis()->SetLabelSize(17);
  //h->GetZaxis()->SetTitleFont(43);    h->GetZaxis()->SetTitleSize(17);
  //h->GetZaxis()->SetLabelFont(43);    h->GetZaxis()->SetLabelSize(15);
  h->GetYaxis()->SetNdivisions(505);
  h->GetXaxis()->SetNdivisions(505);
  int c_mat[]={1,4,2,6,28,36,7,8,9,46,30};
  int m_mat[]={20,25,24,30,27,25,28,33,20,23,29,34};
  h->SetTitle("");

  h->SetMarkerStyle(m_mat[itype1]);  
  h->SetMarkerSize(1.);
  h->SetMarkerColor(c_mat[itype1+1]);
  h->SetLineColor(c_mat[itype1]);

  if(itype1==0){
    h->SetMarkerStyle(24);
    h->SetMarkerColor(4);
    h->SetLineColor(4);
  }
  if(itype1==1){
    h->SetMarkerStyle(25);
    h->SetMarkerColor(2);
    h->SetLineColor(2);
  }
  if(itype1==2){
    h->SetMarkerStyle(26);
    h->SetMarkerColor(1);
    h->SetLineColor(1);
  }
  
  h->SetXTitle("");
  h->SetYTitle("");
  h->GetXaxis()->SetRangeUser(0.5,5.);
  h->GetYaxis()->SetRangeUser(-0.01,0.3);
  //if(itype2>=3) h->GetYaxis()->SetRangeUser(-0.02,0.02);
}

void plot::plot8_v2_pt(){
  const int NCan8=2;
  TCanvas *c8[NCan8];
  TLegend *leg8[NCan8];
  TPad *pad8[NCan8][8];
  for(int i=0;i<NCan8;i++){
    sprintf(name,"c8_%d",i);
    c8[i]=new TCanvas(name,name,1100,700);
    leg8[i]=new TLegend(0.6,0.65,0.88,0.88);
    leg8[i]->SetLineColor(0);
  }
  TLine *line=new TLine(0.5,0.,5.,0.);
  line->SetLineStyle(7);
  line->SetLineWidth(1);
  
  int icent_mat[]={0,1,3,4,5,6};
  int ieta_mat[]={0,2,5,6,9,11,NETA};
  int n1=2;
  int n2=1;
  int ihar=1;
  int iside=2;
  
  int i=0;
  Divide_Pad(c8[i],pad8[i],2,3);
  for(int j=0;j<8;j++){
    c8[i]->cd();
    pad8[i][j]->cd();
    hVn_Zdc_pt_fg[icent_mat[j]][iside][ihar]->Draw();
    setStyle8(hVn_Zdc_pt_fg[icent_mat[j]][iside][ihar],0,j);
    if(j==5) hVn_Zdc_pt_fg[icent_mat[j]][iside][ihar]->SetXTitle("p_{T}");
    if(j==0) hVn_Zdc_pt_fg[icent_mat[j]][iside][ihar]->SetYTitle("v_{2}");
    //hVn_Zdc_pt_fg[icent_mat[j]][iside][ihar]->Draw("same");
    //setStyle8(hVn_Zdc_pt_fg[icent_mat[j]][iside][ihar],1,j);
    //hVn_Zdc_pt_fg[icent_mat[j]][2][ihar]->Draw("same");
    //setStyle8(hVn_Zdc_pt_fg[icent_mat[j]][2][ihar],2,j);
    if(j==0){
      //leg8[i]->AddEntry(hVn_Zdc_pt_fg[icent_mat[j]][iside][ihar],"using #Psi_{1}^{C}","P");
      //leg8[i]->AddEntry(hVn_Zdc_pt_fg[icent_mat[j]][iside][ihar],"using #Psi_{1}^{A}","P");
      //leg8[i]->AddEntry(hVn_Zdc_pt_fg[icent_mat[j]][2][ihar],"using #Psi_{1}^{A+C}","P");
      //leg8[i]->Draw("same");
      //sprintf(name,"%1.1f < #eta <%1.1f",eta_mat,eta_mat[ieta_mat[n2]+1]);
      sprintf(name,"All #eta");
      text.DrawLatex(0.2, 0.7,name);
    }
    sprintf(name,"%s",Cent_mat[icent_mat[j]].c_str());
    text.DrawLatex(0.2, 0.8,name);
    line->Draw("same");
  }

  /*
  i=1;
  Divide_Pad(c8[i],pad8[i],2,3);
  for(int j=0;j<8;j++){
    c8[i]->cd();
    pad8[i][j]->cd();
    hVn_Zdc_pt_fg[icent_mat[n1]][0][ihar][ieta_mat[j]]->Draw();
    setStyle8(hVn_Zdc_pt_fg[icent_mat[n1]][0][ihar][ieta_mat[j]],0,j);
    if(j==5) hVn_Zdc_pt_fg[icent_mat[n1]][0][ihar][ieta_mat[j]]->SetXTitle("p_{T}");
    if(j==0) hVn_Zdc_pt_fg[icent_mat[n1]][0][ihar][ieta_mat[j]]->SetYTitle("v_{2}");
    hVn_Zdc_pt_fg[icent_mat[n1]][1][ihar][ieta_mat[j]]->Draw("same");
    setStyle8(hVn_Zdc_pt_fg[icent_mat[n1]][1][ihar][ieta_mat[j]],1,j);
    hVn_Zdc_pt_fg[icent_mat[n1]][2][ihar][ieta_mat[j]]->Draw("same");
    setStyle8(hVn_Zdc_pt_fg[icent_mat[n1]][2][ihar][ieta_mat[j]],2,j);
    if(j==0){
      leg8[i]->AddEntry(hVn_Zdc_pt_fg[icent_mat[n1]][0][ihar][ieta_mat[j]],"using #Psi_{1}^{C}","P");
      leg8[i]->AddEntry(hVn_Zdc_pt_fg[icent_mat[n1]][1][ihar][ieta_mat[j]],"using #Psi_{1}^{A}","P");
      leg8[i]->AddEntry(hVn_Zdc_pt_fg[icent_mat[n1]][2][ihar][ieta_mat[j]],"using #Psi_{1}^{A+C}","P");
      leg8[i]->Draw("same");
      sprintf(name,"%s",Cent_mat[icent_mat[n1]].c_str());
      text.DrawLatex(0.2, 0.7,name);
    }
    sprintf(name,"%1.1f < #eta <%1.1f",eta_mat[ieta_mat[j]],eta_mat[ieta_mat[j]+1]);
    text.DrawLatex(0.2, 0.8,name);
    line->Draw("same");
  }
  */
  //Save Canvas
  c8[0]->SaveAs("./plots_3PC/V2_pt_allCent.pdf");
  //c8[1]->SaveAs("./plots_3PC/V2_pt_allEta.pdf");
}

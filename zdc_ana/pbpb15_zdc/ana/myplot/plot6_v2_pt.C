void plot::setStyle6(TH1D *h,int itype1,int itype2){
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

void plot::plot6_v2_pt(){
  const int NCan6=2;
  TCanvas *c6[NCan6];
  TLegend *leg6[NCan6];
  TPad *pad6[NCan6][6];
  for(int i=0;i<NCan6;i++){
    sprintf(name,"c6_%d",i);
    c6[i]=new TCanvas(name,name,1100,700);
    leg6[i]=new TLegend(0.6,0.65,0.88,0.88);
    leg6[i]->SetLineColor(0);
  }
  TLine *line=new TLine(0.5,0.,5.,0.);
  line->SetLineStyle(7);
  line->SetLineWidth(1);
  
  int icent_mat[]={0,1,3,4,5,6};
  int ieta_mat[]={0,2,4,7,9,11,NETA,NETA+1};
  int n1=2;
  int n2=6;
  int ihar=1;
  
  int i=0;
  Divide_Pad(c6[i],pad6[i],2,3);
  for(int j=0;j<6;j++){
    c6[i]->cd();
    pad6[i][j]->cd();
    hVn_Zdc_pt_res[icent_mat[j]][0][ihar][ieta_mat[n2]]->Draw();
    setStyle6(hVn_Zdc_pt_res[icent_mat[j]][0][ihar][ieta_mat[n2]],0,j);
    if(j==5) hVn_Zdc_pt_res[icent_mat[j]][0][ihar][ieta_mat[n2]]->SetXTitle("p_{T}");
    if(j==0) hVn_Zdc_pt_res[icent_mat[j]][0][ihar][ieta_mat[n2]]->SetYTitle("v_{2}");
    hVn_Zdc_pt_res[icent_mat[j]][1][ihar][ieta_mat[n2]]->Draw("same");
    setStyle6(hVn_Zdc_pt_res[icent_mat[j]][1][ihar][ieta_mat[n2]],1,j);
    hVn_Zdc_pt_res[icent_mat[j]][2][ihar][ieta_mat[n2]]->Draw("same");
    setStyle6(hVn_Zdc_pt_res[icent_mat[j]][2][ihar][ieta_mat[n2]],2,j);
    if(j==0){
      leg6[i]->AddEntry(hVn_Zdc_pt_res[icent_mat[j]][0][ihar][ieta_mat[n2]],"using #Psi_{1}^{C}","P");
      leg6[i]->AddEntry(hVn_Zdc_pt_res[icent_mat[j]][1][ihar][ieta_mat[n2]],"using #Psi_{1}^{A}","P");
      leg6[i]->AddEntry(hVn_Zdc_pt_res[icent_mat[j]][2][ihar][ieta_mat[n2]],"using #Psi_{1}^{A+C}","P");
      leg6[i]->Draw("same");
      //sprintf(name,"%1.1f < #eta <%1.1f",eta_mat[ieta_mat[n2]],eta_mat[ieta_mat[n2]+1]);
      sprintf(name,"-2.5 < #eta < 2.5");
      text.DrawLatex(0.2, 0.7,name);
    }
    sprintf(name,"%s",Cent_mat[icent_mat[j]].c_str());
    text.DrawLatex(0.2, 0.8,name);
    line->Draw("same");
  }


  i=1;
  Divide_Pad(c6[i],pad6[i],2,3);
  for(int j=0;j<6;j++){
    c6[i]->cd();
    pad6[i][j]->cd();
    hVn_Zdc_pt_res[icent_mat[n1]][0][ihar][ieta_mat[j]]->Draw();
    setStyle6(hVn_Zdc_pt_res[icent_mat[n1]][0][ihar][ieta_mat[j]],0,j);
    if(j==5) hVn_Zdc_pt_res[icent_mat[n1]][0][ihar][ieta_mat[j]]->SetXTitle("p_{T}");
    if(j==0) hVn_Zdc_pt_res[icent_mat[n1]][0][ihar][ieta_mat[j]]->SetYTitle("v_{2}");
    hVn_Zdc_pt_res[icent_mat[n1]][1][ihar][ieta_mat[j]]->Draw("same");
    setStyle6(hVn_Zdc_pt_res[icent_mat[n1]][1][ihar][ieta_mat[j]],1,j);
    hVn_Zdc_pt_res[icent_mat[n1]][2][ihar][ieta_mat[j]]->Draw("same");
    setStyle6(hVn_Zdc_pt_res[icent_mat[n1]][2][ihar][ieta_mat[j]],2,j);
    if(j==0){
      leg6[i]->AddEntry(hVn_Zdc_pt_res[icent_mat[n1]][0][ihar][ieta_mat[j]],"using #Psi_{1}^{C}","P");
      leg6[i]->AddEntry(hVn_Zdc_pt_res[icent_mat[n1]][1][ihar][ieta_mat[j]],"using #Psi_{1}^{A}","P");
      leg6[i]->AddEntry(hVn_Zdc_pt_res[icent_mat[n1]][2][ihar][ieta_mat[j]],"using #Psi_{1}^{A+C}","P");
      leg6[i]->Draw("same");
      sprintf(name,"%s",Cent_mat[icent_mat[n1]].c_str());
      text.DrawLatex(0.2, 0.7,name);
    }
    sprintf(name,"%1.1f < #eta <%1.1f",eta_mat[ieta_mat[j]],eta_mat[ieta_mat[j]+1]);
    text.DrawLatex(0.2, 0.8,name);
    line->Draw("same");
  }
  
  //Save Canvas
  c6[0]->SaveAs("./plots/V2_pt_allCent.pdf");
  c6[1]->SaveAs("./plots/V2_pt_allEta.pdf");
}

void plot::setStyle4(TH1D *h,int itype1,int itype2){
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
  h->GetYaxis()->SetRangeUser(-0.01,0.01);
  //if(itype2>=3) h->GetYaxis()->SetRangeUser(-0.02,0.02);
}

void Find_sqrt_h(TH1D *h){
  int N=h->GetNbinsX();
  for(int i=0;i<N;i++){
    double y=h->GetBinContent(i+1);
    double dy=h->GetBinError(i+1);
    double ycorr=(y/fabs(y)) *sqrt(fabs(y));
    double dycorr=0.5*dy/sqrt(y);
    h->SetBinContent(i+1,ycorr);
    h->SetBinError(i+1,dycorr);
  }
}



void plot::plot4_v1_pt(){
  const int NCan4=4;
  TCanvas *c4[NCan4];
  TLegend *leg4[NCan4];
  TPad *pad4[NCan4][6];
  for(int i=0;i<NCan4;i++){
    sprintf(name,"c4_%d",i);
    c4[i]=new TCanvas(name,name,1100,700);
    leg4[i]=new TLegend(0.6,0.65,0.88,0.88);
    leg4[i]->SetLineColor(0);
  }
  TLine *line=new TLine(0.5,0.,5.,0.);
  line->SetLineStyle(7);
  line->SetLineWidth(1);
  
  int icent_mat[]={0,1,3,4,5,6};
  int ieta_mat[]={0,4,5,6,7,11,NETA,NETA+1};
  int n1=2;
  int n2=7;

  int i=0;
  Divide_Pad(c4[i],pad4[i],2,3);
  for(int j=0;j<6;j++){
    c4[i]->cd();
    pad4[i][j]->cd();
    hVn_Zdc_pt_res[icent_mat[j]][0][0][ieta_mat[n2]]->Draw();
    setStyle4(hVn_Zdc_pt_res[icent_mat[j]][0][0][ieta_mat[n2]],0,j);
    if(j==5) hVn_Zdc_pt_res[icent_mat[j]][0][0][ieta_mat[n2]]->SetXTitle("p_{T}");
    if(j==0) hVn_Zdc_pt_res[icent_mat[j]][0][0][ieta_mat[n2]]->SetYTitle("v_{1}");
    hVn_Zdc_pt_res[icent_mat[j]][1][0][ieta_mat[n2]]->Draw("same");
    setStyle4(hVn_Zdc_pt_res[icent_mat[j]][1][0][ieta_mat[n2]],1,j);
    hVn_Zdc_pt_res[icent_mat[j]][2][0][ieta_mat[n2]]->Draw("same");
    setStyle4(hVn_Zdc_pt_res[icent_mat[j]][2][0][ieta_mat[n2]],2,j);
    if(j==0){
      leg4[i]->AddEntry(hVn_Zdc_pt_res[icent_mat[j]][0][0][ieta_mat[n2]],"using #Psi_{1}^{C}","P");
      leg4[i]->AddEntry(hVn_Zdc_pt_res[icent_mat[j]][1][0][ieta_mat[n2]],"using #Psi_{1}^{A}","P");
      leg4[i]->AddEntry(hVn_Zdc_pt_res[icent_mat[j]][2][0][ieta_mat[n2]],"using #Psi_{1}^{A+C}","P");
      leg4[i]->Draw("same");
      sprintf(name,"-0.8 < #eta < 0.8");
      text.DrawLatex(0.2, 0.7,name);
    }
    sprintf(name,"%s",Cent_mat[icent_mat[j]].c_str());
    text.DrawLatex(0.2, 0.8,name);
    line->Draw("same");
  }


  i=1;
  Divide_Pad(c4[i],pad4[i],2,3);
  for(int j=0;j<6;j++){
    c4[i]->cd();
    pad4[i][j]->cd();
    hVn_Zdc_pt_res[icent_mat[n1]][0][0][ieta_mat[j]]->Draw();
    setStyle4(hVn_Zdc_pt_res[icent_mat[n1]][0][0][ieta_mat[j]],0,j);
    if(j==5) hVn_Zdc_pt_res[icent_mat[n1]][0][0][ieta_mat[j]]->SetXTitle("p_{T}");
    if(j==0) hVn_Zdc_pt_res[icent_mat[n1]][0][0][ieta_mat[j]]->SetYTitle("v_{1}");
    hVn_Zdc_pt_res[icent_mat[n1]][1][0][ieta_mat[j]]->Draw("same");
    setStyle4(hVn_Zdc_pt_res[icent_mat[n1]][1][0][ieta_mat[j]],1,j);
    hVn_Zdc_pt_res[icent_mat[n1]][2][0][ieta_mat[j]]->Draw("same");
    setStyle4(hVn_Zdc_pt_res[icent_mat[n1]][2][0][ieta_mat[j]],2,j);
    if(j==0){
      leg4[i]->AddEntry(hVn_Zdc_pt_res[icent_mat[n1]][0][0][ieta_mat[j]],"using #Psi_{1}^{C}","P");
      leg4[i]->AddEntry(hVn_Zdc_pt_res[icent_mat[n1]][1][0][ieta_mat[j]],"using #Psi_{1}^{A}","P");
      leg4[i]->AddEntry(hVn_Zdc_pt_res[icent_mat[n1]][2][0][ieta_mat[j]],"using #Psi_{1}^{A+C}","P");
      leg4[i]->Draw("same");
      sprintf(name,"%s",Cent_mat[icent_mat[n1]].c_str());
      text.DrawLatex(0.2, 0.7,name);
    }
    sprintf(name,"%1.1f < #eta <%1.1f",eta_mat[ieta_mat[j]],eta_mat[ieta_mat[j]+1]);
    text.DrawLatex(0.2, 0.8,name);
    line->Draw("same");
  }

  i=2;
  Divide_Pad(c4[i],pad4[i],2,3);
  for(int j=0;j<6;j++){
    c4[i]->cd();
    pad4[i][j]->cd();
    hVnEven_pt_res[icent_mat[j]][0][ieta_mat[n2]]->Draw();
    setStyle4(hVnEven_pt_res[icent_mat[j]][0][ieta_mat[n2]],0,j);
    if(j==5) hVnEven_pt_res[icent_mat[j]][0][ieta_mat[n2]]->SetXTitle("p_{T}");
    if(j==0) hVnEven_pt_res[icent_mat[j]][0][ieta_mat[n2]]->SetYTitle("v_{1}");
    hVnOdd_pt_res[icent_mat[j]][0][ieta_mat[n2]]->Draw("same");
    setStyle4(hVnOdd_pt_res[icent_mat[j]][0][ieta_mat[n2]],1,j);
    if(j==0){
      //leg4[i]->AddEntry(hVnEven_pt_res[icent_mat[j]][0][ieta_mat[n2]],"{v_{1}(#Psi_{1}^{A})+v_{1}(#Psi_{1}^{C})}/2","P");
      //leg4[i]->AddEntry(hVnOdd_pt_res[icent_mat[j]][0][ieta_mat[n2]],"{v_{1}(#Psi_{1}^{A})-v_{1}(#Psi_{1}^{C})}/2","P");
      leg4[i]->AddEntry(hVnEven_pt_res[icent_mat[j]][0][ieta_mat[n2]],"v_{1}(Even)","P");
      leg4[i]->AddEntry(hVnOdd_pt_res[icent_mat[j]][0][ieta_mat[n2]],"v_{1}(Odd)","P");
      leg4[i]->Draw("same");
      sprintf(name,"-0.8 < #eta < 0.8");
      text.DrawLatex(0.2, 0.7,name);
    }
    sprintf(name,"%s",Cent_mat[icent_mat[j]].c_str());
    text.DrawLatex(0.2, 0.8,name);
    line->Draw("same");
  }

  i=3;
  Divide_Pad(c4[i],pad4[i],2,3);
  for(int j=0;j<6;j++){
    c4[i]->cd();
    pad4[i][j]->cd();
    hVnEven_pt_res[icent_mat[n1]][0][ieta_mat[j]]->Draw();
    setStyle4(hVnEven_pt_res[icent_mat[n1]][0][ieta_mat[j]],0,j);
    if(j==5) hVnEven_pt_res[icent_mat[n1]][0][ieta_mat[j]]->SetXTitle("p_{T}");
    if(j==0) hVnEven_pt_res[icent_mat[n1]][0][ieta_mat[j]]->SetYTitle("v_{1}");
    hVnOdd_pt_res[icent_mat[n1]][0][ieta_mat[j]]->Draw("same");
    setStyle4(hVnOdd_pt_res[icent_mat[n1]][0][ieta_mat[j]],1,j);
    if(j==0){
      //leg4[i]->AddEntry(hVnEven_pt_res[icent_mat[n1]][0][ieta_mat[j]],"{v_{1}(#Psi_{1}^{A})+v_{1}(#Psi_{1}^{C})}/2","P");
      //leg4[i]->AddEntry(hVnOdd_pt_res[icent_mat[n1]][0][ieta_mat[j]],"{v_{1}(#Psi_{1}^{A})-v_{1}(#Psi_{1}^{C})}/2","P");
      leg4[i]->AddEntry(hVnEven_pt_res[icent_mat[n1]][0][ieta_mat[j]],"v_{1}(Even)","P");
      leg4[i]->AddEntry(hVnOdd_pt_res[icent_mat[n1]][0][ieta_mat[j]],"v_{1}(Odd)","P");
      leg4[i]->Draw("same");
      sprintf(name,"%s",Cent_mat[icent_mat[n1]].c_str());
      text.DrawLatex(0.2, 0.7,name);
    }
    sprintf(name,"%1.1f < #eta <%1.1f",eta_mat[ieta_mat[j]],eta_mat[ieta_mat[j]+1]);
    text.DrawLatex(0.2, 0.8,name);
    line->Draw("same");
  }
  
  //Save Canvas
  c4[0]->SaveAs("./plots/V1_pt_allCent.pdf");
  c4[1]->SaveAs("./plots/V1_pt_allEta.pdf");
  c4[2]->SaveAs("./plots/V1EvenOdd_pt_allCent.pdf");
  c4[3]->SaveAs("./plots/V1EvenOdd_pt_allEta.pdf");
}

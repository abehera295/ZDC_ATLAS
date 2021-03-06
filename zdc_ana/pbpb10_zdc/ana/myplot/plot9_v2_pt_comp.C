void plot::setStyle9(TH1D *h,int itype1,int itype2){
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
    h->SetMarkerStyle(20);
    h->SetMarkerColor(4);
    h->SetLineColor(4);
  }
  if(itype1==1){
    h->SetMarkerStyle(21);
    h->SetMarkerColor(2);
    h->SetLineColor(2);
  }
  if(itype1==2){
    h->SetMarkerStyle(24);
    h->SetMarkerColor(4);
    h->SetLineColor(4);
  }
  if(itype1==3){
    h->SetMarkerStyle(25);
    h->SetMarkerColor(2);
    h->SetLineColor(2);
  }
  
  h->SetXTitle("");
  h->SetYTitle("");
  h->GetXaxis()->SetRangeUser(0.5,5.);
  h->GetYaxis()->SetRangeUser(-0.1,0.3);
  //if(itype2>=3) h->GetYaxis()->SetRangeUser(-0.02,0.02);
}

void plot::plot9_v2_pt_comp(){
  const int NCan9=2;
  TCanvas *c9[NCan9];
  TLegend *leg9[NCan9];
  TPad *pad9[NCan9][6];
  for(int i=0;i<NCan9;i++){
    sprintf(name,"c9_%d",i);
    c9[i]=new TCanvas(name,name,1100,700);
    leg9[i]=new TLegend(0.55,0.5,0.88,0.88);
    leg9[i]->SetLineColor(0);
  }
  TLine *line=new TLine(0.5,0.,5.,0.);
  line->SetLineStyle(7);
  line->SetLineWidth(1);

  TLine *line1=new TLine(0.5,1.,5.,1.);
  line1->SetLineStyle(7);
  line1->SetLineWidth(1);

  int icent_mat[]={0,1,3,4,5,6};
  int ieta_mat[]={0,2,NETA,NETA+1};
  int n1=2;
  int n2=2;
  //FCal has v2,v3,v4,v5
  int ihar=0;//v2
  
  int i=0;
  Divide_Pad(c9[i],pad9[i],2,3);
  for(int j=0;j<6;j++){
    c9[i]->cd();
    pad9[i][j]->cd();
    hVn_FCal_pt_res[icent_mat[j]][2][ihar][ieta_mat[n2]]->Draw();
    setStyle9(hVn_FCal_pt_res[icent_mat[j]][2][ihar][ieta_mat[n2]],0,j);
    if(j==5) hVn_FCal_pt_res[icent_mat[j]][2][ihar][ieta_mat[n2]]->SetXTitle("p_{T}");
    if(j==0) hVn_FCal_pt_res[icent_mat[j]][2][ihar][ieta_mat[n2]]->SetYTitle("v_{n}");
    hVn_FCal_pt_res[icent_mat[j]][2][ihar+1][ieta_mat[n2]]->Draw("same");
    setStyle9(hVn_FCal_pt_res[icent_mat[j]][2][ihar+1][ieta_mat[n2]],1,j);
    hVn_Zdc_pt_res[icent_mat[j]][2][ihar+1][ieta_mat[n2]]->Draw("same");
    setStyle9(hVn_Zdc_pt_res[icent_mat[j]][2][ihar+1][ieta_mat[n2]],2,j);
    hVn_Zdc_pt_res[icent_mat[j]][2][ihar+2][ieta_mat[n2]]->Draw("same");
    setStyle9(hVn_Zdc_pt_res[icent_mat[j]][2][ihar+2][ieta_mat[n2]],3,j);
    if(j==0){
      leg9[i]->AddEntry(hVn_FCal_pt_res[icent_mat[j]][2][ihar][ieta_mat[n2]],"v_{2}(FCal)","P");
      leg9[i]->AddEntry(hVn_FCal_pt_res[icent_mat[j]][2][ihar+1][ieta_mat[n2]],"v_{3}(FCal)","P");
      leg9[i]->AddEntry(hVn_Zdc_pt_res[icent_mat[j]][2][ihar+1][ieta_mat[n2]],"v_{2}(Zdc)","P");
      leg9[i]->AddEntry(hVn_Zdc_pt_res[icent_mat[j]][2][ihar+2][ieta_mat[n2]],"v_{3}(Zdc)","P");
      leg9[i]->Draw("same");
      //sprintf(name,"%1.1f < #eta <%1.1f",eta_mat[ieta_mat[n2]],eta_mat[ieta_mat[n2]+1]);
      sprintf(name,"-2.5 <  #eta < 2.5");
      text.DrawLatex(0.2, 0.7,name);
    }
    sprintf(name,"%s",Cent_mat[icent_mat[j]].c_str());
    text.DrawLatex(0.2, 0.8,name);
    line->Draw("same");
  }

  /*
  i=1;
  Divide_Pad(c9[i],pad9[i],2,3);
  for(int j=0;j<6;j++){
    c9[i]->cd();
    pad9[i][j]->cd();
    hVn_FCal_pt_res[icent_mat[n1]][2][ihar][ieta_mat[j]]->Draw();
    setStyle9(hVn_FCal_pt_res[icent_mat[n1]][2][ihar][ieta_mat[j]],0,j);
    if(j==5) hVn_FCal_pt_res[icent_mat[n1]][2][ihar][ieta_mat[j]]->SetXTitle("p_{T}");
    if(j==0) hVn_FCal_pt_res[icent_mat[n1]][2][ihar][ieta_mat[j]]->SetYTitle("v_{2}");
    hVn_Zdc_pt_res[icent_mat[n1]][2][ihar+1][ieta_mat[j]]->Draw("same");
    setStyle9(hVn_Zdc_pt_res[icent_mat[n1]][2][ihar+1][ieta_mat[j]],1,j);
    if(j==0){
      leg9[i]->AddEntry(hVn_FCal_pt_res[icent_mat[n1]][2][ihar][ieta_mat[j]],"From FCal","P");
      leg9[i]->AddEntry(hVn_Zdc_pt_res[icent_mat[n1]][2][ihar+1][ieta_mat[j]],"From Zdc","P");
      leg9[i]->Draw("same");
      sprintf(name,"%s",Cent_mat[icent_mat[n1]].c_str());
      text.DrawLatex(0.2, 0.7,name);
    }
    sprintf(name,"%1.1f < #eta <%1.1f",eta_mat[ieta_mat[j]],eta_mat[ieta_mat[j]+1]);
    text.DrawLatex(0.2, 0.8,name);
    line->Draw("same");
  }
  */
  //Draw Ratios plot
  TH1D *hR[6][2];
  i=1;
  Divide_Pad(c9[i],pad9[i],2,3);
  for(int j=0;j<6;j++){
    c9[i]->cd();
    pad9[i][j]->cd();
    hR[j][0]=(TH1D*)hVn_Zdc_pt_res[icent_mat[j]][2][ihar+1][ieta_mat[n2]]->Clone();
    hR[j][0]->Divide(hVn_FCal_pt_res[icent_mat[j]][2][ihar][ieta_mat[n2]]);
    hR[j][1]=(TH1D*)hVn_Zdc_pt_res[icent_mat[j]][2][ihar+2][ieta_mat[n2]]->Clone();
    hR[j][1]->Divide(hVn_FCal_pt_res[icent_mat[j]][2][ihar+1][ieta_mat[n2]]);

    hR[j][0]->Draw();
    setStyle9(hR[j][0],0,j);
    setStyle9(hR[j][1],1,j);
    hR[j][0]->GetYaxis()->SetRangeUser(-0.5,1.5);
    if(j==5) hR[j][0]->SetXTitle("p_{T}");
    if(j==0) hR[j][0]->SetYTitle("v_{n}");
    hR[j][1]->Draw("same");
    if(j==0){
      leg9[i]->AddEntry(hR[j][0],"v_{2} Zdc/v_{2} FCal","P");
      leg9[i]->AddEntry(hR[j][1],"v_{3} Zdc/v_{3} FCal","P");
      leg9[i]->Draw("same");
      //sprintf(name,"%1.1f < #eta <%1.1f",eta_mat[ieta_mat[n2]],eta_mat[ieta_mat[n2]+1]);
      sprintf(name,"-2.5 <  #eta < 2.5");
      text.DrawLatex(0.2, 0.7,name);
    }
    sprintf(name,"%s",Cent_mat[icent_mat[j]].c_str());
    text.DrawLatex(0.2, 0.8,name);
    line1->Draw("same");
  }
 
  //Save Canvas
  c9[0]->SaveAs("./plots/FCal_V2_pt_comp_allCent.pdf");
  //c9[1]->SaveAs("./plots/FCal_V2_pt_comp_allEta.pdf");
  c9[1]->SaveAs("./plots/Ratio_V2_pt_comp.pdf");

}

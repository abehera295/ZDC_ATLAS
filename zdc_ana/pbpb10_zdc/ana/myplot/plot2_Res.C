void plot::setStyle2(TH1D *h,int itype1){
  h->GetYaxis()->CenterTitle(0);       h->GetXaxis()->CenterTitle(0);  //h->GetZaxis()->CenterTitle();
  h->GetYaxis()->SetTitleOffset(1.); h->GetXaxis()->SetTitleOffset(1.);
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
  h->SetMarkerColor(c_mat[itype1+1]);
  h->SetLineColor(c_mat[itype1+1]);

  //h->GetYaxis()->SetRangeUser(0.,5.);
}

void plot::plot2_Res(){
  TCanvas *c2[4][2];
  TLegend *leg2[2][2][2];
  for(int i=0;i<4;i++){
    for(int j=0;j<2;j++){
      sprintf(name,"c2_%d_%d",i,j);
      c2[i][j]=new TCanvas(name,name,1000,400);
      c2[i][j]->Divide(2,1);
    }
  }
  
  for(int ihar=0;ihar<2;ihar++){
    //Set-1
    setStyle2(hResCent_flat_FCal[ihar],0);
    setStyle2(hResCent_mix_FCal[ihar],1);
    setStyle2(hResCent_flat_FCal_sub[ihar],0);
    setStyle2(hResCent_mix_FCal_sub[ihar],1);
    //Set-2
    setStyle2(hResCent_flat_Zdc[ihar],0);
    setStyle2(hResCent_mix_Zdc[ihar],1);
    setStyle2(hResCent_flat_Zdc_sub[ihar],0);
    setStyle2(hResCent_mix_Zdc_sub[ihar],1);
    //Set-3
    setStyle2(hResXCent_mix_FCal_ful[ihar],0);
    setStyle2(hResXCent_mix_FCal_sub[ihar],1);
    setStyle2(hResXCent_mix_Zdc_ful[ihar],0);
    setStyle2(hResXCent_mix_Zdc_sub[ihar],1);
    //Set-4
    setStyle2(hResCent_mix_FCal_ful[ihar],0);
    setStyle2(hResCent_mix_FCal_sub[ihar],1);
    setStyle2(hResCent_mix_Zdc_ful[ihar],0);
    setStyle2(hResCent_mix_Zdc_sub[ihar],1);
  }

  for(int ihar=0;ihar<2;ihar++){
    //Type-1
    //For FCal
    c2[0][ihar]->cd(1);
    leg2[0][ihar][0]=new TLegend(0.15,0.25,0.5,0.5);
    leg2[0][ihar][0]->SetLineColor(0);
    leg2[0][ihar][0]->AddEntry(hResCent_flat_FCal[ihar],"Same Ev","P");
    leg2[0][ihar][0]->AddEntry(hResCent_mix_FCal[ihar],"Same Ev - Mixed Ev","P");
    //gPad->SetRightMargin(0.02);
    hResCent_flat_FCal[ihar]->Draw();
    hResCent_flat_FCal[ihar]->SetXTitle("Centrality %");
    sprintf(name,"<cos%d(#Psi_{%d}^{A}-#Psi_{%d}^{C})>",ihar+2,ihar+2,ihar+2);
    hResCent_flat_FCal[ihar]->SetYTitle(name);
    hResCent_mix_FCal[ihar]->Draw("same");
    leg2[0][ihar][0]->Draw("same");
  
    c2[0][ihar]->cd(2);
    leg2[0][ihar][1]=new TLegend(0.15,0.25,0.5,0.5);
    leg2[0][ihar][1]->SetLineColor(0);
    leg2[0][ihar][1]->AddEntry(hResCent_flat_FCal_sub[ihar],"Same Ev","P");
    leg2[0][ihar][1]->AddEntry(hResCent_mix_FCal_sub[ihar],"Same Ev - Mixed Ev","P");
    //gPad->SetRightMargin(0.02);
    hResCent_flat_FCal_sub[ihar]->Draw();
    hResCent_flat_FCal_sub[ihar]->SetXTitle("Centrality %");
    sprintf(name,"Res#left{%d#Psi_{%d}#right}",ihar+2,ihar+2);
    hResCent_flat_FCal_sub[ihar]->SetYTitle(name);
    hResCent_mix_FCal_sub[ihar]->Draw("same");
    leg2[0][ihar][1]->Draw("same");

    //For Zdc
    c2[1][ihar]->cd(1);
    leg2[1][ihar][0]=new TLegend(0.5,0.7,0.88,0.88);
    if(ihar==0) leg2[1][ihar][0]=new TLegend(0.25,0.7,0.7,0.88);
    leg2[1][ihar][0]->SetLineColor(0);
    leg2[1][ihar][0]->AddEntry(hResCent_flat_Zdc[ihar],"Same Ev","P");
    leg2[1][ihar][0]->AddEntry(hResCent_mix_Zdc[ihar],"Same Ev - Mixed Ev","P");
    //gPad->SetRightMargin(0.02);
    hResCent_flat_Zdc[ihar]->Draw();
    if(ihar==1) hResCent_flat_Zdc[ihar]->GetYaxis()->SetRangeUser(0.,0.01);
    hResCent_flat_Zdc[ihar]->SetXTitle("Centrality %");
    sprintf(name,"<cos(#Psi_{1}^{A}-#Psi_{1}^{C})>");
    if(ihar>0) sprintf(name,"<cos%d(#Psi_{1}^{A}-#Psi_{1}^{C})>",ihar+1);
    hResCent_flat_Zdc[ihar]->SetYTitle(name);
    hResCent_mix_Zdc[ihar]->Draw("same");
    leg2[1][ihar][0]->Draw("same");
  
    c2[1][ihar]->cd(2);
    leg2[1][ihar][1]=new TLegend(0.25,0.25,0.7,0.45);
    if(ihar==0) leg2[1][ihar][1]=new TLegend(0.25,0.3,0.7,0.5);
    leg2[1][ihar][1]->SetLineColor(0);
    leg2[1][ihar][1]->AddEntry(hResCent_flat_Zdc_sub[ihar],"Same Ev","P");
    leg2[1][ihar][1]->AddEntry(hResCent_mix_Zdc_sub[ihar],"Same Ev - Mixed Ev","P");
    //gPad->SetRightMargin(0.02);
    hResCent_flat_Zdc_sub[ihar]->Draw();
    if(ihar==1) hResCent_flat_Zdc_sub[ihar]->GetYaxis()->SetRangeUser(0.,0.1);
    hResCent_flat_Zdc_sub[ihar]->SetXTitle("Centrality %");
    sprintf(name,"Res#left{#Psi_{1}#right}");
    if(ihar>0) sprintf(name,"Res#left{%d#Psi_{1}#right}",ihar+1);
    hResCent_flat_Zdc_sub[ihar]->SetYTitle(name);
    hResCent_mix_Zdc_sub[ihar]->Draw("same");
    leg2[1][ihar][1]->Draw("same");

    //Type-2
    //For FCal
    c2[2][ihar]->cd(1);
    leg2[2][ihar][0]=new TLegend(0.2,0.15,0.4,0.3);
    leg2[2][ihar][0]->SetLineColor(0);
    leg2[2][ihar][0]->AddEntry(hResXCent_mix_FCal_sub[ihar],"Sub Detector","P");
    leg2[2][ihar][0]->AddEntry(hResXCent_mix_FCal_ful[ihar],"Full Detector","P");
    //gPad->SetRightMargin(0.02);
    hResXCent_mix_FCal_ful[ihar]->Draw();
    hResXCent_mix_FCal_ful[ihar]->SetXTitle("Centrality %");
    //sprintf(name,"<cos%d(#Psi_{%d}^{A}-#Psi_{%d}^{C})>",ihar+2,ihar+2,ihar+2);
    sprintf(name,"#chi_{%d}",ihar+2);
    hResXCent_mix_FCal_ful[ihar]->SetYTitle(name);
    hResXCent_mix_FCal_sub[ihar]->Draw("same");
    leg2[2][ihar][0]->Draw("same");
  
    c2[2][ihar]->cd(2);
    leg2[2][ihar][1]=new TLegend(0.2,0.15,0.4,0.3);
    leg2[2][ihar][1]->SetLineColor(0);
    leg2[2][ihar][1]->AddEntry(hResCent_mix_FCal_sub[ihar],"Sub Detector","P");
    leg2[2][ihar][1]->AddEntry(hResCent_mix_FCal_ful[ihar],"Full Detector","P");
    //gPad->SetRightMargin(0.02);
    hResCent_mix_FCal_ful[ihar]->Draw();
    hResCent_mix_FCal_ful[ihar]->SetXTitle("Centrality %");
    sprintf(name,"Res#left{%d#Psi_{%d}#right}",ihar+2,ihar+2);
    hResCent_mix_FCal_ful[ihar]->SetYTitle(name);
    hResCent_mix_FCal_sub[ihar]->Draw("same");
    leg2[2][ihar][1]->Draw("same");

    //For Zdc
    c2[3][ihar]->cd(1);
    leg2[3][ihar][0]=new TLegend(0.2,0.15,0.4,0.3);
    leg2[3][ihar][0]->SetLineColor(0);
    leg2[3][ihar][0]->AddEntry(hResXCent_mix_Zdc_sub[ihar],"Sub Detector","P");
    leg2[3][ihar][0]->AddEntry(hResXCent_mix_Zdc_ful[ihar],"Full Detector","P");
    //gPad->SetRightMargin(0.02);
    hResXCent_mix_Zdc_ful[ihar]->Draw();
    hResXCent_mix_Zdc_ful[ihar]->SetMinimum(0.);
    hResXCent_mix_Zdc_ful[ihar]->SetXTitle("Centrality %");
    //sprintf(name,"<cos(#Psi_{1}^{A}-#Psi_{1}^{C})>");
    sprintf(name,"#chi_{%d}",ihar+1);
    hResXCent_mix_Zdc_ful[ihar]->SetYTitle(name);
    hResXCent_mix_Zdc_sub[ihar]->Draw("same");
    leg2[3][ihar][0]->Draw("same");
  
    c2[3][ihar]->cd(2);
    leg2[3][ihar][1]=new TLegend(0.2,0.15,0.4,0.3);
    leg2[3][ihar][1]->SetLineColor(0);
    leg2[3][ihar][1]->AddEntry(hResCent_mix_Zdc_sub[ihar],"Sub Detector","P");
    leg2[3][ihar][1]->AddEntry(hResCent_mix_Zdc_ful[ihar],"Full Detector","P");
    //gPad->SetRightMargin(0.02);
    hResCent_mix_Zdc_ful[ihar]->Draw();
    hResCent_mix_Zdc_ful[ihar]->SetMinimum(0.);
    hResCent_mix_Zdc_ful[ihar]->SetXTitle("Centrality %");
    sprintf(name,"Res#left{#Psi_{1}#right}");
    if(ihar>0) sprintf(name,"Res#left{%d#Psi_{1}#right}",ihar+1);
    hResCent_mix_Zdc_ful[ihar]->SetYTitle(name);
    hResCent_mix_Zdc_sub[ihar]->Draw("same");
    leg2[3][ihar][1]->Draw("same");

    
    //Save hist
    sprintf(name,"./plots/FCal_Res_sub_v%d.pdf",ihar+2);
    c2[0][ihar]->SaveAs(name);
    sprintf(name,"./plots/Zdc_Res_sub_v%d.pdf",ihar+1);
    c2[1][ihar]->SaveAs(name);
    sprintf(name,"./plots/FCal_Res_ful_v%d.pdf",ihar+2);
    c2[2][ihar]->SaveAs(name);
    sprintf(name,"./plots/Zdc_Res_ful_v%d.pdf",ihar+1);
    c2[3][ihar]->SaveAs(name);
  }
}

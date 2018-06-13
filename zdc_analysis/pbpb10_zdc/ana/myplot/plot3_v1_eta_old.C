void plot::setStyle3(TH1D *h,int itype1){
h->GetYaxis()->CenterTitle(0);       h->GetXaxis()->CenterTitle(0);  //h->GetZaxis()->CenterTitle\
();
  h->GetYaxis()->SetTitleOffset(1.); h->GetXaxis()->SetTitleOffset(1.);
  //h->GetZaxis()->SetTitleOffset(1.5);
  h->GetXaxis()->SetTitleFont(43);    h->GetXaxis()->SetTitleSize(20);
  h->GetXaxis()->SetLabelFont(43);    h->GetXaxis()->SetLabelSize(17);
  h->GetYaxis()->SetTitleFont(43);    h->GetYaxis()->SetTitleSize(20);
  h->GetYaxis()->SetLabelFont(43);    h->GetYaxis()->SetLabelSize(17);
  //h->GetZaxis()->SetTitleFont(43);    h->GetZaxis()->SetTitleSize(17);
  //h->GetZaxis()->SetLabelFont(43);    h->GetZaxis()->SetLabelSize(15);
  h->GetYaxis()->SetNdivisions(505);
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
    h->SetMarkerStyle(25);
    h->SetMarkerColor(2);
    h->SetLineColor(2);
  }
  
  h->SetXTitle("#eta");
  h->SetYTitle("v_{1}");
  h->GetYaxis()->SetRangeUser(-0.005,0.005);
}

void plot::plot3_v1_eta(){
  /*const int NCan=4;
  TCanvas *c3[NCan];
  TLegend *leg3[NCan];
  for(int i=0;i<NCan;i++){
    sprintf(name,"c3_%d",i);
    c3[i]=new TCanvas(name,name,1100,700);
    c3[i]->Divide(3,2);
    leg3[i]=new TLegend(0.7,0.75,0.88,0.88);
    leg3[i]->SetLineColor(0);
  }
  const int ipt_mat[]={0,10};
  int icent_mat[]={21,22,23,32};
  for(int i=0;i<NCan;i++){
    c3[i]->cd();
    hV1_eta_res[icent_mat[i]][0]->Draw();
    setStyle3(hV1_eta_res[icent_mat[i]][0],0);
    hV1_eta_res[icent_mat[i]][1]->Draw("same");
    setStyle3(hV1_eta_res[icent_mat[i]][1],1);
    leg3[i]->AddEntry(hV1_eta_res[icent_mat[i]][0],"using #Psi_{1}^{C}","PL");
    leg3[i]->AddEntry(hV1_eta_res[icent_mat[i]][1],"using #Psi_{1}^{A}","PL");
    leg3[i]->Draw("same");
    sprintf(name,"%s",Cent_mat_new[icent_mat[i]].c_str());
    text.DrawLatex(0.2,0.8,name);
    sprintf(name,"V1_eta_res_cent%d.pdf",icent_mat[i]);
    c3[i]->SaveAs(name);
  }
  */

}

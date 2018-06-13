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
  int c_mat[]={1,2,4,6,28,36,7,8,9,46,30};
  int m_mat[]={24,30,27,25,28,33,20,23,29,34};
  h->SetTitle("");

  h->SetMarkerStyle(20);
  h->SetMarkerSize(1.);
  h->SetMarkerColor(c_mat[itype1]);
  h->SetLineColor(c_mat[itype1]);

  //h->GetYaxis()->SetRangeUser(0.,5.);
}

void plot::plot2_Res(){
  TCanvas *c2;
  TLegend *leg2;
  sprintf(name,"c2");
  c2=new TCanvas(name,name,1000,400);
  c2->Divide(2,1);

  TH1D *h;
  h=(TH1D*)hResCent_flat[0]->Clone();
  int N=h->GetNbinsX();
  for(int i=0;i<N;i++){
    double y=hResCent_flat[0]->GetBinContent(i+1);
    double dy=hResCent_flat[0]->GetBinError(i+1);
    h->SetBinContent(i+1,sqrt(fabs(y)));
    h->SetBinError(i+1,0.5*dy/(h->GetBinContent(i+1)));
  }
  setStyle2(hResCent_mix[0],2);
  setStyle2(h,2);
  c2->cd(1);
  gPad->SetRightMargin(0.);
  hResCent_mix[0]->Draw();
  hResCent_mix[0]->SetXTitle("Centrality %");
  hResCent_mix[0]->SetYTitle("<cos(#Psi_{1}^{A}-#Psi_{1}^{C})>");
  c2->cd(2);
  gPad->SetRightMargin(0.);
  h->Draw();
  h->SetXTitle("Centrality %");
  h->SetYTitle("Resolution");
  c2->SaveAs("Res.pdf");

}

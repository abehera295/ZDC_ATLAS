#include <iostream>
#include <fstream>
#include <iomanip>
#include <fstream>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLine.h"
#include "THStack.h"
#include "TLatex.h"
#include "string.h"

enum {
  NDET = 16,
  NSIDE   = 2,
  NMODULE = 2,
  NGAIN   = 2,
  NCH1    = 4,
  NCH2    = 24,
  NHAR    = 4,
  //NCENT   = 20,
  NCENT   = 11,
  NCENT_1p  = 100,
  //NK      = 12,
  NSAMP   = 10,
  //STORE_DEP = 20,
  STORE_DEP = 5,
  NVZ_EP=10,
  NPT_EP=10,
  NCH_EP=2
};

char name[300];
//Calibration hists
TH1D *hPixel_Energy[NSIDE+1][NCENT][NGAIN][3];
TH1D *hStrip_Energy[NSIDE+1][NCENT][NGAIN][3];
TH2F *hStrip_Pixel_Energy[NSIDE+1][NCENT][NGAIN][3];
TH2F *hStrip_Pixel_Energy_allCent[NSIDE+1][NGAIN][3];
TH2F *hGrid_Energy[NSIDE+1][NCENT][NGAIN][3];
TProfile2D *hGrid_Energy2[NSIDE+1][NCENT][NGAIN][3];
TH2F *hPixel_LgHg[NSIDE+1][NCENT][3];
TH2F *hStrip_LgHg[NSIDE+1][NCENT][3];
TLatex tt;
void plot1();
void plot2();
void plot3();

void setStyle(TH1D *h,int itype1){
  h->GetYaxis()->CenterTitle(0);       h->GetXaxis()->CenterTitle(0);  //h->GetZaxis()->CenterTitle();
  h->GetYaxis()->SetTitleOffset(1.); h->GetXaxis()->SetTitleOffset(1.);
  //h->GetZaxis()->SetTitleOffset(1.5);
  h->GetXaxis()->SetTitleFont(43);    h->GetXaxis()->SetTitleSize(25);
  h->GetXaxis()->SetLabelFont(43);    h->GetXaxis()->SetLabelSize(15);
  h->GetYaxis()->SetTitleFont(43);    h->GetYaxis()->SetTitleSize(25);
  h->GetYaxis()->SetLabelFont(43);    h->GetYaxis()->SetLabelSize(15);
  //h->GetZaxis()->SetTitleFont(43);    h->GetZaxis()->SetTitleSize(17);
  //h->GetZaxis()->SetLabelFont(43);    h->GetZaxis()->SetLabelSize(15);
  h->GetYaxis()->SetNdivisions(505);
  h->GetXaxis()->SetNdivisions(505);
  int c_mat[]={4,2,4,2,6,28,36,7,8,9,46,30};
  int m_mat[]={20,20,25,25,25,21,30,27,25,28,33,20,23,29,34};
  h->SetTitle("");

  h->SetMarkerStyle(m_mat[itype1]);
  h->SetMarkerColor(c_mat[itype1]);
  h->SetLineColor(c_mat[itype1]);
  h->SetMarkerSize(1.5);

  //h->GetYaxis()->SetRangeUser(0.,5.);
}


void plot_zdc_channels(){
  
  TFile *fin[3];
  fin[0]=new TFile("./Output_zdc_2010.root","read");
  fin[1]=new TFile("./Output_zdc_2015.root","read");
  fin[2]=new TFile("./Output_zdc_2011.root","read");
  
  for(int i=0;i<3;i++){
    for(int iside=0;iside<NSIDE+1;iside++){
      for(int icent=0;icent<NCENT;icent++){
	for(int ig=0;ig<NGAIN;ig++){
	  sprintf(name,"hPixel_Energy_side%d_cent%d_gain%d",iside,icent,ig);
	  hPixel_Energy[iside][icent][ig][i]=(TH1D*)fin[i]->Get(name);
	  sprintf(name,"hStrip_Energy_side%d_cent%d_gain%d",iside,icent,ig);
	  hStrip_Energy[iside][icent][ig][i]=(TH1D*)fin[i]->Get(name);
	
	  sprintf(name,"hStrip_Pixel_Energy_side%d_cent%d_gain%d",iside,icent,ig);
	  hStrip_Pixel_Energy[iside][icent][ig][i]=(TH2F*)fin[i]->Get(name);
	  sprintf(name,"hGrid_Energy_side%d_cent%d_gain%d",iside,icent,ig);
	  hGrid_Energy[iside][icent][ig][i]=(TH2F*)fin[i]->Get(name);
	  sprintf(name,"hGrid_Energy2_side%d_cent%d_gain%d",iside,icent,ig);
	  hGrid_Energy2[iside][icent][ig][i]=(TProfile2D*)fin[i]->Get(name);
	}
	sprintf(name,"hPixel_LgHg_side%d_cent%d",iside,icent);
	hPixel_LgHg[iside][icent][i]=(TH2F*)fin[i]->Get(name);
	sprintf(name,"hStrip_LgHg_side%d_cent%d",iside,icent);
	hStrip_LgHg[iside][icent][i]=(TH2F*)fin[i]->Get(name);
      }
      for(int ig=0;ig<NGAIN;ig++){
	sprintf(name,"hStrip_Pixel_Energy_allCent_side%d_gain%d",iside,ig);
	hStrip_Pixel_Energy_allCent[iside][ig][i]=(TH2F*)fin[i]->Get(name);
      }
    }
  }

  tt.SetNDC(1);
  tt.SetTextFont(43);
  tt.SetTextSize(25);
  
  plot1();//2D StripvsPixel Energy
  plot2();//2D Grid Energy
  plot3();//1D Pixel Energy
}

void plot1(){
  //Plot
  gStyle->SetOptStat(0);
  TLine *line0=new TLine(0.,0.,5.,0.);
  line0->SetLineStyle(7);
  line0->SetLineWidth(1);

  const int NCan=4;
  TCanvas *c[NGAIN][NCan];
  for(int ig=0;ig<NGAIN;ig++){
    for(int ican=0;ican<NCan;ican++){
      sprintf(name,"c_gain%d_%d",ig,ican);
      c[ig][ican]=new TCanvas(name,name,800,600);
    }
  }
  int ican,id,iside;
  string side_text[2]={"Side C","Side A"};

  for(int ig=0;ig<NGAIN;ig++){
    ican=0,iside=0,id=0;
    c[ig][ican]->cd();
    gPad->SetLogz();
    hStrip_Pixel_Energy_allCent[iside][ig][id]->Draw("colz");
    if(ig==0){
      hStrip_Pixel_Energy_allCent[iside][ig][id]->GetXaxis()->SetRangeUser(0.,4000.);
      hStrip_Pixel_Energy_allCent[iside][ig][id]->GetYaxis()->SetRangeUser(0.,4000.);
    }
    if(ig==1){
      hStrip_Pixel_Energy_allCent[iside][ig][id]->GetXaxis()->SetRangeUser(0.,8000.);
      hStrip_Pixel_Energy_allCent[iside][ig][id]->GetYaxis()->SetRangeUser(0.,24000.);
    }
    hStrip_Pixel_Energy_allCent[iside][ig][id]->GetYaxis()->SetTitleOffset(1.2);
    hStrip_Pixel_Energy_allCent[iside][ig][id]->SetXTitle("Sum of 3 Strip Energies");
    hStrip_Pixel_Energy_allCent[iside][ig][id]->SetYTitle("Sum of Pixel Energies");
    tt.DrawLatex(0.2,0.8,side_text[iside].c_str());
    if(ig==0) tt.DrawLatex(0.6,0.8,"2010 Pb+Pb LG");
    if(ig==1) tt.DrawLatex(0.6,0.8,"2010 Pb+Pb HG");
    
    ican=1,iside=1,id=0;
    c[ig][ican]->cd();
    gPad->SetLogz();
    hStrip_Pixel_Energy_allCent[iside][ig][id]->Draw("colz");
    if(ig==0){
      hStrip_Pixel_Energy_allCent[iside][ig][id]->GetXaxis()->SetRangeUser(0.,4000.);
      hStrip_Pixel_Energy_allCent[iside][ig][id]->GetYaxis()->SetRangeUser(0.,4000.);
    }
    if(ig==1){
      hStrip_Pixel_Energy_allCent[iside][ig][id]->GetXaxis()->SetRangeUser(0.,8000.);
      hStrip_Pixel_Energy_allCent[iside][ig][id]->GetYaxis()->SetRangeUser(0.,24000.);
    }
    hStrip_Pixel_Energy_allCent[iside][ig][id]->GetYaxis()->SetTitleOffset(1.2);
    hStrip_Pixel_Energy_allCent[iside][ig][id]->SetXTitle("Sum of 3 Strip Energies");
    hStrip_Pixel_Energy_allCent[iside][ig][id]->SetYTitle("Sum of Pixel Energies");
    tt.DrawLatex(0.2,0.8,side_text[iside].c_str());
    if(ig==0) tt.DrawLatex(0.6,0.8,"2010 Pb+Pb LG");
    if(ig==1) tt.DrawLatex(0.6,0.8,"2010 Pb+Pb HG");

    
    ican=2,iside=0,id=1;
    c[ig][ican]->cd();
    gPad->SetLogz();
    hStrip_Pixel_Energy_allCent[iside][ig][id]->Draw("colz");
    if(ig==0){
      hStrip_Pixel_Energy_allCent[iside][ig][id]->GetXaxis()->SetRangeUser(0.,4000.);
      hStrip_Pixel_Energy_allCent[iside][ig][id]->GetYaxis()->SetRangeUser(0.,4000.);
    }
    if(ig==1){
      hStrip_Pixel_Energy_allCent[iside][ig][id]->GetXaxis()->SetRangeUser(0.,8000.);
      hStrip_Pixel_Energy_allCent[iside][ig][id]->GetYaxis()->SetRangeUser(0.,24000.);
    }
    hStrip_Pixel_Energy_allCent[iside][ig][id]->GetYaxis()->SetTitleOffset(1.2);
    hStrip_Pixel_Energy_allCent[iside][ig][id]->SetXTitle("Sum of 3 Strip Energies");
    hStrip_Pixel_Energy_allCent[iside][ig][id]->SetYTitle("Sum of Pixel Energies");
    tt.DrawLatex(0.2,0.8,side_text[iside].c_str());
    if(ig==0) tt.DrawLatex(0.6,0.8,"2015 Pb+Pb LG");
    if(ig==1) tt.DrawLatex(0.6,0.8,"2015 Pb+Pb HG");


    ican=3,iside=1,id=1;
    c[ig][ican]->cd();
    gPad->SetLogz();
    hStrip_Pixel_Energy_allCent[iside][ig][id]->Draw("colz");
    if(ig==0){
      hStrip_Pixel_Energy_allCent[iside][ig][id]->GetXaxis()->SetRangeUser(0.,4000.);
      hStrip_Pixel_Energy_allCent[iside][ig][id]->GetYaxis()->SetRangeUser(0.,4000.);
    }
    if(ig==1){
      hStrip_Pixel_Energy_allCent[iside][ig][id]->GetXaxis()->SetRangeUser(0.,8000.);
      hStrip_Pixel_Energy_allCent[iside][ig][id]->GetYaxis()->SetRangeUser(0.,24000.);
    }
    hStrip_Pixel_Energy_allCent[iside][ig][id]->GetYaxis()->SetTitleOffset(1.2);
    hStrip_Pixel_Energy_allCent[iside][ig][id]->SetXTitle("Sum of 3 Strip Energies");
    hStrip_Pixel_Energy_allCent[iside][ig][id]->SetYTitle("Sum of Pixel Energies");
    tt.DrawLatex(0.2,0.8,side_text[iside].c_str());
    if(ig==0) tt.DrawLatex(0.6,0.8,"2015 Pb+Pb LG");
    if(ig==1) tt.DrawLatex(0.6,0.8,"2015 Pb+Pb HG");
  }
  //Save canvas
  c[0][0]->SaveAs("./plots/Strip_Pixel_Energy_sideC_2010_LG.pdf");
  c[0][1]->SaveAs("./plots/Strip_Pixel_Energy_sideA_2010_LG.pdf");
  c[0][2]->SaveAs("./plots/Strip_Pixel_Energy_sideC_2015_LG.pdf");
  c[0][3]->SaveAs("./plots/Strip_Pixel_Energy_sideA_2015_LG.pdf");
  c[1][0]->SaveAs("./plots/Strip_Pixel_Energy_sideC_2010_HG.pdf");
  c[1][1]->SaveAs("./plots/Strip_Pixel_Energy_sideA_2010_HG.pdf");
  c[1][2]->SaveAs("./plots/Strip_Pixel_Energy_sideC_2015_HG.pdf");
  c[1][3]->SaveAs("./plots/Strip_Pixel_Energy_sideA_2015_HG.pdf");

}


void plot2(){
  //Plot
  gStyle->SetOptStat(0);
  TLine *line0=new TLine(0.,0.,5.,0.);
  line0->SetLineStyle(7);
  line0->SetLineWidth(1);

  const int NCan=2;
  TCanvas *c2[NGAIN][NCan];
  for(int ig=0;ig<NGAIN;ig++){
    for(int ican=0;ican<NCan;ican++){
      sprintf(name,"c2_gain%d_%d",ig,ican);
      c2[ig][ican]=new TCanvas(name,name,1200,800);
      c2[ig][ican]->Divide(2,2);
    }
  }
  int ican,id,iside,icent,ip;
  int cent_id[]={3,5,7};
  string side_text[2]={"Side C","Side A"};
  string cent_text[]={"20-30%","40-50%","60-70%"};
  string data_text[2]={"2010 Pb+Pb","2015 Pb+Pb"};
  
  for(int ig=0;ig<NGAIN;ig++){
    for(int ican=0;ican<2;ican++){
      iside=0,ip=0;
      for(int id=0;id<2;id++){
	icent=cent_id[ican];
	c2[ig][ican]->cd(ip+id+1);
	//gPad->SetLogz();
	hGrid_Energy[iside][icent][ig][id]->Draw("colz");
	//hGrid_Energy[iside][icent][ig][id]->GetXaxis()->SetRangeUser(0.,2000.);
	//hGrid_Energy[iside][icent][ig][ig][id]->GetYaxis()->SetRangeUser(0.,2000.);
	hGrid_Energy[iside][icent][ig][id]->GetYaxis()->SetTitleOffset(1.2);
	hGrid_Energy[iside][icent][ig][id]->SetXTitle("Grid Point(X-axis)");
	hGrid_Energy[iside][icent][ig][id]->SetYTitle("Grid Point(Y-axis)");
	tt.DrawLatex(0.2,0.84,side_text[iside].c_str());
	tt.DrawLatex(0.4,0.84,cent_text[ican].c_str());
	tt.DrawLatex(0.6,0.84,data_text[id].c_str());
      }
      iside=1,ip=2;
      for(int id=0;id<2;id++){
        icent=cent_id[ican];
        c2[ig][ican]->cd(ip+id+1);
        //gPad->SetLogz();
        hGrid_Energy[iside][icent][ig][id]->Draw("colz");
        //hGrid_Energy[iside][icent][ig][id]->GetXaxis()->SetRangeUser(0.,2000.);
        //hGrid_Energy[iside][icent][ig][ig][id]->GetYaxis()->SetRangeUser(0.,2000.);
        hGrid_Energy[iside][icent][ig][id]->GetYaxis()->SetTitleOffset(1.2);
      	hGrid_Energy[iside][icent][ig][id]->SetXTitle("Grid Point(X-axis)");
        hGrid_Energy[iside][icent][ig][id]->SetYTitle("Grid Point(Y-axis)");
        tt.DrawLatex(0.2,0.84,side_text[iside].c_str());
        tt.DrawLatex(0.4,0.84,cent_text[ican].c_str());
	tt.DrawLatex(0.6,0.84,data_text[id].c_str());
      }
    }
  }
  //Save canvas
  c2[0][0]->SaveAs("./plots/Grid_Energy_cent_20to30_LG.pdf");
  c2[0][1]->SaveAs("./plots/Grid_Energy_cent_40to50_LG.pdf");
  c2[1][0]->SaveAs("./plots/Grid_Energy_cent_20to30_HG.pdf");
  c2[1][1]->SaveAs("./plots/Grid_Energy_cent_40to50_HG.pdf");
}


void plot3(){
  //Plot
  gStyle->SetOptStat(0);
  TLine *line0=new TLine(0.,0.,5.,0.);
  line0->SetLineStyle(7);
  line0->SetLineWidth(1);

  const int NCan=2;
  TCanvas *c3[NGAIN][NCan];
  TLegend *leg3[NGAIN][NCan];
  for(int ig=0;ig<NGAIN;ig++){
    for(int ican=0;ican<NCan;ican++){
      sprintf(name,"c3_gain%d_%d",ig,ican);
      c3[ig][ican]=new TCanvas(name,name,1200,800);
      c3[ig][ican]->Divide(2,2);
      leg3[ig][ican]=new TLegend(0.5,0.55,0.88,0.88);
      leg3[ig][ican]->SetLineColor(0);
    }
  }
  int ican,id,iside,icent;
  int cent_id[]={0,3,5,7};
  string side_text[2]={"Side C","Side A"};
  string cent_text[]={"0-5%","20-30%","40-50%","60-70%"};

  for(int ig=0;ig<NGAIN;ig++){ 
    ican=0;iside=0;
    for(int i=0;i<4;i++){
      icent=cent_id[i];
      c3[ig][ican]->cd(i+1);
      gPad->SetLogy();
      hPixel_Energy[iside][icent][ig][0]->Draw();
      hPixel_Energy[iside][icent][ig][0]->GetXaxis()->SetRangeUser(0.,1500.);
      if(ig==1) hPixel_Energy[iside][icent][ig][0]->GetXaxis()->SetRangeUser(0.,20000.);
      //hPixel_Energy[iside][icent][ig][0]->GetYaxis()->SetRangeUser(0.,2000.);
      hPixel_Energy[iside][icent][ig][0]->GetYaxis()->SetTitleOffset(1.2);
      hPixel_Energy[iside][icent][ig][0]->SetXTitle("Sum of Pixel Energies");
      hPixel_Energy[iside][icent][ig][1]->Draw("same");
      hPixel_Energy[iside][icent][ig][2]->Draw("same");
      hPixel_Energy[iside][icent][ig][0]->SetLineColor(4);
      hPixel_Energy[iside][icent][ig][1]->SetLineColor(2);
      hPixel_Energy[iside][icent][ig][2]->SetLineColor(6);
      tt.DrawLatex(0.4,0.82,cent_text[i].c_str());
      if(i==0){
	tt.DrawLatex(0.2,0.82,side_text[iside].c_str());
	leg3[ig][0]->AddEntry(hPixel_Energy[iside][icent][ig][0],"2010 Pb+Pb","L");
	leg3[ig][0]->AddEntry(hPixel_Energy[iside][icent][ig][2],"2011 Pb+Pb","L");
	leg3[ig][0]->AddEntry(hPixel_Energy[iside][icent][ig][1],"2015 Pb+Pb","L");
	leg3[ig][0]->Draw("same");
      }
    }

    ican=1,iside=1;
    for(int i=0;i<4;i++){
      icent=cent_id[i];
      c3[ig][ican]->cd(i+1);
      gPad->SetLogy();
      hPixel_Energy[iside][icent][ig][0]->Draw();
      hPixel_Energy[iside][icent][ig][0]->GetXaxis()->SetRangeUser(0.,2500.);
      if(ig==1) hPixel_Energy[iside][icent][ig][0]->GetXaxis()->SetRangeUser(0.,20000.);
      //hPixel_Energy[iside][icent][ig][0]->GetYaxis()->SetRangeUser(0.,2000.);
      hPixel_Energy[iside][icent][ig][0]->GetYaxis()->SetTitleOffset(1.2);
      hPixel_Energy[iside][icent][ig][0]->SetXTitle("Sum of Pixel Energies");
      hPixel_Energy[iside][icent][ig][1]->Draw("same");
      hPixel_Energy[iside][icent][ig][2]->Draw("same");
      hPixel_Energy[iside][icent][ig][0]->SetLineColor(4);
      hPixel_Energy[iside][icent][ig][1]->SetLineColor(2);
      hPixel_Energy[iside][icent][ig][2]->SetLineColor(6);
      tt.DrawLatex(0.4,0.82,cent_text[i].c_str());
      if(i==0){
	tt.DrawLatex(0.2,0.82,side_text[iside].c_str());
	leg3[ig][ican]->AddEntry(hPixel_Energy[iside][icent][ig][0],"2010 Pb+Pb","L");
	leg3[ig][ican]->AddEntry(hPixel_Energy[iside][icent][ig][2],"2011 Pb+Pb","L");
	leg3[ig][ican]->AddEntry(hPixel_Energy[iside][icent][ig][1],"2015 Pb+Pb","L");
	leg3[ig][ican]->Draw("same");
      }
    }
  }
  //Save canvas
  c3[0][0]->SaveAs("./plots/Pixel_Energy_sideC_LG.pdf");
  c3[0][1]->SaveAs("./plots/Pixel_Energy_sideA_LG.pdf");
  c3[1][0]->SaveAs("./plots/Pixel_Energy_sideC_HG.pdf");
  c3[1][1]->SaveAs("./plots/Pixel_Energy_sideA_HG.pdf");

}

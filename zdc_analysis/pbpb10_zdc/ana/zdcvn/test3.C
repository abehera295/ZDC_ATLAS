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

//Variables and Hists
char name[300];
static double PI = acos(-1.0);
const int NCENT=10;
const int NHAR=4;
const int NCH=2;
const int NPT=10;
const int NSIDE=2;
TProfile *hMean_pt[NCENT][NCH+1][NSIDE+1];

void test3(){
  TFile *fin=new TFile("./Output_cum.root","read");
  for (int icent=0; icent<NCENT; icent++) {
    for (int ich=0; ich<NCH+1; ich++) {
      for (int iside=0; iside<NSIDE+1; iside++) {
	sprintf(name, "hMean_pt_cent%d_ch%d_side%d",icent,ich,iside);
	hMean_pt[icent][ich][iside] = (TProfile*)fin->Get(name);
      }
    }
  }

  int icent=0;
  int iside=NSIDE,ich=NCH;
  cout<<"Enter icent :"<<endl;
  cin>>icent;
  int N=hMean_pt[icent][ich][iside]->GetNbinsX();
  for(int i=0;i<N;i++){
    double y,Ny;
    y=hMean_pt[icent][ich][iside]->GetBinContent(i+1);
    Ny=hMean_pt[icent][ich][iside]->GetBinEntries(i+1);
    cout<<i<<"\t\t"<<Ny<<"\t\t"<<y<<endl;
  }
  /*
  TCanvas *c;
  c=new TCanvas("c","c",1000,400);
  c->Divide(2);
  c->cd(1);
  hVn_Zdc_pt_fg[icent][iside][ihar]->Draw();
  c->cd(2);
  hVnc_Zdc_pt_fg[icent][iside][ihar]->Draw();
  //hVn_Zdc_pt_fg[icent][iside][ihar]->GetYaxis()->SetRangeUser(0.,0.5);
  */
}

#define extractor_flat_cxx
#include "./extractor_flat.h"
#include <TH2.h>
#include <TF1.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <stdio.h>
#include <vector>
#include <sstream>
#include <iterator>
#include <algorithm>
#include "TStopwatch.h"

#include "map_zdcPix.C"
//#include "extractor_zdc_psi.C"
//#include "extractor_wei.C"
#include "extractor_res.C"
#include "extractor_vn.C"
//#include "extractor_vn2.C"
//#include "extractor_SP.C"
//#include "extractor_cum.C"
//#include "extractor_3pv.C"
#include "extractor_doAna.C"
//#include "extractor_doAna2.C"

void extractor_flat::set_etaBins() {
  
  float diffEtaArr[] = {-2.5, -2.0, -1.5, -1.0, -0.5, -0.25, 0.0, 0.25, 0.5, 1.0, 1.5, 2.0, 2.5};
  diffEtaBins.assign(diffEtaArr, diffEtaArr + sizeof(diffEtaArr)/sizeof(float));
  //NETA = diffEtaBins.size()-1;
  std::cout << "NDIFFETA  = " << NETA << ": diffEtaBins:: {";
  std::copy(diffEtaBins.begin(), diffEtaBins.end(), std::ostream_iterator<float>(std::cout, ", "));  std::cout << '\b' << '\b' << "}" << std::endl;
  
}


void extractor_flat::set_ptBins() {
  
  float diffPtArr[] = {0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 7.0, 9.0, 12.0};
  diffPtBins.assign(diffPtArr, diffPtArr + sizeof(diffPtArr)/sizeof(float));
  //NPT = diffPtBins.size()-1;
  std::cout << "NPT       = " << NPT << ": diffPtBins :: {";
  std::copy(diffPtBins.begin(), diffPtBins.end(), std::ostream_iterator<float>(std::cout, ", "));  std::cout << '\b' << '\b' << "}" << std::endl;
  
}



extractor_flat::extractor_flat(string filelist,int rno, int iev, int nevts) {
  RNum=rno;
  NEv_id=iev;
  NEv=nevts;
  sprintf(name,"./out/out_%d_%d.root",RNum,NEv_id);
  outname=name;

  if(RNum!=0){
    fChain= new TChain("HeavyIonD3PD","");
    char fname[400];
    ifstream lis(filelist.c_str());
    int cnt=0;
    int cnt2=0;
    cout<<filelist<<endl;
    while(!lis.eof()){
      string filename;
      lis >> filename;
      sprintf(fname,"%s",filename.c_str());
      //cout << fname << endl;
      if(!filename.empty()) {
	fChain->Add(fname);
	cout<<fname<<endl;
	cnt++;
      }
      //cnt++;
      //if(cnt>0) break;
      if(cnt%1==0) cout<<"Added root files : "<<cnt<<endl;
      if(cnt>1000000) {cout<<"Too Many Files"<<endl;break;}
    }
    //Init(chain);
    //Initialise branches of tree
    fChain->SetBranchAddress("RunNumber",   &RunNumber);//run number
    fChain->SetBranchAddress("EventNumber", &EventNumber);//run number
    fChain->SetBranchAddress("lbn", &lbn);//run number
    fChain->SetBranchAddress("vx_n", &vx_n);
    fChain->SetBranchAddress("vx_x", &vx_x);
    fChain->SetBranchAddress("vx_y", &vx_y);
    fChain->SetBranchAddress("vx_z", &vx_z);
    fChain->SetBranchAddress("Fcal_Et", &Fcal_Et);
    fChain->SetBranchAddress("Fcal_Et_p", &Fcal_Et_p);
    fChain->SetBranchAddress("Fcal_Et_n", &Fcal_Et_n);
    fChain->SetBranchAddress("Centrality", &Centrality);
    fChain->SetBranchAddress("nloose", &nloose);
    fChain->SetBranchAddress("ntight", &ntight);
    fChain->SetBranchAddress("trk_n", &trk_n);
    fChain->SetBranchAddress("trk_data", trk_data);
    fChain->SetBranchAddress("Qw", &Qw);
    fChain->SetBranchAddress("Qx", &Qx);
    fChain->SetBranchAddress("Qy", &Qy);
    //ZDC branches
    fChain->SetBranchAddress("Qw_zdc", Qw_zdc);
    fChain->SetBranchAddress("Qx_zdc", Qx_zdc);
    fChain->SetBranchAddress("Qy_zdc", Qy_zdc);
    fChain->SetBranchAddress("Zdc_n", &Zdc_n);
    //fChain->SetBranchAddress("Zdc_Vec_Size", &Zdc_Vec_Size);
    fChain->SetBranchAddress("Zdc_Energy_LG", &Zdc_Energy_LG);
    //fChain->SetBranchAddress("Zdc_Time_LG", &Zdc_Time_LG);
    fChain->SetBranchAddress("Zdc_Energy_HG", &Zdc_Energy_HG);
    //fChain->SetBranchAddress("Zdc_Time_HG", &Zdc_Time_HG);
    fChain->SetBranchAddress("Zdc_Id", &Zdc_Id);
    fChain->SetBranchAddress("Zdc_Side", &Zdc_Side);
    fChain->SetBranchAddress("Zdc_Type", &Zdc_Type);
    fChain->SetBranchAddress("Zdc_Module", &Zdc_Module);
    fChain->SetBranchAddress("Zdc_Channel", &Zdc_Channel);
    //fChain->SetBranchAddress("Zdc_TimeCalib", &Zdc_TimeCalib);
  }
  set_ptBins ();
  set_etaBins();

}

/*
  void extractor_flat::Exec_zdc_psi(int entSel) {
  TStopwatch *tt=new TStopwatch();
  tt->Start();
  if (fChain == 0) return;
  nevents = 0;

  //fOut_zdc_psi = new TFile("Output_zdc_psi.root", "recreate");
  fOut_zdc_psi = new TFile(outname.c_str(),"recreate");

  Init_zdc_psi();
  //FillCalib();
  //FillFlattening();
  Run_zdc_psi();
  //DoAna_vn();
  SaveHistos_zdc_psi();
  
  tt->Stop();
  tt->Print();

  }
*/
Event_EP::Event_EP(int evt){
  memset(Psi_FCal,0,sizeof(Psi_FCal));
  memset(Psi_Zdc,0,sizeof(Psi_Zdc));

  id=evt;
  icent=0;
  ivz=0;
}

Event::Event(int evt){
  memset(Qx_ID,0,sizeof(Qx_ID));
  memset(Qy_ID,0,sizeof(Qy_ID));
  memset(Qw_ID,0,sizeof(Qw_ID));

  memset(Qx_ID_cum,0,sizeof(Qx_ID_cum));
  memset(Qy_ID_cum,0,sizeof(Qy_ID_cum));
  memset(Qw_ID_cum,0,sizeof(Qw_ID_cum));

  memset(Qx_FCal,0,sizeof(Qx_FCal));
  memset(Qy_FCal,0,sizeof(Qy_FCal));
  memset(Qw_FCal,0,sizeof(Qw_FCal));
  memset(Qx_Zdc,0,sizeof(Qx_Zdc));
  memset(Qy_Zdc,0,sizeof(Qy_Zdc));
  memset(Qw_Zdc,0,sizeof(Qw_Zdc));

  id=evt;
  icent=0;
  ivz=0;
}

void extractor_flat::Exec_res(int entSel) {
  TStopwatch *tt=new TStopwatch();
  tt->Start();
  
  if (fChain == 0) return;
  nevents = 0;
  
  //fOut_res = new TFile("Output_res.root", "recreate");
  fOut_res = new TFile(outname.c_str(),"recreate");
  effTool = new Tool();
  ReadHistos_res();
  Init_res();
  Run_res();
  SaveHistos_res();
  
  tt->Stop();
  tt->Print();
}

void extractor_flat::Exec_vn(int entSel) {
  TStopwatch *tt=new TStopwatch();
  tt->Start();
  
  if (fChain == 0) return;
  nevents = 0;
  
  //fOut_vn = new TFile("Output_vn.root", "recreate");
  fOut_vn = new TFile(outname.c_str(),"recreate");
  effTool = new Tool();
  ReadHistos_vn();
  Init_vn();
  Run_vn();
  SaveHistos_vn();
  
  tt->Stop();
  tt->Print();
}

void extractor_flat::Exec_doAna(int entSel) {
  TStopwatch *tt=new TStopwatch();
  tt->Start();

  //fOut_doAna = new TFile("Output_doAna.root", "recreate");
  ReadHistos_doAna();
  Init_doAna();
  Run_doAna();
  SaveHistos_doAna();
  
  tt->Stop();
  tt->Print();
}


/*
void extractor_flat::Exec_vn2(int entSel) {
  TStopwatch *tt=new TStopwatch();
  tt->Start();
  
  if (fChain == 0) return;
  nevents = 0;
  
  //fOut_vn2 = new TFile("Output_vn2.root", "recreate");
  fOut_vn2 = new TFile(outname.c_str(),"recreate");
  effTool = new Tool();
  ReadHistos_vn2();
  Init_vn2();
  Run_vn2();
  SaveHistos_vn2();
  
  tt->Stop();
  tt->Print();
}


void extractor_flat::Exec_SP(int entSel) {
  TStopwatch *tt=new TStopwatch();
  tt->Start();
  
  if (fChain == 0) return;
  nevents = 0;
  
  //fOut_SP = new TFile("Output_SP.root", "recreate");
  fOut_SP = new TFile(outname.c_str(),"recreate");
  effTool = new Tool();
  ReadHistos_SP();
  Init_SP();
  Run_SP();
  SaveHistos_SP();
  
  tt->Stop();
  tt->Print();
}

void extractor_flat::Exec_cum(int entSel) {
  TStopwatch *tt=new TStopwatch();
  tt->Start();
  
  if (fChain == 0) return;
  nevents = 0;
  
  //fOut_cum = new TFile("Output_cum.root", "recreate");
  fOut_cum = new TFile(outname.c_str(),"recreate");
  effTool = new Tool();
  ReadHistos_cum();
  Init_cum();
  Run_cum();
  SaveHistos_cum();
  
  tt->Stop();
  tt->Print();
}
*/


/*
  void extractor_flat::Exec_wei(int entSel) {
  TStopwatch *tt=new TStopwatch();
  tt->Start();
  
  if (fChain == 0) return;
  nevents = 0;
  
  fOut_wei = new TFile("Output_wei.root", "recreate");
  effTool=new Tool();
  Init_wei();
  Run_wei();
  SaveHistos_wei();
  
  tt->Stop();
  tt->Print();
  }
*/
void extractor_flat::get_qVec_zdc(std::vector<std::vector<double> >& Qx, std::vector<std::vector<double> >& Qy, std::vector<double>& Qw) {

  for (int iside=0; iside<NSIDE+1; iside++) {
    std::vector<double> Qx_side; 
    std::vector<double> Qy_side;
    for (int ihar=0; ihar<NHAR; ihar++) {
      Qx_side.push_back(Qx_zdc[iside][ihar]);
      Qy_side.push_back(Qy_zdc[iside][ihar]);
    }
    Qx.push_back(Qx_side);
    Qy.push_back(Qy_side);
    Qw.push_back(Qw_zdc[iside]);
  }
}

//Old version when Q-vectors were not saved in TTree
void extractor_flat::get_qVec(std::vector<std::vector<double> >& Qx, std::vector<std::vector<double> >& Qy, std::vector<double>& Qw) {

  for (int iside=0; iside<NSIDE+1; iside++) {
    std::vector<double> Qx_side(NHAR, 0); 
    std::vector<double> Qy_side(NHAR, 0);
    Qx.push_back(Qx_side);
    Qy.push_back(Qy_side);
    Qw.push_back(0);
  }

  for (int ichn=0; ichn<Zdc_n; ichn++) {
    if (Zdc_Type->at(ichn) == 0 || Zdc_Module->at(ichn) != 1) continue; 

    int side = Zdc_Side->at(ichn)+1;
    if (side==2) side = 1;
    if (!(side ==0 || side == 1)) {std::cout << "Unknown side: chId = " << ichn << ": side = " << side << std::endl; continue;}
    int side_sign = 1;
    if (side == 1) side_sign = -1;

    int chId = ichn;
    if (side == 0 ) chId = ichn - 2;
    if (side == 1 ) chId = ichn - 30;

    double xcord, ycord, r, phi;
    set_xy_ZdcPix_had(chId, xcord, ycord);
    if (side == 1) xcord = -xcord;
    set_rphi_ZdcPix(xcord, ycord, r, phi);

    if (fabs(xcord-1000.0)<1e-6 && fabs(ycord-1000.0)<1e-6) {
      std::cout << "Unknown channel: side = " << side << ": channel = " << ichn << ": type = " << Zdc_Type->at(ichn) << std::endl;
      continue;
    }

    double eAdc[2];
    eAdc[0] = Zdc_Energy_LG->at(ichn);

    for (int ihar=0; ihar<NHAR; ihar++) {
      Qx.at(side).at(ihar)     += eAdc[0]*cos((ihar+1)*phi);
      Qy.at(side).at(ihar)     += eAdc[0]*sin((ihar+1)*phi);
      if (ihar==0) Qw.at(side) += eAdc[0];

      Qx.at(NSIDE).at(ihar)    += side_sign*eAdc[0]*cos((ihar+1)*phi);
      Qy.at(NSIDE).at(ihar)    += side_sign*eAdc[0]*sin((ihar+1)*phi);
      if (ihar==0) Qw.at(NSIDE)+= eAdc[0];
      //if (ihar==0) std::cout << "ch = " << chId << ": x = " << xcord << ": y = " << ycord << ": phi = " << phi*180/3.1416 << ":     eADC =" << eAdc[0] << ": cos = " << eAdc[0]*cos((ihar+1)*phi) << ": sin = " << eAdc[0]*sin((ihar+1)*phi) << ": Qx = " << Qx.at(side).at(ihar) << ": Qy = " << Qy.at(side).at(ihar) << std::endl; 
    }
  }
}



std::vector<double> extractor_flat::discFourrier(const TH1D* h1, int ihar, TFile* f1) {

  std::vector<double> vn;
  if(h1->GetEntries() == 0) {
    //std::cout << "discFourrier ERROR : Empty histogram " << h1->GetName() << std::endl;
    vn.push_back(0.);
    vn.push_back(0.);
    return vn;
  }

  double p = 0, ep = 0;
  double bwid    = h1->GetBinWidth(1)/2;
  int    Nbins   = h1->GetNbinsX();

  double Sum=0,Sum_err=0;//Find vn
  for(int phibin=1; phibin<=Nbins;phibin++){
    float binval = h1->GetBinContent(phibin);
    float phi    = h1->GetBinCenter(phibin);
    Sum         += binval;
    Sum_err     += pow(h1->GetBinError(phibin),2);
    //p += cos( (ihar)*phi) * binval;
    p += cos( (ihar+1)*phi) * binval;    
  }
  double SumAvg = Sum/Nbins;
  //double ErrAvg = sqrt(Sum_err)/Nbins;
  p /=Sum;

  for(int phibin=1; phibin<=Nbins;phibin++){//Find Errors
    float binerr = h1->GetBinError(phibin);
    float phi    = h1->GetBinCenter(phibin);
    //ep  += pow((cos((ihar)*phi)-p)*binerr,2);
    ep  += pow((cos((ihar+1)*phi)-p)*binerr,2);
  }

  //Correct for Bin Width
  double corr = (ihar+1)*bwid/sin((ihar+1)*bwid);
  corr=1.;
  p*=corr;
  ep=sqrt(ep)/Sum*corr;

  vn.push_back(p);
  vn.push_back(ep);

  if (f1 != NULL) {//Write Histogram with Harmonics
    std::string hName = h1->GetName();
    hName += "_fourhar";
    TH1D* hFour = (TH1D*)h1->Clone(hName.c_str());
    hFour->SetMarkerStyle(33);
    hFour->SetMarkerSize (1.3);
    hFour->SetMarkerColor(1);
    hFour->SetLineColor  (1);

    double hLowEdge = h1->GetBinLowEdge(1);
    double hUpEdge  = h1->GetBinLowEdge(h1->GetNbinsX()+1);

    TF1* fHar = (TF1*)NULL;
    std::stringstream      ss_har; ss_har << ihar;
    std::string fName    = "FourrierHar" + ss_har.str();
    std::string fFormula = "[0]*(1+2*[1]*cos(" + ss_har.str() + "*x))";

    fHar = new TF1(fName.c_str(), fFormula.c_str(), hLowEdge, hUpEdge);
    fHar ->SetLineWidth(1);
    fHar ->SetLineStyle(5);
    fHar ->SetLineColor(ihar%8+2);
    hFour->Fit(fHar,"QR");
    fHar ->SetParameter(1, vn.at(0));
    fHar ->SetParError (1, vn.at(1));
    hFour->GetListOfFunctions()->Add(fHar);

    f1->cd(0); if (!(f1->GetListOfKeys()->Contains(hFour->GetName()))) hFour->Write();
    hFour->GetListOfFunctions()->Remove(fHar);

    delete fHar;
  }

  return vn;
}

float extractor_flat::get_pt(int ptbin){
  float ptgev=110;
  if(ptbin<20) ptgev = 0.1*ptbin+0.05;
  else if(ptbin<30) ptgev = 0.2*(ptbin-20)+2.1;
  else if(ptbin<38) ptgev = 0.5*(ptbin-30)+4.2;
  else if(ptbin<54) ptgev = 2*(ptbin-38)+8.8;
  else if(ptbin<60) ptgev = 5*(ptbin-54)+42;
  else if(ptbin<63) ptgev = 10*(ptbin-60)+74;
  return ptgev;
}

int extractor_flat::get_centb(double et) {
  centb = -1;
  int icent;
  for(icent=0;icent<80;icent++){
    if(et>centcuts[icent]) break;
  }

  centb_1p = icent;
  //centb = icent/4.99999999999;
  //centb = icent/9.99999999999;
  for(int i=0;i<NCENT;i++){
    if(centb_1p>=cent_mat[i] && centb_1p<cent_mat[i+1]) {centb=i;break; } 
  }
  return centb;
}

int extractor_flat::get_vzb(double vz){
  float vz_mat[]={-100.,-80.,-60.,-40.,-20.,0.,20.,40.,60.,80.,100.};
  int vzb=-1;
  for(int i=0;i<NVZ_EP;i++){
    if(vz>=vz_mat[i] && vz<vz_mat[i+1]){
      vzb=i;
      break;
    }
  }
  return vzb;
}

int extractor_flat::get_ptBin(double ptVal){
  int bin=-1;
  if (ptVal < diffPtBins.at(0)) return bin;

  for (int ipt=0; ipt<NPT; ipt++) {
    if (ptVal < diffPtBins.at(ipt+1)) {
      bin = ipt;
      break;
    }
  }
  return bin;
}


int extractor_flat::get_etaBin(double etaVal){
  int bin=-1;
  if (etaVal < diffEtaBins.at(0)) return bin;

  for (int ieta=0; ieta<NETA; ieta++) {
    if (etaVal < diffEtaBins.at(ieta+1)) {
      bin = ieta;
      break;
    }
  }

  return bin;
}

int extractor_flat::get_ptBin_TrkWei(int ptb){
  int ptbin=-1;
  if(ptb<7) ptbin=0;//0.5-0.7
  else if(ptb<10)  ptbin=1;//0.7-1
  else if(ptb<15)  ptbin=2;//1-1.5
  else if(ptb<20)  ptbin=3;//1.5-2
  else if(ptb<25)  ptbin=4;//2-3
  else if(ptb<30)  ptbin=5;//3-4
  else if(ptb<34)  ptbin=6;//4-6
  else if(ptb<39)  ptbin=7;//6-10
  else if(ptb<44)  ptbin=8;//10-20

  return ptbin;
}


extractor_flat::~extractor_flat() {
  //if (!fChain) return;
  //delete fChain->GetCurrentFile();
}
/*

Bool_t extractor_flat::Notify() {
  return kTRUE;
}

void extractor_flat::Show(Long64_t entry) {
  if (!fChain) return;
  fChain->Show(entry);
}

Int_t extractor_flat::Cut(Long64_t entry) {
  return 1;
}

Int_t extractor_flat::GetEntry(Long64_t entry) {
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}

Long64_t extractor_flat::LoadTree(Long64_t entry) {
  // Set the environment to read one entry
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (fChain->GetTreeNumber() != fCurrent) {
    fCurrent = fChain->GetTreeNumber();
    Notify();
  }
  return centry;
}
*/

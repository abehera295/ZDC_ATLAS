#include "Tool.h"



Tool::Tool()
{
  /*
  TFile* fInFlat;
  sprintf(name,"/star/data01/pwg/jjia/Run13_AuAu/trkwei1/wei%d.root",runNum);
  fInFlat = TFile::Open(name);
  for(int iZ=0; iZ<nZvtxFlat; iZ++)
    {
      for(int iT=0; iT<nCentFlat; iT++)
	{
	  for(int iC=0; iC<nChgFlat; iC++)
	    {
	      for(int iP=0; iP<nPtFlat+3; iP++)
		{
		  sprintf(name,"hetaphi_z%d_c%d_ch%d_pt%d",iZ,iT,iC,iP);
		  hEtaPhi[iZ][iT][iC][iP] = (TH2D*)fInFlat->Get(name);
		}
	      for(int iP=nPtFlat; iP<nPtFlat+3; iP++)
		{
		  hEtaPhi[iZ][iT][iC][nPtFlat-1]->Add(hEtaPhi[iZ][iT][iC][iP]); // merge the higher pT bins
		}
	      for(int iP=0; iP<nPtFlat; iP++)
		{
		  TH1* proj = (TH1*)hEtaPhi[iZ][iT][iC][iP]->ProjectionY("proj"); proj->Scale(1./nPhiFlat);
		  for(int iy=1; iy<=nEtaFlat; iy++)
		    {
		      double nor = proj->GetBinContent(iy);
		      for(int ix=1; ix<=nPhiFlat; ix++)
			{
			  hEtaPhi[iZ][iT][iC][iP]->SetBinContent(ix,iy,hEtaPhi[iZ][iT][iC][iP]->GetBinContent(ix,iy)/nor);
			  hEtaPhi[iZ][iT][iC][iP]->SetBinError(ix,iy,hEtaPhi[iZ][iT][iC][iP]->GetBinError(ix,iy)/nor);
			}
		    }
		}
	    }
	}
    }
  */
  TFile* fInEff = new TFile("/phenix/plhf/mzhou/AnaCumu_PbPb502/MAIN_binCent/INPUT/trkEff.root","READ");
  allTrue = (TH3D*)fInEff->Get("allTrue");
  allReco = (TH3D*)fInEff->Get("allReco");
  matchedReco = (TH3D*)fInEff->Get("matchedReco");
  unMatchedReco = (TH3D*)fInEff->Get("unMatchedReco");
  matchedReco->Divide(allTrue);
  unMatchedReco->Divide(allReco);

  TFile* fInTrig = new TFile("/phenix/plhf/mzhou/AnaCumu_PbPb502/MAIN_binCent/INPUT/trigWeight.root","READ");
  hTrigAll = (TH1D*)fInTrig->Get("disFCal_Trig0");
  hTrigAll->Rebin(2);
  hTrigMB = (TH1D*)fInTrig->Get("disFCal_Trig1");
  hTrigMB->Rebin(2);
  hTrigAll->Divide(hTrigMB);

  TFile* fInCvt = new TFile("/phenix/plhf/mzhou/AnaCumu_PbPb502/MAIN_binCent/INPUT/hist_cvt.root","READ");
  gCvt_FCal_Cent = (TGraphErrors*)fInCvt->Get("gCvt_FCal_Cent");
}

Tool::~Tool()
{
}

double Tool::detTrigWght(double fcalEt)
{
  if(fcalEt>4.55) return 1./21.8; // no statistics there
  if(fcalEt<4.0 ) return 1.0; // no trigger there

  int iBin = int((fcalEt+1)/hTrigAll->GetBinWidth(1));
  if(hTrigAll->GetBinContent(iBin+1)<=0) return 1;
  else return 1./(hTrigAll->GetBinContent(iBin+1));
}

double Tool::detTrkEff(double fcalEt, double pt, double eta)
{
  if(pt<=0.5) return 0.4; // FIX
  if(pt>=10) return 0.75; // FIX

  int tagFcalEt = 0;
  int tagPt = 0;
  int tagEta = 0;
  for(int i=matchedReco->GetNbinsX(); i>0; i--)
    {
      if(fcalEt>matchedReco->GetXaxis()->GetBinLowEdge(i))
	{
	  tagFcalEt = i-1;
	  break;
	}
    }
  for(int i=matchedReco->GetNbinsY(); i>0; i--)
    {
      if(pt>matchedReco->GetYaxis()->GetBinLowEdge(i))
	{
	  tagPt = i-1;
	  break;
	}
    }
  for(int i=matchedReco->GetNbinsZ(); i>0; i--)
    {
      if(eta>matchedReco->GetZaxis()->GetBinLowEdge(i))
	{
	  tagEta = i-1;
	  break;
	}
    }

  return matchedReco->GetBinContent(tagFcalEt+1,tagPt+1,tagEta+1);
}

double Tool::detTrkFak(double fcalEt, double pt, double eta)
{
  if(pt<=0.5) return 0; // FIX
  if(pt>=10) return 0; // FIX

  int tagFcalEt = 0;
  int tagPt = 0;
  int tagEta = 0;
  for(int i=unMatchedReco->GetNbinsX(); i>0; i--)
    {
      if(fcalEt>unMatchedReco->GetXaxis()->GetBinLowEdge(i))
	{
	  tagFcalEt = i-1;
	  break;
	}
    }
  for(int i=unMatchedReco->GetNbinsY(); i>0; i--)
    {
      if(pt>unMatchedReco->GetYaxis()->GetBinLowEdge(i))
	{
	  tagPt = i-1;
	  break;
	}
    }
  for(int i=unMatchedReco->GetNbinsZ(); i>0; i--)
    {
      if(eta>unMatchedReco->GetZaxis()->GetBinLowEdge(i))
	{
	  tagEta = i-1;
	  break;
	}
    }

  return unMatchedReco->GetBinContent(tagFcalEt+1,tagPt+1,tagEta+1);
}

double Tool::detFlat(double zVtx, int cent, int tagChg, double pt, double eta, double phi)
{
  /*
  double weight = 1;
  int tagZvtx = -1;
  int tagCent = -1;
  int tagPt   = -1;
  int tagEta  = -1;
  int tagPhi  = -1;

  tagZvtx = int((zVtx+100.)/25);
  if(tagZvtx<0 || tagZvtx>=nZvtxFlat) return 1;

  if(cent<1) tagCent=0;
  else if(cent<3)  tagCent=1;
  else if(cent<5)  tagCent=2;
  else if(cent<10) tagCent=3;
  else if(cent<15) tagCent=4;
  else if(cent<20) tagCent=5;
  else if(cent<25) tagCent=6;
  else if(cent<30) tagCent=7;
  else if(cent<35) tagCent=8;
  else if(cent<40) tagCent=9;
  else if(cent<45) tagCent=10;
  else if(cent<50) tagCent=11;
  else if(cent<60) tagCent=12;
  else tagCent=13;

  if(pt>=0.5 && pt<0.7) tagPt = 0;
  else if(pt<1.0) tagPt = 1;
  else if(pt<1.5) tagPt = 2;
  else if(pt<2.0) tagPt = 3;
  else if(pt<3.0) tagPt = 4;
  else tagPt = 5;
  if(tagPt<0) return 1;

  tagEta = int((eta+2.5)/5.*nEtaFlat);
  if(tagEta<0 || tagEta>=nEtaFlat) return 1;
  tagPhi = int(phi/2./TMath::Pi()*nPhiFlat);
  if(tagPhi<0 || tagPhi>=nPhiFlat) return 1;

  weight = hEtaPhi[tagZvtx][tagCent][tagChg][tagPt]->GetBinContent(tagPhi+1,tagEta+1);
  if(weight<0.1) weight=0.1;
  if(weight>2.0) weight=2.0;
  return 1./weight;
  */
}

double Tool::detCent(double fcalEt)
{
  return gCvt_FCal_Cent->Eval(fcalEt);
}

void extractor_flat::Init_wei(){
  set_ptBins ();
  set_etaBins();
  hNTrk=new TH1D("hNTrk","hNTrk",5000,0.,5000.); hNTrk->Sumw2();
  //Track histograms
  hPt = new TH1D("hPt","",1000,0.,25.); hPt->Sumw2();
  hPhi = new TH1D("hPhi","",128,0.,PI2); hPhi->Sumw2();
  hPhi_wei = new TH1D("hPhi_wei","",128,0.,PI2); hPhi_wei->Sumw2();
  hEta = new TH1D("hEta","",100,-2.5,2.5); hEta->Sumw2();
  hCharge = new TH1D("hCharge","",3,-0.5,2.5); hCharge->Sumw2();
  hTrk_w = new TH1D("hTrk_w","",1000,-300.,300.); hTrk_w->Sumw2();
  hTrk_wei = new TH1D("hTrk_wei","",1000,-5.,5.); hTrk_wei->Sumw2();
  hTrk_eff=new TH1D("hTrk_eff","",1000,-5.,5.); hTrk_eff->Sumw2();
  hWeight = new TH1D("hWeight","",100,0.,10.); hWeight->Sumw2();
  //Correlation hists
  for(int icent=0;icent<NCENT;icent++){
    for(int ivz=0;ivz<NVZ_EP;ivz++){
      for(int ipt=0;ipt<NPT_EP;ipt++){
	for(int ich=0;ich<NCH_EP;ich++){
	  sprintf(name,"hEtaPhi_icent%d_ivz%d_ipt%d_ich%d",icent,ivz,ipt,ich);
	  hEtaPhi[icent][ivz][ipt][ich]=new TH2F(name,"",100,-2.5,2.5,128,0.,PI2);
	  sprintf(name,"hEtaPhi2_icent%d_ivz%d_ipt%d_ich%d",icent,ivz,ipt,ich);
	  hEtaPhi2[icent][ivz][ipt][ich]=new TH2F(name,"",100,-2.5,2.5,128,0.,PI2);
	  sprintf(name,"hEtaPhi_wei_icent%d_ivz%d_ipt%d_ich%d",icent,ivz,ipt,ich);
	  hEtaPhi_wei[icent][ivz][ipt][ich]=new TH2F(name,"",100,-2.5,2.5,128,0.,PI2);
	}
      }
    }
  }
  
}



void extractor_flat::Run_wei(){
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;

  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);

    if (nb <=0)        {cout<<"Nbytes became <=0"<<endl;break;}
    if(jentry%20000==0) {cout<<"running event "<<jentry<<endl ;}
    //if (jentry > 200000) break;
    //-------------------------------------------------------------------------------------------------
    centb = get_centb(Fcal_Et);
    vzb = get_vzb(vx_z);
    if (centb < 0) continue;
    if (vzb < 0) continue;
    
    nevents++;

    ntrkQ = 0;

    for(int j=0; j<trk_n; j++){
      //if(trk_Quality2[j] < 3) continue;
      //new data format
      UInt_t trkb = trk_data[j];
      double btmp=0;
      trk_phi0_wrt_PV[j]  = trkb&PHIM;// btmp  = trk_phi0_wrt_PV[j];
      trk_phi0_wrt_PV[j]+=0.5;   trk_phi0_wrt_PV[j] *= PIslice;     trkb >>= PHIS;
      //trk_phibin[j] = btmp/(512)*NPHI; trk_phibin[j]++;

      trk_eta[j]  = trkb&ETAM; btmp  = trk_eta[j];
      trk_eta[j]+=0.5;   trk_eta[j] *= ETAslice;   trk_eta[j] -=ETAoff;  trkb >>= ETAS;
      //trk_etabin[j] = btmp/(500.0)*NETA; trk_etabin[j]++;

      trk_charge[j]= trkb&CHM;       trkb >>= CHS;
      trk_pt[j]=get_pt(trkb&PTM);
      //trk_ptbin[j] = getptbin(trkb&PTM);
      trkb >>= PTS;
      //trk_ptbin[j] = get_ptbin2(trk_pt[j]);
      //trk_Quality2[j]  = trkb&QUALM;
      trk_w[j]  = trkb&QUALM;
      trk_wei[j] =  trk_wei[j] = 1./(trk_w[j]/256.0 * 1.4+0.3);//0.3-1.7;
      
      //double ptgev = trk_pt[j]/1000.0;
      double ptgev=trk_pt[j];
      int ptBin  = get_ptBin(ptgev);       if (ptBin  < 0) continue;
      int etaBin = get_etaBin(trk_eta[j]); if (etaBin < 0) continue;
      int chBin=0;
      if(trk_charge[j]==0) chBin=0;
      if(trk_charge[j]==1) chBin=1;
      int sign_eta = 1;
      if (trk_eta[j] < 0) sign_eta = -1;
      trk_wei[j]=1.;
      trk_eff[j]=1.;
      trk_eff[j]=1./effTool->detTrkEff(Centrality,trk_eta[j],trk_pt[j]);//Trk Efficiency from Tool.cxx
      trk_phi[j]=trk_phi0_wrt_PV[j];
      hEta->Fill(trk_eta[j]);
      hPhi->Fill(trk_phi0_wrt_PV[j]);
      hPhi_wei->Fill(trk_phi0_wrt_PV[j],trk_eff[j]);
      hPt->Fill(ptgev);
      hCharge->Fill(trk_charge[j]);
      hTrk_w->Fill(trk_w[j]);
      hTrk_wei->Fill(trk_wei[j]);
      hTrk_eff->Fill(trk_eff[j]);
      hEtaPhi[centb][vzb][ptBin][chBin]->Fill(trk_eta[j],trk_phi[j]);
      ntrkQ++;
    }
    hNTrk->Fill(ntrkQ);
  }

  for(int icent=0;icent<NCENT;icent++){
    for(int ivz=0;ivz<NVZ_EP;ivz++){
      for(int ipt=0;ipt<NPT_EP;ipt++){
	for(int ich=0;ich<NCH_EP;ich++){
	  for(int ieta=0;ieta<hEtaPhi[icent][ivz][ipt][ich]->GetNbinsX();ieta++){
	    double avg=0.,avg_e=0.;
	    int cnt=0;
	    for(int iphi=0;iphi<hEtaPhi[icent][ivz][ipt][ich]->GetNbinsY();iphi++){
	      avg+=hEtaPhi[icent][ivz][ipt][ich]->GetBinContent(ieta+1,iphi+1);
	      avg_e+=pow(hEtaPhi[icent][ivz][ipt][ich]->GetBinError(ieta+1,iphi+1),2);
	      cnt++;
	    }
	    avg/=cnt;
	    avg_e=sqrt(avg_e/cnt);
	    for(int iphi=0;iphi<hEtaPhi[icent][ivz][ipt][ich]->GetNbinsY();iphi++){
	      double etaphi=hEtaPhi[icent][ivz][ipt][ich]->GetBinContent(ieta+1,iphi+1);
	      double etaphi_e=hEtaPhi[icent][ivz][ipt][ich]->GetBinError(ieta+1,iphi+1);
	      double wei=0.,wei_e=0.;
	      if(etaphi!=0.) wei=avg/etaphi;
	      if(wei<0.3) wei=0.3;
	      if(wei>1.7) wei=1.7;
	      if(avg!=0. && etaphi!=0.) wei_e=wei*sqrt(pow(avg_e/avg,2)+pow(etaphi_e/etaphi,2));;
	      hWeight->Fill(wei);
	      hEtaPhi_wei[icent][ivz][ipt][ich]->SetBinContent(ieta+1,iphi+1,wei);
	      hEtaPhi_wei[icent][ivz][ipt][ich]->SetBinError(ieta+1,iphi+1,fabs(wei_e));
	      hEtaPhi2[icent][ivz][ipt][ich]->SetBinContent(ieta+1,iphi+1,etaphi*wei);
	      hEtaPhi2[icent][ivz][ipt][ich]->SetBinError(ieta+1,iphi+1,etaphi_e*fabs(wei));
	    }
	  }
	}//ich
      }//ipt
    }//ivz
  }//icent
}


void extractor_flat::SaveHistos_wei(){
  fOut_wei->cd();
  hNTrk->Write();
  hPt->Write();
  hPhi->Write();
  hPhi_wei->Write();
  hEta->Write();
  hCharge->Write();
  hWeight->Write();
  hTrk_w->Write();
  hTrk_wei->Write();
  hTrk_eff->Write();
  for(int icent=0;icent<NCENT;icent++){
    for(int ivz=0;ivz<NVZ_EP;ivz++){
      for(int ipt=0;ipt<NPT_EP;ipt++){
	for(int ich=0;ich<NCH_EP;ich++){
	  hEtaPhi[icent][ivz][ipt][ich]->Write();
	  hEtaPhi2[icent][ivz][ipt][ich]->Write();
	  hEtaPhi_wei[icent][ivz][ipt][ich]->Write();
	}
      }
    }
  }
  fOut_wei->Close();
}

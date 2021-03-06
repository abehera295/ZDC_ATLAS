//#define withBG

void extractor_flat::ReadHistos_mu() {
  sprintf(name,"../do_calib/EPCalib_hists/EP_%d_flat.root",RNum);
  TFile *fin1=new TFile(name,"read");
  //For FCal Recenter
  sprintf(name,"hevts");
  hevts=(TH2F*)fin1->Get(name);
  for(int iz=0;iz<NZ;iz++){//FCal
    sprintf(name,"hEt_vz_z%d",iz);
    hEt_vz[iz]=(TProfile*)fin1->Get(name);
  }
  for(int id=0;id<NDET;id++){
    for(int iz=0;iz<NZ;iz++){
      sprintf(name,"hQw_FCal_det%d_z%d",id,iz);
      hQw_FCal[id][iz]=(TProfile*)fin1->Get(name);
      for(int ihar=0;ihar<NHAR;ihar++){//2,3,4,5
	sprintf(name,"hQx_FCal_det%d_z%d_h%d",id,iz,ihar);
	hQx_FCal[id][iz][ihar]=(TProfile*)fin1->Get(name);
	sprintf(name,"hQy_FCal_det%d_z%d_h%d",id,iz,ihar);
	hQy_FCal[id][iz][ihar]=(TProfile*)fin1->Get(name);
      }
    }
  }
  //For Zdc Recenter
  for (int iside=0; iside<NSIDE+1; iside++) {
    for (int ihar=0; ihar<NHAR; ihar++) {
      sprintf(name, "hQx_Zdc_side%d_har%d", iside, ihar);
      hQx_Zdc[iside][ihar] = (TProfile*)fin1->Get(name);
      sprintf(name, "hQy_Zdc_side%d_har%d", iside, ihar);
      hQy_Zdc[iside][ihar] = (TProfile*)fin1->Get(name);
    }
  }
  //For Zdc Flattening
  for (int icent=0; icent<NCENT_1p; icent++) {
    for (int iside=0; iside<NSIDE+1; iside++) {
      for (int ihar=0; ihar<NHAR; ihar++) {
	sprintf(name, "hFlatCos_Zdc_cent%d_side%d_har%d", icent, iside, ihar);
	hFlatCos_Zdc[icent][iside][ihar] = (TProfile*)fin1->Get(name);
	sprintf(name, "hFlatSin_Zdc_cent%d_side%d_har%d", icent, iside, ihar);
	hFlatSin_Zdc[icent][iside][ihar] = (TProfile*)fin1->Get(name);
      }
    }
  }
  /*
  //For Track Weights
  sprintf(name,"../do_calib/TrkWei_hists/Wei%d.root",RNum);
  TFile *fin2=new TFile(name,"read");
  TH1D* hproj;
  for(int iz=0;iz<NZ;iz++){
    for(int ic=0;ic<NCENT_TrkWei;ic++){
      for(int ich=0;ich<NCH;ich++){
        for(int ipt=0;ipt<NPT_TrkWei;ipt++){
          sprintf(name,"hetaphi_z%d_c%d_ch%d_pt%d",iz,ic,ich,ipt);
          hetaphi[iz][ic][ich][ipt]= (TH2F*)fin2->Get(name);
        }
        for(int ipt=MPT+1;ipt<NPT_TrkWei;ipt++){//highest pT bin is 3-4 GeV
          //hetaphi[iz][ic][ich][MPT] ->Add(hetaphi[iz][ic][ich][ipt]);
        }
	for(int ipt=0;ipt<NPT_TrkWei;ipt++){//normalize
          TH2F *htmp =(TH2F*)hetaphi[iz][ic][ich][ipt]->Clone();
          hproj = (TH1D*) htmp->ProjectionY("proj"); hproj->Scale(1./NPHI_TrkWei);
          for(int iy=1;iy<=NETA_TrkWei;iy++){
            double nor = hproj->GetBinContent(iy);
            for(int ix=1;ix<=NPHI_TrkWei;ix++){
              if(nor!=0.) hetaphi[iz][ic][ich][ipt]->SetBinContent(ix,iy,htmp->GetBinContent(ix,iy)/nor);
              if(nor!=0.) hetaphi[iz][ic][ich][ipt]->SetBinError(ix,iy,htmp->GetBinError(ix,iy)/nor);
            }
          }
        }
      }
    }
  }
  */
  
}

void extractor_flat::Run_mu() {
  int evts_min=NEv_id*NEv;
  int evts_max=(NEv_id+1)*NEv;

  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  cout<<"Run events from "<<evts_min<<" to "<<evts_max<<endl;
  //FCalEt Shift
  for(int ic=0;ic<NCENT_1p;ic++){
    for(int iz=0;iz<NZ;iz++){
      evt[iz][ic]=hevts->GetBinContent(ic+1,iz+1);
    }
  }
  for(int iz=0;iz<NZ;iz++){
    for(int ic=0;ic<NCENT_1p;ic++){
      FCal_Et_m[iz][ic] = hEt_vz[iz]->GetBinContent(ic+1);//fcal mean
    }
  }
  //offset correction
  for(int id=0;id<NDET;id++){
    for(int iz=0;iz<NZ;iz++){
      for(int ic=0;ic<NCENT_1p;ic++){
	Qw_m_FCal[id][iz][ic] = hQw_FCal[id][iz]->GetBinContent(ic+1);
      }
      for(int ihar=0;ihar<NHAR;ihar++){//1,2,3,4
	//hqx_FCal[id][iz][ihar]->Divide(hQx_FCal[id][iz][ihar],hQw_FCal[id][iz]);
	//hqy_FCal[id][iz][ihar]->Divide(hQy_FCal[id][iz][ihar],hQw_FCal[id][iz]);
	for(int ic=0;ic<NCENT_1p;ic++){
	  Qx_m_FCal[id][iz][ic][ihar] = hQx_FCal[id][iz][ihar]->GetBinContent(ic+1);
	  Qy_m_FCal[id][iz][ic][ihar] = hQy_FCal[id][iz][ihar]->GetBinContent(ic+1);
	  //qx_m_FCal[id][iz][ic][ihar] = hqx[id][iz][ihar]->GetBinContent(ic+1);
	  //qy_m_FCal[id][iz][ic][ihar] = hqy[id][iz][ihar]->GetBinContent(ic+1);
	}
      }
    }
  }

  //Event Loop starts
  int cnt_ev=0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    //Long64_t ientry = LoadTree(jentry);
    //if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);
    if (nb <=0)        {cout<<"Nbytes became <=0"<<endl;break;}
    
    if(jentry<evts_min) continue;
    if(jentry>=evts_max) break;
    if(jentry%10000==0) {cout<<"running event "<<jentry<<endl ;}
    //if(jentry>10) break;

    //-------------------------------------------------------------------------------------------------
    centb = get_centb(FCal_Et);
    vzb=get_vzb(vx_z);
    if (centb < 0) continue;
    if (vzb < 0) continue;
    nevents++;

    FillGlobal_mu();

    //Do FCal Recenter Correction
    int iz = (vx_z+100.)/25;if(iz<0) iz=0; if(iz>NZ-1) iz=NZ-1;
    int ic = Centrality;
    int N = evt[iz][ic];
    if(ic<100&&N>5){//only calib events for bins that have more than 5 events.
      double scal = FCal_Et/FCal_Et_m[iz][ic];
      for(int id=0;id<NDET;id++){
	for(int ihar=0;ihar<NHAR;ihar++){//2,3,4,5
	  //subtract the contribution from present event, avoid self-bias.
	  double Qxm = (Qx_m_FCal[id][iz][ic][ihar]*N-Qx_FCal[id][ihar]);
	  double Qym = (Qy_m_FCal[id][iz][ic][ihar]*N-Qy_FCal[id][ihar]);
	  double Qwm = (Qw_m_FCal[id][iz][ic]*N-Qw_FCal[id]);
	 
	  double shiftx =  Qxm/Qwm*Qw_FCal[id]*scal;
	  double shifty =  Qym/Qwm*Qw_FCal[id]*scal;
	  Qx_FCal[id][ihar] -= shiftx;
	  Qy_FCal[id][ihar] -= shifty;
	  hQxS_FCal[id][iz][ihar]->Fill(ic,shiftx);
	  hQyS_FCal[id][iz][ihar]->Fill(ic,shifty);
	  //hqxS_FCal[id][iz][ihar]->Fill(ic,Qxm/Qwm);
	  //hqyS_FCal[id][iz][ihar]->Fill(ic,Qym/Qwm);
	}
      }
    }

    
    //Do ZDC Recenter and Flattening Correction
    for (int iside=0; iside<NSIDE+1; iside++) {
      for (int ihar=0; ihar<NHAR; ihar++) {
	double PsiZ_raw = atan2(Qy_Zdc[iside][ihar], Qx_Zdc[iside][ihar])/double(ihar+1);
	Qx_Zdc[iside][ihar]  = ( Qx_Zdc[iside][ihar] - hQx_Zdc[iside][ihar]->GetBinContent(ic+1) );
        Qy_Zdc[iside][ihar]  = ( Qy_Zdc[iside][ihar] - hQy_Zdc[iside][ihar]->GetBinContent(ic+1) );	
        double PsiZ_rec = atan2(Qy_Zdc[iside][ihar], Qx_Zdc[iside][ihar])/double(ihar+1);
        double psi_rec2pi = (ihar+1)*PsiZ_rec;
	double deltaPsi = 0;
	
	for (int ik=0; ik<NK; ik++) {
          double flatcos = hFlatCos_Zdc[ic][iside][ihar]->GetBinContent(ik+1);
          double flatsin = hFlatSin_Zdc[ic][iside][ihar]->GetBinContent(ik+1);
          double cosPsi_rec = cos( (ik+1)*psi_rec2pi );
          double sinPsi_rec = sin( (ik+1)*psi_rec2pi );
          deltaPsi += (-flatsin*cosPsi_rec + flatcos*sinPsi_rec)*2/(ik+1);
        }
	Qx_Zdc[iside][ihar]=cos(psi_rec2pi+deltaPsi);
	Qy_Zdc[iside][ihar]=sin(psi_rec2pi+deltaPsi);
	//double PsiF_Zdc = atan2( sin(psiC2pi+deltaPsi), cos(psiC2pi+deltaPsi) )/double(ihar+1);
	double PsiZ_flat = atan2( Qy_Zdc[iside][ihar], Qx_Zdc[iside][ihar])/double(ihar+1);
	hPsiZ_raw[centb][iside][ihar]->Fill(PsiZ_raw);
	hPsiZ_rec[centb][iside][ihar]->Fill(PsiZ_rec);
	hPsiZ_flat[centb][iside][ihar]->Fill(PsiZ_flat);
      }
    }

    //float Psi_FCal[NSIDE+1][NHAR];
    for (int iside=0; iside<NSIDE+1; iside++) {
      for (int ihar=0; ihar<NHAR; ihar++) {
	int ieta;
	if(iside==0) Psi_FCal[iside][ihar] = atan2(Qy_FCal[0][ihar],Qx_FCal[0][ihar])/double(ihar+1);
	if(iside==1) Psi_FCal[iside][ihar] = atan2(Qy_FCal[1][ihar],Qx_FCal[1][ihar])/double(ihar+1);
	if(iside==2) Psi_FCal[iside][ihar] = atan2(Qy_FCal[0][ihar]+Qy_FCal[1][ihar],Qx_FCal[0][ihar]+Qx_FCal[1][ihar])/double(ihar+1);
	hPsi_FCal[centb][iside][ihar]->Fill(Psi_FCal[iside][ihar]);
      }
    }
    //float Psi_Zdc[NSIDE+1][NHAR];
    for (int iside=0; iside<NSIDE+1; iside++) {
      for (int ihar=0; ihar<NHAR; ihar++) {
	//int side_sign=1;
	//if(iside==1) side_sign=-1;
	//Psi_Zdc[iside][ihar] = atan2(Qy_Zdc[iside][ihar],side_sign*Qx_Zdc[iside][ihar])/double(ihar+1);
	Psi_Zdc[iside][ihar] = atan2(Qy_Zdc[iside][ihar],Qx_Zdc[iside][ihar])/double(ihar+1);
	hPsi_Zdc[centb][iside][ihar]->Fill(Psi_Zdc[iside][ihar]);
      }
    }

    EVENT_PTR_EP ev = new Event_EP(jentry);
    ev->icent = centb;
    ev->ivz = vzb;
    for (int iside=0; iside<NSIDE+1; iside++) {
      for (int ihar=0; ihar<NHAR; ihar++) {
	ev->Psi_FCal[iside][ihar]=Psi_FCal[iside][ihar];
	ev->Psi_Zdc[iside][ihar]=Psi_Zdc[iside][ihar];
      }
    }
    double dPsi_FCal,dPsi_Zdc;
    for (int ihar=0; ihar<NHAR; ihar++) {
      //hRes     [centb][ihar]->Fill(cos( (ihar+1)* (Psi.at (0).at(ihar) - Psi.at (1).at(ihar)) ));
      //hRes_recn[centb][ihar]->Fill(cos( (ihar+1)* (PsiC.at(0).at(ihar) - PsiC.at(1).at(ihar)) ));
      //hRes_flat[centb][ihar]->Fill(cos( (ihar+1)* (Psi_Zdc.at(0).at(ihar) - Psi_Zdc.at(1).at(ihar)) ));
      hRes_flat_FCal[centb][ihar]->Fill(cos( (ihar+1)*(Psi_FCal[0][ihar] - Psi_FCal[1][ihar]) ));
      hRes_flat_Zdc[centb][ihar]->Fill(cos( (ihar+1)*(Psi_Zdc[0][0] - Psi_Zdc[1][0]) ));

      dPsi_FCal = (ihar+1)*(Psi_FCal[1][ihar]-Psi_FCal[0][ihar] );
      hResCent_FCal_fg[ihar]->Fill(centb,cos(dPsi_FCal));
      int nn;
      //Wrap to -2pi to 2pi
      nn=fabs(dPsi_FCal)/(2*PI);
      if(dPsi_FCal < 0.) dPsi_FCal += nn*2*PI;
      if(dPsi_FCal > 0.) dPsi_FCal -= nn*2*PI;
      //Wrap to -pi to pi
      if (dPsi_FCal < -PI) dPsi_FCal += 2*PI;
      if (dPsi_FCal >  PI) dPsi_FCal -= 2*PI;
      dPsi_Zdc = (ihar+1)*(Psi_Zdc[1][0]-Psi_Zdc[0][0] );
      //dPsi_Zdc = (ihar+1)*(Psi_Zdc[1][ihar]-Psi_Zdc[0][ihar] );
      hResCent_Zdc_fg[ihar]->Fill(centb,cos(dPsi_Zdc));
      //Wrap to -2pi to 2pi
      nn=fabs(dPsi_Zdc)/(2*PI);
      if(dPsi_Zdc < 0.) dPsi_Zdc += nn*2*PI;
      if(dPsi_Zdc > 0.) dPsi_Zdc -= nn*2*PI;
      //Wrap to -pi to pi
      if (dPsi_Zdc < -PI) dPsi_Zdc += 2*PI;
      if (dPsi_Zdc >  PI) dPsi_Zdc -= 2*PI;
      hDphi_FCal_fg[centb][ihar]->Fill(dPsi_FCal);
      hDphi_Zdc_fg[centb][ihar]->Fill(dPsi_Zdc);
#ifdef WithBG
      //Background calculation 
      Pool_size =Pool_EP[centb][vzb].size();
      if(Pool_size>0){
	for (int isize=0; isize<Pool_size; isize++) {
	  dPsi_FCal = (ihar+1)*( Psi_FCal[1][ihar] - (Pool_EP[centb][vzb].at(isize)->Psi_FCal[0][ihar]) );
	  hResCent_FCal_bg[ihar]->Fill(centb,cos(dPsi_FCal));
	  //Wrap to -2pi to 2pi
	  nn=fabs(dPsi_FCal)/(2*PI);
	  if(dPsi_FCal < 0.) dPsi_FCal += nn*2*PI;
	  if(dPsi_FCal > 0.) dPsi_FCal -= nn*2*PI;
	  //Wrap to -pi to pi
	  if (dPsi_FCal < -PI) dPsi_FCal += 2*PI;
	  if (dPsi_FCal >  PI) dPsi_FCal -= 2*PI;
	  hDphi_FCal_bg[centb][ihar]->Fill(dPsi_FCal);
	  
	  dPsi_FCal = (ihar+1)*( Psi_FCal[0][ihar] - (Pool_EP[centb][vzb].at(isize)->Psi_FCal[1][ihar]) );
	  hResCent_FCal_bg[ihar]->Fill(centb,cos(dPsi_FCal));
	  //Wrap to -2pi to 2pi
          nn=fabs(dPsi_FCal)/(2*PI);
          if(dPsi_FCal < 0.) dPsi_FCal += nn*2*PI;
          if(dPsi_FCal > 0.) dPsi_FCal -= nn*2*PI;
          //Wrap to -pi to pi
          if (dPsi_FCal < -PI) dPsi_FCal += 2*PI;
          if (dPsi_FCal >  PI) dPsi_FCal -= 2*PI;
	  hDphi_FCal_bg[centb][ihar]->Fill(-dPsi_FCal);

	  dPsi_Zdc = (ihar+1)*( Psi_Zdc[1][0] - (Pool_EP[centb][vzb].at(isize)->Psi_Zdc[0][0]) );
	  //dPsi_Zdc = (ihar+1)*( Psi_Zdc[1][ihar] - (Pool_EP[centb][vzb].at(isize)->Psi_Zdc[0][ihar]) );
	  hResCent_Zdc_bg[ihar]->Fill(centb,cos(dPsi_Zdc));
	  //Wrap to -2pi to 2pi
	  nn=fabs(dPsi_Zdc)/(2*PI);
	  if(dPsi_Zdc < 0.) dPsi_Zdc += nn*2*PI;
	  if(dPsi_Zdc > 0.) dPsi_Zdc -= nn*2*PI;
	  //Wrap to -pi to pi
	  if (dPsi_Zdc < -PI) dPsi_Zdc += 2*PI;
	  if (dPsi_Zdc >  PI) dPsi_Zdc -= 2*PI;
	  hDphi_Zdc_bg[centb][ihar]->Fill(dPsi_Zdc);

	  dPsi_Zdc = (ihar+1)*( Psi_Zdc[0][0] - (Pool_EP[centb][vzb].at(isize)->Psi_Zdc[1][0]) );
	  //dPsi_Zdc = (ihar+1)*( Psi_Zdc[0][ihar] - (Pool_EP[centb][vzb].at(isize)->Psi_Zdc[1][ihar]) );
	  hResCent_Zdc_bg[ihar]->Fill(centb,cos(dPsi_Zdc));
	  //Wrap to -2pi to 2pi
	  nn=fabs(dPsi_Zdc)/(2*PI);
	  if(dPsi_Zdc < 0.) dPsi_Zdc += nn*2*PI;
	  if(dPsi_Zdc > 0.) dPsi_Zdc -= nn*2*PI;
	  //Wrap to -pi to pi
	  if (dPsi_Zdc < -PI) dPsi_Zdc += 2*PI;
	  if (dPsi_Zdc >  PI) dPsi_Zdc -= 2*PI;
	  hDphi_Zdc_bg[centb][ihar]->Fill(-dPsi_Zdc);
	}
      }
#endif
    }

    if(Centrality<60) FillTrack_mu(Psi_FCal,Psi_Zdc);

#ifdef withBG    
    //Mixed Event Pool
    Pool_size =Pool_EP[centb][vzb].size();
    if(Pool_size>=STORE_DEP){
      delete Pool_EP[centb][vzb].at(0);
      Pool_EP[centb][vzb].erase( (Pool_EP[centb][vzb]).begin() );
    }
    Pool_EP[centb][vzb].push_back( ev );
    //delete ev;
#endif
    cnt_ev++;
  }//Event Loop ends
  cout<<"Finished running all events."<<endl;
}


void extractor_flat::FillGlobal_mu() {

  double eTot[NSIDE][NGAIN] = {{0}};
  double ePix[NSIDE][NGAIN] = {{0}};

  for (int ichn=0; ichn<Zdc_n; ichn++) {

    int side = Zdc_Side->at(ichn)+1;
    if (side==2) side = 1;
    if (!(side ==0 || side == 1)) {std::cout << "Unknown side: chId = " << ichn << ": side = " << side << std::endl; continue;}

    int module = Zdc_Module->at(ichn);
    if (module < 0 || module > 3) {std::cout << "Unknown module: chId = " << ichn << ": module = " << module << std::endl; continue;}

    int chId = ichn;
    if (side == 0 && module != 0) chId = ichn - 2;
    if (side == 1 && module != 0) chId = ichn - 30;
    //-------------------------------------------------------------------------------------------------

    if (Zdc_Type->at(ichn) == 0) {
      eTot[side][0] += Zdc_Energy_LG->at(ichn);
      eTot[side][1] += Zdc_Energy_HG->at(ichn);
    }
    if (Zdc_Type->at(ichn) == 1 && module !=0) {
      ePix[side][0] += Zdc_Energy_LG->at(ichn);
      ePix[side][1] += Zdc_Energy_HG->at(ichn);
    }
    //-------------------------------------------------------------------------------------------------
  }

  //-------------------------------------------------------------------------------------------------
  for (int iside=0; iside<NSIDE; iside++ ) {
    h_eTot           [iside]->Fill(eTot[iside][0]);
    h_ePix           [iside]->Fill(ePix[iside][0]);
    hEt_eTot         [iside]->Fill(FCal_Et*1000, eTot[iside][0]);
    hEt_ePix         [iside]->Fill(FCal_Et*1000, ePix[iside][0]);
    hNtrk_eTot       [iside]->Fill(ntrkQ, eTot[iside][0]);
    hNtrk_ePix       [iside]->Fill(ntrkQ, ePix[iside][0]);
  }

  hEt  ->Fill(FCal_Et*1000);
  hNtrk->Fill(ntrkQ);
  hCent->Fill(centb);
  hEt_ntrk->Fill(FCal_Et*1000, ntrkQ);
  //-------------------------------------------------------------------------------------------------
}

int get_dPhiPsibin(float dphi){
  float dphi_min=PI/4.;
  int idphi=-1;
  if(dphi>=0. && dphi<dphi_min) idphi=0;
  else if(dphi>=dphi_min && dphi<2.*dphi_min) idphi=1;
  else if(dphi>=2.*dphi_min && dphi<3.*dphi_min) idphi=2;
  else if(dphi>=3.*dphi_min && dphi<=4.*dphi_min) idphi=3;
  else idphi=-1;
  return idphi;
}

void extractor_flat::FillTrack_mu(float Psi_FCal[3][4],float Psi_Zdc[3][4]) {

  hMu_n_raw->Fill(Mu_n);
  int ntrk_mu=0;
  //Run track loop
  for(int itrk=0; itrk<Mu_n; itrk++){
    hMu_eta_raw->Fill(Mu_eta->at(itrk));
    hMu_phi_raw->Fill(Mu_phi->at(itrk));
    hMu_pt_raw->Fill(Mu_pt->at(itrk));
    hMu_eloss_all_raw->Fill(Mu_eloss->at(itrk));
    hMu_ms_phi_raw->Fill(Mu_ms_phi->at(itrk));
    hMu_ms_theta_raw->Fill(Mu_ms_theta->at(itrk));
    hMu_ms_qoverp_raw->Fill(fabs(1./Mu_ms_qoverp->at(itrk)/1000.));
    if(Mu_pt->at(itrk)<4. || Mu_pt->at(itrk)>14.) continue;
    if(fabs(Mu_eta->at(itrk))>2.) continue;
    int ipt=get_ptBin_Mu(Mu_pt->at(itrk));
    int ieta=get_etaBin_Mu(Mu_eta->at(itrk));
    if(ipt<0) continue;
    if(ieta<0) continue;
    ntrk_mu++;
    //Fill track hists after selection cuts
    hMu_eta->Fill(Mu_eta->at(itrk));
    hMu_phi->Fill(Mu_phi->at(itrk));
    hMu_pt->Fill(Mu_pt->at(itrk));
    hMu_eloss_all->Fill(Mu_eloss->at(itrk));
    hMu_ms_phi->Fill(Mu_ms_phi->at(itrk));
    hMu_ms_theta->Fill(Mu_ms_theta->at(itrk));
    hMu_ms_qoverp->Fill(fabs(1./Mu_ms_qoverp->at(itrk)/1000.));

    //Fill phi-Psi hists
    int nn;
    float dphi_FCal,dphi_Zdc;
    for (int iside=0; iside<NSIDE+1; iside++){
      for (int ihar=0; ihar<NHAR; ihar++){
	//FCal
	dphi_FCal = (ihar+1)*(Mu_phi->at(itrk)-Psi_FCal[iside][ihar]);
	hMu_Vn_FCal_pt[centb][iside][ihar]->Fill(Mu_pt->at(itrk),cos(dphi_FCal));
	hMu_Vn_FCal_pt[NCENT][iside][ihar]->Fill(Mu_pt->at(itrk),cos(dphi_FCal));
	hMu_Vn_FCal_cent[ipt][iside][ihar]->Fill(Centrality,cos(dphi_FCal));
	//Wrap to -2pi to 2pi
	nn=fabs(dphi_FCal)/(2*PI);
	if(dphi_FCal < 0.) dphi_FCal += nn*2*PI;
	if(dphi_FCal > 0.) dphi_FCal -= nn*2*PI;
	//Wrap to -pi to pi
	if (dphi_FCal < -PI) dphi_FCal += 2*PI;
	if (dphi_FCal >  PI) dphi_FCal -= 2*PI;
	//int idphi_FCal=get_dPhiPsibin(fabs(dphi_FCal));
	hMu_dphi_yield_FCal[centb][iside][ipt][ihar]->Fill(fabs(dphi_FCal));
	hMu_dphi_eloss_FCal[centb][iside][ipt][ihar]->Fill(fabs(dphi_FCal),Mu_eloss->at(itrk));
	
	//Zdc
	dphi_Zdc = (ihar+1)*(Mu_phi->at(itrk)-Psi_Zdc[iside][0]);
	hMu_Vn_Zdc_eta[centb][iside][ihar]->Fill(Mu_eta->at(itrk),cos(dphi_Zdc));
	hMu_Vn_Zdc_eta[NCENT][iside][ihar]->Fill(Mu_eta->at(itrk),cos(dphi_Zdc));
	hMu_Vn_FCal_cent[ieta][iside][ihar]->Fill(Centrality,cos(dphi_Zdc));
	//Wrap to -2pi to 2pi
	nn=fabs(dphi_Zdc)/(2*PI);
	if(dphi_Zdc < 0.) dphi_Zdc += nn*2*PI;
	if(dphi_Zdc > 0.) dphi_Zdc -= nn*2*PI;
	//Wrap to -pi to pi
	if (dphi_Zdc < -PI) dphi_Zdc += 2*PI;
	if (dphi_Zdc >  PI) dphi_Zdc -= 2*PI;
	//int idphi_Zdc=get_dPhiPsibin(fabs(dphi_Zdc));
	hMu_dphi_yield_Zdc[centb][iside][ieta][ihar]->Fill(fabs(dphi_Zdc));
	hMu_dphi_eloss_Zdc[centb][iside][ieta][ihar]->Fill(fabs(dphi_Zdc),Mu_eloss->at(itrk));
      }
    }

  }//Track loop ends
  hMu_n->Fill(ntrk_mu);

}


void extractor_flat::Init_mu() {
  //set_ptBins ();
  //set_etaBins();
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  for (int iside=0; iside<NSIDE; iside++) {
    sprintf(name, "h_eTot_side%d", iside);
    h_eTot[iside] = new TH1D(name, "", 4000, -0.5, 3999.5);
    sprintf(name, "h_ePix_side%d", iside);
    h_ePix[iside] = new TH1D(name, "", 4000, -0.5, 3999.5);
   
    sprintf(name, "h_et_eTot_side%d", iside);
    hEt_eTot[iside] = new TH2F(name, "", 500, 0, 5000, 400, 0, 4000);
    sprintf(name, "h_et_ePix_side%d", iside);
    hEt_ePix[iside] = new TH2F(name, "", 500, 0, 5000, 200, 0, 2000);
      
    sprintf(name, "h_ntrk_eTot_side%d", iside);
    hNtrk_eTot[iside] = new TH2F(name, "", 400, 0, 4000, 400, 0, 4000);
    sprintf(name, "h_ntrk_ePix_side%d", iside);
    hNtrk_ePix[iside] = new TH2F(name, "", 400, 0, 4000, 200, 0, 2000);
  }
  
  hEt      = new TH1D("hEt", "", 5000, 0, 5000);
  hNtrk    = new TH1D("hNtrk", "", 5000, 0, 5000);
  hCent    = new TH1D("hCent", "", NCENT, -0.5, NCENT-0.5);
  hEt_ntrk = new TH2F("h_et_ntrk", "", 400, 0, 4000, 400, 0, 4000);

  for(int id=0;id<NDET;id++){
    for(int iz=0;iz<NZ;iz++){
      for(int ihar=0;ihar<NHAR;ihar++){//1,2,3,4,
	sprintf(name,"hQxS_FCal_det%d_z%d_h%d",id,iz,ihar);
	hQxS_FCal[id][iz][ihar]=new TProfile(name,"",NCENT,-0.5,NCENT-0.5,"s");
	sprintf(name,"hQyS_FCal_det%d_z%d_h%d",id,iz,ihar);
	hQyS_FCal[id][iz][ihar]=new TProfile(name,"",NCENT,-0.5,NCENT-0.5,"s");
      }
    }
  }
  //ZDC Correction
  for (int icent=0; icent<NCENT; icent++) {
    for (int iside=0; iside<NSIDE+1; iside++) {
      for (int ihar=0; ihar<NHAR; ihar++) {
	sprintf(name, "hPsiZ_raw_cent%d_side%d_har%d", icent, iside, ihar);
	hPsiZ_raw[icent][iside][ihar] = new TH1D(name, "", 100, -2*PI, 2*PI); 

	sprintf(name, "hPsiZ_rec_cent%d_side%d_har%d", icent, iside, ihar);
	hPsiZ_rec[icent][iside][ihar] = new TH1D(name, "", 100, -2*PI, 2*PI);

	sprintf(name, "hPsiZ_flat_cent%d_side%d_har%d", icent, iside, ihar);
	hPsiZ_flat[icent][iside][ihar] = new TH1D(name, "", 100, -2*PI, 2*PI);
      }
    }
  }
  for (int icent=0; icent<NCENT; icent++) {
    for (int iside=0; iside<NSIDE+1; iside++) {
      for (int ihar=0; ihar<NHAR; ihar++) {
	sprintf(name, "hPsi_FCal_cent%d_side%d_har%d", icent, iside, ihar);
	hPsi_FCal[icent][iside][ihar] = new TH1D(name, "", 100, -2*PI, 2*PI);
	sprintf(name, "hPsi_Zdc_cent%d_side%d_har%d", icent, iside, ihar);
	hPsi_Zdc[icent][iside][ihar] = new TH1D(name, "", 100, -2*PI, 2*PI);
      }
    }
  }
   

  //Mean pt and eta distributions
  for(int icent=0;icent<NCENT;icent++){
    for(int ich=0;ich<NCH+1;ich++){
      for(int iside=0;iside<NSIDE+1;iside++){
        sprintf(name,"hMean_pt_cent%d_ch%d_side%d",icent,ich,iside);
        hMean_pt[icent][ich][iside]=new TProfile(name,"",NPT,0.-0.5,NPT-0.5);
        sprintf(name,"hMean_eta_cent%d_ch%d_side%d",icent,ich,iside);
        hMean_eta[icent][ich][iside]=new TProfile(name,"",NETA,0.-0.5,NETA-0.5);
      }
    }
  }

  double x1=3.2;
  //int n1=640;
  int n1=64;
  int n2=2;
  for (int icent=0; icent<NCENT; icent++) {
    for (int ihar=0; ihar<NHAR; ihar++) {
      sprintf(name, "hDphi_FCal_fg_cent%d_har%d", icent, ihar);
      hDphi_FCal_fg[icent][ihar] = new TH1D(name, "", n2*n1,-n2*x1,n2*x1);
      sprintf(name, "hDphi_FCal_bg_cent%d_har%d", icent, ihar);
      hDphi_FCal_bg[icent][ihar] = new TH1D(name, "", n2*n1,-n2*x1,n2*x1);
      sprintf(name, "hDphi_Zdc_fg_cent%d_har%d", icent, ihar);
      hDphi_Zdc_fg[icent][ihar] = new TH1D(name, "", n2*n1,-n2*x1,n2*x1);
      sprintf(name, "hDphi_Zdc_bg_cent%d_har%d", icent, ihar);
      hDphi_Zdc_bg[icent][ihar] = new TH1D(name, "", n2*n1,-n2*x1,n2*x1);
    }
  }
  
  for (int icent=0; icent<NCENT; icent++) {
    for (int ihar=0; ihar<NHAR; ihar++) {

      sprintf(name, "hRes_cent%d_har%d", icent, ihar);
      hRes[icent][ihar] = new TH1D(name, "", 200, -1.0, 1.0);

      sprintf(name, "hRes_recn_cent%d_har%d", icent, ihar);
      hRes_recn[icent][ihar] = new TH1D(name, "", 200, -1.0, 1.0);

      sprintf(name, "hRes_flat_FCal_cent%d_har%d", icent, ihar);
      hRes_flat_FCal[icent][ihar] = new TH1D(name, "", 200, -1.0, 1.0);

      sprintf(name, "hRes_flat_Zdc_cent%d_har%d", icent, ihar);
      hRes_flat_Zdc[icent][ihar] = new TH1D(name, "", 200, -1.0, 1.0);
    }
  }

  for (int ihar=0; ihar<NHAR; ihar++) {
    sprintf(name, "hResCent_FCal_fg_har%d", ihar);
    hResCent_FCal_fg[ihar] = new TProfile(name, "", NCENT, 0, NCENT);
    sprintf(name, "hResCent_FCal_bg_har%d", ihar);
    hResCent_FCal_bg[ihar] = new TProfile(name, "", NCENT, 0, NCENT);

    sprintf(name, "hResCent_Zdc_fg_har%d", ihar);
    hResCent_Zdc_fg[ihar] = new TProfile(name, "", NCENT, 0, NCENT);
    sprintf(name, "hResCent_Zdc_bg_har%d", ihar);
    hResCent_Zdc_bg[ihar] = new TProfile(name, "", NCENT, 0, NCENT);
  }
  //Muon hists
  for(int icent=0;icent<NCENT;icent++){
    for (int iside=0; iside<NSIDE+1; iside++){
      for (int ihar=0; ihar<NHAR; ihar++){
	for(int ipt=0;ipt<NPT_Mu;ipt++){
	  //Yield for n(\phi-\Psi)
	  sprintf(name,"hMu_dphi_yield_FCal_cent%d_side%d_pt%d_har%d",icent,iside,ipt,ihar);
	  hMu_dphi_yield_FCal[icent][iside][ipt][ihar]=new TH1D(name,"",4,0.,PI);
	  //Momentum Imbalance for n(\phi-\Psi)
	  sprintf(name,"hMu_dphi_eloss_FCal_cent%d_side%d_pt%d_har%d",icent,iside,ipt,ihar);
	  hMu_dphi_eloss_FCal[icent][iside][ipt][ihar]=new TH2F(name,"",4,0.,PI,100,-1.,1.);
	}
	for(int ieta=0;ieta<NETA_Mu;ieta++){
	  //Yield for n(\phi-\Psi)
	  sprintf(name,"hMu_dphi_yield_Zdc_cent%d_side%d_eta%d_har%d",icent,iside,ieta,ihar);
	  hMu_dphi_yield_Zdc[icent][iside][ieta][ihar]=new TH1D(name,"",4,0.,PI);
	  //Momentum Imbalance for n(\phi-\Psi)
	  sprintf(name,"hMu_dphi_eloss_Zdc_cent%d_side%d_eta%d_har%d",icent,iside,ieta,ihar);
	  hMu_dphi_eloss_Zdc[icent][iside][ieta][ihar]=new TH2F(name,"",4,0.,PI,100,-1.,1.);
	}
      }
    }
  }
  //Muon FLow hists
  for (int iside=0; iside<NSIDE+1; iside++){
    for (int ihar=0; ihar<NHAR; ihar++){
      for(int icent=0;icent<NCENT+1;icent++){
	sprintf(name,"hMu_Vn_FCal_pt_cent%d_side%d_har%d",icent,iside,ihar);
	hMu_Vn_FCal_pt[icent][iside][ihar]=new TProfile(name,"",NPT_Mu,pt_mat_Mu);
	sprintf(name,"hMu_Vn_Zdc_eta_cent%d_side%d_har%d",icent,iside,ihar);
	hMu_Vn_Zdc_eta[icent][iside][ihar]=new TProfile(name,"",NETA_Mu,eta_mat_Mu);
      }
      for(int ipt=0;ipt<NPT_Mu;ipt++){
	sprintf(name,"hMu_Vn_FCal_cent_pt%d_side%d_har%d",ipt,iside,ihar);
        hMu_Vn_FCal_cent[ipt][iside][ihar]=new TProfile(name,"",NCENT,cent_mat_double);
      }
      for(int ieta=0;ieta<NETA_Mu;ieta++){
	sprintf(name,"hMu_Vn_Zdc_cent_eta%d_side%d_har%d",ieta,iside,ihar);
        hMu_Vn_Zdc_cent[ieta][iside][ihar]=new TProfile(name,"",NCENT,cent_mat_double);
      }
    }
  }
  //Muon Track hists
  //Raw hists before selection
  hMu_n_raw = new TH1D("hMu_n_raw", "", 50, 0., 50.);
  hMu_pt_raw = new TH1D("hMu_pt_raw","",1000,0.,25.);
  hMu_phi_raw = new TH1D("hMu_phi_raw","",50,-PI2,PI2);
  hMu_eta_raw = new TH1D("hMu_eta_raw","",100,-2.5,2.5);
  hMu_eloss_all_raw = new TH1D("hMu_eloss_all_raw","",100,-1.,1.);
  hMu_ms_phi_raw = new TH1D("hMu_ms_phi_raw","",50,-PI2,PI2);
  hMu_ms_theta_raw = new TH1D("hMu_ms_theta_raw","",25,0.,PI);
  hMu_ms_qoverp_raw = new TH1D("hMu_ms_qoverp_raw","",250,0.,25.);
  //After track selection cuts
  hMu_n = new TH1D("hMu_n", "", 50, 0., 50.);
  hMu_pt = new TH1D("hMu_pt","",1000,0.,25.);
  hMu_phi = new TH1D("hMu_phi","",50,-PI2,PI2);
  hMu_eta = new TH1D("hMu_eta","",100,-2.5,2.5);
  hMu_eloss_all = new TH1D("hMu_eloss_all","",100,-1.,1.);
  hMu_ms_phi = new TH1D("hMu_ms_phi","",50,-PI2,PI2);
  hMu_ms_theta = new TH1D("hMu_ms_theta","",25,0.,PI);
  hMu_ms_qoverp = new TH1D("hMu_ms_qoverp","",250,0.,25.);
    
  for(int icent=0;icent<NCENT;icent++){
    for(int ipt=0;ipt<NPT;ipt++){
      sprintf(name,"hTrk_eff2_cent%d_pt%d",icent,ipt);
      hTrk_eff2[icent][ipt]=new TProfile(name,"",50,-2.5,2.5);
    }
  }
  hTrk_eff=new TH1D("hTrk_eff","",200,-2.,2.);
  //Correlation hists
  //hEtaPhi=new TH2F("hEtaPhi","",100,-2.5,2.5,50,0.,7.); hEtaPhi->Sumw2();
  //hEtaPhi_wei=new TH2F("hEtaPhi_wei","",100,-2.5,2.5,50,0.,7.); hEtaPhi_wei->Sumw2();

}


void extractor_flat::SaveHistos_mu() {
  fOut_mu->cd();

  for (int iside=0; iside<NSIDE; iside++) {
    h_eTot[iside]->Write();
    h_ePix[iside]->Write();
    hEt_eTot[iside]->Write();
    hEt_ePix[iside]->Write();
    hNtrk_eTot[iside]->Write();
    hNtrk_ePix[iside]->Write();
  }
  //FCal Correction
  for(int id=0;id<NDET;id++){
    for(int iz=0;iz<NZ;iz++){
      for(int ihar=0;ihar<NHAR;ihar++){//1,2,3,4,
	hQxS_FCal[id][iz][ihar]->Write();
	hQyS_FCal[id][iz][ihar]->Write();
      }
    }
  }
  //ZDC Correction
  for (int icent=0; icent<NCENT; icent++) {
    for (int iside=0; iside<NSIDE+1; iside++) {
      for (int ihar=0; ihar<NHAR; ihar++) {
	hPsiZ_raw[icent][iside][ihar]->Write();
	hPsiZ_rec[icent][iside][ihar]->Write();
	hPsiZ_flat[icent][iside][ihar]->Write();
      }
    }
  }

  for (int icent=0; icent<NCENT; icent++) {
    for (int iside=0; iside<NSIDE+1; iside++) {
      for (int ihar=0; ihar<NHAR; ihar++) {
	hPsi_FCal[icent][iside][ihar]->Write();
	hPsi_Zdc[icent][iside][ihar]->Write();
      }
    }
  }
  for (int icent=0; icent<NCENT; icent++) {
    for (int ihar=0; ihar<NHAR; ihar++) {
      hDphi_FCal_fg[icent][ihar]->Write();
      hDphi_FCal_bg[icent][ihar]->Write();
      hRes_flat_FCal[icent][ihar]->Write();
      
      hDphi_Zdc_fg[icent][ihar]->Write();
      hDphi_Zdc_bg[icent][ihar]->Write();
      hRes_flat_Zdc[icent][ihar]->Write();
    }
  }
  for(int icent=0;icent<NCENT;icent++){
    for(int ich=0;ich<NCH+1;ich++){
      for(int iside=0;iside<NSIDE+1;iside++){
        hMean_pt[icent][ich][iside]->Write();
        hMean_eta[icent][ich][iside]->Write();
      }
    }
  }

  hEt->Write();
  hNtrk->Write();
  hCent->Write();
  hEt_ntrk->Write();

  //Muon hists
  for(int icent=0;icent<NCENT;icent++){
    for (int iside=0; iside<NSIDE+1; iside++){
      for (int ihar=0; ihar<NHAR; ihar++){
	for(int ipt=0;ipt<NPT_Mu;ipt++){
	  hMu_dphi_yield_FCal[icent][iside][ipt][ihar]->Write();
	  hMu_dphi_eloss_FCal[icent][iside][ipt][ihar]->Write();	
	}
	for(int ieta=0;ieta<NETA_Mu;ieta++){
	  hMu_dphi_yield_Zdc[icent][iside][ieta][ihar]->Write();
	  hMu_dphi_eloss_Zdc[icent][iside][ieta][ihar]->Write();
	}
      }
    }
  }
  //Muon Flow hists
  for (int iside=0; iside<NSIDE+1; iside++){
    for (int ihar=0; ihar<NHAR; ihar++){
      for(int icent=0;icent<NCENT+1;icent++){
        hMu_Vn_FCal_pt[icent][iside][ihar]->Write();
        hMu_Vn_Zdc_eta[icent][iside][ihar]->Write();
      }
      for(int ipt=0;ipt<NPT_Mu;ipt++){
        hMu_Vn_FCal_cent[ipt][iside][ihar]->Write();
      }
      for(int ieta=0;ieta<NETA_Mu;ieta++){
        hMu_Vn_Zdc_cent[ieta][iside][ihar]->Write();
      }
    }
  }
  //Muon Track hists
  //Raw hists Before selection
  hMu_n_raw->Write();
  hMu_pt_raw->Write();
  hMu_phi_raw->Write();
  hMu_eta_raw->Write();
  hMu_eloss_all_raw->Write();
  hMu_ms_phi_raw->Write();
  hMu_ms_theta_raw->Write();
  hMu_ms_qoverp_raw->Write();
  //After selection cuts
  hMu_n->Write();
  hMu_pt->Write();
  hMu_phi->Write();
  hMu_eta->Write();
  hMu_eloss_all->Write();
  hMu_ms_phi->Write();
  hMu_ms_theta->Write();
  hMu_ms_qoverp->Write();

  fOut_mu->Close();
}

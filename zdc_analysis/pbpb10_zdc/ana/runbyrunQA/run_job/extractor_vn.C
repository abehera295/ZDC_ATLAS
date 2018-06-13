void extractor_flat::ReadHistos_vn() {
  /*
    TFile *fIn_vn=new TFile("./Output_zdc_psi.root","read");
    for (int iside=0; iside<NSIDE+1; iside++) {
    for (int ihar=0; ihar<NHAR; ihar++) {
    sprintf(name, "hProf_sumx_side%d_har%d", iside, ihar);
    h_prof_sumx[iside][ihar] = (TProfile*)fIn_vn->Get(name);
    sprintf(name, "hProf_sumy_side%d_har%d", iside, ihar);
    h_prof_sumy[iside][ihar] = (TProfile*)fIn_vn->Get(name);
    }
    }
    for (int icent=0; icent<NCENT; icent++) {
    for (int iside=0; iside<NSIDE+1; iside++) {
    for (int ihar=0; ihar<NHAR; ihar++) {
    sprintf(name, "hProf_flatcos_cent%d_side%d_har%d", icent, iside, ihar);
    h_prof_flatCos[icent][iside][ihar] = (TProfile*)fIn_vn->Get(name);
    sprintf(name, "hProf_flatsin_cent%d_side%d_har%d", icent, iside, ihar);
    h_prof_flatSin[icent][iside][ihar] = (TProfile*)fIn_vn->Get(name);
    }
    }
    }

    for (int icent=0; icent<NCENT; icent++) {
    for (int ihar=0; ihar<NHAR; ihar++) {
    sprintf(name, "h_Dphi_fg_cent%d_har%d", icent, ihar);
    hDphi_fg[icent][ihar] = (TH1D*)fIn_vn->Get(name);
    sprintf(name, "h_Dphi_bg_cent%d_har%d", icent, ihar);
    hDphi_bg[icent][ihar] = (TH1D*)fIn_vn->Get(name);

    sprintf(name, "h_Res_flat_cent%d_har%d", icent, ihar);
    hRes_flat[icent][ihar] = (TH1D*)fIn_vn->Get(name);
    }
    }
  */
}

void extractor_flat::Run_vn() {

  int evts_min=NEv_id*NEv;
  int evts_max=(NEv_id+1)*NEv;

  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;

  /*
    for (int icent=0; icent<NCENT_1p; icent++) {
    for (int iside=0; iside<NSIDE; iside++) {
    for (int ihar=0; ihar<NHAR; ihar++) {
    (Psi_FCal_store[icent][iside][ihar]).assign(0,0);
    (Psi_Zdc_store[icent][iside][ihar]).assign(0,0);
    }
    }
    }
  */
  cout<<"Run events from "<<evts_min<<" to "<<evts_max<<endl;

  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);

    if (nb <=0)        {cout<<"Nbytes became <=0"<<endl;break;}
    if(jentry<evts_min || jentry>=evts_max) continue;
    if(jentry%5000==0) {cout<<"running event "<<jentry<<endl ;}
    //if (jentry > 80000) break;
    //-------------------------------------------------------------------------------------------------
    centb = get_centb(Fcal_Et);
    vzb=get_vzb(vx_z);
    if (centb < 0) continue;
    if (vzb < 0) continue;
    nevents++;

    FillGlobal_vn();
    
    //float Psi_FCal[NSIDE+1][NHAR];
    for (int iside=0; iside<NSIDE+1; iside++) {
      for (int ihar=0; ihar<NHAR; ihar++) {
	int ieta;
	if(iside==0) Psi_FCal[iside][ihar] = atan2(Qy[2][ihar],Qx[2][ihar])/double(ihar+2);
	if(iside==1) Psi_FCal[iside][ihar] = atan2(Qy[3][ihar],Qx[3][ihar])/double(ihar+2);
	if(iside==2) Psi_FCal[iside][ihar] = atan2(Qy[2][ihar]+Qy[3][ihar],Qx[2][ihar]+Qx[3][ihar])/double(ihar+2);
	hPsi_FCal[centb][iside][ihar]->Fill(Psi_FCal[iside][ihar]);
      }
    }
    //float Psi_Zdc[NSIDE+1][NHAR];
    for (int iside=0; iside<NSIDE+1; iside++) {
      for (int ihar=0; ihar<NHAR; ihar++) {
	//int side_sign=1;
	//if(iside==1) side_sign=-1;
	//Psi_Zdc[iside][ihar] = atan2(Qy_zdc[iside][ihar],side_sign*Qx_zdc[iside][ihar])/double(ihar+1);
	Psi_Zdc[iside][ihar] = atan2(Qy_zdc[iside][ihar],Qx_zdc[iside][ihar])/double(ihar+1);
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
      hRes_flat_FCal[centb][ihar]->Fill(cos( (ihar+2)*(Psi_FCal[0][ihar] - Psi_FCal[1][ihar]) ));
      hRes_flat_Zdc[centb][ihar]->Fill(cos( (ihar+1)*(Psi_Zdc[0][0] - Psi_Zdc[1][0]) ));
      dPsi_FCal = (ihar+2)*(Psi_FCal[1][ihar]-Psi_FCal[0][ihar] );
      hResCent_FCal_fg[ihar]->Fill(centb,cos(dPsi_FCal));
      int nn;
      //Wrap to -2pi to 2pi
      nn=fabs(dPsi_FCal)/(2*PI);
      if(dPsi_FCal < 0.) dPsi_FCal += nn*2*PI;
      if(dPsi_FCal > 0.) dPsi_FCal -= nn*2*PI;
      //Wrap to -pi to pi
      if (dPsi_FCal < -PI) dPsi_FCal += 2*PI;
      if (dPsi_FCal >  PI) dPsi_FCal -= 2*PI;
      /*
      for(int ihar2=0;ihar2<NHAR;ihar2++){
	if(dPsi_FCal < -(2.*ihar2+1)*PI && dPsi_FCal >= -(2.*ihar2+3)*PI) dPsi_FCal+=2.*(ihar2+1)*PI; 
	if(dPsi_FCal > (2.*ihar2+1)*PI && dPsi_FCal <= (2.*ihar2+3)*PI) dPsi_FCal-=2.*(ihar2+1)*PI; 
      }
      */
      dPsi_Zdc = (ihar+1)*(Psi_Zdc[1][0]-Psi_Zdc[0][0] );
      hResCent_Zdc_fg[ihar]->Fill(centb,cos(dPsi_Zdc));
      //Wrap to -2pi to 2pi
      nn=fabs(dPsi_Zdc)/(2*PI);
      if(dPsi_Zdc < 0.) dPsi_Zdc += nn*2*PI;
      if(dPsi_Zdc > 0.) dPsi_Zdc -= nn*2*PI;
      //Wrap to -pi to pi
      if (dPsi_Zdc < -PI) dPsi_Zdc += 2*PI;
      if (dPsi_Zdc >  PI) dPsi_Zdc -= 2*PI;
      /*
      for(int ihar2=0;ihar2<NHAR;ihar2++){
	if(dPsi_Zdc < -(2.*ihar2+1)*PI && dPsi_Zdc >= -(2.*ihar2+3)*PI) dPsi_Zdc+=2.*(ihar2+1)*PI; 
	if(dPsi_Zdc > (2.*ihar2+1)*PI && dPsi_Zdc <= (2.*ihar2+3)*PI) dPsi_Zdc-=2.*(ihar2+1)*PI; 
      }
      */
      hDphi_FCal_fg[centb][ihar]->Fill(dPsi_FCal);
      hDphi_Zdc_fg[centb][ihar]->Fill(dPsi_Zdc);

      //Background calculation 
      Pool_size =Pool_EP[centb][vzb].size();
      if(Pool_size>0){
	for (int isize=0; isize<Pool_size; isize++) {
	  dPsi_FCal = (ihar+2)*( Psi_FCal[1][ihar] - (Pool_EP[centb][vzb].at(isize)->Psi_FCal[0][ihar]) );
	  hResCent_FCal_bg[ihar]->Fill(centb,cos(dPsi_FCal));
	  //Wrap to -2pi to 2pi
	  nn=fabs(dPsi_FCal)/(2*PI);
	  if(dPsi_FCal < 0.) dPsi_FCal += nn*2*PI;
	  if(dPsi_FCal > 0.) dPsi_FCal -= nn*2*PI;
	  //Wrap to -pi to pi
	  if (dPsi_FCal < -PI) dPsi_FCal += 2*PI;
	  if (dPsi_FCal >  PI) dPsi_FCal -= 2*PI;
	  /*
	  for(int ihar2=0;ihar2<NHAR;ihar2++){
	    if(dPsi_FCal < -(2.*ihar2+1)*PI && dPsi_FCal >= -(2.*ihar2+3)*PI) dPsi_FCal+=2.*(ihar2+1)*PI;
	    if(dPsi_FCal > (2.*ihar2+1)*PI && dPsi_FCal <= (2.*ihar2+3)*PI) dPsi_FCal-=2.*(ihar2+1)*PI;
	  }
	  */
	  hDphi_FCal_bg[centb][ihar]->Fill(dPsi_FCal);
	  
	  dPsi_FCal = (ihar+2)*( Psi_FCal[0][ihar] - (Pool_EP[centb][vzb].at(isize)->Psi_FCal[1][ihar]) );
	  hResCent_FCal_bg[ihar]->Fill(centb,cos(dPsi_FCal));
	  //Wrap to -2pi to 2pi
          nn=fabs(dPsi_FCal)/(2*PI);
          if(dPsi_FCal < 0.) dPsi_FCal += nn*2*PI;
          if(dPsi_FCal > 0.) dPsi_FCal -= nn*2*PI;
          //Wrap to -pi to pi
          if (dPsi_FCal < -PI) dPsi_FCal += 2*PI;
          if (dPsi_FCal >  PI) dPsi_FCal -= 2*PI;
	  /*
	  for(int ihar2=0;ihar2<NHAR;ihar2++){
	    if(dPsi_FCal < -(2.*ihar2+1)*PI && dPsi_FCal >= -(2.*ihar2+3)*PI) dPsi_FCal+=2.*(ihar2+1)*PI;
	    if(dPsi_FCal > (2.*ihar2+1)*PI && dPsi_FCal <= (2.*ihar2+3)*PI) dPsi_FCal-=2.*(ihar2+1)*PI;
	  }
	  */
	  hDphi_FCal_bg[centb][ihar]->Fill(-dPsi_FCal);

	  dPsi_Zdc = (ihar+1)*( Psi_Zdc[1][0] - (Pool_EP[centb][vzb].at(isize)->Psi_Zdc[0][0]) );
	  hResCent_Zdc_bg[ihar]->Fill(centb,cos(dPsi_Zdc));
	  //Wrap to -2pi to 2pi
	  nn=fabs(dPsi_Zdc)/(2*PI);
	  if(dPsi_Zdc < 0.) dPsi_Zdc += nn*2*PI;
	  if(dPsi_Zdc > 0.) dPsi_Zdc -= nn*2*PI;
	  //Wrap to -pi to pi
	  if (dPsi_Zdc < -PI) dPsi_Zdc += 2*PI;
	  if (dPsi_Zdc >  PI) dPsi_Zdc -= 2*PI;
	  /*
	  for(int ihar2=0;ihar2<NHAR;ihar2++){
	    if(dPsi_Zdc < -(2.*ihar2+1)*PI && dPsi_Zdc >= -(2.*ihar2+3)*PI) dPsi_Zdc+=2.*(ihar2+1)*PI;
	    if(dPsi_Zdc > (2.*ihar2+1)*PI && dPsi_Zdc <= (2.*ihar2+3)*PI) dPsi_Zdc-=2.*(ihar2+1)*PI;
	  }
	  */
	  hDphi_Zdc_bg[centb][ihar]->Fill(dPsi_Zdc);

	  dPsi_Zdc = (ihar+1)*( Psi_Zdc[0][0] - (Pool_EP[centb][vzb].at(isize)->Psi_Zdc[1][0]) );
	  hResCent_Zdc_bg[ihar]->Fill(centb,cos(dPsi_Zdc));
	  //Wrap to -2pi to 2pi
	  nn=fabs(dPsi_Zdc)/(2*PI);
	  if(dPsi_Zdc < 0.) dPsi_Zdc += nn*2*PI;
	  if(dPsi_Zdc > 0.) dPsi_Zdc -= nn*2*PI;
	  //Wrap to -pi to pi
	  if (dPsi_Zdc < -PI) dPsi_Zdc += 2*PI;
	  if (dPsi_Zdc >  PI) dPsi_Zdc -= 2*PI;
	  /*
	  for(int ihar2=0;ihar2<NHAR;ihar2++){
	    if(dPsi_Zdc < -(2.*ihar2+1)*PI && dPsi_Zdc >= -(2.*ihar2+3)*PI) dPsi_Zdc+=2.*(ihar2+1)*PI;
	    if(dPsi_Zdc > (2.*ihar2+1)*PI && dPsi_Zdc <= (2.*ihar2+3)*PI) dPsi_Zdc-=2.*(ihar2+1)*PI;
	  }
	  */
	  hDphi_Zdc_bg[centb][ihar]->Fill(-dPsi_Zdc);

	}
      }


      /*
	//Old code
      //Psi Store creation 
      int str_size0 = (Psi_FCal_store[centb_1p][0][ihar]).size();
      int str_size1 = (Psi_FCal_store[centb_1p][1][ihar]).size();
      if(str_size0>0){
      for (int isiz=0; isiz<str_size0; isiz++) {
      //dPsi = (ihar+1)*( Psi_FCal.at(0).at(ihar) - (Psi_store[centb_1p][1][ihar]).at(isiz) );
      dPsi_FCal = (ihar+2)*( Psi_FCal[1][ihar] - (Psi_FCal_store[centb_1p][0][ihar]).at(isiz) );
      if (dPsi_FCal < -PI) dPsi_FCal += 2*PI;
      if (dPsi_FCal >  PI) dPsi_FCal -= 2*PI;
      hDphi_FCal_bg[centb][ihar]->Fill(dPsi_FCal);
      }
      }
      if(str_size1>0){
      for (int isiz=0; isiz<str_size1; isiz++) {
      dPsi_FCal = (ihar+2)*( Psi_FCal[0][ihar] - (Psi_FCal_store[centb_1p][1][ihar]).at(isiz) );
      if (dPsi_FCal < -PI) dPsi_FCal += 2*PI;
      if (dPsi_FCal >  PI) dPsi_FCal -= 2*PI;
      hDphi_FCal_bg[centb][ihar]->Fill(-dPsi_FCal);
      }
      }
      int str_size2 = (Psi_Zdc_store[centb_1p][0][ihar]).size();
      int str_size3 = (Psi_Zdc_store[centb_1p][1][ihar]).size();
      if(str_size2>0){
      for (int isiz=0; isiz<str_size2; isiz++) {
      //dPsi = (ihar+1)*( Psi_Zdc.at(0).at(ihar) - (Psi_store[centb_1p][1][ihar]).at(isiz) );
      dPsi_Zdc = (ihar+1)*( Psi_Zdc[1][ihar] - (Psi_Zdc_store[centb_1p][0][ihar]).at(isiz) );
      if (dPsi_Zdc < -PI) dPsi_Zdc += 2*PI;
      if (dPsi_Zdc >  PI) dPsi_Zdc -= 2*PI;
      hDphi_Zdc_bg[centb][ihar]->Fill(dPsi_Zdc);
      }
      }
      if(str_size3>0){
      for (int isiz=0; isiz<str_size3; isiz++) {
      dPsi_Zdc = (ihar+1)*( Psi_Zdc[0][ihar] - (Psi_Zdc_store[centb_1p][1][ihar]).at(isiz) );
      if (dPsi_Zdc < -PI) dPsi_Zdc += 2*PI;
      if (dPsi_Zdc >  PI) dPsi_Zdc -= 2*PI;
      hDphi_Zdc_bg[centb][ihar]->Fill(-dPsi_Zdc);
      }
      }

      */
    }
    
    FillTrack_vn(Psi_FCal,Psi_Zdc);
    
    /*
    //Old code
    //Fill store for Psi background
    for (int iside=0; iside<NSIDE+1; iside++) {
    for (int ihar=0; ihar<NHAR; ihar++) {
    (Psi_FCal_store[centb_1p][iside][ihar]).push_back(Psi_FCal[iside][ihar]);
    (Psi_Zdc_store[centb_1p][iside][ihar]).push_back(Psi_Zdc[iside][ihar]);
    if ((Psi_FCal_store[centb_1p][iside][ihar]).size() >= STORE_DEP) (Psi_FCal_store[centb_1p][iside][ihar]).erase((Psi_FCal_store[centb_1p][iside][ihar]).begin());
    if ((Psi_Zdc_store[centb_1p][iside][ihar]).size() >= STORE_DEP) (Psi_Zdc_store[centb_1p][iside][ihar]).erase((Psi_Zdc_store[centb_1p][iside][ihar]).begin());
    }
    }
    */
    
    
    //Mixed Event Pool
    Pool_size =Pool_EP[centb][vzb].size();
    if(Pool_size>=STORE_DEP){
      delete Pool_EP[centb][vzb].at(0);
      Pool_EP[centb][vzb].erase( (Pool_EP[centb][vzb]).begin() );
    }
    Pool_EP[centb][vzb].push_back( ev );
    //delete ev;
    
  }//Event Loop ends
  cout<<"Finished running all events."<<endl;
}

void extractor_flat::FillGlobal_vn() {

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
    hEt_eTot         [iside]->Fill(Fcal_Et*1000, eTot[iside][0]);
    hEt_ePix         [iside]->Fill(Fcal_Et*1000, ePix[iside][0]);
    hNtrk_eTot       [iside]->Fill(ntrkQ, eTot[iside][0]);
    hNtrk_ePix       [iside]->Fill(ntrkQ, ePix[iside][0]);
  }

  hEt  ->Fill(Fcal_Et*1000);
  hNtrk->Fill(ntrkQ);
  hCent->Fill(centb);
  hEt_ntrk->Fill(Fcal_Et*1000, ntrkQ);
  //-------------------------------------------------------------------------------------------------
}


void extractor_flat::FillTrack_vn(float Psi_FCal[3][4],float Psi_Zdc[3][4]) {
  ntrkQ = 0;
  for (int icent=0; icent<NCENT; icent++) {
    for (int iside=0; iside<NSIDE+1; iside++) {
      for (int ihar=0; ihar<NHAR; ihar++) {
	hVn_FCal_tmp_fg[icent][iside][ihar]->Reset();
	hVn_Zdc_tmp_fg[icent][iside][ihar]->Reset();
	/*
	hVnS_FCal_tmp_fg[icent][iside][ihar]->Reset();
	hVnS_Zdc_tmp_fg[icent][iside][ihar]->Reset();
	*/
	hVnW_tmp_fg[icent][iside][ihar]->Reset();
	for(int isize=0;isize<5;isize++){
	  hVn_FCal_tmp_bg[icent][iside][ihar][isize]->Reset();
	  hVn_Zdc_tmp_bg[icent][iside][ihar][isize]->Reset();
	  /*
	  hVnS_FCal_tmp_bg[icent][iside][ihar][isize]->Reset();
	  hVnS_Zdc_tmp_bg[icent][iside][ihar][isize]->Reset();
	  */
	  hVnW_tmp_bg[icent][iside][ihar][isize]->Reset();
	}
      }
    }
  }

  //Run track loop
  for(int itrk=0; itrk<trk_n; itrk++){
    //if(trk_Quality2[itrk] < 3) continue;
  
      
    //new data format
    UInt_t trkb = trk_data[itrk];
    double btmp=0;
    trk_phi0_wrt_PV[itrk]  = trkb&PHIM;// btmp  = trk_phi0_wrt_PV[itrk];
    trk_phi0_wrt_PV[itrk]+=0.5;   trk_phi0_wrt_PV[itrk] *= PIslice;     trkb >>= PHIS;
    //trk_phibin[itrk] = btmp/(512)*NPHI; trk_phibin[itrk]++;

    trk_eta[itrk]  = trkb&ETAM; btmp  = trk_eta[itrk];
    trk_eta[itrk]+=0.5;   trk_eta[itrk] *= ETAslice;   trk_eta[itrk] -=ETAoff;  trkb >>= ETAS;
    //trk_etabin[itrk] = btmp/(500.0)*NETA; trk_etabin[itrk]++;

    trk_charge[itrk]= trkb&CHM;       trkb >>= CHS;
    trk_pt[itrk]=get_pt(trkb&PTM);
    //trk_ptbin[itrk] = getptbin(trkb&PTM);
    trkb >>= PTS;
    //trk_ptbin[itrk] = get_ptbin2(trk_pt[itrk]);
    trk_Quality2[itrk]  = trkb&QUALM;
    //trk_w[itrk]  = trkb&QUALM;

    //double ptgev = trk_pt[itrk]/1000.0;
    double ptgev=trk_pt[itrk];
    int ptBin  = get_ptBin(ptgev);       if (ptBin  < 0) continue;
    int etaBin = get_etaBin(trk_eta[itrk]); if (etaBin < 0) continue;
    int sign_eta = 1;
    if (trk_eta[itrk] < 0) sign_eta = -1;

    trk_wei[itrk] = 1.;
    trk_eff[itrk]=1.;
    trk_wei[itrk]= 1./(fabs(trk_w[itrk])/256.0 * 1.4+0.3);//0.3-1.7;
    //trk_wei[itrk] = 1.;
    trk_eff[itrk]=1./effTool->detTrkEff(Centrality,trk_eta[itrk],trk_pt[itrk]);//Trk Efficiency from Tool.cxx
    trk_phi[itrk]=trk_phi0_wrt_PV[itrk];
    hEta->Fill(trk_eta[itrk]);
    hPhi->Fill(trk_phi[itrk]);
    hPhi_wei->Fill(trk_phi[itrk],trk_wei[itrk]);
    hPt->Fill(ptgev);
    hTrk_w->Fill(trk_w[itrk]);
    hTrk_wei->Fill(trk_wei[itrk]);
    hTrk_eff2[centb][ptBin]->Fill(trk_eta[itrk],trk_eff[itrk]);
    hTrk_eff->Fill(trk_eff[itrk]);
    //hEtaPhi->Fill(trk_eta[itrk],trk_phi[itrk]);
    //hEtaPhi_wei->Fill(trk_eta[itrk],trk_phi[itrk],trk_eff[itrk]);

    //Apply track efficiency
    trk_wei[itrk]/=trk_eff[itrk];
    
    double corr1_FCal_fg,corr2_FCal_fg,corr1_Zdc_fg,corr2_Zdc_fg;
    double corr1_FCal_bg,corr2_FCal_bg,corr1_Zdc_bg,corr2_Zdc_bg;
    for (int iside=0; iside<NSIDE+1; iside++) {
      for (int ihar=0; ihar<NHAR; ihar++) {
	//Fill foreground
	corr1_FCal_fg = cos( (ihar+2)*(trk_phi[itrk]-Psi_FCal[iside][ihar]) )*trk_wei[itrk];
	corr2_FCal_fg = sin( (ihar+2)*(trk_phi[itrk]-Psi_FCal[iside][ihar]) )*trk_wei[itrk];
	corr1_Zdc_fg = cos( (ihar+1)*(trk_phi[itrk]-Psi_Zdc[iside][0]) )*trk_wei[itrk];
        corr2_Zdc_fg = sin( (ihar+1)*(trk_phi[itrk]-Psi_Zdc[iside][0]) )*trk_wei[itrk];
	hVn_FCal_tmp_fg[centb][iside][ihar]->Fill(ptBin,etaBin,corr1_FCal_fg);
	hVn_FCal_tmp_fg[centb][iside][ihar]->Fill(NPT,etaBin,corr1_FCal_fg);
	hVn_FCal_tmp_fg[centb][iside][ihar]->Fill(ptBin,NETA,corr1_FCal_fg);
	hVn_FCal_tmp_fg[centb][iside][ihar]->Fill(ptBin,NETA+1,corr1_FCal_fg*sign_eta);
	hVn_FCal_tmp_fg[centb][iside][ihar]->Fill(NPT,NETA,corr1_FCal_fg);
	hVn_FCal_tmp_fg[centb][iside][ihar]->Fill(NPT,NETA+1,corr1_FCal_fg*sign_eta);
	//hVn_FCal_tmp_fg[centb][iside][ihar]->Fill(NPT+1,etaBin,corr1_FCal_fg*sign_eta);
	if(fabs(trk_eta[itrk])<0.8){
	  hVn_FCal_tmp_fg[centb][iside][ihar]->Fill(ptBin,NETA+2,corr1_FCal_fg);
	  hVn_FCal_tmp_fg[centb][iside][ihar]->Fill(ptBin,NETA+3,corr1_FCal_fg*sign_eta);
	  hVn_FCal_tmp_fg[centb][iside][ihar]->Fill(NPT,NETA+2,corr1_FCal_fg);
	  hVn_FCal_tmp_fg[centb][iside][ihar]->Fill(NPT,NETA+3,corr1_FCal_fg*sign_eta);
	}
	/*
	hVnS_FCal_tmp_fg[centb][iside][ihar]->Fill(ptBin,etaBin,corr2_FCal_fg);
        hVnS_FCal_tmp_fg[centb][iside][ihar]->Fill(NPT,etaBin,corr2_FCal_fg);
	hVnS_FCal_tmp_fg[centb][iside][ihar]->Fill(ptBin,NETA,corr2_FCal_fg);
	hVnS_FCal_tmp_fg[centb][iside][ihar]->Fill(ptBin,NETA+1,corr2_FCal_fg*sign_eta);
        hVnS_FCal_tmp_fg[centb][iside][ihar]->Fill(NPT,NETA,corr2_FCal_fg);
        hVnS_FCal_tmp_fg[centb][iside][ihar]->Fill(NPT,NETA+1,corr2_FCal_fg*sign_eta);
	//hVnS_FCal_tmp_fg[centb][iside][ihar]->Fill(NPT+1,etaBin,corr2_FCal_fg*sign_eta);
        if(fabs(trk_eta[itrk])<0.8){
          hVnS_FCal_tmp_fg[centb][iside][ihar]->Fill(ptBin,NETA+2,corr2_FCal_fg);
          hVnS_FCal_tmp_fg[centb][iside][ihar]->Fill(ptBin,NETA+3,corr2_FCal_fg*sign_eta);
          hVnS_FCal_tmp_fg[centb][iside][ihar]->Fill(NPT,NETA+2,corr2_FCal_fg);
          hVnS_FCal_tmp_fg[centb][iside][ihar]->Fill(NPT,NETA+3,corr2_FCal_fg*sign_eta);
        }
	*/

	hVn_Zdc_tmp_fg[centb][iside][ihar]->Fill(ptBin,etaBin,corr1_Zdc_fg);
        hVn_Zdc_tmp_fg[centb][iside][ihar]->Fill(NPT,etaBin,corr1_Zdc_fg);
        hVn_Zdc_tmp_fg[centb][iside][ihar]->Fill(ptBin,NETA,corr1_Zdc_fg);
	hVn_Zdc_tmp_fg[centb][iside][ihar]->Fill(ptBin,NETA+1,corr1_Zdc_fg*sign_eta);
        hVn_Zdc_tmp_fg[centb][iside][ihar]->Fill(NPT,NETA,corr1_Zdc_fg);
        hVn_Zdc_tmp_fg[centb][iside][ihar]->Fill(NPT,NETA+1,corr1_Zdc_fg*sign_eta);
	//hVn_Zdc_tmp_fg[centb][iside][ihar]->Fill(NPT+1,etaBin,corr1_Zdc_fg*sign_eta);
        if(fabs(trk_eta[itrk])<0.8){
          hVn_Zdc_tmp_fg[centb][iside][ihar]->Fill(ptBin,NETA+2,corr1_Zdc_fg);
          hVn_Zdc_tmp_fg[centb][iside][ihar]->Fill(ptBin,NETA+3,corr1_Zdc_fg*sign_eta);
          hVn_Zdc_tmp_fg[centb][iside][ihar]->Fill(NPT,NETA+2,corr1_Zdc_fg);
          hVn_Zdc_tmp_fg[centb][iside][ihar]->Fill(NPT,NETA+3,corr1_Zdc_fg*sign_eta);
        }
	/*
        hVnS_Zdc_tmp_fg[centb][iside][ihar]->Fill(ptBin,etaBin,corr2_Zdc_fg);
        hVnS_Zdc_tmp_fg[centb][iside][ihar]->Fill(NPT,etaBin,corr2_Zdc_fg);
        hVnS_Zdc_tmp_fg[centb][iside][ihar]->Fill(ptBin,NETA,corr2_Zdc_fg);
	hVnS_Zdc_tmp_fg[centb][iside][ihar]->Fill(ptBin,NETA+1,corr2_Zdc_fg*sign_eta);
        hVnS_Zdc_tmp_fg[centb][iside][ihar]->Fill(NPT,NETA,corr2_Zdc_fg);
        hVnS_Zdc_tmp_fg[centb][iside][ihar]->Fill(NPT,NETA+1,corr2_Zdc_fg*sign_eta);
	//hVnS_Zdc_tmp_fg[centb][iside][ihar]->Fill(NPT+1,etaBin,corr2_Zdc_fg*sign_eta);
        if(fabs(trk_eta[itrk])<0.8){
          hVnS_Zdc_tmp_fg[centb][iside][ihar]->Fill(ptBin,NETA+2,corr2_Zdc_fg);
          hVnS_Zdc_tmp_fg[centb][iside][ihar]->Fill(ptBin,NETA+3,corr2_Zdc_fg*sign_eta);
          hVnS_Zdc_tmp_fg[centb][iside][ihar]->Fill(NPT,NETA+2,corr2_Zdc_fg);
          hVnS_Zdc_tmp_fg[centb][iside][ihar]->Fill(NPT,NETA+3,corr2_Zdc_fg*sign_eta);
        }
	*/
	hVnW_tmp_fg[centb][iside][ihar]->Fill(ptBin,etaBin,trk_wei[itrk]);
        hVnW_tmp_fg[centb][iside][ihar]->Fill(NPT,etaBin,trk_wei[itrk]);
        hVnW_tmp_fg[centb][iside][ihar]->Fill(ptBin,NETA,trk_wei[itrk]);
        hVnW_tmp_fg[centb][iside][ihar]->Fill(ptBin,NETA+1,trk_wei[itrk]);
        hVnW_tmp_fg[centb][iside][ihar]->Fill(NPT,NETA,trk_wei[itrk]);
	hVnW_tmp_fg[centb][iside][ihar]->Fill(NPT,NETA+1,trk_wei[itrk]);
	//hVnW_tmp_fg[centb][iside][ihar]->Fill(NPT+1,etaBin,trk_wei[itrk]);
        if(fabs(trk_eta[itrk])<0.8){
	  hVnW_tmp_fg[centb][iside][ihar]->Fill(ptBin,NETA+2,trk_wei[itrk]);
	  hVnW_tmp_fg[centb][iside][ihar]->Fill(ptBin,NETA+3,trk_wei[itrk]);
	  hVnW_tmp_fg[centb][iside][ihar]->Fill(NPT,NETA+2,trk_wei[itrk]);
	  hVnW_tmp_fg[centb][iside][ihar]->Fill(NPT,NETA+3,trk_wei[itrk]);
	}


	
	//Fill background
	Pool_size =Pool_EP[centb][vzb].size();
	if(Pool_size>0){
	  for (int isize=0; isize<Pool_size; isize++) {
	    float Psi_FCal_temp=Pool_EP[centb][vzb].at(isize)->Psi_FCal[iside][ihar];
	    corr1_FCal_bg = cos( (ihar+2)*(trk_phi[itrk]-Psi_FCal_temp) )*trk_wei[itrk];
	    corr2_FCal_bg = sin( (ihar+2)*(trk_phi[itrk]-Psi_FCal_temp) )*trk_wei[itrk];
	    hVn_FCal_tmp_bg[centb][iside][ihar][isize]->Fill(ptBin,etaBin,corr1_FCal_bg);
	    hVn_FCal_tmp_bg[centb][iside][ihar][isize]->Fill(NPT,etaBin,corr1_FCal_bg);
	    hVn_FCal_tmp_bg[centb][iside][ihar][isize]->Fill(ptBin,NETA,corr1_FCal_bg);
	    hVn_FCal_tmp_bg[centb][iside][ihar][isize]->Fill(ptBin,NETA+1,corr1_FCal_bg*sign_eta);
	    hVn_FCal_tmp_bg[centb][iside][ihar][isize]->Fill(NPT,NETA,corr1_FCal_bg);
	    hVn_FCal_tmp_bg[centb][iside][ihar][isize]->Fill(NPT,NETA+1,corr1_FCal_bg*sign_eta);
	    //hVn_FCal_tmp_bg[centb][iside][ihar][isize]->Fill(NPT+1,etaBin,corr1_FCal_bg*sign_eta);
	    if(fabs(trk_eta[itrk])<0.8){
	      hVn_FCal_tmp_bg[centb][iside][ihar][isize]->Fill(ptBin,NETA+2,corr1_FCal_bg);
	      hVn_FCal_tmp_bg[centb][iside][ihar][isize]->Fill(ptBin,NETA+3,corr1_FCal_bg*sign_eta);
	      hVn_FCal_tmp_bg[centb][iside][ihar][isize]->Fill(NPT,NETA+2,corr1_FCal_bg);
	      hVn_FCal_tmp_bg[centb][iside][ihar][isize]->Fill(NPT,NETA+3,corr1_FCal_bg*sign_eta);
	    }
	    /*
	    hVnS_FCal_tmp_bg[centb][iside][ihar][isize]->Fill(ptBin,etaBin,corr2_FCal_bg);
	    hVnS_FCal_tmp_bg[centb][iside][ihar][isize]->Fill(NPT,etaBin,corr2_FCal_bg);
	    hVnS_FCal_tmp_bg[centb][iside][ihar][isize]->Fill(ptBin,NETA,corr2_FCal_bg);
	    hVnS_FCal_tmp_bg[centb][iside][ihar][isize]->Fill(ptBin,NETA+1,corr2_FCal_bg*sign_eta);
	    hVnS_FCal_tmp_bg[centb][iside][ihar][isize]->Fill(NPT,NETA,corr2_FCal_bg);
	    hVnS_FCal_tmp_bg[centb][iside][ihar][isize]->Fill(NPT,NETA+1,corr2_FCal_bg*sign_eta);
	    //hVnS_FCal_tmp_bg[centb][iside][ihar][isize]->Fill(NPT+1,etaBin,corr2_FCal_bg*sign_eta);
	    if(fabs(trk_eta[itrk])<0.8){
	      hVnS_FCal_tmp_bg[centb][iside][ihar][isize]->Fill(ptBin,NETA+2,corr2_FCal_bg);
	      hVnS_FCal_tmp_bg[centb][iside][ihar][isize]->Fill(ptBin,NETA+3,corr2_FCal_bg*sign_eta);
	      hVnS_FCal_tmp_bg[centb][iside][ihar][isize]->Fill(NPT,NETA+2,corr2_FCal_bg);
	      hVnS_FCal_tmp_bg[centb][iside][ihar][isize]->Fill(NPT,NETA+3,corr2_FCal_bg*sign_eta);
	    }
	    */

	    float Psi_Zdc_temp=Pool_EP[centb][vzb].at(isize)->Psi_Zdc[iside][0];
	    corr1_Zdc_bg = cos( (ihar+1)*(trk_phi[itrk]-Psi_Zdc_temp) )*trk_wei[itrk];
	    corr2_Zdc_bg = sin( (ihar+1)*(trk_phi[itrk]-Psi_Zdc_temp) )*trk_wei[itrk];
	    hVn_Zdc_tmp_bg[centb][iside][ihar][isize]->Fill(ptBin,etaBin,corr1_Zdc_bg);
	    hVn_Zdc_tmp_bg[centb][iside][ihar][isize]->Fill(NPT,etaBin,corr1_Zdc_bg);
	    hVn_Zdc_tmp_bg[centb][iside][ihar][isize]->Fill(ptBin,NETA,corr1_Zdc_bg);
	    hVn_Zdc_tmp_bg[centb][iside][ihar][isize]->Fill(ptBin,NETA+1,corr1_Zdc_bg*sign_eta);
	    hVn_Zdc_tmp_bg[centb][iside][ihar][isize]->Fill(NPT,NETA,corr1_Zdc_bg);
	    hVn_Zdc_tmp_bg[centb][iside][ihar][isize]->Fill(NPT,NETA+1,corr1_Zdc_bg*sign_eta);
	    //hVn_Zdc_tmp_bg[centb][iside][ihar][isize]->Fill(NPT+1,etaBin,corr1_Zdc_bg*sign_eta);
	    if(fabs(trk_eta[itrk])<0.8){
	      hVn_Zdc_tmp_bg[centb][iside][ihar][isize]->Fill(ptBin,NETA+2,corr1_Zdc_bg);
	      hVn_Zdc_tmp_bg[centb][iside][ihar][isize]->Fill(ptBin,NETA+3,corr1_Zdc_bg*sign_eta);
	      hVn_Zdc_tmp_bg[centb][iside][ihar][isize]->Fill(NPT,NETA+2,corr1_Zdc_bg);
	      hVn_Zdc_tmp_bg[centb][iside][ihar][isize]->Fill(NPT,NETA+3,corr1_Zdc_bg*sign_eta);
	    }
	    /*
	    hVnS_Zdc_tmp_bg[centb][iside][ihar][isize]->Fill(ptBin,etaBin,corr2_Zdc_bg);
	    hVnS_Zdc_tmp_bg[centb][iside][ihar][isize]->Fill(NPT,etaBin,corr2_Zdc_bg);
	    hVnS_Zdc_tmp_bg[centb][iside][ihar][isize]->Fill(ptBin,NETA,corr2_Zdc_bg);
	    hVnS_Zdc_tmp_bg[centb][iside][ihar][isize]->Fill(ptBin,NETA+1,corr2_Zdc_bg*sign_eta);
	    hVnS_Zdc_tmp_bg[centb][iside][ihar][isize]->Fill(NPT,NETA,corr2_Zdc_bg);
	    hVnS_Zdc_tmp_bg[centb][iside][ihar][isize]->Fill(NPT,NETA+1,corr2_Zdc_bg*sign_eta);
	    //hVnS_Zdc_tmp_bg[centb][iside][ihar][isize]->Fill(NPT+1,etaBin,corr2_Zdc_bg*sign_eta);
	    if(fabs(trk_eta[itrk])<0.8){
	      hVnS_Zdc_tmp_bg[centb][iside][ihar][isize]->Fill(ptBin,NETA+2,corr2_Zdc_bg);
	      hVnS_Zdc_tmp_bg[centb][iside][ihar][isize]->Fill(ptBin,NETA+3,corr2_Zdc_bg*sign_eta);
	      hVnS_Zdc_tmp_bg[centb][iside][ihar][isize]->Fill(NPT,NETA+2,corr2_Zdc_bg);
	      hVnS_Zdc_tmp_bg[centb][iside][ihar][isize]->Fill(NPT,NETA+3,corr2_Zdc_bg*sign_eta);
	    }
	    */
	
	    hVnW_tmp_bg[centb][iside][ihar][isize]->Fill(ptBin,etaBin,trk_wei[itrk]);
	    hVnW_tmp_bg[centb][iside][ihar][isize]->Fill(NPT,etaBin,trk_wei[itrk]);
	    hVnW_tmp_bg[centb][iside][ihar][isize]->Fill(ptBin,NETA,trk_wei[itrk]);
	    hVnW_tmp_bg[centb][iside][ihar][isize]->Fill(ptBin,NETA+1,trk_wei[itrk]);
	    hVnW_tmp_bg[centb][iside][ihar][isize]->Fill(NPT,NETA,trk_wei[itrk]);
	    hVnW_tmp_bg[centb][iside][ihar][isize]->Fill(NPT,NETA+1,trk_wei[itrk]);
	    //hVnW_tmp_bg[centb][iside][ihar][isize]->Fill(NPT+1,etaBin,trk_wei[itrk]);
	    if(fabs(trk_eta[itrk])<0.8){
	      hVnW_tmp_bg[centb][iside][ihar][isize]->Fill(ptBin,NETA+2,trk_wei[itrk]);
	      hVnW_tmp_bg[centb][iside][ihar][isize]->Fill(ptBin,NETA+3,trk_wei[itrk]);
	      hVnW_tmp_bg[centb][iside][ihar][isize]->Fill(NPT,NETA+2,trk_wei[itrk]);
	      hVnW_tmp_bg[centb][iside][ihar][isize]->Fill(NPT,NETA+3,trk_wei[itrk]);
	    }
	  }
	}
	
      }
    }

    ntrkQ++;
  }
  //Fill Flow histograms
  double Vn_FCal,VnS_FCal,Vn_Zdc,VnS_Zdc,Vnw;
  for (int iside=0; iside<NSIDE+1; iside++) {
    for (int ihar=0; ihar<NHAR; ihar++) {
      //Foreground
      for (int ipt=0; ipt<NPT+1;ipt++){
	for (int ieta=0; ieta<NETA+4;ieta++){
	  Vn_FCal=hVn_FCal_tmp_fg[centb][iside][ihar]->GetBinContent(ipt+1,ieta+1);
	  Vn_Zdc=hVn_Zdc_tmp_fg[centb][iside][ihar]->GetBinContent(ipt+1,ieta+1);
	  /*
	  VnS_FCal=hVnS_FCal_tmp_fg[centb][iside][ihar]->GetBinContent(ipt+1,ieta+1);
	  VnS_Zdc=hVnS_Zdc_tmp_fg[centb][iside][ihar]->GetBinContent(ipt+1,ieta+1);
	  */
	  Vnw=hVnW_tmp_fg[centb][iside][ihar]->GetBinContent(ipt+1,ieta+1);
	  if(Vnw!=0.){
	    Vn_FCal/=Vnw;
	    Vn_Zdc/=Vnw;
	    hVn_FCal_fg[centb][iside][ihar]->Fill(ipt,ieta,Vn_FCal,Vnw);
	    hVn_Zdc_fg[centb][iside][ihar]->Fill(ipt,ieta,Vn_Zdc,Vnw);
	    /*
	    VnS_FCal/=Vnw;
	    VnS_Zdc/=Vnw;
	    hVnS_FCal_fg[centb][iside][ihar]->Fill(ipt,ieta,VnS_FCal,Vnw);
	    hVnS_Zdc_fg[centb][iside][ihar]->Fill(ipt,ieta,VnS_Zdc,Vnw);
	    */
	  }
	}
      }
      //For Background
      Pool_size =Pool_EP[centb][vzb].size();
      for(int isize=0;isize<Pool_size;isize++){
	for (int ipt=0; ipt<NPT+1;ipt++){
	  for (int ieta=0; ieta<NETA+4;ieta++){
	    Vn_FCal=hVn_FCal_tmp_bg[centb][iside][ihar][isize]->GetBinContent(ipt+1,ieta+1);
	    Vn_Zdc=hVn_Zdc_tmp_bg[centb][iside][ihar][isize]->GetBinContent(ipt+1,ieta+1);
	    /*
	    VnS_FCal=hVnS_FCal_tmp_bg[centb][iside][ihar][isize]->GetBinContent(ipt+1,ieta+1);
	    VnS_Zdc=hVnS_Zdc_tmp_bg[centb][iside][ihar][isize]->GetBinContent(ipt+1,ieta+1);
	    */
	    Vnw=hVnW_tmp_bg[centb][iside][ihar][isize]->GetBinContent(ipt+1,ieta+1);
	    if(Vnw!=0.){
	      Vn_FCal/=Vnw;
	      Vn_Zdc/=Vnw;	      
	      hVn_FCal_bg[centb][iside][ihar]->Fill(ipt,ieta,Vn_FCal,Vnw);
	      hVn_Zdc_bg[centb][iside][ihar]->Fill(ipt,ieta,Vn_Zdc,Vnw);
	      /*
	      VnS_FCal/=Vnw;
	      VnS_Zdc/=Vnw;
	      hVnS_FCal_bg[centb][iside][ihar]->Fill(ipt,ieta,VnS_FCal,Vnw);
	      hVnS_Zdc_bg[centb][iside][ihar]->Fill(ipt,ieta,VnS_Zdc,Vnw);
	      */
	    }
	  }
	}
      }
    }//ihar
  }//iside
  hNTrk->Fill(ntrkQ);
  /*
  //Delete pointers to temp hists
  for (int icent=0; icent<NCENT; icent++) {
    for (int iside=0; iside<NSIDE+1; iside++) {
      for (int ihar=0; ihar<NHAR; ihar++) {
	delete hVn_FCal_tmp_fg[icent][iside][ihar];
	delete hVnS_FCal_tmp_fg[icent][iside][ihar];
	delete hVn_Zdc_tmp_fg[icent][iside][ihar];
	delete hVnS_Zdc_tmp_fg[icent][iside][ihar];
        delete  hVnW_tmp_fg[icent][iside][ihar];
	for(int isize=0;isize<5;isize++){
	  delete hVn_FCal_tmp_bg[icent][iside][ihar][isize];
 	  delete hVnS_FCal_tmp_bg[icent][iside][ihar][isize];
	  delete hVn_Zdc_tmp_bg[icent][iside][ihar][isize];
	  delete hVnS_Zdc_tmp_bg[icent][iside][ihar][isize];
	  delete hVnW_tmp_bg[icent][iside][ihar][isize];
	}
      }
    }
  }
  */

}


void extractor_flat::Init_vn() {
  //set_ptBins ();
  //set_etaBins();

  for (int iside=0; iside<NSIDE; iside++) {
    sprintf(name, "h_eTot_side%d", iside);
    h_eTot[iside] = new TH1D(name, "", 4000, -0.5, 3999.5); h_eTot[iside]->Sumw2();
    sprintf(name, "h_ePix_side%d", iside);
    h_ePix[iside] = new TH1D(name, "", 4000, -0.5, 3999.5); h_ePix[iside]->Sumw2();
   
    sprintf(name, "h_et_eTot_side%d", iside);
    hEt_eTot[iside] = new TH2F(name, "", 500, 0, 5000, 400, 0, 4000); hEt_eTot[iside]->Sumw2();
    sprintf(name, "h_et_ePix_side%d", iside);
    hEt_ePix[iside] = new TH2F(name, "", 500, 0, 5000, 200, 0, 2000); hEt_ePix[iside]->Sumw2();
      
    sprintf(name, "h_ntrk_eTot_side%d", iside);
    hNtrk_eTot[iside] = new TH2F(name, "", 400, 0, 4000, 400, 0, 4000); hNtrk_eTot[iside]->Sumw2();
    sprintf(name, "h_ntrk_ePix_side%d", iside);
    hNtrk_ePix[iside] = new TH2F(name, "", 400, 0, 4000, 200, 0, 2000); hNtrk_ePix[iside]->Sumw2();
  }
  
  hEt      = new TH1D("hEt", "", 5000, 0, 5000); hEt->Sumw2();
  hNtrk    = new TH1D("hNtrk", "", 5000, 0, 5000); hNtrk->Sumw2();
  hCent    = new TH1D("hCent", "", NCENT, -0.5, NCENT-0.5); hCent->Sumw2();
  hEt_ntrk = new TH2F("h_et_ntrk", "", 400, 0, 4000, 400, 0, 4000); hEt_ntrk->Sumw2();

  for (int iside=0; iside<NSIDE+1; iside++) {
    for (int ihar=0; ihar<NHAR; ihar++) {
      sprintf(name, "hProf_sumx_side%d_har%d", iside, ihar);
      h_prof_sumx[iside][ihar] = new TProfile(name," ",100,-0.5,99.5,-3000,3000,"s");
      sprintf(name, "hProf_sumy_side%d_har%d", iside, ihar);
      h_prof_sumy[iside][ihar] = new TProfile(name," ",100,-0.5,99.5,-3000,3000,"s");
      /*
      for (int icent=0; icent<NCENT; icent++) {
	sprintf(name, "h_Q2D_cent%d_side%d_har%d", icent, iside, ihar);
	h_Q2D[icent][iside][ihar] = new TH2F(name, "", 200, -1.0, 1.0, 100, -1.0, 1.0);
	h_Q2D[icent][iside][ihar]->Sumw2();
      }
      */
    }
  }

  for (int icent=0; icent<NCENT; icent++) {
    for (int iside=0; iside<NSIDE+1; iside++) {
      for (int ihar=0; ihar<NHAR; ihar++) {
	sprintf(name, "hProf_flatcos_cent%d_side%d_har%d", icent, iside, ihar);
	h_prof_flatCos[icent][iside][ihar] = new TProfile(name, "", NK, -0.5, NK-0.5, -1.1, 1.1);
	sprintf(name, "hProf_flatsin_cent%d_side%d_har%d", icent, iside, ihar);
	h_prof_flatSin[icent][iside][ihar] = new TProfile(name, "", NK, -0.5, NK-0.5, -1.1, 1.1);
      }
    }
  }

  for (int icent=0; icent<NCENT; icent++) {
    for (int iside=0; iside<NSIDE+1; iside++) {
      for (int ihar=0; ihar<NHAR; ihar++) {

	sprintf(name, "hPsi_cent%d_side%d_har%d", icent, iside, ihar);
	hPsi[icent][iside][ihar] = new TH1D(name, "", 100, -2*PI, 2*PI); 
	hPsi[icent][iside][ihar]->Sumw2();

	sprintf(name, "hPsiC_cent%d_side%d_har%d", icent, iside, ihar);
	hPsiC[icent][iside][ihar] = new TH1D(name, "", 100, -2*PI, 2*PI);
	hPsiC[icent][iside][ihar]->Sumw2();
            
	sprintf(name, "hPsi_FCal_cent%d_side%d_har%d", icent, iside, ihar);
	hPsi_FCal[icent][iside][ihar] = new TH1D(name, "", 100, -2*PI, 2*PI);
	hPsi_FCal[icent][iside][ihar]->Sumw2();
	sprintf(name, "hPsi_Zdc_cent%d_side%d_har%d", icent, iside, ihar);
	hPsi_Zdc[icent][iside][ihar] = new TH1D(name, "", 100, -2*PI, 2*PI);
	hPsi_Zdc[icent][iside][ihar]->Sumw2();

      }
    }
  }
   
  for (int icent=0; icent<NCENT; icent++) {
    for (int iside=0; iside<NSIDE+1; iside++) {
      for (int ihar=0; ihar<NHAR; ihar++) {

	sprintf(name, "h_Qx_cent%d_side%d_har%d", icent, iside, ihar);
	h_Qx[icent][iside][ihar] = new TH1D(name, "", 2000, -1000, 1000); h_Qx[icent][iside][ihar]->Sumw2();
	sprintf(name, "h_Qy_cent%d_side%d_har%d", icent, iside, ihar);
	h_Qy[icent][iside][ihar] = new TH1D(name, "", 2000, -1000, 1000); h_Qy[icent][iside][ihar]->Sumw2();
            
	sprintf(name, "h_QxC_cent%d_side%d_har%d", icent, iside, ihar);
	h_QxC[icent][iside][ihar] = new TH1D(name, "", 2000, -1000, 1000); h_QxC[icent][iside][ihar]->Sumw2();
	sprintf(name, "h_QyC_cent%d_side%d_har%d", icent, iside, ihar);
	h_QyC[icent][iside][ihar] = new TH1D(name, "", 2000, -1000, 1000); h_QyC[icent][iside][ihar]->Sumw2();
      }
    }
  }
  /*
  for (int icent=0; icent<NCENT; icent++) {
    for (int ihar=0; ihar<NHAR; ihar++) {
      sprintf(name, "hDphi_FCal_fg_cent%d_har%d", icent, ihar);
      hDphi_FCal_fg[icent][ihar] = new TH1D(name, "", 640, -3.2, 3.2); hDphi_FCal_fg[icent][ihar]->Sumw2();
      sprintf(name, "hDphi_FCal_bg_cent%d_har%d", icent, ihar);
      hDphi_FCal_bg[icent][ihar] = new TH1D(name, "", 640, -3.2, 3.2); hDphi_FCal_bg[icent][ihar]->Sumw2();
      sprintf(name, "hDphi_Zdc_fg_cent%d_har%d", icent, ihar);
      hDphi_Zdc_fg[icent][ihar] = new TH1D(name, "", 640, -3.2, 3.2); hDphi_Zdc_fg[icent][ihar]->Sumw2();
      sprintf(name, "hDphi_Zdc_bg_cent%d_har%d", icent, ihar);
      hDphi_Zdc_bg[icent][ihar] = new TH1D(name, "", 640, -3.2, 3.2); hDphi_Zdc_bg[icent][ihar]->Sumw2();
    }
  }
  */
  double x1=3.2;
  //int n1=640;
  int n1=64;
  int n2=2;
  for (int icent=0; icent<NCENT; icent++) {
    for (int ihar=0; ihar<NHAR; ihar++) {
      sprintf(name, "hDphi_FCal_fg_cent%d_har%d", icent, ihar);
      hDphi_FCal_fg[icent][ihar] = new TH1D(name, "", n2*n1,-n2*x1,n2*x1); hDphi_FCal_fg[icent][ihar]->Sumw2();
      sprintf(name, "hDphi_FCal_bg_cent%d_har%d", icent, ihar);
      hDphi_FCal_bg[icent][ihar] = new TH1D(name, "", n2*n1,-n2*x1,n2*x1); hDphi_FCal_bg[icent][ihar]->Sumw2();
      sprintf(name, "hDphi_Zdc_fg_cent%d_har%d", icent, ihar);
      hDphi_Zdc_fg[icent][ihar] = new TH1D(name, "", n2*n1,-n2*x1,n2*x1); hDphi_Zdc_fg[icent][ihar]->Sumw2();
      sprintf(name, "hDphi_Zdc_bg_cent%d_har%d", icent, ihar);
      hDphi_Zdc_bg[icent][ihar] = new TH1D(name, "", n2*n1,-n2*x1,n2*x1); hDphi_Zdc_bg[icent][ihar]->Sumw2();
    }
  }
  
  for (int icent=0; icent<NCENT; icent++) {
    for (int ihar=0; ihar<NHAR; ihar++) {

      sprintf(name, "hRes_cent%d_har%d", icent, ihar);
      hRes[icent][ihar] = new TH1D(name, "", 200, -1.0, 1.0);
      hRes[icent][ihar]->Sumw2();

      sprintf(name, "hRes_recn_cent%d_har%d", icent, ihar);
      hRes_recn[icent][ihar] = new TH1D(name, "", 200, -1.0, 1.0);
      hRes_recn[icent][ihar]->Sumw2();

      sprintf(name, "hRes_flat_FCal_cent%d_har%d", icent, ihar);
      hRes_flat_FCal[icent][ihar] = new TH1D(name, "", 200, -1.0, 1.0);
      hRes_flat_FCal[icent][ihar]->Sumw2();
      sprintf(name, "hRes_flat_Zdc_cent%d_har%d", icent, ihar);
      hRes_flat_Zdc[icent][ihar] = new TH1D(name, "", 200, -1.0, 1.0);
      hRes_flat_Zdc[icent][ihar]->Sumw2();
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

  for (int icent=0; icent<NCENT; icent++) {
    for (int iside=0; iside<NSIDE+1; iside++) {
      for (int ihar=0; ihar<NHAR; ihar++) {
	//Vn Cos
	sprintf(name,"hV%d_FCal_fg_cent%d_side%d",ihar+1,icent,iside);
	hVn_FCal_fg[icent][iside][ihar]  = new TProfile2D(name, "",NPT+1,0.-0.5,NPT+1-0.5,NETA+4,0.-0.5,NETA+4-0.5); hVn_FCal_fg[icent][iside][ihar]->Sumw2();
	sprintf(name,"hV%d_FCal_bg_cent%d_side%d",ihar+1,icent,iside);
	hVn_FCal_bg[icent][iside][ihar]  = new TProfile2D(name, "",NPT+1,0.-0.5,NPT+1-0.5,NETA+4,0.-0.5,NETA+4-0.5); hVn_FCal_bg[icent][iside][ihar]->Sumw2();
	sprintf(name,"hV%d_Zdc_fg_cent%d_side%d",ihar+1,icent,iside);
	hVn_Zdc_fg[icent][iside][ihar]  = new TProfile2D(name, "",NPT+1,0.-0.5,NPT+1-0.5,NETA+4,0.-0.5,NETA+4-0.5); hVn_Zdc_fg[icent][iside][ihar]->Sumw2();
	sprintf(name,"hV%d_Zdc_bg_cent%d_side%d",ihar+1,icent,iside);
	hVn_Zdc_bg[icent][iside][ihar]  = new TProfile2D(name, "",NPT+1,0.-0.5,NPT+1-0.5,NETA+4,0.-0.5,NETA+4-0.5); hVn_Zdc_bg[icent][iside][ihar]->Sumw2();
	/*
	//Vn Sin
	sprintf(name,"hV%dS_FCal_fg_cent%d_side%d",ihar+1,icent,iside);
	hVnS_FCal_fg[icent][iside][ihar]  = new TProfile2D(name, "",NPT+1,0.-0.5,NPT+1-0.5,NETA+4,0.-0.5,NETA+4-0.5); hVnS_FCal_fg[icent][iside][ihar]->Sumw2();
	sprintf(name,"hV%dS_FCal_bg_cent%d_side%d",ihar+1,icent,iside);
	hVnS_FCal_bg[icent][iside][ihar]  = new TProfile2D(name, "",NPT+1,0.-0.5,NPT+1-0.5,NETA+4,0.-0.5,NETA+4-0.5); hVnS_FCal_bg[icent][iside][ihar]->Sumw2();
	sprintf(name,"hV%dS_Zdc_fg_cent%d_side%d",ihar+1,icent,iside);
	hVnS_Zdc_fg[icent][iside][ihar]  = new TProfile2D(name, "",NPT+1,0.-0.5,NPT+1-0.5,NETA+4,0.-0.5,NETA+4-0.5); hVnS_Zdc_fg[icent][iside][ihar]->Sumw2();
	sprintf(name,"hV%dS_Zdc_bg_cent%d_side%d",ihar+1,icent,iside);
	hVnS_Zdc_bg[icent][iside][ihar]  = new TProfile2D(name, "",NPT+1,0.-0.5,NPT+1-0.5,NETA+4,0.-0.5,NETA+4-0.5); hVnS_Zdc_bg[icent][iside][ihar]->Sumw2();
	*/
      }
    }
  }

  //Temp hists 
  for (int icent=0; icent<NCENT; icent++) {
    for (int iside=0; iside<NSIDE+1; iside++) {
      for (int ihar=0; ihar<NHAR; ihar++) {
	//Vn Cos
	sprintf(name,"hV%d_FCal_tmp_fg_cent%d_side%d",ihar+1,icent,iside);
	hVn_FCal_tmp_fg[icent][iside][ihar]=new TH2F(name, "",NPT+1,0.-0.5,NPT+1-0.5,NETA+4,0.-0.5,NETA+4-0.5);
	hVn_FCal_tmp_fg[icent][iside][ihar]->Sumw2();
	sprintf(name,"hV%d_Zdc_tmp_fg_cent%d_side%d",ihar+1,icent,iside);
	hVn_Zdc_tmp_fg[icent][iside][ihar]=new TH2F(name, "",NPT+1,0.-0.5,NPT+1-0.5,NETA+4,0.-0.5,NETA+4-0.5);
	hVn_Zdc_tmp_fg[icent][iside][ihar]->Sumw2();
	/*
	//Vn Sin
	sprintf(name,"hV%dS_FCal_tmp_fg_cent%d_side%d",ihar+1,icent,iside);
	hVnS_FCal_tmp_fg[icent][iside][ihar]=new TH2F(name, "",NPT+1,0.-0.5,NPT+1-0.5,NETA+4,0.-0.5,NETA+4-0.5);
	hVnS_FCal_tmp_fg[icent][iside][ihar]->Sumw2();
	sprintf(name,"hV%dS_Zdc_tmp_fg_cent%d_side%d",ihar+1,icent,iside);
	hVnS_Zdc_tmp_fg[icent][iside][ihar]=new TH2F(name, "",NPT+1,0.-0.5,NPT+1-0.5,NETA+4,0.-0.5,NETA+4-0.5);
	hVnS_Zdc_tmp_fg[icent][iside][ihar]->Sumw2();
	*/
	sprintf(name,"hV%dW_tmp_fg_cent%d_side%d",ihar+1,icent,iside);
        hVnW_tmp_fg[icent][iside][ihar]=new TH2F(name, "",NPT+1,0.-0.5,NPT+1-0.5,NETA+4,0.-0.5,NETA+4-0.5);
        hVnW_tmp_fg[icent][iside][ihar]->Sumw2();
	for(int isize=0;isize<5;isize++){
	  //Vn Cos
	  sprintf(name,"hV%d_FCal_tmp_bg_cent%d_side%d_size%d",ihar+1,icent,iside,isize);
	  hVn_FCal_tmp_bg[icent][iside][ihar][isize]=new TH2F(name, "",NPT+1,0.-0.5,NPT+1-0.5,NETA+4,0.-0.5,NETA+4-0.5);
	  hVn_FCal_tmp_bg[icent][iside][ihar][isize]->Sumw2();
	  sprintf(name,"hV%d_Zdc_tmp_bg_cent%d_side%d_size%d",ihar+1,icent,iside,isize);
	  hVn_Zdc_tmp_bg[icent][iside][ihar][isize]=new TH2F(name, "",NPT+1,0.-0.5,NPT+1-0.5,NETA+4,0.-0.5,NETA+4-0.5);
	  hVn_Zdc_tmp_bg[icent][iside][ihar][isize]->Sumw2();
	  /*
	  //Vn Sin
	  sprintf(name,"hV%dS_FCal_tmp_bg_cent%d_side%d_size%d",ihar+1,icent,iside,isize);
	  hVnS_FCal_tmp_bg[icent][iside][ihar][isize]=new TH2F(name, "",NPT+1,0.-0.5,NPT+1-0.5,NETA+4,0.-0.5,NETA+4-0.5);
	  hVnS_FCal_tmp_bg[icent][iside][ihar][isize]->Sumw2();
	  sprintf(name,"hV%dS_Zdc_tmp_bg_cent%d_side%d_size%d",ihar+1,icent,iside,isize);
	  hVnS_Zdc_tmp_bg[icent][iside][ihar][isize]=new TH2F(name, "",NPT+1,0.-0.5,NPT+1-0.5,NETA+4,0.-0.5,NETA+4-0.5);
	  hVnS_Zdc_tmp_bg[icent][iside][ihar][isize]->Sumw2();
	  */
	  sprintf(name,"hV%dW_tmp_bg_cent%d_side%d_size%d",ihar+1,icent,iside,isize);
	  hVnW_tmp_bg[icent][iside][ihar][isize]=new TH2F(name, "",NPT+1,0.-0.5,NPT+1-0.5,NETA+4,0.-0.5,NETA+4-0.5);
	  hVnW_tmp_bg[icent][iside][ihar][isize]->Sumw2();
	}
      }
    }
  }

  /*
    for (int icent=0; icent<NCENT; icent++) {
    for (int iside=0; iside<NSIDE+1; iside++) {
    for (int ihar=0; ihar<NHAR; ihar++) {
    for (int ipt=0; ipt<NPT+1; ipt++) {
    std::vector<TH1D*> hTemp_vec1_fg(NETA+4, (TH1D*)0);
    std::vector<TH1D*> hTemp_vec2_fg(NETA+4, (TH1D*)0);
    std::vector<TH1D*> hTemp_vec1_bg(NETA+4, (TH1D*)0);
    std::vector<TH1D*> hTemp_vec2_bg(NETA+4, (TH1D*)0);
    for (int ieta=0; ieta<NETA+4; ieta++) {
    //if (ipt==NPT && ieta >= NETA) continue;
    sprintf(name, "hV%d_fg_cent%d_side%d_pt%d_eta%d", ihar+1, icent, iside, ipt, ieta);
    hTemp_vec1_fg.at(ieta) = new TH1D(name, "", 220, -1.1, 1.1); hTemp_vec1_fg.at(ieta)->Sumw2();
    sprintf(name, "hV%dS_fg_cent%d_side%d_pt%d_eta%d", ihar+1, icent, iside, ipt, ieta); 
    hTemp_vec2_fg.at(ieta) = new TH1D(name, "", 220, -1.1, 1.1); hTemp_vec2_fg.at(ieta)->Sumw2();
    sprintf(name, "hV%d_bg_cent%d_side%d_pt%d_eta%d", ihar+1, icent, iside, ipt, ieta);
    hTemp_vec1_bg.at(ieta) = new TH1D(name, "", 220, -1.1, 1.1); hTemp_vec1_bg.at(ieta)->Sumw2();
    sprintf(name, "hV%dS_bg_cent%d_side%d_pt%d_eta%d", ihar+1, icent, iside, ipt, ieta); 
    hTemp_vec2_bg.at(ieta) = new TH1D(name, "", 220, -1.1, 1.1); hTemp_vec2_bg.at(ieta)->Sumw2();
    }
    (hVn_fg[icent][iside][ihar]).push_back( hTemp_vec1_fg);
    (hVnS_fg[icent][iside][ihar]).push_back( hTemp_vec2_fg);
    (hVn_bg[icent][iside][ihar]).push_back( hTemp_vec1_bg);
    (hVnS_bg[icent][iside][ihar]).push_back( hTemp_vec2_bg);
    }
    }
    }
    }
  */
  //Track histograms
  hNTrk    = new TH1D("hNTrk", "", 5000, 0, 5000); hNtrk->Sumw2();
  hPt = new TH1D("hPt","",1000,0.,25.); hPt->Sumw2();
  hPhi = new TH1D("hPhi","",50,0.,7); hPhi->Sumw2();
  hPhi_wei = new TH1D("hPhi_wei","",50,0.,7.); hPhi_wei->Sumw2();
  hEta = new TH1D("hEta","",100,-2.5,2.5); hEta->Sumw2();
  hTrk_w = new TH1D("hTrk_w","",1000,-300.,300.); hTrk_w->Sumw2();
  hTrk_wei = new TH1D("hTrk_wei","",1000,-5.,5.); hTrk_wei->Sumw2();
  for(int icent=0;icent<NCENT;icent++){
    for(int ipt=0;ipt<NPT;ipt++){
      sprintf(name,"hTrk_eff2_cent%d_pt%d",icent,ipt);
      hTrk_eff2[icent][ipt]=new TProfile(name,"",50,-2.5,2.5);
      hTrk_eff2[icent][ipt]->Sumw2();
    }
  }
  hTrk_eff=new TH1D("hTrk_eff","",200,-2.,2.);
  //Correlation hists
  //hEtaPhi=new TH2F("hEtaPhi","",100,-2.5,2.5,50,0.,7.); hEtaPhi->Sumw2();
  //hEtaPhi_wei=new TH2F("hEtaPhi_wei","",100,-2.5,2.5,50,0.,7.); hEtaPhi_wei->Sumw2();

}


void extractor_flat::SaveHistos_vn() {
  fOut_vn->cd();

  for (int iside=0; iside<NSIDE; iside++) {
    h_eTot[iside]->Write();
    h_ePix[iside]->Write();
    hEt_eTot[iside]->Write();
    hEt_ePix[iside]->Write();
    hNtrk_eTot[iside]->Write();
    hNtrk_ePix[iside]->Write();
  }
  for (int icent=0; icent<NCENT; icent++) {
    for (int iside=0; iside<NSIDE+1; iside++) {
      for (int ihar=0; ihar<NHAR; ihar++) {
	//hPsi[icent][iside][ihar]->Write();
	//hPsiC[icent][iside][ihar]->Write();
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
  for (int icent=0; icent<NCENT; icent++) {
    for (int iside=0; iside<NSIDE+1; iside++) {
      for (int ihar=0; ihar<NHAR; ihar++) {
	//Vn Cos
	hVn_FCal_fg[icent][iside][ihar]->Write();
	hVn_FCal_bg[icent][iside][ihar]->Write();
	hVn_Zdc_fg[icent][iside][ihar]->Write();
	hVn_Zdc_bg[icent][iside][ihar]->Write();
	/*
	//Vn Sin
	hVnS_FCal_fg[icent][iside][ihar]->Write();
	hVnS_FCal_bg[icent][iside][ihar]->Write();
	hVnS_Zdc_fg[icent][iside][ihar]->Write();
	hVnS_Zdc_bg[icent][iside][ihar]->Write();
	*/
      }
    }
  }
  
  
  /*
    for (int icent=0; icent<NCENT; icent++) {
    for (int iside=0; iside<NSIDE+1; iside++) {
    for (int ihar=0; ihar<NHAR; ihar++) {
    for (int ipt=0; ipt<NPT+1; ipt++) {
    for (int ieta=0; ieta<NETA+4; ieta++) {
    //if (ipt==NPT && ieta >= NETA) continue;
    hVn_fg[icent][iside][ihar].at(ipt).at(ieta)->Write();
    hVnS_fg[icent][iside][ihar].at(ipt).at(ieta)->Write();
    hVn_bg[icent][iside][ihar].at(ipt).at(ieta)->Write();
    hVnS_bg[icent][iside][ihar].at(ipt).at(ieta)->Write();
    }
    }
    }
    }
    }
  */
  //Event histograms
  for(int icent=0;icent<NCENT;icent++){
    for(int ipt=0;ipt<NPT;ipt++){
      hTrk_eff2[icent][ipt]->Write();
    }
  }
  hEt->Write();
  hNtrk->Write();
  hCent->Write();
  hEt_ntrk->Write();
  //Track histograms
  hNTrk->Write();
  hPt->Write();
  hPhi->Write();
  hPhi_wei->Write();
  hEta->Write();
  hTrk_w->Write();
  hTrk_wei->Write();
  hTrk_eff->Write();
  //hEtaPhi->Write();
  //hEtaPhi_wei->Write();

  fOut_vn->Close();
}

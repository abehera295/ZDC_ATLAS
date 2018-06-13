void extractor_flat::ReadHistos_res() {
  sprintf(name,"../do_calib/EPCalib_hists/EP_%d_recenter.root",RNum);
  TFile *fin1=new TFile(name,"read");
  //For FCal Recenter
  sprintf(name,"hevts");
  hevts=(TH2F*)fin1->Get(name);
  for(int iz=0;iz<NZ;iz++){//FCal
    sprintf(name,"hEt_z%d",iz);
    hEt_Vz[iz]=(TProfile*)fin1->Get(name);
  }
  for(int id=0;id<NDET;id++){
    for(int iz=0;iz<NZ;iz++){
      sprintf(name,"hQw_d%d_z%d",id,iz);
      hQw[id][iz]=(TProfile*)fin1->Get(name);
      for(int ihar=0;ihar<NHAR;ihar++){//2,3,4,5
	sprintf(name,"hQx_d%d_z%d_h%d",id,iz,ihar);
	hQx[id][iz][ihar]=(TProfile*)fin1->Get(name);
	sprintf(name,"hQy_d%d_z%d_h%d",id,iz,ihar);
	hQy[id][iz][ihar]=(TProfile*)fin1->Get(name);
      }
    }
  }
  //For Zdc Recenter
  for (int iside=0; iside<NSIDE+1; iside++) {
    for (int ihar=0; ihar<NHAR; ihar++) {
      sprintf(name, "hProf_sumx_side%d_har%d", iside, ihar);
      h_prof_sumx[iside][ihar] = (TProfile*)fin1->Get(name);
      sprintf(name, "hProf_sumy_side%d_har%d", iside, ihar);
      h_prof_sumy[iside][ihar] = (TProfile*)fin1->Get(name);
    }
  }

  sprintf(name,"../do_calib/EPCalib_hists/EP_%d_flat.root",RNum);
  TFile *fin2=new TFile(name,"read");
  for (int icent=0; icent<NCENT_1p; icent++) {
    for (int iside=0; iside<NSIDE+1; iside++) {
      for (int ihar=0; ihar<NHAR; ihar++) {
	sprintf(name, "hProf_flatcos_cent%d_side%d_har%d", icent, iside, ihar);
	h_prof_flatCos[icent][iside][ihar] = (TProfile*)fin2->Get(name);
	sprintf(name, "hProf_flatsin_cent%d_side%d_har%d", icent, iside, ihar);
	h_prof_flatSin[icent][iside][ihar] = (TProfile*)fin2->Get(name);
      }
    }
  }
  /*  
  for (int icent=0; icent<NCENT; icent++) {
    for (int ihar=0; ihar<NHAR; ihar++) {
      sprintf(name, "h_Dphi_fg_cent%d_har%d", icent, ihar);
      hDphi_fg[icent][ihar] = (TH1D*)fin1->Get(name);
      sprintf(name, "h_Dphi_bg_cent%d_har%d", icent, ihar);
      hDphi_bg[icent][ihar] = (TH1D*)fin1->Get(name);
      
      sprintf(name, "h_Res_flat_cent%d_har%d", icent, ihar);
      hRes_flat[icent][ihar] = (TH1D*)fin1->Get(name);
    }
  }
  */
  cout<<"Finished reading histograms"<<endl;
}

void extractor_flat::Run_res() {

  int evts_min=NEv_id*NEv;
  int evts_max=(NEv_id+1)*NEv;

  int nentries = fChain->GetEntries();
  Long64_t nbytes = 0, nb=0;

  cout<<"Run events from "<<evts_min<<" to "<<evts_max<<endl;
  //FCalEt Shift
  for(int ic=0;ic<NCENT_1p;ic++){
    for(int iz=0;iz<NZ;iz++){
      evt[iz][ic]=hevts->GetBinContent(ic+1,iz+1);
    }
  }
  for(int iz=0;iz<NZ;iz++){
    for(int ic=0;ic<NCENT_1p;ic++){
      qmean[iz][ic] = hEt_Vz[iz]->GetBinContent(ic+1);//fcal mean
    }
  }

  //offset correction
  for(int id=0;id<NDET;id++){
    for(int iz=0;iz<NZ;iz++){
      for(int ic=0;ic<NCENT_1p;ic++){
	Qwmean[id][iz][ic] = hQw[id][iz]->GetBinContent(ic+1);
      }
      for(int ihar=0;ihar<NHAR;ihar++){//2,3,4,5
	//hqx[id][iz][ihar]->Divide(hQx[id][iz][ihar],hQw[id][iz]);
	//hqy[id][iz][ihar]->Divide(hQy[id][iz][ihar],hQw[id][iz]);
	for(int ic=0;ic<NCENT_1p;ic++){
	  Qxmean[id][iz][ic][ihar] = hQx[id][iz][ihar]->GetBinContent(ic+1);
	  Qymean[id][iz][ic][ihar] = hQy[id][iz][ihar]->GetBinContent(ic+1);
	  //qxmean[id][iz][ic][ihar] = hqx[id][iz][ihar]->GetBinContent(ic+1);
	  //qymean[id][iz][ic][ihar] = hqy[id][iz][ihar]->GetBinContent(ic+1);
	}
      }
    }
  }


  //Event Loop starts
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    //Long64_t ientry = LoadTree(jentry);
    //if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);
    if (nb <=0)        {cout<<"Nbytes became <=0"<<endl;break;}
    if(jentry<evts_min || jentry>=evts_max) continue;
    if(jentry%10000==0) {cout<<"running event "<<jentry<<endl ;}
    //-------------------------------------------------------------------------------------------------
    centb = get_centb(Fcal_Et);
    vzb=get_vzb(vx_z);
    if (centb < 0) continue;
    if (vzb < 0) continue;
    nevents++;

    FillGlobal_res();

    //Do FCal Recenter Correction
    int iz = (vx_z+100.)/25;if(iz<0) iz=0; if(iz>NZ-1) iz=NZ-1;
    int ic = Centrality;
    int N = evt[iz][ic];
    if(ic<100&&N>5){//only calib events for bins that have more than 5 events.
      double scal = Fcal_Et/qmean[iz][ic];
      for(int id=0;id<NDET;id++){
        for(int ihar=0;ihar<NHAR;ihar++){//2,3,4,5
          //subtract the contribution from present event, avoid self-bias.
          double Qxm = (Qxmean[id][iz][ic][ihar]*N-Qx[id][ihar]);
          double Qym = (Qymean[id][iz][ic][ihar]*N-Qy[id][ihar]);
          double Qwm = (Qwmean[id][iz][ic]*N-Qw[id]);

          double shiftx =  Qxm/Qwm*Qw[id]*scal;
          double shifty =  Qym/Qwm*Qw[id]*scal;
          Qx[id][ihar] -= shiftx;
          Qy[id][ihar] -= shifty;
          hQxS[id][iz][ihar]->Fill(ic,shiftx);
          hQyS[id][iz][ihar]->Fill(ic,shifty);
          //hqxS[id][iz][ihar]->Fill(ic,Qxm/Qwm);
          //hqyS[id][iz][ihar]->Fill(ic,Qym/Qwm);
        }
      }
    }
    
    //Do ZDC Recenter and Flattening Correction
    for (int iside=0; iside<NSIDE+1; iside++) {
      for (int ihar=0; ihar<NHAR; ihar++) {
	double PsiZ_raw = atan2(Qy_zdc[iside][ihar], Qx_zdc[iside][ihar])/double(ihar+1);
	Qx_zdc[iside][ihar]  = ( Qx_zdc[iside][ihar] - h_prof_sumx[iside][ihar]->GetBinContent(ic+1) );
        Qy_zdc[iside][ihar]  = ( Qy_zdc[iside][ihar] - h_prof_sumy[iside][ihar]->GetBinContent(ic+1) );	
        double PsiZ_rec = atan2(Qy_zdc[iside][ihar], Qx_zdc[iside][ihar])/double(ihar+1);
        double psi_rec2pi = (ihar+1)*PsiZ_rec;
	double deltaPsi = 0;
	
	for (int ik=0; ik<NK; ik++) {
          double flatcos = h_prof_flatCos[ic][iside][ihar]->GetBinContent(ik+1);
          double flatsin = h_prof_flatSin[ic][iside][ihar]->GetBinContent(ik+1);
          double cosPsi_rec = cos( (ik+1)*psi_rec2pi );
          double sinPsi_rec = sin( (ik+1)*psi_rec2pi );
          deltaPsi += (-flatsin*cosPsi_rec + flatcos*sinPsi_rec)*2/(ik+1);
        }
	Qx_zdc[iside][ihar]=cos(psi_rec2pi+deltaPsi);
	Qy_zdc[iside][ihar]=sin(psi_rec2pi+deltaPsi);
	//double PsiF_zdc = atan2( sin(psiC2pi+deltaPsi), cos(psiC2pi+deltaPsi) )/double(ihar+1);
	double PsiZ_flat = atan2( Qy_zdc[iside][ihar], Qx_zdc[iside][ihar])/double(ihar+1);
	hPsiZ_raw[centb][iside][ihar]->Fill(PsiZ_raw);
	hPsiZ_rec[centb][iside][ihar]->Fill(PsiZ_rec);
	hPsiZ_flat[centb][iside][ihar]->Fill(PsiZ_flat);
      }
    }

    //float Psi_FCal[NSIDE+1][NHAR];
    for (int iside=0; iside<NSIDE+1; iside++) {
      for (int ihar=0; ihar<NHAR; ihar++) {
	int ieta;
	if(iside==0) Psi_FCal[iside][ihar] = atan2(Qy[0][ihar],Qx[0][ihar])/double(ihar+1);
	if(iside==1) Psi_FCal[iside][ihar] = atan2(Qy[1][ihar],Qx[1][ihar])/double(ihar+1);
	if(iside==2) Psi_FCal[iside][ihar] = atan2(Qy[0][ihar]+Qy[1][ihar],Qx[0][ihar]+Qx[1][ihar])/double(ihar+1);
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


    }
    
    //FillTrack_res(Psi_FCal,Psi_Zdc);
        
    //Mixed Event Pool
    Pool_size =Pool_EP[centb][vzb].size();
    if(Pool_size>=STORE_DEP){
      delete Pool_EP[centb][vzb].at(0);
      Pool_EP[centb][vzb].erase( (Pool_EP[centb][vzb]).begin() );
    }
    Pool_EP[centb][vzb].push_back( ev );
    //delete ev;

  }//Event Loop ends
  //cout<<"cnt_ev = "<<cnt_ev<<endl;
  cout<<"Finished running all events."<<endl;
}

void extractor_flat::FillGlobal_res() {

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


void extractor_flat::Init_res() {
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

  //FCal Correction
  /*
  TH1D *htmp;
  for(int id=0;id<NDET;id++){
    for(int iz=0;iz<NZ;iz++){
      for(int ihar=0;ihar<NHAR;ihar++){//2,3,4,5
	sprintf(name,"hqx_d%d_z%d_h%d",id,iz,ihar);
	htmp = new TH1D(name,"",NCENT_1p,-0.5,NCENT_1p-0.5);  htmp->Sumw2();  //setstyle(htmp); 
	htmp->GetXaxis()->SetTitle("centrality");    htmp->GetYaxis()->SetTitle("q_{x} [GeV]");
	hqx[id][iz][ihar] =htmp;
	sprintf(name,"hqy_d%d_z%d_h%d",id,iz,ihar);
	htmp = new TH1D(name,"",NCENT_1p,-0.5,NCENT_1p-0.5);  htmp->Sumw2();  //setstyle(htmp); 
	htmp->GetXaxis()->SetTitle("centrality");    htmp->GetYaxis()->SetTitle("q_{y} [GeV]");
	hqy[id][iz][ihar] =htmp;
      }
    }
  }
  */

  for(int id=0;id<NDET;id++){
    for(int iz=0;iz<NZ;iz++){
      for(int ihar=0;ihar<NHAR;ihar++){//1,2,3,4,
	sprintf(name,"hQxS_d%d_z%d_h%d",id,iz,ihar);
	hQxS[id][iz][ihar]=new TProfile(name,"",NCENT,-0.5,NCENT-0.5,"s");
	hQxS[id][iz][ihar]->Sumw2();
	sprintf(name,"hQyS_d%d_z%d_h%d",id,iz,ihar);
	hQyS[id][iz][ihar]=new TProfile(name,"",NCENT,-0.5,NCENT-0.5,"s");
	hQxS[id][iz][ihar]->Sumw2();
      }
    }
  }
  
  //ZDC Correction
  for (int icent=0; icent<NCENT; icent++) {
    for (int iside=0; iside<NSIDE+1; iside++) {
      for (int ihar=0; ihar<NHAR; ihar++) {
	sprintf(name, "hPsiZ_raw_cent%d_side%d_har%d", icent, iside, ihar);

	hPsiZ_raw[icent][iside][ihar] = new TH1D(name, "", 100, -2*PI, 2*PI); 
	hPsiZ_raw[icent][iside][ihar]->Sumw2();

	sprintf(name, "hPsiZ_rec_cent%d_side%d_har%d", icent, iside, ihar);
	hPsiZ_rec[icent][iside][ihar] = new TH1D(name, "", 100, -2*PI, 2*PI);
	hPsiZ_rec[icent][iside][ihar]->Sumw2();

	sprintf(name, "hPsiZ_flat_cent%d_side%d_har%d", icent, iside, ihar);
	hPsiZ_flat[icent][iside][ihar] = new TH1D(name, "", 100, -2*PI, 2*PI);
	hPsiZ_flat[icent][iside][ihar]->Sumw2();
      }
    }
  }
  for (int icent=0; icent<NCENT; icent++) {
    for (int iside=0; iside<NSIDE+1; iside++) {
      for (int ihar=0; ihar<NHAR; ihar++) {
	sprintf(name, "hPsi_FCal_cent%d_side%d_har%d", icent, iside, ihar);
	hPsi_FCal[icent][iside][ihar] = new TH1D(name, "", 100, -2*PI, 2*PI);
	hPsi_FCal[icent][iside][ihar]->Sumw2();
	sprintf(name, "hPsi_Zdc_cent%d_side%d_har%d", icent, iside, ihar);
	hPsi_Zdc[icent][iside][ihar] = new TH1D(name, "", 100, -2*PI, 2*PI);
	hPsi_Zdc[icent][iside][ihar]->Sumw2();
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

}


void extractor_flat::SaveHistos_res() {
  fOut_res->cd();
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
	hQxS[id][iz][ihar]->Write();
	hQyS[id][iz][ihar]->Write();
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

  for (int ihar=0; ihar<NHAR; ihar++) {
    hResCent_FCal_fg[ihar]->Write();
    hResCent_FCal_bg[ihar]->Write();
    hResCent_Zdc_fg[ihar]->Write();
    hResCent_Zdc_bg[ihar]->Write();
  }

  hEt->Write();
  hNtrk->Write();
  hCent->Write();
  hEt_ntrk->Write();

  fOut_res->Close();
}

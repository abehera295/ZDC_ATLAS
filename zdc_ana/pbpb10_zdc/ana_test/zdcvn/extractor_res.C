void extractor_flat::ReadHistos_res() {
  /*
    TFile *fIn_res=new TFile("./Output_zdc_psi.root","read");
    for (int iside=0; iside<NSIDE+1; iside++) {
    for (int ihar=0; ihar<NHAR; ihar++) {
    sprintf(name, "hProf_sumx_side%d_har%d", iside, ihar);
    h_prof_sumx[iside][ihar] = (TProfile*)fIn_res->Get(name);
    sprintf(name, "hProf_sumy_side%d_har%d", iside, ihar);
    h_prof_sumy[iside][ihar] = (TProfile*)fIn_res->Get(name);
    }
    }
    for (int icent=0; icent<NCENT; icent++) {
    for (int iside=0; iside<NSIDE+1; iside++) {
    for (int ihar=0; ihar<NHAR; ihar++) {
    sprintf(name, "hProf_flatcos_cent%d_side%d_har%d", icent, iside, ihar);
    h_prof_flatCos[icent][iside][ihar] = (TProfile*)fIn_res->Get(name);
    sprintf(name, "hProf_flatsin_cent%d_side%d_har%d", icent, iside, ihar);
    h_prof_flatSin[icent][iside][ihar] = (TProfile*)fIn_res->Get(name);
    }
    }
    }

    for (int icent=0; icent<NCENT; icent++) {
    for (int ihar=0; ihar<NHAR; ihar++) {
    sprintf(name, "h_Dphi_fg_cent%d_har%d", icent, ihar);
    hDphi_fg[icent][ihar] = (TH1D*)fIn_res->Get(name);
    sprintf(name, "h_Dphi_bg_cent%d_har%d", icent, ihar);
    hDphi_bg[icent][ihar] = (TH1D*)fIn_res->Get(name);

    sprintf(name, "h_Res_flat_cent%d_har%d", icent, ihar);
    hRes_flat[icent][ihar] = (TH1D*)fIn_res->Get(name);
    }
    }
  */
}

void extractor_flat::Run_res() {

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
  double Psi_ECal[4][NSIDE+1][NHAR];
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);

    if (nb <=0)        {cout<<"Nbytes became <=0"<<endl;break;}
    if(jentry<evts_min || jentry>=evts_max) continue;
    if(jentry%10000==0) {cout<<"running event "<<jentry<<endl ;}
    //if (jentry > 10000) break;
    //-------------------------------------------------------------------------------------------------
    centb = get_centb(Fcal_Et);
    vzb=get_vzb(vx_z);
    if (centb < 0) continue;
    if (vzb < 0) continue;
    nevents++;

    FillGlobal_res();
    
    //float Psi_FCal[NSIDE+1][NHAR];
    for (int iside=0; iside<NSIDE+1; iside++) {
      for (int ihar=0; ihar<NHAR; ihar++) {
	//FCal
	if(iside==0) Psi_FCal[iside][ihar] = atan2(Qy[0][ihar],Qx[0][ihar])/double(ihar+1);
	if(iside==1) Psi_FCal[iside][ihar] = atan2(Qy[1][ihar],Qx[1][ihar])/double(ihar+1);
	if(iside==2) Psi_FCal[iside][ihar] = atan2(Qy[0][ihar]+Qy[1][ihar],Qx[0][ihar]+Qx[1][ihar])/double(ihar+1);
	hPsi_FCal[centb][iside][ihar]->Fill(Psi_FCal[iside][ihar]);
	//EMCal
	//2,3(0<eta<0.5) , 4,5(0.5<eta<1.5) , 6,7(1.5<eta<2.1) , 8,9(2.1<eta<2.7)
	//10,11(2.7<eta<3.2) , 12,13(3.2<eta<4) , 14,15(4<eta<5)
	int ieta1;
	for(int iem=0;iem<4;iem++){
	  ieta1=2*(iem+1);
	  if(iside==0) Psi_ECal[iem][iside][ihar] = atan2(Qy[ieta1][ihar],Qx[ieta1][ihar])/double(ihar+1);
	  if(iside==1) Psi_ECal[iem][iside][ihar] = atan2(Qy[ieta1+1][ihar],Qx[ieta1+1][ihar])/double(ihar+1);
	  if(iside==2) Psi_ECal[iem][iside][ihar] = atan2(Qy[ieta1][ihar]+Qy[ieta1+1][ihar],Qx[ieta1][ihar]+Qx[ieta1+1][ihar])/double(ihar+1);
	}
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
      /*
      for(int ihar2=0;ihar2<NHAR;ihar2++){
	if(dPsi_FCal < -(2.*ihar2+1)*PI && dPsi_FCal >= -(2.*ihar2+3)*PI) dPsi_FCal+=2.*(ihar2+1)*PI; 
	if(dPsi_FCal > (2.*ihar2+1)*PI && dPsi_FCal <= (2.*ihar2+3)*PI) dPsi_FCal-=2.*(ihar2+1)*PI; 
      }
      */
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
      /*
      for(int ihar2=0;ihar2<NHAR;ihar2++){
	if(dPsi_Zdc < -(2.*ihar2+1)*PI && dPsi_Zdc >= -(2.*ihar2+3)*PI) dPsi_Zdc+=2.*(ihar2+1)*PI; 
	if(dPsi_Zdc > (2.*ihar2+1)*PI && dPsi_Zdc <= (2.*ihar2+3)*PI) dPsi_Zdc-=2.*(ihar2+1)*PI; 
      }
      */
      hDphi_FCal_fg[centb][ihar]->Fill(dPsi_FCal);
      hDphi_Zdc_fg[centb][ihar]->Fill(dPsi_Zdc);
      
      double dPsi_EF[4][2][2],dPsi_EZ[4][2][2],dPsi_FZ[2][2];
      int cnt=0;
      for(int i=0;i<2;i++){
	for(int j=0;j<2;j++){
	  for(int iem=0;iem<4;iem++){
	    dPsi_EF[iem][i][j]=(ihar+1)*(Psi_ECal[iem][i][ihar]-Psi_FCal[j][ihar]);
	    dPsi_EZ[iem][i][j]=(ihar+1)*(Psi_ECal[iem][i][ihar]-Psi_Zdc[j][0]);
	    hCorrCent_EF_fg[ihar][iem][i][j]->Fill(centb,cos(dPsi_EF[iem][i][j]));
	    hCorrCent_EZ_fg[ihar][iem][i][j]->Fill(centb,cos(dPsi_EZ[iem][i][j]));
	  }
	  dPsi_FZ[i][j]=(ihar+1)*(Psi_FCal[i][ihar]-Psi_Zdc[j][0]);
	  hCorrCent_FZ_fg[ihar][i][j]->Fill(centb,cos(dPsi_FZ[i][j]));
	}
      }
      for(int iem=0;iem<4;iem++){
	double dPsi_EE=(ihar+1)*(Psi_ECal[iem][0][ihar]-Psi_ECal[iem][1][ihar]);
	hCorrCent_EE_fg[ihar][iem]->Fill(centb,cos(dPsi_EE));
      }
      double dPsi_FF=(ihar+1)*(Psi_FCal[0][ihar]-Psi_FCal[1][ihar]);
      double dPsi_ZZ=(ihar+1)*(Psi_Zdc[0][0]-Psi_Zdc[1][0]);
      hCorrCent_FF_fg[ihar]->Fill(centb,cos(dPsi_FF));
      hCorrCent_ZZ_fg[ihar]->Fill(centb,cos(dPsi_ZZ));
      
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
	  /*
	  for(int ihar2=0;ihar2<NHAR;ihar2++){
	    if(dPsi_FCal < -(2.*ihar2+1)*PI && dPsi_FCal >= -(2.*ihar2+3)*PI) dPsi_FCal+=2.*(ihar2+1)*PI;
	    if(dPsi_FCal > (2.*ihar2+1)*PI && dPsi_FCal <= (2.*ihar2+3)*PI) dPsi_FCal-=2.*(ihar2+1)*PI;
	  }
	  */
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
	  /*
	  for(int ihar2=0;ihar2<NHAR;ihar2++){
	    if(dPsi_FCal < -(2.*ihar2+1)*PI && dPsi_FCal >= -(2.*ihar2+3)*PI) dPsi_FCal+=2.*(ihar2+1)*PI;
	    if(dPsi_FCal > (2.*ihar2+1)*PI && dPsi_FCal <= (2.*ihar2+3)*PI) dPsi_FCal-=2.*(ihar2+1)*PI;
	  }
	  */
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
	  /*
	  for(int ihar2=0;ihar2<NHAR;ihar2++){
	    if(dPsi_Zdc < -(2.*ihar2+1)*PI && dPsi_Zdc >= -(2.*ihar2+3)*PI) dPsi_Zdc+=2.*(ihar2+1)*PI;
	    if(dPsi_Zdc > (2.*ihar2+1)*PI && dPsi_Zdc <= (2.*ihar2+3)*PI) dPsi_Zdc-=2.*(ihar2+1)*PI;
	  }
	  */
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
      dPsi_FCal = (ihar+1)*( Psi_FCal[1][ihar] - (Psi_FCal_store[centb_1p][0][ihar]).at(isiz) );
      if (dPsi_FCal < -PI) dPsi_FCal += 2*PI;
      if (dPsi_FCal >  PI) dPsi_FCal -= 2*PI;
      hDphi_FCal_bg[centb][ihar]->Fill(dPsi_FCal);
      }
      }
      if(str_size1>0){
      for (int isiz=0; isiz<str_size1; isiz++) {
      dPsi_FCal = (ihar+1)*( Psi_FCal[0][ihar] - (Psi_FCal_store[centb_1p][1][ihar]).at(isiz) );
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
    
    //FillTrack_res(Psi_FCal,Psi_Zdc);
    
    /*
    //Old code
    //Fill store for Psi background
    for (int iside=0; iside<NSIDE+1; iside++) {
    for (int ihar=0; ihar<NHAR; ihar++) {
    (Psi_FCal_store[centb_1p][iside][ihar]).push_back(Psi_FCal[iside][ihar]);
    (Psi_Zdc_store[centb_1p][iside][ihar]).push_back(Psi_Zdc[iside][ihar]);
    if ((Psi_FCal_store[centb_1p][iside][ihar]).size() >= STORE_DEP) (Psi_FCal_store[centb_1p][iside][ihar]).erase((Psi_FCal_store[centb_1p][iside][ihar]).begin());
    if ((Psi_Zdc_store[centb_1p][iside][ihar]).size() >= STORE_DEP) (Psi_Zdc_store[centb_1p][iside][ihar]).erase((Psi_Zdc_store[centb_1p][isid
e][ihar]).begin());
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
  
  //TProfiles hResCent
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
  //3-SubEvent hRes
  for (int ihar=0; ihar<NHAR; ihar++) {
    for(int i=0;i<2;i++){
      for(int j=0;j<2;j++){
	for(int iem=0;iem<4;iem++){
	  sprintf(name, "hCorrCent_EF_fg_har%d_EMCal%d_%d_%d",ihar,iem,i,j);
	  hCorrCent_EF_fg[ihar][iem][i][j]= new TProfile(name, "", NCENT, 0, NCENT);
	  sprintf(name, "hCorrCent_EZ_fg_har%d_EMCal%d_%d_%d",ihar,iem,i,j);
	  hCorrCent_EZ_fg[ihar][iem][i][j]= new TProfile(name, "", NCENT, 0, NCENT);
	}
	sprintf(name, "hCorrCent_FZ_fg_har%d_%d_%d",ihar,i,j);
	hCorrCent_FZ_fg[ihar][i][j]= new TProfile(name, "", NCENT, 0, NCENT);
      }
    }
    for(int iem=0;iem<4;iem++){
      sprintf(name, "hCorrCent_EE_fg_har%d_EMCal%d",ihar,iem);
      hCorrCent_EE_fg[ihar][iem]= new TProfile(name, "", NCENT, 0, NCENT);
    }
    sprintf(name, "hCorrCent_FF_fg_har%d",ihar);
    hCorrCent_FF_fg[ihar]= new TProfile(name, "", NCENT, 0, NCENT);
    sprintf(name, "hCorrCent_ZZ_fg_har%d",ihar);
    hCorrCent_ZZ_fg[ihar]= new TProfile(name, "", NCENT, 0, NCENT);
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
  
  //2-SubEvent hRes
  for (int ihar=0; ihar<NHAR; ihar++) {
    hResCent_FCal_fg[ihar]->Write();
    hResCent_FCal_bg[ihar]->Write();
    hResCent_Zdc_fg[ihar]->Write();
    hResCent_Zdc_bg[ihar]->Write();
  }

  //3-SubEvent hRes
  for (int ihar=0; ihar<NHAR; ihar++) {
    for(int i=0;i<2;i++){
      for(int j=0;j<2;j++){
	for(int iem=0;iem<4;iem++){
	  hCorrCent_EF_fg[ihar][iem][i][j]->Write();
	  hCorrCent_EZ_fg[ihar][iem][i][j]->Write();
	}
	hCorrCent_FZ_fg[ihar][i][j]->Write();
      }
    }
    for(int iem=0;iem<4;iem++){
      hCorrCent_EE_fg[ihar][iem]->Write();
    }
    hCorrCent_FF_fg[ihar]->Write();
    hCorrCent_ZZ_fg[ihar]->Write();
  }
  
  hEt->Write();
  hNtrk->Write();
  hCent->Write();
  hEt_ntrk->Write();

  fOut_res->Close();
}

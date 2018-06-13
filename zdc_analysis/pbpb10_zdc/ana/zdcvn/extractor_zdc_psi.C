
void extractor_flat::Run_zdc_psi() {
  int evts_min=NEv_id*NEv;
  int evts_max=(NEv_id+1)*NEv;

  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;

  for (int icent=0; icent<NCENT_1p; icent++) {
    for (int iside=0; iside<NSIDE; iside++) {
      for (int ihar=0; ihar<NHAR; ihar++) {
	(Psi_store[icent][iside][ihar]).assign(0,0);
      }
    }
  }

  for (Long64_t jentry=0; jentry<nentries;jentry++) {
     
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);

    if (nb <=0)        {cout<<"Nbytes became <=0"<<endl;break;}
    if(jentry<evts_min || jentry>=evts_max) continue;
    if(jentry%10000==0) {cout<<"running event "<<jentry<<endl ;}
    //if (ientry > 20000) break;
    //-------------------------------------------------------------------------------------------------
    centb = get_centb(Fcal_Et);
    if (centb < 0) continue;
    nevents++;
    FillGlobal();
    /*
    std::vector<std::vector<double> > Qx;
    std::vector<std::vector<double> > Qy;
    std::vector<double> Qw;
    get_qVec_zdc(Qx, Qy, Qw);

    std::vector<std::vector<double> > QxC  = Qx;
    std::vector<std::vector<double> > QyC  = Qy;
    std::vector<std::vector<double> > Psi  = Qy;
    std::vector<std::vector<double> > PsiC = Qy;
    std::vector<std::vector<double> > PsiF = Qy;
    */
    double PsiF[NSIDE+1][NHAR];
    for (int iside=0; iside<NSIDE+1; iside++) {
      for (int ihar=0; ihar<NHAR; ihar++) {
	/*
	QxC.at(iside).at(ihar)  = ( Qx.at(iside).at(ihar) - h_prof_sumx[iside][ihar]->GetBinContent(centb+1) );
	QyC.at(iside).at(ihar)  = ( Qy.at(iside).at(ihar) - h_prof_sumy[iside][ihar]->GetBinContent(centb+1) );
	Psi.at(iside).at(ihar)  = atan2(Qy.at(iside).at(ihar) , Qx.at(iside).at(ihar) )/double(ihar+1);
	PsiC.at(iside).at(ihar) = atan2(QyC.at(iside).at(ihar), QxC.at(iside).at(ihar))/double(ihar+1);

	double psiC2pi = (ihar+1)*PsiC.at(iside).at(ihar);
	double deltaPsi = 0;
	for (int ik=0; ik<NK; ik++) {
	  double flatcos = h_prof_flatCos[centb][iside][ihar]->GetBinContent(ik+1);
	  double flatsin = h_prof_flatSin[centb][iside][ihar]->GetBinContent(ik+1);
	  double cosPsiC = cos( (ik+1)*psiC2pi );
	  double sinPsiC = sin( (ik+1)*psiC2pi );

	  deltaPsi += (-flatsin*cosPsiC + flatcos*sinPsiC)*2/(ik+1);
	}
	PsiF.at(iside).at(ihar) = atan2( sin(psiC2pi+deltaPsi), cos(psiC2pi+deltaPsi) )/double(ihar+1); 
	*/
	/*
	h_Qx [centb][iside][ihar]->Fill(Qx.at(iside).at(ihar));
	h_Qy [centb][iside][ihar]->Fill(Qy.at(iside).at(ihar));
	h_QxC[centb][iside][ihar]->Fill(QxC.at(iside).at(ihar));
	h_QyC[centb][iside][ihar]->Fill(QyC.at(iside).at(ihar));

	hPsi [centb][iside][ihar]->Fill(Psi .at(iside).at(ihar));
	hPsiC[centb][iside][ihar]->Fill(PsiC.at(iside).at(ihar));
	*/
	//PsiF.at(iside).at(ihar) = atan2(Qy.at(iside).at(ihar) , Qx.at(iside).at(ihar) )/double(ihar+1);
	PsiF[iside][ihar]=atan2(Qy_zdc[iside][ihar] , Qx_zdc[iside][ihar] )/double(ihar+1);
	hPsiF[centb][iside][ihar]->Fill(PsiF[iside][ihar]);

	(Psi_store[centb_1p][iside][ihar]).push_back(PsiF[iside][ihar]);
	if ((Psi_store[centb_1p][iside][ihar]).size() >= STORE_DEP) (Psi_store[centb_1p][iside][ihar]).erase((Psi_store[centb_1p][iside][ihar]).begin());
      }
    }

    for (int ihar=0; ihar<NHAR; ihar++) {
      //hRes     [centb][ihar]->Fill(cos( (ihar+1)* (Psi.at (0).at(ihar) - Psi.at (1).at(ihar)) ));
      //hRes_recn[centb][ihar]->Fill(cos( (ihar+1)* (PsiC.at(0).at(ihar) - PsiC.at(1).at(ihar)) ));
      //hRes_flat[centb][ihar]->Fill(cos( (ihar+1)* (PsiF.at(0).at(ihar) - PsiF.at(1).at(ihar)) ));
      hRes_flat[centb][ihar]->Fill(cos( (ihar+1)* (PsiF[0][ihar] - PsiF[1][ihar]) ));

      //double dPsi = (ihar+1)*( PsiF.at(0).at(ihar) -  PsiF.at(1).at(ihar) );
      //double dPsi = (ihar+1)*( PsiF.at(1).at(ihar) -  PsiF.at(0).at(ihar) );
      double dPsi = (ihar+1)*( PsiF[1][ihar] -  PsiF[0][ihar] );
      if (dPsi < -PI) dPsi += 2*PI;
      if (dPsi >  PI) dPsi -= 2*PI;
      hDphi_fg[centb][ihar]->Fill(dPsi);

      //int str_dep = (Psi_store[centb_1p][0][ihar]).size();
      int str_dep = (Psi_store[centb_1p][1][ihar]).size();
      for (int isiz=0; isiz<str_dep; isiz++) {
	//dPsi = (ihar+1)*( PsiF.at(0).at(ihar) - (Psi_store[centb_1p][1][ihar]).at(isiz) );
	dPsi = (ihar+1)*( PsiF[1][ihar] - (Psi_store[centb_1p][0][ihar]).at(isiz) );
	if (dPsi < -PI) dPsi += 2*PI;
	if (dPsi >  PI) dPsi -= 2*PI;
	hDphi_bg[centb][ihar]->Fill( dPsi);
      }
    }

    //FillTrackVn(PsiF);
  }

}


void extractor_flat::FillFlattening(){

  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;

  for (Long64_t jentry=0; jentry<nentries;jentry++) {

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);

    if (nb <=0)        {cout<<"Nbytes became <=0"<<endl;break;}
    if(jentry%100000==0) {cout<<"In FillFlattening:: event "<<jentry<<endl ;}
    //if (ientry > 20000) break;
    //-------------------------------------------------------------------------------------------------
    centb = get_centb(Fcal_Et);
    if (centb < 0) continue;
    nevents++;

    std::vector<std::vector<double> > Qx;
    std::vector<std::vector<double> > Qy;
    std::vector<double> Qw;
    get_qVec_zdc(Qx, Qy, Qw);

    std::vector<std::vector<double> > QxC  = Qx;
    std::vector<std::vector<double> > QyC  = Qy;
    std::vector<std::vector<double> > PsiC = Qy;

    for (int iside=0; iside<NSIDE+1; iside++) {
      for (int ihar=0; ihar<NHAR; ihar++) {

	QxC.at(iside).at(ihar)  = ( Qx.at(iside).at(ihar) - h_prof_sumx[iside][ihar]->GetBinContent(centb+1) );
	QyC.at(iside).at(ihar)  = ( Qy.at(iside).at(ihar) - h_prof_sumy[iside][ihar]->GetBinContent(centb+1) );
	PsiC.at(iside).at(ihar) = atan2(QyC.at(iside).at(ihar), QxC.at(iside).at(ihar))/double(ihar+1);

	double psi2pi = PsiC.at(iside).at(ihar)*(ihar+1);
	for (int ik=0; ik<NK; ik++) {
	  h_prof_flatCos[centb][iside][ihar]->Fill(ik, cos( (ik+1)*psi2pi) );
	  h_prof_flatSin[centb][iside][ihar]->Fill(ik, sin( (ik+1)*psi2pi) );
	}
      }
    }
  }
}


void extractor_flat::FillCalib() {

  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  sampSel.assign(NCENT, 0);

  double QxMax=0, QyMax=0, QwMax=0;

  for (Long64_t jentry=0; jentry<nentries;jentry++) {

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   

    if (nb <=0)        {cout<<"Nbytes became <=0"<<endl;break;}
    if(jentry%100000==0) {cout<<"In FillCalib::  event "<<jentry<<endl ;}
    //if (ientry > 20000) break;
    //-------------------------------------------------------------------------------------------------
    centb = get_centb(Fcal_Et);
    if (centb < 0) continue;

    sampSel[centb]++;
    if (sampSel[centb]>NSAMP) sampSel[centb] = NSAMP;

    std::vector<std::vector<double> > Qx;
    std::vector<std::vector<double> > Qy;
    std::vector<double> Qw;
    get_qVec_zdc(Qx, Qy, Qw);
    if (Qx.at(0).at(0) > QxMax) QxMax = Qx.at(0).at(0);
    if (Qy.at(0).at(0) > QyMax) QyMax = Qy.at(0).at(0);
    if (Qw.at(0) > QwMax) QwMax = Qw.at(0);

    for (int iside=0; iside<NSIDE+1; iside++) {
      for (int ihar=0; ihar<NHAR; ihar++) {
	h_prof_sumx [iside][ihar]->Fill(centb, Qx.at(iside).at(ihar));
	h_prof_sumy [iside][ihar]->Fill(centb, Qy.at(iside).at(ihar));
	h_Q2D[centb][iside][ihar]->Fill(Qx.at(iside).at(ihar)/Qw.at(iside), Qy.at(iside).at(ihar)/Qw.at(iside));
      }
    }

    FillGlobal();
  }

  std::cout << "QxMax = " << QxMax << ": QyMax = " << QyMax << ": QwMax = " << QwMax << std::endl;
}


void extractor_flat::FillGlobal() {

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

void extractor_flat::Init_zdc_psi() {

  set_ptBins ();
  set_etaBins();

  //char name[500];

  for (int iside=0; iside<NSIDE; iside++) {
    sprintf(name, "h_eTot_side%d", iside);
    h_eTot[iside] = new TH1D(name, "", 4000, -0.5, 3999.5); h_eTot[iside]->Sumw2();
    sprintf(name, "h_ePix_side%d", iside);
    h_ePix[iside] = new TH1D(name, "", 4000, -0.5, 3999.5); h_ePix[iside]->Sumw2();
   
    sprintf(name, "h_et_eTot_side%d", iside);
    hEt_eTot[iside] = new TH2D(name, "", 500, 0, 5000, 400, 0, 4000); hEt_eTot[iside]->Sumw2();
    sprintf(name, "h_et_ePix_side%d", iside);
    hEt_ePix[iside] = new TH2D(name, "", 500, 0, 5000, 200, 0, 2000); hEt_ePix[iside]->Sumw2();
      
    sprintf(name, "h_ntrk_eTot_side%d", iside);
    hNtrk_eTot[iside] = new TH2D(name, "", 400, 0, 4000, 400, 0, 4000); hNtrk_eTot[iside]->Sumw2();
    sprintf(name, "h_ntrk_ePix_side%d", iside);
    hNtrk_ePix[iside] = new TH2D(name, "", 400, 0, 4000, 200, 0, 2000); hNtrk_ePix[iside]->Sumw2();
  }
  
  hEt      = new TH1D("hEt", "", 5000, 0, 5000); hEt->Sumw2();
  hNtrk    = new TH1D("hNtrk", "", 5000, 0, 5000); hNtrk->Sumw2();
  hCent    = new TH1D("hCent", "", NCENT, -0.5, NCENT-0.5); hCent->Sumw2();
  hEt_ntrk = new TH2D("h_et_ntrk", "", 400, 0, 4000, 400, 0, 4000); hEt_ntrk->Sumw2();

  for (int iside=0; iside<NSIDE+1; iside++) {
    for (int ihar=0; ihar<NHAR; ihar++) {
      sprintf(name, "hProf_sumx_side%d_har%d", iside, ihar);
      h_prof_sumx[iside][ihar] = new TProfile(name," ",100,-0.5,99.5,-3000,3000,"s");
      sprintf(name, "hProf_sumy_side%d_har%d", iside, ihar);
      h_prof_sumy[iside][ihar] = new TProfile(name," ",100,-0.5,99.5,-3000,3000,"s");

      for (int icent=0; icent<NCENT; icent++) {
	sprintf(name, "h_Q2D_cent%d_side%d_har%d", icent, iside, ihar);
	h_Q2D[icent][iside][ihar] = new TH2D(name, "", 200, -1.0, 1.0, 100, -1.0, 1.0);
	h_Q2D[icent][iside][ihar]->Sumw2();
      }
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

	sprintf(name, "h_Psi_cent%d_side%d_har%d", icent, iside, ihar);
	hPsi[icent][iside][ihar] = new TH1D(name, "", 100, -PI, PI); 
	hPsi[icent][iside][ihar]->Sumw2();

	sprintf(name, "h_PsiC_cent%d_side%d_har%d", icent, iside, ihar);
	hPsiC[icent][iside][ihar] = new TH1D(name, "", 100, -PI, PI);
	hPsiC[icent][iside][ihar]->Sumw2();
            
	sprintf(name, "h_PsiF_cent%d_side%d_har%d", icent, iside, ihar);
	hPsiF[icent][iside][ihar] = new TH1D(name, "", 100, -PI, PI);
	hPsiF[icent][iside][ihar]->Sumw2();
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

  for (int icent=0; icent<NCENT; icent++) {
    for (int ihar=0; ihar<NHAR; ihar++) {
      sprintf(name, "h_Dphi_fg_cent%d_har%d", icent, ihar);
      hDphi_fg[icent][ihar] = new TH1D(name, "", 640, -3.2, 3.2); hDphi_fg[icent][ihar]->Sumw2();
      sprintf(name, "h_Dphi_bg_cent%d_har%d", icent, ihar);
      hDphi_bg[icent][ihar] = new TH1D(name, "", 640, -3.2, 3.2); hDphi_bg[icent][ihar]->Sumw2();
    }
  }

  
  for (int icent=0; icent<NCENT; icent++) {
    for (int ihar=0; ihar<NHAR; ihar++) {

      sprintf(name, "h_Res_cent%d_har%d", icent, ihar);
      hRes[icent][ihar] = new TH1D(name, "", 200, -1.0, 1.0);
      hRes[icent][ihar]->Sumw2();

      sprintf(name, "h_Res_recn_cent%d_har%d", icent, ihar);
      hRes_recn[icent][ihar] = new TH1D(name, "", 200, -1.0, 1.0);
      hRes_recn[icent][ihar]->Sumw2();

      sprintf(name, "h_Res_flat_cent%d_har%d", icent, ihar);
      hRes_flat[icent][ihar] = new TH1D(name, "", 200, -1.0, 1.0);
      hRes_flat[icent][ihar]->Sumw2();
    }
  }

  for (int ihar=0; ihar<NHAR; ihar++) {
    sprintf(name, "h_ResCent_har%d", ihar);
    hResCent[ihar] = new TH1D(name, "", NCENT, 0, 100);
    hResCent[ihar]->Sumw2();

    sprintf(name, "h_ResCent_recn_har%d", ihar);
    hResCent_recn[ihar] = new TH1D(name, "", NCENT, 0, 100);
    hResCent_recn[ihar]->Sumw2();
      
    sprintf(name, "h_ResCent_flat_har%d", ihar);
    hResCent_flat[ihar] = new TH1D(name, "", NCENT, 0, 100);
    hResCent_flat[ihar]->Sumw2();
      
    sprintf(name, "h_ResCent_mix_har%d", ihar);
    hResCent_mix [ihar] = new TH1D(name, "", NCENT, 0, 100);
    hResCent_mix [ihar]->Sumw2();
  }
}


void extractor_flat::SaveHistos_zdc_psi() {
  std::cout << "Saving" << std::endl;

  fOut_zdc_psi->cd();
  /*
  for (int iside=0; iside<NSIDE; iside++) {
    h_eTot[iside]->Write();
    h_ePix[iside]->Write();
    hEt_eTot[iside]->Write();
    hEt_ePix[iside]->Write();
    hNtrk_eTot[iside]->Write();
    hNtrk_ePix[iside]->Write();
  }
  
  
    for (int iside=0; iside<NSIDE+1; iside++) {
    for (int ihar=0; ihar<NHAR; ihar++) {
      h_prof_sumx[iside][ihar]->Write();
      h_prof_sumy[iside][ihar]->Write();
      for (int icent=0; icent<NCENT; icent++) {
	h_Q2D[icent][iside][ihar]->Write();
      }
    }
  }

  for (int icent=0; icent<NCENT; icent++) {
    for (int iside=0; iside<NSIDE+1; iside++) {
      for (int ihar=0; ihar<NHAR; ihar++) {
	h_prof_flatCos[icent][iside][ihar]->Write();
	h_prof_flatSin[icent][iside][ihar]->Write();
      }
    }
  }
  
  for (int icent=0; icent<NCENT; icent++) {
    for (int iside=0; iside<NSIDE+1; iside++) {
      for (int ihar=0; ihar<NHAR; ihar++) {
	hPsi[icent][iside][ihar]->Write();
	hPsiC[icent][iside][ihar]->Write();
	hPsiF[icent][iside][ihar]->Write();
      }
    }
  }
  
  for (int icent=0; icent<NCENT; icent++) {
    for (int iside=0; iside<NSIDE+1; iside++) {
      for (int ihar=0; ihar<NHAR; ihar++) {
	h_Qx[icent][iside][ihar]->Write();
	h_Qy[icent][iside][ihar]->Write();
	h_QxC[icent][iside][ihar]->Write();
	h_QyC[icent][iside][ihar]->Write();
      }
    }
  }

  for (int icent=0; icent<NCENT; icent++) {
    for (int ihar=0; ihar<NHAR; ihar++) {
      hDphi_fg[icent][ihar]->Write();
      hDphi_bg[icent][ihar]->Write();
      hDphi_rat[icent][ihar]->Write();
    }
  }

  
  for (int icent=0; icent<NCENT; icent++) {
    for (int ihar=0; ihar<NHAR; ihar++) {
      hRes[icent][ihar]->Write();
      hRes_recn[icent][ihar]->Write();
      hRes_flat[icent][ihar]->Write();
    }
  }

  for (int ihar=0; ihar<NHAR; ihar++) {
    hResCent[ihar]->Write();
    hResCent_recn[ihar]->Write();
    hResCent_flat[ihar]->Write();
    hResCent_mix [ihar]->Write();
  }
  hEt->Write();
  hNtrk->Write();
  hCent->Write();
  hEt_ntrk->Write();
*/

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
        hPsiF[icent][iside][ihar]->Write();
      }
    }
  }
  for (int icent=0; icent<NCENT; icent++) {
    for (int ihar=0; ihar<NHAR; ihar++) {
      hDphi_fg[icent][ihar]->Write();
      hDphi_bg[icent][ihar]->Write();
      //hDphi_rat[icent][ihar]->Write();
    }
  }


  for (int icent=0; icent<NCENT; icent++) {
    for (int ihar=0; ihar<NHAR; ihar++) {
      //hRes[icent][ihar]->Write();
      //hRes_recn[icent][ihar]->Write();
      hRes_flat[icent][ihar]->Write();
    }
  }
  hEt->Write();
  hNtrk->Write();
  hCent->Write();
  hEt_ntrk->Write();
  fOut_zdc_psi->Close();
  //fOut_zdc_psi->Write();
  std::cout << "Saving Finished" << std::endl;
  std::cout << " total events: " << nevents << std::endl;

}

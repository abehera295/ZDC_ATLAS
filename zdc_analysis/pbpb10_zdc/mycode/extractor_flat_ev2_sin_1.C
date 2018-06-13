#define extractor_flat_cxx
#include "extractor_flat.h"
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

#include "map_zdcPix.C"


void extractor_flat::set_etaBins() {

   float diffEtaArr[] = {-2.5, -2.0, -1.5, -1.0, -0.5, -0.25, 0.0, 0.25, 0.5, 1.0, 1.5, 2.0, 2.5};
   diffEtaBins.assign(diffEtaArr, diffEtaArr + sizeof(diffEtaArr)/sizeof(float));
   NETA = diffEtaBins.size()-1;
   std::cout << "NDIFFETA  = " << NETA << ": diffEtaBins:: {";
   std::copy(diffEtaBins.begin(), diffEtaBins.end(), std::ostream_iterator<float>(std::cout, ", "));  std::cout << '\b' << '\b' << "}" << std::endl;;
}


void extractor_flat::set_ptBins() {

   float diffPtArr[] = {0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 7.0, 9.0, 12.0};
   diffPtBins.assign(diffPtArr, diffPtArr + sizeof(diffPtArr)/sizeof(float));
   NPT = diffPtBins.size()-1;
   std::cout << "NPT       = " << NPT << ": diffPtBins :: {";
   std::copy(diffPtBins.begin(), diffPtBins.end(), std::ostream_iterator<float>(std::cout, ", "));  std::cout << '\b' << '\b' << "}" << std::endl;;
}


void extractor_flat::exec(int entSel) {
   if (fChain == 0) return;
   nevents = 0;

   fOut = new TFile("Output_zdcPix_flat_wtrk_local_wv2_sin.root", "recreate");

   SetUpAll();
   FillCalib();
   FillFlattening();
   run();
   doAna();

   SaveHistos();
}


void extractor_flat::run() {

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
      if(jentry%100000==0) {cout<<"running event "<<jentry<<endl ;}

      //-------------------------------------------------------------------------------------------------
      centb = get_centb(Fcal_Et);
      if (centb < 0) continue;
      nevents++;

      std::vector<std::vector<double> > Qx;
      std::vector<std::vector<double> > Qy;
      std::vector<double> Qw;
      get_qVec(Qx, Qy, Qw);

      std::vector<std::vector<double> > QxC  = Qx;
      std::vector<std::vector<double> > QyC  = Qy;
      std::vector<std::vector<double> > Psi  = Qy;
      std::vector<std::vector<double> > PsiC = Qy;
      std::vector<std::vector<double> > PsiF = Qy;

      for (int iside=0; iside<NSIDE+1; iside++) {
         for (int ihar=0; ihar<NHAR; ihar++) {

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

            h_Qx [centb][iside][ihar]->Fill(Qx.at(iside).at(ihar));
            h_Qy [centb][iside][ihar]->Fill(Qy.at(iside).at(ihar));
            h_QxC[centb][iside][ihar]->Fill(QxC.at(iside).at(ihar));
            h_QyC[centb][iside][ihar]->Fill(QyC.at(iside).at(ihar));

            hPsi [centb][iside][ihar]->Fill(Psi .at(iside).at(ihar));
            hPsiC[centb][iside][ihar]->Fill(PsiC.at(iside).at(ihar));
            hPsiF[centb][iside][ihar]->Fill(PsiF.at(iside).at(ihar));

            (Psi_store[centb_1p][iside][ihar]).push_back(PsiF.at(iside).at(ihar));
            if ((Psi_store[centb_1p][iside][ihar]).size() >= STORE_DEP) (Psi_store[centb_1p][iside][ihar]).erase((Psi_store[centb_1p][iside][ihar]).begin());
         }
      }

      for (int ihar=0; ihar<NHAR; ihar++) {
         hRes     [centb][ihar]->Fill(cos( (ihar+1)* (Psi.at (0).at(ihar) - Psi.at (1).at(ihar)) ));
         hRes_recn[centb][ihar]->Fill(cos( (ihar+1)* (PsiC.at(0).at(ihar) - PsiC.at(1).at(ihar)) ));
         hRes_flat[centb][ihar]->Fill(cos( (ihar+1)* (PsiF.at(0).at(ihar) - PsiF.at(1).at(ihar)) ));

         //double dPsi = (ihar+1)*( PsiF.at(0).at(ihar) -  PsiF.at(1).at(ihar) );
         double dPsi = (ihar+1)*( PsiF.at(1).at(ihar) -  PsiF.at(0).at(ihar) );
         if (dPsi < -PI) dPsi += 2*PI;
         if (dPsi >  PI) dPsi -= 2*PI;
         hDphi_fg[centb][ihar]->Fill(dPsi);

         //int str_dep = (Psi_store[centb_1p][0][ihar]).size();
         int str_dep = (Psi_store[centb_1p][1][ihar]).size();
         for (int isiz=0; isiz<str_dep; isiz++) {
            //dPsi = (ihar+1)*( PsiF.at(0).at(ihar) - (Psi_store[centb_1p][1][ihar]).at(isiz) );
            dPsi = (ihar+1)*( PsiF.at(1).at(ihar) - (Psi_store[centb_1p][0][ihar]).at(isiz) );
            if (dPsi < -PI) dPsi += 2*PI;
            if (dPsi >  PI) dPsi -= 2*PI;
            hDphi_bg[centb][ihar]->Fill( dPsi);
         }
      }

      FillTrackVn(PsiF);
   }
}


void extractor_flat::doAna() {

   for (int ihar=0; ihar<NHAR; ihar++) {
      for (int icent=0; icent<NCENT; icent++) {
         hResCent     [ihar]->SetBinContent(icent+1, hRes     [icent][ihar]->GetMean());
         hResCent     [ihar]->SetBinError  (icent+1, hRes     [icent][ihar]->GetMeanError());
         hResCent_recn[ihar]->SetBinContent(icent+1, hRes_recn[icent][ihar]->GetMean());
         hResCent_recn[ihar]->SetBinError  (icent+1, hRes_recn[icent][ihar]->GetMeanError());
         hResCent_flat[ihar]->SetBinContent(icent+1, hRes_flat[icent][ihar]->GetMean());
         hResCent_flat[ihar]->SetBinError  (icent+1, hRes_flat[icent][ihar]->GetMeanError());

         char name[500];
         sprintf(name, "h_Dphi_rat_cent%d_har%d", icent, ihar);
         hDphi_rat[icent][ihar] = (TH1D*)hDphi_fg[icent][ihar]->Clone(name);
         hDphi_rat[icent][ihar]->Divide(hDphi_bg[icent][ihar]);
         double norm_sc = hDphi_bg[icent][ihar]->Integral()/(2*PI);
         hDphi_rat[icent][ihar]->Scale(norm_sc);

         std::vector<double> res_vec = discFourrier(hDphi_rat[icent][ihar], 1);
         hResCent_mix[ihar]->SetBinContent(icent+1, res_vec.at(0) );
         hResCent_mix[ihar]->SetBinError  (icent+1, res_vec.at(1) );
      }
   }

   for (int icent=0; icent<NCENT; icent++) {
      for (int iside=0; iside<NSIDE+1; iside++) {
         for (int ihar=0; ihar<NHAR; ihar++) {

            double reso_2sub   = hResCent_mix[0]->GetBinContent(icent+1);
            double reso_2sub_e = hResCent_mix[0]->GetBinError  (icent+1);
            if (iside == NSIDE) {
               reso_2sub   = sqrt(2*reso_2sub);
               reso_2sub_e = reso_2sub_e/reso_2sub;
            }

            for (int ieta=0; ieta<NETA+2; ieta++) {
               for (int ipt=0; ipt<NPT; ipt++) {
                  double vn     = (hVn [icent][iside][ihar]).at(ipt).at(ieta)->GetMean();
                  double vn_e   = (hVn [icent][iside][ihar]).at(ipt).at(ieta)->GetMeanError();
                  (hVn_pt_raw[icent][iside][ihar]).at(ieta)->SetBinContent(ipt+1, vn  );
                  (hVn_pt_raw[icent][iside][ihar]).at(ieta)->SetBinError  (ipt+1, vn_e);
                  (hVn_pt_res[icent][iside][ihar]).at(ieta)->SetBinContent(ipt+1, vn/reso_2sub);
                  (hVn_pt_res[icent][iside][ihar]).at(ieta)->SetBinError  (ipt+1, vn/reso_2sub*sqrt( pow(vn_e/vn,2) + pow(reso_2sub_e/reso_2sub,2)));
                  
                  double vns    = (hVnS[icent][iside][ihar]).at(ipt).at(ieta)->GetMean();
                  double vns_e  = (hVnS[icent][iside][ihar]).at(ipt).at(ieta)->GetMeanError();
                  (hVnS_pt_raw[icent][iside][ihar]).at(ieta)->SetBinContent(ipt+1, vns  );
                  (hVnS_pt_raw[icent][iside][ihar]).at(ieta)->SetBinError  (ipt+1, vns_e);
                  (hVnS_pt_res[icent][iside][ihar]).at(ieta)->SetBinContent(ipt+1, vns/reso_2sub);
                  (hVnS_pt_res[icent][iside][ihar]).at(ieta)->SetBinError  (ipt+1, vns/reso_2sub*sqrt( pow(vns_e/vns,2) + pow(reso_2sub_e/reso_2sub,2)));
               }
            }

            for (int ipt=0; ipt<NPT+1; ipt++) {
               for (int ieta=0; ieta<NETA; ieta++) {
                  double vn    = (hVn [icent][iside][ihar]).at(ipt).at(ieta)->GetMean();
                  double vn_e  = (hVn [icent][iside][ihar]).at(ipt).at(ieta)->GetMeanError();
                  (hVn_eta_raw[icent][iside][ihar]).at(ipt)->SetBinContent(ieta+1, vn  );
                  (hVn_eta_raw[icent][iside][ihar]).at(ipt)->SetBinError  (ieta+1, vn_e);
                  (hVn_eta_res[icent][iside][ihar]).at(ipt)->SetBinContent(ieta+1, vn/reso_2sub);
                  (hVn_eta_res[icent][iside][ihar]).at(ipt)->SetBinError  (ieta+1, vn/reso_2sub*sqrt( pow(vn_e/vn,2) + pow(reso_2sub_e/reso_2sub,2)));
                  
                  double vns   = (hVnS[icent][iside][ihar]).at(ipt).at(ieta)->GetMean();
                  double vns_e = (hVnS[icent][iside][ihar]).at(ipt).at(ieta)->GetMeanError();
                  (hVnS_eta_raw[icent][iside][ihar]).at(ipt)->SetBinContent(ieta+1, vns  );
                  (hVnS_eta_raw[icent][iside][ihar]).at(ipt)->SetBinError  (ieta+1, vns_e);
                  (hVnS_eta_res[icent][iside][ihar]).at(ipt)->SetBinContent(ieta+1, vns/reso_2sub);
                  (hVnS_eta_res[icent][iside][ihar]).at(ipt)->SetBinError  (ieta+1, vns/reso_2sub*sqrt( pow(vns_e/vns,2) + pow(reso_2sub_e/reso_2sub,2)));
               }
            }

            for (int ieta=0; ieta<NETA+2; ieta++) {
               for (int ipt=0; ipt<NPT+1; ipt++) {
                  if (ieta>NETA-1 && ipt == NPT) continue;
                  double vn   = (hVn[icent][iside][ihar]).at(ipt).at(ieta)->GetMean();
                  double vn_e = (hVn[icent][iside][ihar]).at(ipt).at(ieta)->GetMeanError();
                  (hVn_cent_raw[iside][ihar]).at(ipt).at(ieta)->SetBinContent(icent+1, vn);      
                  (hVn_cent_raw[iside][ihar]).at(ipt).at(ieta)->SetBinError  (icent+1, vn_e);      
                  (hVn_cent_res[iside][ihar]).at(ipt).at(ieta)->SetBinContent(icent+1, vn/reso_2sub);
                  (hVn_cent_res[iside][ihar]).at(ipt).at(ieta)->SetBinError  (icent+1, vn/reso_2sub*sqrt( pow(vn_e/vn,2) + pow(reso_2sub_e/reso_2sub,2)));
                  
                  double vns   = (hVnS[icent][iside][ihar]).at(ipt).at(ieta)->GetMean();
                  double vns_e = (hVnS[icent][iside][ihar]).at(ipt).at(ieta)->GetMeanError();
                  (hVnS_cent_raw[iside][ihar]).at(ipt).at(ieta)->SetBinContent(icent+1, vns);      
                  (hVnS_cent_raw[iside][ihar]).at(ipt).at(ieta)->SetBinError  (icent+1, vns_e);      
                  (hVnS_cent_res[iside][ihar]).at(ipt).at(ieta)->SetBinContent(icent+1, vns/reso_2sub);
                  (hVnS_cent_res[iside][ihar]).at(ipt).at(ieta)->SetBinError  (icent+1, vns/reso_2sub*sqrt( pow(vns_e/vns,2) + pow(reso_2sub_e/reso_2sub,2)));
               }
            }
         
         }
      }
   }

}  


void extractor_flat::FillTrackVn(std::vector<std::vector<double> > PsiF) {

   ntrkQ = 0;
   for(int j=0; j<trk_n; j++){
      if(trk_Quality2[j] < 3) continue;
      ntrkQ++;

      double ptgev = trk_pt[j]/1000.0;
      int ptBin  = get_ptBin(ptgev);       if (ptBin  < 0) continue;
      int etaBin = get_etaBin(trk_eta[j]); if (etaBin < 0) continue;
      int sign_eta = 1;
      if (trk_eta[j] < 0) sign_eta = -1;

      for (int iside=0; iside<NSIDE+1; iside++) {
         for (int ihar=0; ihar<NHAR; ihar++) {
            double corr1 = cos( (ihar+1)*(trk_phi0_wrt_PV[j]-PsiF.at(iside).at(0)) );
            double corr2 = sin( (ihar+1)*(trk_phi0_wrt_PV[j]-PsiF.at(iside).at(0)) );
            
            (hVn[centb][iside][ihar]).at(ptBin).at(etaBin)->Fill(corr1);
            (hVn[centb][iside][ihar]).at(NPT  ).at(etaBin)->Fill(corr1);
            (hVn[centb][iside][ihar]).at(ptBin).at(NETA  )->Fill(corr1*sign_eta);
            (hVn[centb][iside][ihar]).at(ptBin).at(NETA+1)->Fill(corr1);
            
            (hVnS[centb][iside][ihar]).at(ptBin).at(etaBin)->Fill(corr2);
            (hVnS[centb][iside][ihar]).at(NPT  ).at(etaBin)->Fill(corr2);
            (hVnS[centb][iside][ihar]).at(ptBin).at(NETA  )->Fill(corr2*sign_eta);
            (hVnS[centb][iside][ihar]).at(ptBin).at(NETA+1)->Fill(corr2);
         }
      }
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

      //-------------------------------------------------------------------------------------------------
      centb = get_centb(Fcal_Et);
      if (centb < 0) continue;
      nevents++;

      std::vector<std::vector<double> > Qx;
      std::vector<std::vector<double> > Qy;
      std::vector<double> Qw;
      get_qVec(Qx, Qy, Qw);

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

      //-------------------------------------------------------------------------------------------------
      centb = get_centb(Fcal_Et);
      if (centb < 0) continue;

      sampSel[centb]++;
      if (sampSel[centb]>NSAMP) sampSel[centb] = NSAMP;

      std::vector<std::vector<double> > Qx;
      std::vector<std::vector<double> > Qy;
      std::vector<double> Qw;
      get_qVec(Qx, Qy, Qw);
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
      std::cout << "discFourrier ERROR : Empty histogram " << h1->GetName() << std::endl;
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
      p += cos( (ihar)*phi) * binval;
   }
   double SumAvg = Sum/Nbins;
   //double ErrAvg = sqrt(Sum_err)/Nbins;
   p /=Sum;

   for(int phibin=1; phibin<=Nbins;phibin++){//Find Errors
      float binerr = h1->GetBinError(phibin);
      float phi    = h1->GetBinCenter(phibin);
      ep  += pow((cos((ihar)*phi)-p)*binerr,2);
   }

   //Correct for Bin Width
   double corr = (ihar)*bwid/sin((ihar)*bwid);
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


int extractor_flat::get_centb(double et) {

   centb = -1;
   int icent;
   for(icent=0;icent<97;icent++){
      if(et>centcuts[icent]) break;
   }

   centb_1p = icent;
   centb = icent/4.99999999999;
   return centb;
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


void extractor_flat::SetUpAll() {

   set_ptBins ();
   set_etaBins();

   char name[500];

   for (int iside=0; iside<NSIDE; iside++) {
      sprintf(name, "h_eTot_side%d", iside);
      h_eTot[iside] = new TH1D(name, "", 4000, -0.5, 3999.5); h_eTot[iside]->Sumw2();
      sprintf(name, "h_ePix_side%d", iside);
      h_ePix[iside] = new TH1D(name, "", 4000, -0.5, 3999.5); h_ePix[iside]->Sumw2();
   
      sprintf(name, "h_et_eTot_side%d", iside);
      hEt_eTot[iside] = new TH2D(name, "", 500, 0, 5000, 400, 0, 4000);
      sprintf(name, "h_et_ePix_side%d", iside);
      hEt_ePix[iside] = new TH2D(name, "", 500, 0, 5000, 200, 0, 2000);
      
      sprintf(name, "h_ntrk_eTot_side%d", iside);
      hNtrk_eTot[iside] = new TH2D(name, "", 400, 0, 4000, 400, 0, 4000);
      sprintf(name, "h_ntrk_ePix_side%d", iside);
      hNtrk_ePix[iside] = new TH2D(name, "", 400, 0, 4000, 200, 0, 2000);
   }
  
   hEt      = new TH1D("hEt", "", 5000, 0, 5000);
   hNtrk    = new TH1D("hNtrk", "", 5000, 0, 5000);
   hCent    = new TH1D("hCent", "", NCENT, -0.5, NCENT-0.5);
   hEt_ntrk = new TH2D("h_et_ntrk", "", 400, 0, 4000, 400, 0, 4000);

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
      for (int iside=0; iside<NSIDE+1; iside++) {
         for (int ihar=0; ihar<NHAR; ihar++) {
            for (int ipt=0; ipt<NPT+1; ipt++) {
               std::vector<TH1D*> hTemp_vec1(NETA+2, (TH1D*)0);
               std::vector<TH1D*> hTemp_vec2(NETA+2, (TH1D*)0);

               for (int ieta=0; ieta<NETA+2; ieta++) {
                  if (ipt==NPT && ieta >= NETA) continue;
                  sprintf(name, "hV%d_cent%d_side%d_pt%d_eta%d", ihar+1, icent, iside, ipt, ieta);
                  hTemp_vec1.at(ieta) = new TH1D(name, "", 220, -1.1, 1.1);
                  sprintf(name, "hV%dS_cent%d_side%d_pt%d_eta%d", ihar+1, icent, iside, ipt, ieta);
                  hTemp_vec2.at(ieta) = new TH1D(name, "", 220, -1.1, 1.1);
               }
               (hVn [icent][iside][ihar]).push_back( hTemp_vec1);
               (hVnS[icent][iside][ihar]).push_back( hTemp_vec2);
            }
         }
      }
   }

   for (int icent=0; icent<NCENT; icent++) {
      for (int iside=0; iside<NSIDE+1; iside++) {
         for (int ihar=0; ihar<NHAR; ihar++) {
            for (int ipt=0; ipt<NPT+1; ipt++) {
               sprintf(name, "hV%d_eta_raw_cent%d_side%d_pt%d", ihar+1, icent, iside, ipt);
               (hVn_eta_raw[icent][iside][ihar]).push_back(new TH1D(name, "", NETA, &diffEtaBins[0]));
               sprintf(name, "hV%d_eta_res_cent%d_side%d_pt%d", ihar+1, icent, iside, ipt);
               (hVn_eta_res[icent][iside][ihar]).push_back(new TH1D(name, "", NETA, &diffEtaBins[0]));
               
               sprintf(name, "hV%dS_eta_raw_cent%d_side%d_pt%d", ihar+1, icent, iside, ipt);
               (hVnS_eta_raw[icent][iside][ihar]).push_back(new TH1D(name, "", NETA, &diffEtaBins[0]));
               sprintf(name, "hV%dS_eta_res_cent%d_side%d_pt%d", ihar+1, icent, iside, ipt);
               (hVnS_eta_res[icent][iside][ihar]).push_back(new TH1D(name, "", NETA, &diffEtaBins[0]));
            }
         }
      }
   }

   for (int icent=0; icent<NCENT; icent++) {
      for (int iside=0; iside<NSIDE+1; iside++) {
         for (int ihar=0; ihar<NHAR; ihar++) {
            for (int ieta=0; ieta<NETA+2; ieta++) {
               sprintf(name, "hV%d_pt_raw_cent%d_side%d_eta%d", ihar+1, icent, iside, ieta);
               (hVn_pt_raw[icent][iside][ihar]).push_back( new TH1D(name, "", NPT, &diffPtBins[0]) );
               sprintf(name, "hV%d_pt_res_cent%d_side%d_eta%d", ihar+1, icent, iside, ieta);
               (hVn_pt_res[icent][iside][ihar]).push_back( new TH1D(name, "", NPT, &diffPtBins[0]) );
               
               sprintf(name, "hV%dS_pt_raw_cent%d_side%d_eta%d", ihar+1, icent, iside, ieta);
               (hVnS_pt_raw[icent][iside][ihar]).push_back( new TH1D(name, "", NPT, &diffPtBins[0]) );
               sprintf(name, "hV%dS_pt_res_cent%d_side%d_eta%d", ihar+1, icent, iside, ieta);
               (hVnS_pt_res[icent][iside][ihar]).push_back( new TH1D(name, "", NPT, &diffPtBins[0]) );
            }
         }
      }
   }

   for (int iside=0; iside<NSIDE+1; iside++) {
      for (int ihar=0; ihar<NHAR; ihar++) {
         for (int ipt=0; ipt<NPT+1; ipt++) {

            std::vector<TH1D*> hTemp_vec1(NETA+2, (TH1D*)0);
            std::vector<TH1D*> hTemp_vec2(NETA+2, (TH1D*)0);
            std::vector<TH1D*> hTemp_vec3(NETA+2, (TH1D*)0);
            std::vector<TH1D*> hTemp_vec4(NETA+2, (TH1D*)0);

            for (int ieta=0; ieta<NETA+2; ieta++) {
               if (ipt==NPT && ieta >= NETA) continue;
               sprintf(name, "hV%d_cent_raw_side%d_pt%d_eta%d", ihar+1, iside, ipt, ieta);
               hTemp_vec1.at(ieta) = new TH1D(name, "", NCENT, 0, 100);
               sprintf(name, "hV%d_cent_res_side%d_pt%d_eta%d", ihar+1, iside, ipt, ieta);
               hTemp_vec2.at(ieta) = new TH1D(name, "", NCENT, 0, 100);
               
               sprintf(name, "hV%dS_cent_raw_side%d_pt%d_eta%d", ihar+1, iside, ipt, ieta);
               hTemp_vec3.at(ieta) = new TH1D(name, "", NCENT, 0, 100);
               sprintf(name, "hV%dS_cent_res_side%d_pt%d_eta%d", ihar+1, iside, ipt, ieta);
               hTemp_vec4.at(ieta) = new TH1D(name, "", NCENT, 0, 100);
            }
            (hVn_cent_raw [iside][ihar]).push_back(hTemp_vec1);
            (hVn_cent_res [iside][ihar]).push_back(hTemp_vec2);
            (hVnS_cent_raw[iside][ihar]).push_back(hTemp_vec3);
            (hVnS_cent_res[iside][ihar]).push_back(hTemp_vec4);
         }
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


void extractor_flat::SaveHistos() {
   std::cout << "Saving" << std::endl;

   fOut->cd();
   fOut->Write();
   std::cout << "Saving Finished" << std::endl;
   std::cout << " total events: " << nevents << std::endl;

}

extractor_flat::extractor_flat(string filelist) {
   TChain* chain     = new TChain("HeavyIonD3PD","");

   char fname[400];
   ifstream lis(filelist.c_str());
   int cnt=0;
   while(!lis.eof()){
      string filename;
      lis >> filename;
      sprintf(fname,"%s",filename.c_str());
      cout << fname << endl;
      if(!filename.empty()) {
         chain    ->Add(fname);
      }
      cnt++;
      if(cnt>1000000) {cout<<"Too Many Files"<<endl;break;}
   }

   Init(chain);
}

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

extractor_flat::~extractor_flat() {
   //if (!fChain) return;
   //delete fChain->GetCurrentFile();
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

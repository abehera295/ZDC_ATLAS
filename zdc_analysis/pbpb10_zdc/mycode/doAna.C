#include <iostream>
#include <vector>
#include <sstream>
#include <iterator>
#include <algorithm>
#include <TFile.h>
#include <TH2.h>
#include <TH2.h>
#include <TF1.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TCanvas.h>

std::vector<double> discFourrier(const TH1D* h1, int ihar, TFile* f1=NULL);
static double PI;

void doAna() {

   PI = acos(-1.0);

   static const int NSIDE = 2;
   static const int NCENT = 20;
   static const int NHAR  = 3;
   int NPT ;
   int NETA;

   std::vector<float> diffPtBins;
   std::vector<float> diffEtaBins;

   float diffPtArr[] = {0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 7.0, 9.0, 12.0};
   diffPtBins.assign(diffPtArr, diffPtArr + sizeof(diffPtArr)/sizeof(float));
   NPT = diffPtBins.size()-1;

   float diffEtaArr[] = {-2.5, -2.0, -1.5, -1.0, -0.5, -0.25, 0.0, 0.25, 0.5, 1.0, 1.5, 2.0, 2.5};
   diffEtaBins.assign(diffEtaArr, diffEtaArr + sizeof(diffEtaArr)/sizeof(float));
   NETA = diffEtaBins.size()-1;

   char name [500];
   char name1[500];

   TFile* fIn  = new TFile("Output_zdcPix_flat_wtrk.root");

   //---------------------------------------------------------------------------------------------------------------------------------//
   TH1D* hDphi_fg [NCENT][NHAR];
   TH1D* hDphi_bg [NCENT][NHAR];
   TH1D* hDphi_rat[NCENT][NHAR];

   TH1D* hRes      [NCENT][NHAR];
   TH1D* hRes_recn [NCENT][NHAR];
   TH1D* hRes_flat [NCENT][NHAR];

   TH1D* hResCent      [NHAR];
   TH1D* hResCent_recn [NHAR];
   TH1D* hResCent_flat [NHAR];
   TH1D* hResCent_mix  [NHAR];

   std::vector<std::vector<TH1D*> > hVn         [NCENT][NSIDE+1][NHAR];
   std::vector<TH1D*>               hVn_eta_raw [NCENT][NSIDE+1][NHAR];
   std::vector<TH1D*>               hVn_eta_res [NCENT][NSIDE+1][NHAR];
   std::vector<TH1D*>               hVn_pt_raw  [NCENT][NSIDE+1][NHAR];
   std::vector<TH1D*>               hVn_pt_res  [NCENT][NSIDE+1][NHAR];
   std::vector<std::vector<TH1D*> > hVn_cent_raw[NSIDE+1][NHAR];
   std::vector<std::vector<TH1D*> > hVn_cent_res[NSIDE+1][NHAR];
   //---------------------------------------------------------------------------------------------------------------------------------//

   //---------------------------------------------------------------------------------------------------------------------------------//
   for (int icent=0; icent<NCENT; icent++) {
      for (int ihar=0; ihar<NHAR; ihar++) {

         sprintf(name , "h_Res_cent%d_har%d", icent, ihar);
         hRes[icent][ihar] = (TH1D*)fIn->Get(name);

         sprintf(name, "h_Res_recn_cent%d_har%d", icent, ihar);
         hRes_recn[icent][ihar] = (TH1D*)fIn->Get(name);

         sprintf(name, "h_Res_flat_cent%d_har%d", icent, ihar);
         hRes_flat[icent][ihar] = (TH1D*)fIn->Get(name);
      }
   }

   for (int icent=0; icent<NCENT; icent++) {
      for (int ihar=0; ihar<NHAR; ihar++) {
         sprintf(name, "h_Dphi_fg_cent%d_har%d", icent, ihar);
         hDphi_fg[icent][ihar] = (TH1D*)fIn->Get(name);
         sprintf(name, "h_Dphi_bg_cent%d_har%d", icent, ihar);
         hDphi_bg[icent][ihar] = (TH1D*)fIn->Get(name);
      }
   }

   for (int icent=0; icent<NCENT; icent++) {
      for (int iside=0; iside<NSIDE+1; iside++) {
         for (int ihar=0; ihar<NHAR; ihar++) {
            for (int ipt=0; ipt<NPT+1; ipt++) {
               std::vector<TH1D*> hTemp_vec1(NETA+2, (TH1D*)0);

               for (int ieta=0; ieta<NETA+2; ieta++) {
                  if (ipt==NPT && ieta >= NETA) continue;
                  sprintf(name, "hV%d_cent%d_side%d_pt%d_eta%d", ihar+1, icent, iside, ipt, ieta);
                  hTemp_vec1.at(ieta) = (TH1D*)fIn->Get(name);
                  hTemp_vec1.at(ieta)->Sumw2();
               }
               (hVn[icent][iside][ihar]).push_back( hTemp_vec1);
            }
         }
      }
   }
   //---------------------------------------------------------------------------------------------------------------------------------//


   TFile* fOut = new TFile("Output_zdcPix_flat_wtrk_doAna.root", "recreate");
   //---------------------------------------------------------------------------------------------------------------------------------//
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

   for (int icent=0; icent<NCENT; icent++) {
      for (int iside=0; iside<NSIDE+1; iside++) {
         for (int ihar=0; ihar<NHAR; ihar++) {
            for (int ipt=0; ipt<NPT+1; ipt++) {
               sprintf(name, "hV%d_eta_raw_cent%d_side%d_pt%d", ihar+1, icent, iside, ipt);
               (hVn_eta_raw[icent][iside][ihar]).push_back(new TH1D(name, "", NETA, &diffEtaBins[0]));
               sprintf(name, "hV%d_eta_res_cent%d_side%d_pt%d", ihar+1, icent, iside, ipt);
               (hVn_eta_res[icent][iside][ihar]).push_back(new TH1D(name, "", NETA, &diffEtaBins[0]));
            }
         }
      }
   }

   for (int icent=0; icent<NCENT; icent++) {
      for (int iside=0; iside<NSIDE+1; iside++) {
         for (int ihar=0; ihar<NHAR; ihar++) {
            for (int ieta=0; ieta<NETA+2; ieta++) {
               sprintf(name, "hV%d_pt_raw_cent%d_side%d_pt%d", ihar+1, icent, iside, ieta);
               (hVn_pt_raw[icent][iside][ihar]).push_back( new TH1D(name, "", NPT, &diffPtBins[0]) );
               sprintf(name, "hV%d_pt_res_cent%d_side%d_pt%d", ihar+1, icent, iside, ieta);
               (hVn_pt_res[icent][iside][ihar]).push_back( new TH1D(name, "", NPT, &diffPtBins[0]) );
            }
         }
      }
   }

   for (int iside=0; iside<NSIDE+1; iside++) {
      for (int ihar=0; ihar<NHAR; ihar++) {
         for (int ipt=0; ipt<NPT+1; ipt++) {
            std::vector<TH1D*> hTemp_vec1(NETA+2, (TH1D*)0);
            std::vector<TH1D*> hTemp_vec2(NETA+2, (TH1D*)0);

            for (int ieta=0; ieta<NETA+2; ieta++) {
               if (ipt==NPT && ieta >= NETA) continue;
               sprintf(name, "hV%d_cent_raw_side%d_pt%d_eta%d", ihar+1, iside, ipt, ieta);
               hTemp_vec1.at(ieta) = new TH1D(name, "", NCENT, 0, 100);
               sprintf(name, "hV%d_cent_res_side%d_pt%d_eta%d", ihar+1, iside, ipt, ieta);
               hTemp_vec2.at(ieta) = new TH1D(name, "", NCENT, 0, 100);
            }
            (hVn_cent_raw[iside][ihar]).push_back(hTemp_vec1);
            (hVn_cent_res[iside][ihar]).push_back(hTemp_vec2);
         }
      }
   }
   //---------------------------------------------------------------------------------------------------------------------------------//

   //---------------------------------------------------------------------------------------------------------------------------------//
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

            double reso_2sub   = hResCent_mix[ihar]->GetBinContent(icent+1);
            double reso_2sub_e = hResCent_mix[ihar]->GetBinError  (icent+1);
            if (iside == NSIDE) {
               reso_2sub   = sqrt(2*reso_2sub);
               reso_2sub_e = reso_2sub_e/reso_2sub;
            }

            for (int ieta=0; ieta<NETA+2; ieta++) {
               for (int ipt=0; ipt<NPT; ipt++) {
                  double vn   = (hVn[icent][iside][ihar]).at(ipt).at(ieta)->GetMean();
                  double vn_e = (hVn[icent][iside][ihar]).at(ipt).at(ieta)->GetMeanError();
                  (hVn_pt_raw[icent][iside][ihar]).at(ieta)->SetBinContent(ipt+1, vn  );
                  (hVn_pt_raw[icent][iside][ihar]).at(ieta)->SetBinError  (ipt+1, vn_e);
                  (hVn_pt_res[icent][iside][ihar]).at(ieta)->SetBinContent(ipt+1, vn/reso_2sub);
                  (hVn_pt_res[icent][iside][ihar]).at(ieta)->SetBinError  (ipt+1, vn/reso_2sub*sqrt( pow(vn_e/vn,2) + pow(reso_2sub_e/reso_2sub,2)));
               }
            }

            for (int ipt=0; ipt<NPT+1; ipt++) {
               for (int ieta=0; ieta<NETA; ieta++) {
                  double vn   = (hVn[icent][iside][ihar]).at(ipt).at(ieta)->GetMean();
                  double vn_e = (hVn[icent][iside][ihar]).at(ipt).at(ieta)->GetMeanError();
                  (hVn_eta_raw[icent][iside][ihar]).at(ipt)->SetBinContent(ieta+1, vn  );
                  (hVn_eta_raw[icent][iside][ihar]).at(ipt)->SetBinError  (ieta+1, vn_e);
                  (hVn_eta_res[icent][iside][ihar]).at(ipt)->SetBinContent(ieta+1, vn/reso_2sub);
                  (hVn_eta_res[icent][iside][ihar]).at(ipt)->SetBinError  (ieta+1, vn/reso_2sub*sqrt( pow(vn_e/vn,2) + pow(reso_2sub_e/reso_2sub,2)));
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
               }
            }

         }
      }
   }

   fOut->cd(0);
   fOut->Write();
}


std::vector<double> discFourrier(const TH1D* h1, int ihar, TFile* f1) {

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

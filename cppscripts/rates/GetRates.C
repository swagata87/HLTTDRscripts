#define GetRates_cxx
#include "GetRates.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <math.h>       /* sqrt */
#include <algorithm>    // std::sort
#include <vector>       // std::vector

void GetRates::SingleEle(float given_pt,float given_absEta_low=0.0,float given_absEta_high=4.0)
{
  ///
   std::cout << "Running SingleEle" << std::endl; 
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   float myweight=0;
   float sum2_myweight=0;
   float nEvt_passed=0.;
   float nEvt_passed_wt=0.;

   bool WP_tight = 1;
   std::cout << "WP_tight " << WP_tight << std::endl;
   if (WP_tight) {
     std::cout << "running WP tight, 70% eff" << std::endl;
   }
   else {
     std::cout << "running WP loose, 80% eff" << std::endl;
   }
   float pt_cut= given_pt; //32.0;
   std::cout << "pt cut " << pt_cut << std::endl;

   float abs_eta_cut_high = given_absEta_high;
   float abs_eta_cut_low = given_absEta_low;
   std::cout << "eta cut is:"  <<  abs_eta_cut_low << " < |eta| < " << abs_eta_cut_high << std::endl;

   Long64_t nbytes = 0, nb = 0;

   /////////
   ////EVENT LOOP
   /////////
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     myweight=weightV2;

     ////
     /////LOOP ele
     ////
      int nEle_passed=0;
      for (int i=0; i<nrEgs; i++) {
	//// Barrel cuts
	float ecaliso_cut_EB= (WP_tight==1) ? 4.8 : 7.5; 
	float hcaliso_cut_EB= (WP_tight==1) ? 9.5 : 13.0; 
	if ( ( fabs(eg_eta[i]) > 0.8 ) && ( fabs(eg_eta[i]) < 1.479 )  ) {
	  hcaliso_cut_EB = (WP_tight==1) ? 15.0 : 18.0; 
	}
	float pms2_cut_EB= (WP_tight==1) ? 42.0 : 55.0; 
	float sieie_cut_EB= (WP_tight==1) ? 0.0115 : 0.0129; 
	float hoe_EB = (eg_hcalHForHoverE[i])/(eg_energy[i]);
	float hoe_cut_EB = (WP_tight==1) ? 0.17 : 0.175; 
	float ooemoop_cut_EB = (WP_tight==1) ? 0.035 : 0.04;
	if ( ( fabs(eg_eta[i]) > 0.8 ) && ( fabs(eg_eta[i]) < 1.479 )  ) {
	  ooemoop_cut_EB = (WP_tight==1) ? 0.08: 0.08;
	}
	float deta_cut_EB = (WP_tight==1) ? 0.003 : 0.003;
	if ( ( fabs(eg_eta[i]) > 0.8 ) && ( fabs(eg_eta[i]) < 1.479 )  ) {
	  deta_cut_EB = (WP_tight==1) ? 0.009 : 0.009;
	}
	float dphi_cut_EB = (WP_tight==1) ? 0.02 : 0.02;
	if ( ( fabs(eg_eta[i]) > 0.8 ) && ( fabs(eg_eta[i]) < 1.479 )  ) {
	  dphi_cut_EB = (WP_tight==1) ? 0.09 : 0.09;
	}
	float npix_cut_EB = (WP_tight==1) ? 2 : 2;
	float chi2_cut_EB = (WP_tight==1) ? 50.0 : 50.0;
	float trkisohlt_cut_EB = (WP_tight==1) ? 2.0 : 2.5;
	float trkisol1_cut_EB = (WP_tight==1) ? 4.0 : 5.5; 
	if ( ( fabs(eg_eta[i]) > 0.8 ) && ( fabs(eg_eta[i]) < 1.479 )  ) {
	  trkisol1_cut_EB = (WP_tight==1) ? 8.0 : 8.0;
	}
	////end of barrel cuts
	
	////endcap cuts
	float hoe_EE = (eg_hgcalHForHoverE[i])/(eg_energy[i]);
	float hoe_cut_EE = (WP_tight==1) ? 0.15 + (5.0/(eg_energy[i])) : 0.15 + (5.0/(eg_energy[i]));
	float vv_cut_EE = (WP_tight==1) ? 0.8*0.8 : 0.85*0.85; 
	float ww_cut_EE = (WP_tight==1) ? 8*8 : 8.5*8.5; 
	float hgcaliso_cut_EE = (WP_tight==1) ? 130.0 : 150.0;
	if ( ( fabs(eg_eta[i]) > 2.0 ) ) {
	  hgcaliso_cut_EE = (WP_tight==1) ? 340.0 : 350.0;
	}
	float pms2_cut_EE= (WP_tight==1) ? 65.0 : 75.0;
	float ooemoop_cut_EE = (WP_tight==1) ? 0.01 : 0.04;
	float deta_cut_EE = (WP_tight==1) ? 0.003 : 0.004;
	float dphi_cut_EE = (WP_tight==1) ? 0.02 : 0.04;
	float npix_cut_EE = (WP_tight==1) ? 2 : 2;
	float chi2_cut_EE = (WP_tight==1) ? 50.0 : 50.0;
	float trkisohlt_cut_EE = (WP_tight==1) ? 1.5 : 2.2; 
	float trkisol1_cut_EE = (WP_tight==1) ? 5.5 : 5.5;
	/// end of endcap cuts
	
	/// for variables common for barrel and endcap, rename in a generic way that works for both barrel+endcap
	float pms2_cut = 9999;
	float hoe = 0;
	float hoe_cut = 9999;
	float ooemoop_cut = 9999;
	float deta_cut = 9999;
	float dphi_cut = 9999;
	float npix_cut = 0;
	float chi2_cut = 9999;
	float trkisohlt_cut = 9999;
	float trkisol1_cut = 9999;
	

	if ( fabs(eg_eta[i]) < 1.479 ) {
	  pms2_cut = pms2_cut_EB;
	  hoe = hoe_EB;
	  hoe_cut = hoe_cut_EB;
	  ooemoop_cut = ooemoop_cut_EB;
	  deta_cut = deta_cut_EB;
	  dphi_cut = dphi_cut_EB;
	  npix_cut = npix_cut_EB;
	  chi2_cut = chi2_cut_EB;
	  trkisohlt_cut = trkisohlt_cut_EB;
	  trkisol1_cut = trkisol1_cut_EB;
	 
	}
	else  {
	  pms2_cut = pms2_cut_EE;
	  hoe = hoe_EE;
	  hoe_cut = hoe_cut_EE;
	  ooemoop_cut = ooemoop_cut_EE;
	  deta_cut = deta_cut_EE;
	  dphi_cut = dphi_cut_EE;
	  npix_cut = npix_cut_EE;
	  chi2_cut = chi2_cut_EE;
	  trkisohlt_cut = trkisohlt_cut_EE;
	  trkisol1_cut = trkisol1_cut_EE;
	
	}

	////define L1 pass/fail
	bool passL1 = false;
	if (  fabs(eg_eta[i]) > 2.4 ) {
	  passL1=true;
	}
	if ( (fabs(eg_eta[i])<2.4) &&  (eg_l1pho_et[i]>0) && (eg_l1pho_etThres[i]>51.) && (eg_l1pho_passQual[i]) ) {
	  passL1=true;
	}
	if ( (fabs(eg_eta[i])<2.4) && (eg_l1ele_et[i]>0) && (eg_l1ele_etThresNonIso[i]>36.) && (eg_l1ele_passQual[i]) ) {
	  passL1=true;
	}
	if ( (fabs(eg_eta[i])<2.4) && (eg_l1ele_et[i]>0) && (eg_l1ele_etThresIso[i]>28.) && (eg_l1ele_passQual[i]) && (eg_l1ele_passIsol[i]) ) {
	  passL1=true;
	}

	if ( 
	    (myweight>0.) && 
	    ( fabs(eg_eta[i]) >= abs_eta_cut_low ) && ( fabs(eg_eta[i]) < abs_eta_cut_high ) && 
	    (eg_et[i]>pt_cut) &&  
	    (passL1) && 
	    (hoe<hoe_cut) &&  
	    (eg_sigma2vv[i]<vv_cut_EE) &&  
	    (eg_sigma2ww[i]<ww_cut_EE) &&
	    (eg_hgcaliso_layerclus[i]<hgcaliso_cut_EE) && 
	    (eg_pms2[i]<pms2_cut) && 
	    (eg_invEInvP[i]<ooemoop_cut) &&
	    (eg_trkDEtaSeed[i]<deta_cut) && 
	    (eg_trkDPhi[i]<dphi_cut) && 
	    (eg_nLayerIT[i]>npix_cut) && 
	    (eg_normChi2[i]<chi2_cut) &&  
	    (eg_hltisov6[i]<trkisohlt_cut) && 
	    (eg_l1iso[i]<trkisol1_cut) &&
	    (eg_ecaliso[i]<ecaliso_cut_EB ) && 
	    (eg_hcalPFIsol_default[i]<hcaliso_cut_EB) && 
	    (eg_sigmaIEtaIEta[i]<sieie_cut_EB) 
	   
	     ) {
	  nEle_passed=nEle_passed+1;
	}


      } ////end of ele loop
      
      if ( nEle_passed>0 ) {  
	nEvt_passed = nEvt_passed+1;
	nEvt_passed_wt = nEvt_passed_wt+myweight;
	sum2_myweight=sum2_myweight+(myweight*myweight);

      }
   }

   std::cout << "nEvt_passed " << nEvt_passed << std::endl;
   std::cout << "rate " << nEvt_passed_wt ; 
   std::cout << " +/- " << sqrt(sum2_myweight) << std::endl;

}

/////////////////////  Single Photon nonIsolated /////////////
void GetRates::SinglePhoNonIso(float given_pt)
{
  std::cout << "running Single Photon NonIsolated" << std::endl;
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();

  float myweight=0;
  float sum2_myweight=0;
  float nEvt_passed=0.;
  float nEvt_passed_wt=0.;

  float pt_cut=given_pt;
  std::cout << "pt cut " << pt_cut << std::endl;
  Long64_t nbytes = 0, nb = 0;

  /////////
  ////EVENT LOOP
  /////////
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    myweight=weightV2;

    ////
    /////LOOP SC
    ////
    int nPho_Passed=0;
    for (int i=0; i<nrEgs; i++) {

      //// Barrel cuts
      // float sieie_cut_EB= 0.013;
      float hoe_EB = (eg_hcalHForHoverE[i])/(eg_energy[i]);
      float hoe_cut_EB = 0.30;
       ////end of barrel cuts
      
      ////endcap cuts
      float hoe_EE = (eg_hgcalHForHoverE[i])/(eg_energy[i]);
      float hoe_cut_EE = 0.30;
      //float vv_cut_EE = 0.9*0.9;
      // float ww_cut_EE = 9*9;
      /// end of endcap cuts
      
      /// for variables common for barrel and endcap, rename in a generic way that works for both barrel+endcap
      float hoe = 0;
      float hoe_cut = 9999;
      if ( fabs(eg_eta[i]) < 1.479 ) {
	hoe = hoe_EB;
	hoe_cut = hoe_cut_EB;
      }
      else  {
	hoe = hoe_EE;
	hoe_cut = hoe_cut_EE;
      }

      ////define L1 pass/fail
      bool passL1 = false;
      if (  fabs(eg_eta[i]) > 2.4 ) {
	passL1=true;
      }
      if ( (fabs(eg_eta[i])<2.4) &&  (eg_l1pho_et[i]>0) && (eg_l1pho_etThres[i]>51.) && (eg_l1pho_passQual[i]) ) {
	passL1=true;
      }

      if ( 
	  (myweight>0.) && 
	  (eg_et[i]>pt_cut) &&  
	  (passL1) && 
	  (hoe<hoe_cut) // &&  
	  // (eg_sigma2vv[i]<vv_cut_EE) &&  
	  // (eg_sigma2ww[i]<ww_cut_EE) &&
	  // (eg_sigmaIEtaIEta[i]<sieie_cut_EB) 
	     
	   ) {
	nPho_Passed=nPho_Passed+1;
      }


    } ////end of ele loop
      
    if ( nPho_Passed>0 ) {  
      nEvt_passed = nEvt_passed+1;
      nEvt_passed_wt = nEvt_passed_wt+myweight;
      sum2_myweight=sum2_myweight+(myweight*myweight);

    }
  }
  
  std::cout << "Photon Trigger NonIso" << std::endl;
  std::cout << "nEvt_passed " << nEvt_passed << std::endl;
  std::cout << "rate " << nEvt_passed_wt ; 
  std::cout << " +/- " << sqrt(sum2_myweight) << std::endl;

}

/////////////////Single photon isolated barrel-only /////////////////
void GetRates::SinglePhoIsoEBonly(float given_pt)
{
  std::cout << "running Single Photon Isolated EBonly" << std::endl;
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();

  float myweight=0;
  float sum2_myweight=0;
  float nEvt_passed=0.;
  float nEvt_passed_wt=0.;

  float pt_cut=given_pt;
  std::cout << "pt cut " << pt_cut << std::endl;
  Long64_t nbytes = 0, nb = 0;

  /////////
  ////EVENT LOOP
  /////////
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    myweight=weightV2;

    ////
    /////LOOP SC
    ////
    int nPho_Passed=0;
    for (int i=0; i<nrEgs; i++) {

      //// Barrel cuts
      float sieie_cut_EB= 0.01;
      float ecaliso_cut_EB= 3.0 + (0.02*eg_et[i]);
      float hcaliso_cut_EB= 5.3 + (0.02*eg_et[i]);
      float hoe_EB = (eg_hcalHForHoverE[i])/(eg_energy[i]);
      float hoe_cut_EB = 0.05;
       ////end of barrel cuts

      bool passL1 = false;

      if ( (eg_l1pho_et[i]>0) && (eg_l1pho_etThres[i]>51.) && (eg_l1pho_passQual[i]) ) {
	passL1=true;
      }

      if ( (fabs(eg_eta[i])<2.4) && (eg_l1pho_et[i]>0) && (eg_l1pho_etThresIso[i]>36.) && (eg_l1pho_passQual[i]) && (eg_l1pho_passIsol[i]) ) {
	passL1=true;
      }

      if ( 
	  (myweight>0.) && 
	  (eg_et[i]>pt_cut) &&
	  ( fabs(eg_eta[i]) < 1.479) &&
	  (passL1) && 
	  (hoe_EB<hoe_cut_EB) &&	  
	  (eg_ecaliso[i]<ecaliso_cut_EB ) && 
	  (eg_hcalPFIsol_default[i]<hcaliso_cut_EB) && 
	  (eg_sigmaIEtaIEta[i]<sieie_cut_EB) 
	   ) {
	nPho_Passed=nPho_Passed+1;
      }


    } ////end of ele loop
      
    if ( nPho_Passed>0 ) {  
      nEvt_passed = nEvt_passed+1;
      nEvt_passed_wt = nEvt_passed_wt+myweight;
      sum2_myweight=sum2_myweight+(myweight*myweight);

    }
  }
  
  std::cout << "Photon Trigger Iso EB-only" << std::endl;
  std::cout << "nEvt_passed " << nEvt_passed << std::endl;
  std::cout << "rate " << nEvt_passed_wt ; 
  std::cout << " +/- " << sqrt(sum2_myweight) << std::endl;

}
////

/////// double Ele ////////
void GetRates::DoubleEle(float given_pt=25.0)
{
  ///
   std::cout << "Running DoubleEle" << std::endl; 
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   float myweight=0;
   float sum2_myweight=0;
   float nEvt_passed=0.;
   float nEvt_passed_wt=0.;

   float pt_cut= given_pt; //25.0;
   std::cout << "pt cut " << pt_cut << std::endl;
   Long64_t nbytes = 0, nb = 0;

   /////////
   ////EVENT LOOP
   /////////
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     // std::cout << "\nNEW EVT" << std::endl;
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     myweight=weightV2;

      int nEle_passed_leg_1=0;
      int nEle_passed_leg_2=0;
      std::vector<int> indx_leg_1;
      std::vector<int> indx_leg_2;
      std::vector<int> indx_leg_12;
      indx_leg_12.clear();
      indx_leg_1.clear();
      indx_leg_2.clear();


     ////
     /////LOOP ele
     ////
      for (int i=0; i<nrEgs; i++) {

	//// Barrel cuts // DoubleEle HLT
	float pms2_cut_EB= 55.0; 
	float sieie_cut_EB= 0.013; 
	float hoe_EB = (eg_hcalHForHoverE[i])/(eg_energy[i]);
	float hoe_cut_EB = 0.19; 
	////end of barrel cuts
	
	////endcap cuts // DoubleEle HLT  
	float hoe_EE = (eg_hgcalHForHoverE[i])/(eg_energy[i]);
	float hoe_cut_EE = 0.19; 
	float vv_cut_EE = 0.8*0.8; 
	float pms2_cut_EE= 75.0;
	/// end of endcap cuts
	
	/// for variables common for barrel and endcap, rename in a generic way that works for both barrel+endcap
	float pms2_cut = 9999;
	float hoe = 0;
	float hoe_cut = 9999;

	if ( fabs(eg_eta[i]) < 1.479 ) {
	  pms2_cut = pms2_cut_EB;
	  hoe = hoe_EB;
	  hoe_cut = hoe_cut_EB;
	}
	else  {
	  pms2_cut = pms2_cut_EE;
	  hoe = hoe_EE;
	  hoe_cut = hoe_cut_EE;
	}

	///
	//// Using Double TkElectron 25, 12 & Double StaEG 37,24
	///
	////define L1 pass/fail, high pT leg
	bool passL1_highpt = false;
	if (  fabs(eg_eta[i]) > 2.4 ) {
	  passL1_highpt=true;
	}
	if ( (fabs(eg_eta[i])<2.4) &&  (eg_l1pho_et[i]>0) && (eg_l1pho_etThres[i]>37.) && (eg_l1pho_passQual[i]) ) {
	  passL1_highpt=true;
	}
	if ( (fabs(eg_eta[i])<2.4) && (eg_l1ele_et[i]>0) && (eg_l1ele_etThresNonIso[i]>25.) && (eg_l1ele_passQual[i]) ) {
	  passL1_highpt=true;
	}
	/////////

	////define L1 pass/fail, low pT leg
	bool passL1_lowpt = false;
	if (  fabs(eg_eta[i]) > 2.4 ) {
	  passL1_lowpt=true;
	}
	if ( (fabs(eg_eta[i])<2.4) &&  (eg_l1pho_et[i]>0) && (eg_l1pho_etThres[i]>24.) && (eg_l1pho_passQual[i]) ) {
	  passL1_lowpt=true;
	}
	if ( (fabs(eg_eta[i])<2.4) && (eg_l1ele_et[i]>0) && (eg_l1ele_etThresNonIso[i]>12.) && (eg_l1ele_passQual[i]) ) {
	  passL1_lowpt=true;
	}

	////leg_1
	if ( 
	    (myweight>0.) && 
	    (eg_et[i]>pt_cut) &&  
	    (passL1_highpt) && 
	    (hoe<hoe_cut) &&  
	    (eg_sigma2vv[i]<vv_cut_EE) &&  
	    (eg_pms2[i]<pms2_cut) && 
	    (eg_sigmaIEtaIEta[i]<sieie_cut_EB) 
	   
	     ) {
	  nEle_passed_leg_1=nEle_passed_leg_1+1; // how many ele passing leg1 cuts
	  indx_leg_1.push_back(i); // what are their index
	    }

	////leg_2
	if ( 
	    (myweight>0.) && 
	    (eg_et[i]>pt_cut) &&  
	    (passL1_lowpt) && 
	    (hoe<hoe_cut) &&  
	    (eg_sigma2vv[i]<vv_cut_EE) &&  
	    (eg_pms2[i]<pms2_cut) && 
	    (eg_sigmaIEtaIEta[i]<sieie_cut_EB) 
	   
	     ) {
	  nEle_passed_leg_2=nEle_passed_leg_2+1; //how many ele passing leg2 cuts 
	  indx_leg_2.push_back(i); // what are their index 
	}
      } ////end of ele loop


      /// merge vectors indx_leg_1 & indx_leg_2 into indx_leg_12
      indx_leg_12.insert(indx_leg_12.begin(), indx_leg_1.begin(), indx_leg_1.end());
      indx_leg_12.insert(indx_leg_12.end(), indx_leg_2.begin(), indx_leg_2.end());

      ///sort indx_leg_12
      std::sort(indx_leg_12.begin(), indx_leg_12.end()); 

      //// if an ele passed both leg1 cuts and leg2 cuts then their will be duplicate entries in indx_leg_12
      ///remove duplicate entries from indx_leg_12
      auto last = std::unique(indx_leg_12.begin(), indx_leg_12.end());
      indx_leg_12.erase(last, indx_leg_12.end());

      //there should be at least 1 ele passing leg1 cuts, at least 1 ele passing leg2 cuts, and these 2 ele must not be the same ele// 
      if ( nEle_passed_leg_1>0 && nEle_passed_leg_2>0 && indx_leg_12.size()>1 ) {  
	nEvt_passed = nEvt_passed+1;
	nEvt_passed_wt = nEvt_passed_wt+myweight;
	sum2_myweight=sum2_myweight+(myweight*myweight);
      }
   }

   std::cout << "nEvt_passed " << nEvt_passed << std::endl;
   std::cout << "rate " << nEvt_passed_wt ; 
   std::cout << " +/- " << sqrt(sum2_myweight) << std::endl;

}

///// diphoton trigger
void GetRates::DoublePho(float given_pt1=30.0,float given_pt2=18.0)
{
  ///
   std::cout << "Running DoublePho" << std::endl; 
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   float myweight=0;
   float sum2_myweight=0;
   float nEvt_passed=0.;
   float nEvt_passed_wt=0.;

   float pt_cut1= given_pt1; 
   float pt_cut2= given_pt2; 
   std::cout << "pt cuts " << pt_cut1 << " , " << pt_cut2  << std::endl;
   Long64_t nbytes = 0, nb = 0;

   /////////
   ////EVENT LOOP
   /////////
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     // std::cout << "\nNEW EVT" << std::endl;
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     myweight=weightV2;

      int nPho_passed_leg_1=0;
      int nPho_passed_leg_2=0;
      std::vector<int> indx_leg_1;
      std::vector<int> indx_leg_2;
      std::vector<int> indx_leg_12;
      indx_leg_12.clear();
      indx_leg_1.clear();
      indx_leg_2.clear();

     ////
     /////LOOP pho
     ////
      for (int i=0; i<nrEgs; i++) {
	
	//// Barrel cuts // DoublePho HLT
	float sieie_cut_EB= 0.0113; 
	float hoe_EB = (eg_hcalHForHoverE[i])/(eg_energy[i]);
	float hoe_cut_EB = 0.18; 
	float ecaliso_cut_EB=  8.5;
	float hcaliso_cut_EB= 11;
	if ( ( fabs(eg_eta[i]) > 0.8 ) && ( fabs(eg_eta[i]) < 1.479 )  ) {
	  hcaliso_cut_EB = 12; 
	}
	////end of barrel cuts
	
	////endcap cuts // DoublePho HLT  
	float hoe_EE = (eg_hgcalHForHoverE[i])/(eg_energy[i]);
	float hoe_cut_EE = 0.125; 
	float vv_cut_EE =0.8*0.8; 
	float ww_cut_EE = 8*8;
	float hgcaliso_cut_EE =140; 
	if ( ( fabs(eg_eta[i]) > 2.0 ) ) {
	  hgcaliso_cut_EE = 370; 
	}
	/// end of endcap cuts
	
	/// for variables common for barrel and endcap, rename in a generic way that works for both barrel+endcap
	float hoe = 0;
	float hoe_cut = 9999;

	if ( fabs(eg_eta[i]) < 1.479 ) {
	  hoe = hoe_EB;
	  hoe_cut = hoe_cut_EB;
	}
	else  {
	  hoe = hoe_EE;
	  hoe_cut = hoe_cut_EE;
	}

	//	Double StaEG 37,24 
	//      Double TkIsoPhoton 22, 12
	///
	////define L1 pass/fail, high pT leg
	bool passL1_highpt = false;
	if (  fabs(eg_eta[i]) > 2.4 ) {
	  passL1_highpt=true;
	}
	if ( (fabs(eg_eta[i])<2.4) &&  (eg_l1pho_et[i]>0) && (eg_l1pho_etThres[i]>37.) && (eg_l1pho_passQual[i]) ) {
	  passL1_highpt=true;
	}
	if ( (fabs(eg_eta[i])<2.4) && (eg_l1pho_et[i]>0) && (eg_l1pho_etThresIso[i]>22.) && (eg_l1pho_passQual[i]) && (eg_l1pho_passIsol[i]) ) {
	  passL1_highpt=true;
	}
	/////////

	////define L1 pass/fail, low pT leg
	bool passL1_lowpt = false;
	if (  fabs(eg_eta[i]) > 2.4 ) {
	  passL1_lowpt=true;
	}
	if ( (fabs(eg_eta[i])<2.4) &&  (eg_l1pho_et[i]>0) && (eg_l1pho_etThres[i]>24.) && (eg_l1pho_passQual[i]) ) {
	  passL1_lowpt=true;
	}
	if ( (fabs(eg_eta[i])<2.4) && (eg_l1pho_et[i]>0) && (eg_l1pho_etThresIso[i]>12.) && (eg_l1pho_passQual[i]) && (eg_l1pho_passIsol[i]) ) {
	  passL1_lowpt=true;
	}

	////leg_1
	if ( 
	    (myweight>0.) && 
	    (eg_et[i]>pt_cut1) &&  
	    (passL1_highpt) && 
	    (hoe<hoe_cut) &&  
	    (eg_sigma2ww[i]<ww_cut_EE) &&
	    (eg_hgcaliso_layerclus[i]<hgcaliso_cut_EE) && 
	    (eg_ecaliso[i]<ecaliso_cut_EB ) && 
	    (eg_hcalPFIsol_default[i]<hcaliso_cut_EB) && 
	    (eg_sigma2vv[i]<vv_cut_EE) &&  
	    (eg_sigmaIEtaIEta[i]<sieie_cut_EB) 
	   
	     ) {
	  nPho_passed_leg_1=nPho_passed_leg_1+1; // how many pho passing leg1 cuts
	  indx_leg_1.push_back(i); // what are their index
	    }

	////leg_2
	if ( 
	    (myweight>0.) && 
	    (eg_et[i]>pt_cut2) &&  
	    (passL1_lowpt) && 
	    (hoe<hoe_cut) &&  
	    (eg_sigma2ww[i]<ww_cut_EE) &&
            (eg_hgcaliso_layerclus[i]<hgcaliso_cut_EE) &&
            (eg_ecaliso[i]<ecaliso_cut_EB ) &&
            (eg_hcalPFIsol_default[i]<hcaliso_cut_EB) &&
	    (eg_sigma2vv[i]<vv_cut_EE) &&  
	    (eg_sigmaIEtaIEta[i]<sieie_cut_EB) 
	   
	     ) {
	  nPho_passed_leg_2=nPho_passed_leg_2+1; //how many pho passing leg2 cuts 
	  indx_leg_2.push_back(i); // what are their index 
	}
      } ////end of pho loop

      /// merge vectors indx_leg_1 & indx_leg_2 into indx_leg_12
      indx_leg_12.insert(indx_leg_12.begin(), indx_leg_1.begin(), indx_leg_1.end());
      indx_leg_12.insert(indx_leg_12.end(), indx_leg_2.begin(), indx_leg_2.end());

      ///sort indx_leg_12
      std::sort(indx_leg_12.begin(), indx_leg_12.end()); 

      //// if an pho passed both leg1 cuts and leg2 cuts then their will be duplicate entries in indx_leg_12
      ///remove duplicate entries from indx_leg_12
      auto last = std::unique(indx_leg_12.begin(), indx_leg_12.end());
      indx_leg_12.erase(last, indx_leg_12.end());

      //there should be at least 1 pho passing leg1 cuts, at least 1 pho passing leg2 cuts, and these 2 pho must not be the same pho// 
      if ( nPho_passed_leg_1>0 && nPho_passed_leg_2>0 && indx_leg_12.size()>1 ) {  
	nEvt_passed = nEvt_passed+1;
	nEvt_passed_wt = nEvt_passed_wt+myweight;
	sum2_myweight=sum2_myweight+(myweight*myweight);
      }
   }

   std::cout << "nEvt_passed " << nEvt_passed << std::endl;
   std::cout << "rate " << nEvt_passed_wt ; 
   std::cout << " +/- " << sqrt(sum2_myweight) << std::endl;

}

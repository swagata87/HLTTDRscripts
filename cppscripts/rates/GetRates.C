#define GetRates_cxx
#include "GetRates.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <math.h>       /* sqrt */

void GetRates::Loop()
{
  ///
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   float myweight=0;
   float sum2_myweight=0;
   float nEvt_passed=0.;
   float nEvt_passed_wt=0.;

   bool WP_tight = 0;
   std::cout << "WP_tight " << WP_tight << std::endl;
   if (WP_tight) {
     std::cout << "running WP tight, 70% eff" << std::endl;
   }
   else {
     std::cout << "running WP loose, 80% eff" << std::endl;
   }
   float pt_cut=32.0;
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

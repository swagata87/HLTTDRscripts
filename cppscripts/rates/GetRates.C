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
   float sum_myweight=0;
   float nEvt_passed=0.;
   float nEvt_passed_wt=0.;
   
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
	float ecaliso_cut_EB= 6.0; //4.8;  //6.0;
	float hcaliso_cut_EB= 13.0; //9.5; //13;
	if ( ( fabs(eg_eta[i]) > 0.8 ) && ( fabs(eg_eta[i]) < 1.479 )  ) {
	  hcaliso_cut_EB =18; //15.0;  //18;
	}
	float pms2_cut_EB= 55.0; //42.0; //55.0;
	float sieie_cut_EB=0.013; //0.012; //0.013;
	float hoe_EB = (eg_hcalHForHoverE[i])/(eg_energy[i]);
	float hoe_cut_EB = 0.25; //0.17;  //0.25;
	float ooemoop_cut_EB = 0.04; //0.035;  //0.04;
	if ( ( fabs(eg_eta[i]) > 0.8 ) && ( fabs(eg_eta[i]) < 1.479 )  ) {
	  ooemoop_cut_EB = 0.08;
	}
	float deta_cut_EB = 0.003;
	if ( ( fabs(eg_eta[i]) > 0.8 ) && ( fabs(eg_eta[i]) < 1.479 )  ) {
	  deta_cut_EB = 0.009;
	}
	float dphi_cut_EB = 0.02;
	if ( ( fabs(eg_eta[i]) > 0.8 ) && ( fabs(eg_eta[i]) < 1.479 )  ) {
	  dphi_cut_EB = 0.09;
	}
	float npix_cut_EB = 2;
	float chi2_cut_EB = 50.0;
	float trkisohlt_cut_EB =3.5; //2.4; //3.5;
	float trkisol1_cut_EB =5.5; // 4.0; //5.5;
	if ( ( fabs(eg_eta[i]) > 0.8 ) && ( fabs(eg_eta[i]) < 1.479 )  ) {
	  trkisol1_cut_EB = 8.0;
	}
	////end of barrel cuts

	////endcap cuts
	float hoe_EE = (eg_hgcalHForHoverE[i])/(eg_energy[i]);
	float hoe_cut_EE = 0.15 + (5.0/(eg_energy[i]));
	float vv_cut_EE = 0.9*0.9; //0.8*0.8;  //0.9*0.9;
	float ww_cut_EE = 9*9; //8*8;  //9*9;
	float hgcaliso_cut_EE = 150.0; //130.0; //150.0;
	if ( ( fabs(eg_eta[i]) > 2.0 ) ) {
	  hgcaliso_cut_EE =350; // 340.0;  //350;
	}
	float pms2_cut_EE= 75.0; //65.0; //75.0;
	float ooemoop_cut_EE = 0.04; //0.01; //0.04;
	float deta_cut_EE =0.004;  //0.003; //0.004;
	float dphi_cut_EE = 0.04;  //0.02; //0.04;
	float npix_cut_EE = 2;
	float chi2_cut_EE = 50.0;
	float trkisohlt_cut_EE = 3.0; //2.0; //3.0;
	float trkisol1_cut_EE = 5.5;
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
	     (eg_hcaliso[i]<hcaliso_cut_EB) && 
	     (eg_sigmaIEtaIEta[i]<sieie_cut_EB) 
	     ) {
	  nEle_passed=nEle_passed+1;
	}


      } ////end of ele loop
      
      if ( nEle_passed>0 ) {  
	nEvt_passed = nEvt_passed+1;
	nEvt_passed_wt = nEvt_passed_wt+myweight;
	sum_myweight=sum_myweight+myweight;

      }
   }

   //std::cout << "sum_myweight " << sum_myweight << std::endl;
   std::cout << "nEvt_passed " << nEvt_passed << std::endl;
   float avg_myweight = sum_myweight/nEvt_passed;
   // std::cout << "avg_myweight = sum_myweight/nEvt_passed = " << avg_myweight << std::endl;

   std::cout << "rate " << nEvt_passed_wt << std::endl;   
   std::cout << "rate error " << avg_myweight * sqrt(nEvt_passed) << std::endl;

}

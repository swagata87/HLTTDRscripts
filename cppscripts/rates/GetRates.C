#define GetRates_cxx
#include "GetRates.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <math.h>       /* sqrt */

void GetRates::Loop()
{
//   In a ROOT session, you can do:
//      root> .L GetRates.C
//      root> GetRates t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   int reject_EMenriched=1; // this should be 1 only for QCD STD, otherwise 0
   std::cout << "reject_EMenriched " << reject_EMenriched << std::endl;
   float myweight=0;

   float nEvt_passed=0.;
   
   float pt_cut=32.0;
   std::cout << "pt cut " << pt_cut << std::endl;
   Long64_t nbytes = 0, nb = 0;

   /////////
   ////LOOP
   /////////
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      //std::cout << "\neventnr " << eventnr << std::endl;
      if ( (reject_EMenriched==1) && (path_Gen_QCDEmEnrichingNoBCToEFilter==true) ) {
	//std::cout << "passed path_Gen_QCDEmEnrichingNoBCToEFilter, reject event " << std::endl;
	continue;
      } 
      //std::cout << "weight " << weight << std::endl;
      //std::cout << "nrEgs " << nrEgs << std::endl;
      //std::cout << "eg_energy " << *eg_energy << std::endl;      

      if (weight>0) {
	myweight=weight;
	}
      //      std::cout << "start egamma loop " << std::endl;
      ////
      /////LOOP ele
      ////
      int nEle_passed=0;
      for (int i=0; i<nrEgs; i++) {
	//std::cout << i << std::endl;
	//std::cout << eg_energy[i] << std::endl;

	//// Barrel cuts
	float ecaliso_cut_EB=4.8;  //6.0;
	float hcaliso_cut_EB=9.5; //13;
	if ( ( fabs(eg_eta[i]) > 0.8 ) && ( fabs(eg_eta[i]) < 1.479 )  ) {
	  hcaliso_cut_EB =15.0;  //18;
	}
	float pms2_cut_EB=42.0; //55.0;
	float sieie_cut_EB=0.012; //0.013;
	float hoe_EB = (eg_hcalHForHoverE[i])/(eg_energy[i]);
	float hoe_cut_EB =0.17;  //0.25;
	float ooemoop_cut_EB =0.035;  //0.04;
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
	float trkisohlt_cut_EB =2.4; //3.5;
	float trkisol1_cut_EB = 4.0; //5.5;
	if ( ( fabs(eg_eta[i]) > 0.8 ) && ( fabs(eg_eta[i]) < 1.479 )  ) {
	  trkisol1_cut_EB = 8.0;
	}
	////end of barrel cuts

	////endcap cuts
	float hoe_EE = (eg_hgcalHForHoverE[i])/(eg_energy[i]);
	float hoe_cut_EE = 0.15 + (5.0/(eg_energy[i]));
	float vv_cut_EE = 0.8*0.8;  //0.9*0.9;
	float ww_cut_EE = 8*8;  //9*9;
	float hgcaliso_cut_EE = 130.0; //150.0;
	if ( ( fabs(eg_eta[i]) > 2.0 ) ) {
	  hgcaliso_cut_EE = 340.0;  //350;
	}
	float pms2_cut_EE= 65.0; //75.0;
	float ooemoop_cut_EE = 0.01; //0.04;
	float deta_cut_EE = 0.003; //0.004;
	float dphi_cut_EE = 0.02; //0.04;
	float npix_cut_EE = 2;
	float chi2_cut_EE = 50.0;
	float trkisohlt_cut_EE = 2.0; //3.0;
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

	if ( (weight>0.) && (eg_et[i]>pt_cut) &&  
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
	     ///(eg_hltisov6[i]<trkisohlt_cut) && 
	     (eg_l1iso[i]<trkisol1_cut) &&
	     (eg_ecaliso[i]<ecaliso_cut_EB ) && 
	     (eg_hcaliso[i]<hcaliso_cut_EB) && 
	     (eg_sigmaIEtaIEta[i]<sieie_cut_EB) 
	     ) {
	  nEle_passed=nEle_passed+1;
	}


      } ////end of ele loop
      
      if ( nEle_passed>0 ) {  
	//	nEvt_passed = nEvt_passed+1;
	nEvt_passed = nEvt_passed+weight;
      }
   }

   std::cout << "myweight " << myweight << std::endl;
   std::cout << "nEvt_passed " << nEvt_passed << std::endl;
   //std::cout << "rate " << myweight*nEvt_passed << std::endl;
   std::cout << "rate " << nEvt_passed << std::endl;   
   //std::cout << "rate error " << myweight * sqrt(nEvt_passed) << std::endl;
   if (nEvt_passed==0) {
     std::cout << "correct rate error for 0 event = " << myweight*1.84;
   }
}

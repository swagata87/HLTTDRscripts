#define CutOpt_cxx
#include "CutOpt.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void CutOpt::Loop()
{
//   In a ROOT session, you can do:
//      root> .L CutOpt.C
//      root> CutOpt t
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
  
  TFile* outputFile = new TFile("hist_DY_WP80.root","RECREATE");

  TH1D*  num_ele_geneta_hoe = new TH1D("num_ele_geneta_hoe", "num_ele_geneta_hoe", 30, -3, 3); 
  TH1D*  num_ele_geneta_vv = new TH1D("num_ele_geneta_vv", "num_ele_geneta_vv", 30, -3, 3); 
  TH1D*  num_ele_geneta_ww = new TH1D("num_ele_geneta_ww", "num_ele_geneta_ww", 30, -3, 3); 
  TH1D*  num_ele_geneta_hgcaliso = new TH1D("num_ele_geneta_hgcaliso", "num_ele_geneta_hgcaliso", 30, -3, 3); 
  TH1D*  num_ele_geneta_pms2 = new TH1D("num_ele_geneta_pms2", "num_ele_geneta_pms2", 30, -3, 3); 
  TH1D*  num_ele_geneta_ooemoop = new TH1D("num_ele_geneta_ooemoop", "num_ele_geneta_ooemoop", 30, -3, 3); 
  TH1D*  num_ele_geneta_deta = new TH1D("num_ele_geneta_deta", "num_ele_geneta_deta", 30, -3, 3); 
  TH1D*  num_ele_geneta_dphi = new TH1D("num_ele_geneta_dphi", "num_ele_geneta_dphi", 30, -3, 3); 
  TH1D*  num_ele_geneta_npix = new TH1D("num_ele_geneta_npix", "num_ele_geneta_npix", 30, -3, 3); 
  TH1D*  num_ele_geneta_chi2 = new TH1D("num_ele_geneta_chi2", "num_ele_geneta_chi2", 30, -3, 3); 
  TH1D*  num_ele_geneta_trkisohlt = new TH1D("num_ele_geneta_trkisohlt", "num_ele_geneta_trkisohlt", 30, -3, 3);   
  TH1D*  num_ele_geneta_trkisol1 = new TH1D("num_ele_geneta_trkisol1", "num_ele_geneta_trkisol1", 30, -3, 3);   
  TH1D*  num_ele_geneta_ecaliso = new TH1D("num_ele_geneta_ecaliso", "num_ele_geneta_ecaliso", 30, -3, 3);   
  TH1D*  num_ele_geneta_hcaliso = new TH1D("num_ele_geneta_hcaliso", "num_ele_geneta_hcaliso", 30, -3, 3);   
  TH1D*  num_ele_geneta_sieie = new TH1D("num_ele_geneta_sieie", "num_ele_geneta_sieie", 30, -3, 3);

  TH1D*  den_ele_genpt_EB = new TH1D("den_ele_genpt_EB", "den_ele_genpt_EB", 40, 0, 120); 
  TH1D*  den_trkele_genpt_EB = new TH1D("den_trkele_genpt_EB", "den_trkele_genpt_EB", 40, 0, 120); 
  TH1D*  den_trkele_genpt_EE = new TH1D("den_trkele_genpt_EE", "den_trkele_genpt_EE", 40, 0, 120); 

  TH1D*  den_trkele_geneta = new TH1D("den_trkele_geneta", "den_trkele_geneta", 30, -3, 3); 

  TH1D*  den_ele_genpt_EE = new TH1D("den_ele_genpt_EE", "den_ele_genpt_EE", 40, 0, 120); 
  TH1D*  den_ele_geneta = new TH1D("den_ele_geneta", "den_ele_geneta", 30, -3, 3); 

  TH1D*  num_ele_geneta_all = new TH1D("num_ele_geneta_all", "num_ele_geneta_all", 30, -3, 3); 
  TH1D*  num_ele_geneta_passL1 = new TH1D("num_ele_geneta_passL1", "num_ele_geneta_passL1", 30, -3, 3); 

  TH1D*  num_ele_genpt_ecaliso_EB = new TH1D("num_ele_genpt_ecaliso_EB", "num_ele_genpt_ecaliso_EB", 40, 0, 120);
  TH1D*  num_ele_genpt_hcaliso_EB = new TH1D("num_ele_genpt_hcaliso_EB", "num_ele_genpt_hcaliso_EB", 40, 0, 120);

  TH1D*  num_ele_genpt_pms2_EB = new TH1D("num_ele_genpt_pms2_EB", "num_ele_genpt_pms2_EB", 40, 0, 120);
  TH1D*  num_ele_genpt_pms2_EE = new TH1D("num_ele_genpt_pms2_EE", "num_ele_genpt_pms2_EE", 40, 0, 120);

  TH1D*  num_ele_genpt_sieie_EB = new TH1D("num_ele_genpt_sieie_EB", "num_ele_genpt_sieie_EB", 40, 0, 120);

  TH1D*  num_ele_genpt_hoe_EB = new TH1D("num_ele_genpt_hoe_EB", "num_ele_genpt_hoe_EB", 40, 0, 120);
  TH1D*  num_ele_genpt_hoe_EE = new TH1D("num_ele_genpt_hoe_EE", "num_ele_genpt_hoe_EE", 40, 0, 120);

  TH1D*  num_ele_genpt_ooemoop_EB = new TH1D("num_ele_genpt_ooemoop_EB", "num_ele_genpt_ooemoop_EB", 40, 0, 120);
  TH1D*  num_ele_genpt_ooemoop_EE = new TH1D("num_ele_genpt_ooemoop_EE", "num_ele_genpt_ooemoop_EE", 40, 0, 120);

  TH1D*  num_ele_genpt_deta_EB = new TH1D("num_ele_genpt_deta_EB", "num_ele_genpt_deta_EB", 40, 0, 120);
  TH1D*  num_ele_genpt_deta_EE = new TH1D("num_ele_genpt_deta_EE", "num_ele_genpt_deta_EE", 40, 0, 120);

  TH1D*  num_ele_genpt_dphi_EB = new TH1D("num_ele_genpt_dphi_EB", "num_ele_genpt_dphi_EB", 40, 0, 120);
  TH1D*  num_ele_genpt_dphi_EE = new TH1D("num_ele_genpt_dphi_EE", "num_ele_genpt_dphi_EE", 40, 0, 120);

  TH1D*  num_ele_genpt_npix_EB = new TH1D("num_ele_genpt_npix_EB", "num_ele_genpt_npix_EB", 40, 0, 120);
  TH1D*  num_ele_genpt_chi2_EB = new TH1D("num_ele_genpt_chi2_EB", "num_ele_genpt_chi2_EB", 40, 0, 120);
  TH1D*  num_ele_genpt_trkisohlt_EB = new TH1D("num_ele_genpt_trkisohlt_EB", "num_ele_genpt_trkisohlt_EB", 40, 0, 120);
  TH1D*  num_ele_genpt_trkisol1_EB = new TH1D("num_ele_genpt_trkisol1_EB", "num_ele_genpt_trkisol1_EB", 40, 0, 120);

  TH1D*  num_ele_genpt_npix_EE = new TH1D("num_ele_genpt_npix_EE", "num_ele_genpt_npix_EE", 40, 0, 120);
  TH1D*  num_ele_genpt_chi2_EE = new TH1D("num_ele_genpt_chi2_EE", "num_ele_genpt_chi2_EE", 40, 0, 120);
  TH1D*  num_ele_genpt_trkisohlt_EE = new TH1D("num_ele_genpt_trkisohlt_EE", "num_ele_genpt_trkisohlt_EE", 40, 0, 120);
  TH1D*  num_ele_genpt_trkisol1_EE = new TH1D("num_ele_genpt_trkisol1_EE", "num_ele_genpt_trkisol1_EE", 40, 0, 120);

  TH1D*  num_ele_genpt_passL1_EB = new TH1D("num_ele_genpt_passL1_EB", "num_ele_genpt_passL1_EB", 40, 0, 120);
  TH1D*  num_ele_genpt_passL1_EE = new TH1D("num_ele_genpt_passL1_EE", "num_ele_genpt_passL1_EE", 40, 0, 120);
  TH1D*  num_ele_genpt_vv_EE = new TH1D("num_ele_genpt_vv_EE", "num_ele_genpt_vv_EE", 40, 0, 120);
  TH1D*  num_ele_genpt_ww_EE = new TH1D("num_ele_genpt_ww_EE", "num_ele_genpt_ww_EE", 40, 0, 120);
  TH1D*  num_ele_genpt_hgcaliso_EE = new TH1D("num_ele_genpt_hgcaliso_EE", "num_ele_genpt_hgcaliso_EE", 40, 0, 120);

  TH1D*  num_ele_genpt_all_EB = new TH1D("num_ele_genpt_all_EB", "num_ele_genpt_all_EB", 40, 0, 120);
  TH1D*  num_ele_genpt_all_EE = new TH1D("num_ele_genpt_all_EE", "num_ele_genpt_all_EE", 40, 0, 120);

  //num_ele_genpt_vv_EE
  
  Long64_t nentries = fChain->GetEntriesFast();
  
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
    for (int i=0; i<nrEgs; i++) {

      float ecaliso_cut_EB=6.0;
      float hcaliso_cut_EB=13;
      if ( ( fabs(eg_eta[i]) > 0.8 ) && ( fabs(eg_eta[i]) < 1.479 )  ) {
	hcaliso_cut_EB =18.0;
      }
      float pms2_cut_EB=55.0;
      float sieie_cut_EB=0.013;
      float hoe_EB = (eg_hcalHForHoverE[i])/(eg_energy[i]);
      float hoe_cut_EB = 0.25;
      float ooemoop_cut_EB = 0.04;
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
      float trkisohlt_cut_EB = 3.5;
      float trkisol1_cut_EB = 5.5;
      if ( ( fabs(eg_eta[i]) > 0.8 ) && ( fabs(eg_eta[i]) < 1.479 )  ) {
	trkisol1_cut_EB = 8.0;
      }
      
      //////////////////////////////////////////////////////
      float hoe_EE = (eg_hgcalHForHoverE[i])/(eg_energy[i]);
      float hoe_cut_EE = 0.15 + (5.0/(eg_energy[i]));
      float vv_cut_EE = 0.9*0.9;
      float ww_cut_EE = 9*9;
      float hgcaliso_cut_EE = 150.0;
      if ( ( fabs(eg_eta[i]) > 2.0 ) ) {
	hgcaliso_cut_EE = 350;
      } 
      float pms2_cut_EE=75.0;
      float ooemoop_cut_EE = 0.04;
      float deta_cut_EE = 0.004;
      float dphi_cut_EE = 0.04;
      float npix_cut_EE = 2;
      float chi2_cut_EE = 50.0;
      float trkisohlt_cut_EE = 3.0;
      float trkisol1_cut_EE = 5.5;
      
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
      
      ///endcap start    
      if ( (eg_gen_et[i]>10.0) && ((fabs(eg_eta[i]))>1.56) && ((fabs(eg_eta[i]))<2.80) ) {
	den_ele_genpt_EE->Fill(eg_gen_et[i]);
	
	if (eg_nGsf[i]>0) {
	  den_trkele_genpt_EE->Fill(eg_gen_et[i]);
	}
	if (passL1) {
	  num_ele_genpt_passL1_EE->Fill(eg_gen_et[i]);
	}
	if (hoe_EE<hoe_cut_EE) {
	  num_ele_genpt_hoe_EE->Fill(eg_gen_et[i]);
	}
	//
	if (eg_sigma2vv[i]<vv_cut_EE) {
        num_ele_genpt_vv_EE->Fill(eg_gen_et[i]);
	}
	//
	if (eg_sigma2ww[i]<ww_cut_EE) {
        num_ele_genpt_ww_EE->Fill(eg_gen_et[i]);
	}
	//
	if (eg_hgcaliso_layerclus[i]<hgcaliso_cut_EE) {
	  num_ele_genpt_hgcaliso_EE->Fill(eg_gen_et[i]);
	}
	//
	if (eg_pms2[i]<pms2_cut_EE) {
	  num_ele_genpt_pms2_EE->Fill(eg_gen_et[i]);
	}
	//
	if (eg_invEInvP[i]<ooemoop_cut_EE) {
	  num_ele_genpt_ooemoop_EE->Fill(eg_gen_et[i]);
	}
	//
	if (eg_trkDEtaSeed[i]<deta_cut_EE) {
	  num_ele_genpt_deta_EE->Fill(eg_gen_et[i]);
	}
	//
	if (eg_trkDPhi[i]<dphi_cut_EE) {
	  num_ele_genpt_dphi_EE->Fill(eg_gen_et[i]);
	}
	//
	//
	if (eg_nLayerIT[i]>npix_cut_EE) {
	  num_ele_genpt_npix_EE->Fill(eg_gen_et[i]);
	}
	//
	if (eg_normChi2[i]<chi2_cut_EE) {
	  num_ele_genpt_chi2_EE->Fill(eg_gen_et[i]);
	}
	//
	if (eg_hltisov6[i]<trkisohlt_cut_EE) {
	num_ele_genpt_trkisohlt_EE->Fill(eg_gen_et[i]);
	}
	//
	if (eg_l1iso[i]<trkisol1_cut_EE) {
	num_ele_genpt_trkisol1_EE->Fill(eg_gen_et[i]);
	}
	//
	//all cuts                                                                         
	if ( passL1 && (hoe_EE<hoe_cut_EE) &&  (eg_sigma2vv[i]<vv_cut_EE) &&  (eg_sigma2ww[i]<ww_cut_EE) && (eg_hgcaliso_layerclus[i]<hgcaliso_cut_EE) && 
	     (eg_pms2[i]<pms2_cut_EE) && (eg_invEInvP[i]<ooemoop_cut_EE) && 
	     (eg_trkDEtaSeed[i]<deta_cut_EE) && (eg_trkDPhi[i]<dphi_cut_EE) && (eg_nLayerIT[i]>npix_cut_EE) && (eg_normChi2[i]<chi2_cut_EE) &&  
	     (eg_hltisov6[i]<trkisohlt_cut_EE) && (eg_l1iso[i]<trkisol1_cut_EE)  ) {
	  num_ele_genpt_all_EE->Fill(eg_gen_et[i]);
	}
	
      } // endcap end
      
      if ( (eg_gen_et[i]>10.0) && ((fabs(eg_eta[i]))<1.44) ) {
	den_ele_genpt_EB->Fill(eg_gen_et[i]);
	
	if (eg_nGsf[i]>0) {
	  den_trkele_genpt_EB->Fill(eg_gen_et[i]);
	}
	
	if (eg_ecaliso[i]<ecaliso_cut_EB) {
	  num_ele_genpt_ecaliso_EB->Fill(eg_gen_et[i]);
	}
	//
	if (eg_hcaliso[i]<hcaliso_cut_EB) {
	  num_ele_genpt_hcaliso_EB->Fill(eg_gen_et[i]);
	}
	//
	if (eg_pms2[i]<pms2_cut_EB) {
	  num_ele_genpt_pms2_EB->Fill(eg_gen_et[i]);
	}
	//
	if (eg_sigmaIEtaIEta[i]<sieie_cut_EB) {
	  num_ele_genpt_sieie_EB->Fill(eg_gen_et[i]);
	}
	//
	if (hoe_EB<hoe_cut_EB) {
	  num_ele_genpt_hoe_EB->Fill(eg_gen_et[i]);
	}
	//
	if (eg_invEInvP[i]<ooemoop_cut_EB) {
	  num_ele_genpt_ooemoop_EB->Fill(eg_gen_et[i]);
	}
	//
	if (eg_trkDEtaSeed[i]<deta_cut_EB) {
	  num_ele_genpt_deta_EB->Fill(eg_gen_et[i]);
	}
	//
	if (eg_trkDPhi[i]<dphi_cut_EB) {
	  num_ele_genpt_dphi_EB->Fill(eg_gen_et[i]);
	}
	//
	if (eg_nLayerIT[i]>npix_cut_EB) {
	  num_ele_genpt_npix_EB->Fill(eg_gen_et[i]);
	}
	//
	if (eg_normChi2[i]<chi2_cut_EB) {
	  num_ele_genpt_chi2_EB->Fill(eg_gen_et[i]);
	}
	//
	if (eg_hltisov6[i]<trkisohlt_cut_EB) {
	  num_ele_genpt_trkisohlt_EB->Fill(eg_gen_et[i]);
	}
	//
	if (eg_l1iso[i]<trkisol1_cut_EB) {
	  num_ele_genpt_trkisol1_EB->Fill(eg_gen_et[i]);
	}
	//
	//
	if (passL1) {
	  num_ele_genpt_passL1_EB->Fill(eg_gen_et[i]);
	}
	//
	//all cuts
	if ( ( eg_ecaliso[i]<ecaliso_cut_EB  ) && (eg_hcaliso[i]<hcaliso_cut_EB) && 
	     (eg_pms2[i]<pms2_cut_EB) && (eg_sigmaIEtaIEta[i]<sieie_cut_EB) && (hoe_EB<hoe_cut_EB) && (eg_invEInvP[i]<ooemoop_cut_EB) && 
	     (eg_trkDEtaSeed[i]<deta_cut_EB) && (eg_trkDPhi[i]<dphi_cut_EB) && (eg_nLayerIT[i]>npix_cut_EB) && (eg_normChi2[i]<chi2_cut_EB) && (eg_hltisov6[i]<trkisohlt_cut_EB) && 
	     (eg_l1iso[i]<trkisol1_cut_EB) &&  passL1) {
	  num_ele_genpt_all_EB->Fill(eg_gen_et[i]);
	} 
	
      } //barrel
      
      ///////////////////////////////////////////////
      /////// Barrel + Endcap //////////////////////
      //////////////////////////////////////////////
      if ( (eg_gen_et[i]>0.) &&  (eg_et[i]>40.0) ) {
	den_ele_geneta->Fill(eg_gen_eta[i]);
	//
	if (eg_nGsf[i]>0) {
	  den_trkele_geneta->Fill(eg_gen_eta[i]);
	}
	//
	if ( passL1 ) {
	  num_ele_geneta_passL1->Fill(eg_gen_eta[i]);
	}
	//
	if (hoe<hoe_cut) {
	  num_ele_geneta_hoe->Fill(eg_gen_eta[i]);
	}
	//
	if (eg_sigma2vv[i]<vv_cut_EE) {
	  num_ele_geneta_vv->Fill(eg_gen_eta[i]);
	}
	//
	if (eg_sigma2ww[i]<ww_cut_EE) {
	  num_ele_geneta_ww->Fill(eg_gen_eta[i]);
	}
	//
	if (eg_hgcaliso_layerclus[i]<hgcaliso_cut_EE) {
	  num_ele_geneta_hgcaliso->Fill(eg_gen_eta[i]);
	}
	//
	if (eg_pms2[i]<pms2_cut) {
	  num_ele_geneta_pms2->Fill(eg_gen_eta[i]);
	}
	//
	if (eg_invEInvP[i]<ooemoop_cut) {
	  num_ele_geneta_ooemoop->Fill(eg_gen_eta[i]);
	}
	//
	if (eg_trkDEtaSeed[i]<deta_cut) {
	  num_ele_geneta_deta->Fill(eg_gen_eta[i]);
	}
	//
	if (eg_trkDPhi[i]<dphi_cut) {
	  num_ele_geneta_dphi->Fill(eg_gen_eta[i]);
	}
	//
	if (eg_nLayerIT[i]>npix_cut) {
	  num_ele_geneta_npix->Fill(eg_gen_eta[i]);
	}
	//
	if (eg_normChi2[i]<chi2_cut) {
	  num_ele_geneta_chi2->Fill(eg_gen_eta[i]);
	}
	//
	if (eg_hltisov6[i]<trkisohlt_cut) {
	  num_ele_geneta_trkisohlt->Fill(eg_gen_eta[i]);
	}
	//
	if (eg_l1iso[i]<trkisol1_cut) {
	  num_ele_geneta_trkisol1->Fill(eg_gen_eta[i]);
	}
	//
	if ( eg_ecaliso[i]<ecaliso_cut_EB  ) {
	  num_ele_geneta_ecaliso->Fill(eg_gen_eta[i]);
	}
	//
	if (eg_hcaliso[i]<hcaliso_cut_EB) {
	  num_ele_geneta_hcaliso->Fill(eg_gen_eta[i]);
	}
	//
	if (eg_sigmaIEtaIEta[i]<sieie_cut_EB) {
	  num_ele_geneta_sieie->Fill(eg_gen_eta[i]);
	}
	//
	//all cuts
	if ( passL1 && (hoe<hoe_cut) &&  (eg_sigma2vv[i]<vv_cut_EE) &&  (eg_sigma2ww[i]<ww_cut_EE) && (eg_hgcaliso_layerclus[i]<hgcaliso_cut_EE) && (eg_pms2[i]<pms2_cut) && (eg_invEInvP[i]<ooemoop_cut) &&
	     (eg_trkDEtaSeed[i]<deta_cut) && (eg_trkDPhi[i]<dphi_cut) && (eg_nLayerIT[i]>npix_cut) && (eg_normChi2[i]<chi2_cut) &&  (eg_hltisov6[i]<trkisohlt_cut) && (eg_l1iso[i]<trkisol1_cut) && 
	     ( eg_ecaliso[i]<ecaliso_cut_EB  ) && (eg_hcaliso[i]<hcaliso_cut_EB) && (eg_sigmaIEtaIEta[i]<sieie_cut_EB) ) {
	  num_ele_geneta_all->Fill(eg_gen_eta[i]);
	}
	
      } //
    }
    
  }
  den_ele_genpt_EE->Write();
  den_ele_genpt_EB->Write();
  den_trkele_genpt_EB->Write();
  den_trkele_genpt_EE->Write();
  den_ele_geneta->Write();
  num_ele_geneta_all->Write();

  num_ele_genpt_ecaliso_EB->Write();
  num_ele_genpt_hcaliso_EB->Write();
  num_ele_genpt_pms2_EB->Write();
  num_ele_genpt_pms2_EE->Write();

  num_ele_genpt_sieie_EB->Write();

  num_ele_genpt_hoe_EB->Write();
  num_ele_genpt_hoe_EE->Write();

  num_ele_genpt_ooemoop_EB->Write();
  num_ele_genpt_ooemoop_EE->Write();

  num_ele_genpt_deta_EB->Write();
  num_ele_genpt_deta_EE->Write();

  num_ele_genpt_dphi_EB->Write();
  num_ele_genpt_dphi_EE->Write();

  num_ele_genpt_npix_EB->Write();
  num_ele_genpt_chi2_EB->Write();
  num_ele_genpt_trkisohlt_EB->Write();
  num_ele_genpt_trkisol1_EB->Write();

  num_ele_genpt_npix_EE->Write();
  num_ele_genpt_chi2_EE->Write();
  num_ele_genpt_trkisohlt_EE->Write();
  num_ele_genpt_trkisol1_EE->Write();

  num_ele_genpt_passL1_EB->Write();
  num_ele_genpt_passL1_EE->Write();
  num_ele_genpt_vv_EE->Write();
  num_ele_genpt_ww_EE->Write();
  num_ele_genpt_all_EB->Write();
  num_ele_genpt_all_EE->Write();
  num_ele_genpt_hgcaliso_EE->Write();

  num_ele_geneta_passL1->Write();
  den_trkele_geneta->Write();
  num_ele_geneta_hoe->Write();
  num_ele_geneta_vv->Write();
  num_ele_geneta_ww->Write();
  num_ele_geneta_hgcaliso->Write();
  num_ele_geneta_pms2->Write();
  num_ele_geneta_ooemoop->Write();
  num_ele_geneta_deta->Write();
  num_ele_geneta_dphi->Write();
  num_ele_geneta_npix->Write();
  num_ele_geneta_chi2->Write();
  num_ele_geneta_trkisohlt->Write();
  num_ele_geneta_trkisol1->Write();
  num_ele_geneta_ecaliso->Write();
  num_ele_geneta_hcaliso->Write();
  num_ele_geneta_sieie->Write();
}

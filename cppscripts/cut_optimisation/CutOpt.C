#define CutOpt_cxx
#include "CutOpt.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>

void CutOpt::SingleEle(bool run_tight_wp=1)
{
  if (fChain == 0) return;
  
  TFile* outputFile; // = new TFile("hist_DY_WP70.root","RECREATE");

  bool WP_tight = run_tight_wp;
  std::cout << "WP_tight " << WP_tight << std::endl;
  if (WP_tight) {
    std::cout << "running WP tight, 70% eff" << std::endl;
    outputFile = new TFile("hist_Zprime_WP70.root","RECREATE");

  }
  else {
    std::cout << "running WP loose, 80% eff" << std::endl;
    outputFile = new TFile("hist_Zprime_WP80.root","RECREATE");
  }

  const Int_t NBINS_eta = 17;
  Double_t edges_eta[NBINS_eta + 1] = {-3.0,-2.7,-2.4,-2.0,-1.56,-1.44,-1.0,-0.6,-0.2,0.2,0.6,1.0,1.44,1.56,2.0,2.4,2.7,3.0};

  const Int_t NBINS_pt = 16;
  //  Double_t edges_pt[NBINS_pt + 1] = {20,24,28,32,36,40,44,48,52,56,60,65,70,80,150};
  Double_t edges_pt[NBINS_pt + 1] = {20,30,40,50,60,70,80,90,100,150,200,250,300,350,400,450,500};

  const Int_t NBINS_pu = 8;
  Double_t edges_pu[NBINS_pu + 1] = {150,170,180,190,200,210,220,230,250};

  TH1D*  den_ele_PU_EB = new TH1D("den_ele_PU_EB", "den_ele_PU_EB", NBINS_pu, edges_pu); 
  TH1D*  den_ele_PU_EE = new TH1D("den_ele_PU_EE", "den_ele_PU_EE", NBINS_pu, edges_pu); 

  TH1D*  num_ele_PU_EB_all = new TH1D("num_ele_PU_EB_all", "num_ele_PU_EB_all", NBINS_pu, edges_pu); 
  TH1D*  num_ele_PU_EE_all = new TH1D("num_ele_PU_EE_all", "num_ele_PU_EE_all", NBINS_pu, edges_pu); 

  TH1D*  num_ele_geneta_hoe = new TH1D("num_ele_geneta_hoe", "num_ele_geneta_hoe", NBINS_eta, edges_eta); 
  TH1D*  num_ele_geneta_vv = new TH1D("num_ele_geneta_vv", "num_ele_geneta_vv", NBINS_eta, edges_eta); 
  TH1D*  num_ele_geneta_ww = new TH1D("num_ele_geneta_ww", "num_ele_geneta_ww", NBINS_eta, edges_eta); 
  TH1D*  num_ele_geneta_hgcaliso = new TH1D("num_ele_geneta_hgcaliso", "num_ele_geneta_hgcaliso", NBINS_eta, edges_eta); 
  TH1D*  num_ele_geneta_pms2 = new TH1D("num_ele_geneta_pms2", "num_ele_geneta_pms2", NBINS_eta, edges_eta); 
  TH1D*  num_ele_geneta_ooemoop = new TH1D("num_ele_geneta_ooemoop", "num_ele_geneta_ooemoop", NBINS_eta, edges_eta); 
  TH1D*  num_ele_geneta_deta = new TH1D("num_ele_geneta_deta", "num_ele_geneta_deta", NBINS_eta, edges_eta); 
  TH1D*  num_ele_geneta_dphi = new TH1D("num_ele_geneta_dphi", "num_ele_geneta_dphi", NBINS_eta, edges_eta); 
  TH1D*  num_ele_geneta_npix = new TH1D("num_ele_geneta_npix", "num_ele_geneta_npix", NBINS_eta, edges_eta); 
  TH1D*  num_ele_geneta_chi2 = new TH1D("num_ele_geneta_chi2", "num_ele_geneta_chi2", NBINS_eta, edges_eta); 
  TH1D*  num_ele_geneta_trkisohlt = new TH1D("num_ele_geneta_trkisohlt", "num_ele_geneta_trkisohlt", NBINS_eta, edges_eta);   
  TH1D*  num_ele_geneta_trkisol1 = new TH1D("num_ele_geneta_trkisol1", "num_ele_geneta_trkisol1", NBINS_eta, edges_eta);   
  TH1D*  num_ele_geneta_ecaliso = new TH1D("num_ele_geneta_ecaliso", "num_ele_geneta_ecaliso", NBINS_eta, edges_eta);   
  TH1D*  num_ele_geneta_hcaliso = new TH1D("num_ele_geneta_hcaliso", "num_ele_geneta_hcaliso", NBINS_eta, edges_eta);   
  TH1D*  num_ele_geneta_sieie = new TH1D("num_ele_geneta_sieie", "num_ele_geneta_sieie", NBINS_eta, edges_eta);

  TH1D*  den_ele_genpt_EB = new TH1D("den_ele_genpt_EB", "den_ele_genpt_EB", NBINS_pt, edges_pt); 
  TH1D*  den_trkele_genpt_EB = new TH1D("den_trkele_genpt_EB", "den_trkele_genpt_EB", NBINS_pt, edges_pt); 
  TH1D*  den_trkele_genpt_EE = new TH1D("den_trkele_genpt_EE", "den_trkele_genpt_EE", NBINS_pt, edges_pt); 

  TH1D*  den_trkele_geneta = new TH1D("den_trkele_geneta", "den_trkele_geneta", NBINS_eta, edges_eta); 

  TH1D*  den_ele_genpt_EE = new TH1D("den_ele_genpt_EE", "den_ele_genpt_EE", NBINS_pt, edges_pt); 
  TH1D*  den_ele_geneta = new TH1D("den_ele_geneta", "den_ele_geneta", NBINS_eta, edges_eta); 

  TH1D*  num_ele_geneta_all = new TH1D("num_ele_geneta_all", "num_ele_geneta_all", NBINS_eta, edges_eta); 
  TH1D*  num_ele_geneta_passL1 = new TH1D("num_ele_geneta_passL1", "num_ele_geneta_passL1", NBINS_eta, edges_eta); 

  TH1D*  num_ele_genpt_ecaliso_EB = new TH1D("num_ele_genpt_ecaliso_EB", "num_ele_genpt_ecaliso_EB", NBINS_pt, edges_pt);
  TH1D*  num_ele_genpt_hcaliso_EB = new TH1D("num_ele_genpt_hcaliso_EB", "num_ele_genpt_hcaliso_EB", NBINS_pt, edges_pt);

  TH1D*  num_ele_genpt_pms2_EB = new TH1D("num_ele_genpt_pms2_EB", "num_ele_genpt_pms2_EB", NBINS_pt, edges_pt);
  TH1D*  num_ele_genpt_pms2_EE = new TH1D("num_ele_genpt_pms2_EE", "num_ele_genpt_pms2_EE", NBINS_pt, edges_pt);

  TH1D*  num_ele_genpt_sieie_EB = new TH1D("num_ele_genpt_sieie_EB", "num_ele_genpt_sieie_EB", NBINS_pt, edges_pt);

  TH1D*  num_ele_genpt_hoe_EB = new TH1D("num_ele_genpt_hoe_EB", "num_ele_genpt_hoe_EB", NBINS_pt, edges_pt);
  TH1D*  num_ele_genpt_hoe_EE = new TH1D("num_ele_genpt_hoe_EE", "num_ele_genpt_hoe_EE", NBINS_pt, edges_pt);

  TH1D*  num_ele_genpt_ooemoop_EB = new TH1D("num_ele_genpt_ooemoop_EB", "num_ele_genpt_ooemoop_EB", NBINS_pt, edges_pt);
  TH1D*  num_ele_genpt_ooemoop_EE = new TH1D("num_ele_genpt_ooemoop_EE", "num_ele_genpt_ooemoop_EE", NBINS_pt, edges_pt);

  TH1D*  num_ele_genpt_deta_EB = new TH1D("num_ele_genpt_deta_EB", "num_ele_genpt_deta_EB", NBINS_pt, edges_pt);
  TH1D*  num_ele_genpt_deta_EE = new TH1D("num_ele_genpt_deta_EE", "num_ele_genpt_deta_EE", NBINS_pt, edges_pt);

  TH1D*  num_ele_genpt_dphi_EB = new TH1D("num_ele_genpt_dphi_EB", "num_ele_genpt_dphi_EB", NBINS_pt, edges_pt);
  TH1D*  num_ele_genpt_dphi_EE = new TH1D("num_ele_genpt_dphi_EE", "num_ele_genpt_dphi_EE", NBINS_pt, edges_pt);

  TH1D*  num_ele_genpt_npix_EB = new TH1D("num_ele_genpt_npix_EB", "num_ele_genpt_npix_EB", NBINS_pt, edges_pt);
  TH1D*  num_ele_genpt_chi2_EB = new TH1D("num_ele_genpt_chi2_EB", "num_ele_genpt_chi2_EB", NBINS_pt, edges_pt);
  TH1D*  num_ele_genpt_trkisohlt_EB = new TH1D("num_ele_genpt_trkisohlt_EB", "num_ele_genpt_trkisohlt_EB", NBINS_pt, edges_pt);
  TH1D*  num_ele_genpt_trkisol1_EB = new TH1D("num_ele_genpt_trkisol1_EB", "num_ele_genpt_trkisol1_EB", NBINS_pt, edges_pt);

  TH1D*  num_ele_genpt_npix_EE = new TH1D("num_ele_genpt_npix_EE", "num_ele_genpt_npix_EE", NBINS_pt, edges_pt);
  TH1D*  num_ele_genpt_chi2_EE = new TH1D("num_ele_genpt_chi2_EE", "num_ele_genpt_chi2_EE", NBINS_pt, edges_pt);
  TH1D*  num_ele_genpt_trkisohlt_EE = new TH1D("num_ele_genpt_trkisohlt_EE", "num_ele_genpt_trkisohlt_EE", NBINS_pt, edges_pt);
  TH1D*  num_ele_genpt_trkisol1_EE = new TH1D("num_ele_genpt_trkisol1_EE", "num_ele_genpt_trkisol1_EE", NBINS_pt, edges_pt);

  TH1D*  num_ele_genpt_passL1_EB = new TH1D("num_ele_genpt_passL1_EB", "num_ele_genpt_passL1_EB", NBINS_pt, edges_pt);
  TH1D*  num_ele_genpt_passL1_EE = new TH1D("num_ele_genpt_passL1_EE", "num_ele_genpt_passL1_EE", NBINS_pt, edges_pt);
  TH1D*  num_ele_genpt_vv_EE = new TH1D("num_ele_genpt_vv_EE", "num_ele_genpt_vv_EE", NBINS_pt, edges_pt);
  TH1D*  num_ele_genpt_ww_EE = new TH1D("num_ele_genpt_ww_EE", "num_ele_genpt_ww_EE", NBINS_pt, edges_pt);
  TH1D*  num_ele_genpt_hgcaliso_EE = new TH1D("num_ele_genpt_hgcaliso_EE", "num_ele_genpt_hgcaliso_EE", NBINS_pt, edges_pt);

  TH1D*  num_ele_genpt_all_EB = new TH1D("num_ele_genpt_all_EB", "num_ele_genpt_all_EB", NBINS_pt, edges_pt);
  TH1D*  num_ele_genpt_all_EE = new TH1D("num_ele_genpt_all_EE", "num_ele_genpt_all_EE", NBINS_pt, edges_pt);

  //num_ele_genpt_vv_EE
  
  Long64_t nentries = fChain->GetEntriesFast();

  
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
    for (int i=0; i<nrEgs; i++) {

      //// Barrel cuts
      float ecaliso_cut_EB= (WP_tight==1) ? 4.4 + (0.02*eg_et[i]) : 9.0 + (0.02*eg_et[i]); 
      float hcaliso_cut_EB= (WP_tight==1) ? 12+(0.02*eg_et[i]) : 19+(0.02*eg_et[i]); 
      //      if ( ( fabs(eg_eta[i]) > 0.8 ) && ( fabs(eg_eta[i]) < 1.479 )  ) {
      //hcaliso_cut_EB = (WP_tight==1) ? 15.0 + (0.02*eg_et[i]) : 18.0 + (0.02*eg_et[i]); 
      //}
      float pms2_cut_EB= (WP_tight==1) ? 42.0 : 55.0; 
      float sieie_cut_EB= (WP_tight==1) ? 0.012 : 0.013; 
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
      float vv_cut_EE = (WP_tight==1) ? (0.8*0.8)+(0.0008*eg_et[i]) : (0.85*0.85)+(0.0008*eg_et[i]); 
      float ww_cut_EE = (WP_tight==1) ? (8*8)+(0.04*eg_et[i]) : (8.5*8.5)+(0.04*eg_et[i]); 
      float hgcaliso_cut_EE = (WP_tight==1) ? 130.0 + (0.05*eg_energy[i]) : 150.0 + (0.05*eg_energy[i]);
      if ( ( fabs(eg_eta[i]) > 2.0 ) ) {
	hgcaliso_cut_EE = (WP_tight==1) ? 340.0+ (0.05*eg_energy[i]) : 350.0+ (0.05*eg_energy[i]);
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
      
      bool passL1 = false;
      if (  fabs(eg_eta[i]) > 2.8 ) {
	passL1=true;
      }
      if ( (fabs(eg_eta[i])<2.8) &&  (eg_l1pho_et[i]>0) && (eg_l1pho_etThres[i]>51.) && (eg_l1pho_passQual[i]) ) {
	passL1=true;
      }
      if ( (fabs(eg_eta[i])<2.4) && (eg_l1ele_et[i]>0) && (eg_l1ele_etThresNonIso[i]>36.) && (eg_l1ele_passQual[i]) ) {
	passL1=true;
      }
      if ( (fabs(eg_eta[i])<2.4) && (eg_l1ele_et[i]>0) && (eg_l1ele_etThresIso[i]>28.) && (eg_l1ele_passQual[i]) && (eg_l1ele_passIsol[i]) ) {
	passL1=true;
      }
      
      ///endcap start    
      if ( (eg_gen_et[i]>10.0) && ((fabs(eg_eta[i]))>1.56) && ((fabs(eg_eta[i]))<2.70) ) {
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
    
      //PU , endcap
      if ( (eg_gen_et[i]>40.0) && ((fabs(eg_eta[i]))>1.56) && ((fabs(eg_eta[i]))<2.70) ) {
	//std::cout << "nrPtHats " << nrPtHats << std::endl;
	den_ele_PU_EE->Fill(nrPtHats);
	if ( passL1 && (hoe_EE<hoe_cut_EE) &&  (eg_sigma2vv[i]<vv_cut_EE) &&  (eg_sigma2ww[i]<ww_cut_EE) && (eg_hgcaliso_layerclus[i]<hgcaliso_cut_EE) &&
             (eg_pms2[i]<pms2_cut_EE) && (eg_invEInvP[i]<ooemoop_cut_EE) &&
             (eg_trkDEtaSeed[i]<deta_cut_EE) && (eg_trkDPhi[i]<dphi_cut_EE) && (eg_nLayerIT[i]>npix_cut_EE) && (eg_normChi2[i]<chi2_cut_EE) &&
             (eg_hltisov6[i]<trkisohlt_cut_EE) && (eg_l1iso[i]<trkisol1_cut_EE)  ) {
          num_ele_PU_EE_all->Fill(nrPtHats);
        }
      }

      if ( (eg_gen_et[i]>10.0) && ((fabs(eg_eta[i]))<1.44) ) {
	den_ele_genpt_EB->Fill(eg_gen_et[i]);
	
	if (eg_nGsf[i]>0) {
	  den_trkele_genpt_EB->Fill(eg_gen_et[i]);
	}
	
	if (eg_ecaliso[i]<ecaliso_cut_EB) {
	  num_ele_genpt_ecaliso_EB->Fill(eg_gen_et[i]);
	}
	//
	if (eg_hcalPFIsol_default[i] <hcaliso_cut_EB) {
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
	if ( ( eg_ecaliso[i]<ecaliso_cut_EB  ) && (eg_hcalPFIsol_default[i]<hcaliso_cut_EB) && 
	     (eg_pms2[i]<pms2_cut_EB) && (eg_sigmaIEtaIEta[i]<sieie_cut_EB) && (hoe_EB<hoe_cut_EB) && (eg_invEInvP[i]<ooemoop_cut_EB) && 
	     (eg_trkDEtaSeed[i]<deta_cut_EB) && (eg_trkDPhi[i]<dphi_cut_EB) && (eg_nLayerIT[i]>npix_cut_EB) && (eg_normChi2[i]<chi2_cut_EB) && (eg_hltisov6[i]<trkisohlt_cut_EB) && 
	     (eg_l1iso[i]<trkisol1_cut_EB) &&  passL1) {
	  num_ele_genpt_all_EB->Fill(eg_gen_et[i]);
	} 
	
      } //barrel

      //PU, barrel
      if ( (eg_gen_et[i]>40.0) && ((fabs(eg_eta[i]))<1.44) ) {
	den_ele_PU_EB->Fill(nrPtHats);
	if ( ( eg_ecaliso[i]<ecaliso_cut_EB  ) && (eg_hcalPFIsol_default[i]<hcaliso_cut_EB) &&
             (eg_pms2[i]<pms2_cut_EB) && (eg_sigmaIEtaIEta[i]<sieie_cut_EB) && (hoe_EB<hoe_cut_EB) && (eg_invEInvP[i]<ooemoop_cut_EB) &&
             (eg_trkDEtaSeed[i]<deta_cut_EB) && (eg_trkDPhi[i]<dphi_cut_EB) && (eg_nLayerIT[i]>npix_cut_EB) && (eg_normChi2[i]<chi2_cut_EB) && (eg_hltisov6[i]<trkisohlt_cut_EB) &&
             (eg_l1iso[i]<trkisol1_cut_EB) &&  passL1) {
          num_ele_PU_EB_all->Fill(nrPtHats);
        }
      }
      
      ///////////////////////////////////////////////
      /////// Barrel + Endcap //////////////////////
      //////////////////////////////////////////////
      if ( (eg_gen_et[i]>40.) &&  (eg_gen_et[i]<500.0) ) {
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
	if (eg_hcalPFIsol_default[i]<hcaliso_cut_EB) {
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
	     ( eg_ecaliso[i]<ecaliso_cut_EB  ) && (eg_hcalPFIsol_default[i]<hcaliso_cut_EB) && (eg_sigmaIEtaIEta[i]<sieie_cut_EB) ) {
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

  den_ele_PU_EE->Write();
  num_ele_PU_EE_all->Write();

  den_ele_PU_EB->Write();
  num_ele_PU_EB_all->Write();

}

///////////////////////////////////////////////////////////
/////////////// Single Pho Non Isolated //////////////////
////////////////////////////////////////////////////////////
void CutOpt::SinglePhoNonIso()
{
  if (fChain == 0) return;
  
  TFile* outputFile = new TFile("histSinglePhoNoIso_Sig_Zprime.root","RECREATE");

  const Int_t NBINS_eta = 17;
  Double_t edges_eta[NBINS_eta + 1] = {-3.0,-2.7,-2.4,-2.0,-1.56,-1.44,-1.0,-0.6,-0.2,0.2,0.6,1.0,1.44,1.56,2.0,2.4,2.7,3.0};

  //  const Int_t NBINS_pt = 6;
  //  Double_t edges_pt[NBINS_pt + 1] = {200,250,300,400,600,800,1000};

  const Int_t NBINS_pt = 16;
  //  Double_t edges_pt[NBINS_pt + 1] = {0,20,40,60,80,100,150,200,300,400,500,600,700,800,900,1000,1200,1400,1600,1800,2000};
  Double_t edges_pt[NBINS_pt + 1] = {20,30,40,50,60,70,80,90,100,150,200,250,300,350,400,450,500};

  TH1D*  den_ele_geneta = new TH1D("den_ele_geneta", "den_ele_geneta", NBINS_eta, edges_eta); 
  TH1D*  num_ele_geneta_all = new TH1D("num_ele_geneta_all", "num_ele_geneta_all", NBINS_eta, edges_eta); 

  TH1D*  den_ele_genpt_EE = new TH1D("den_ele_genpt_EE", "den_ele_genpt_EE", NBINS_pt, edges_pt); 
  TH1D*  num_ele_genpt_all_EE = new TH1D("num_ele_genpt_all_EE", "num_ele_genpt_all_EE", NBINS_pt, edges_pt);

  TH1D*  den_ele_genpt_EB = new TH1D("den_ele_genpt_EB", "den_ele_genpt_EB", NBINS_pt, edges_pt); 
  TH1D*  num_ele_genpt_all_EB = new TH1D("num_ele_genpt_all_EB", "num_ele_genpt_all_EB", NBINS_pt, edges_pt);

  Long64_t nentries = fChain->GetEntriesFast();
  
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
    for (int i=0; i<nrEgs; i++) {

      //// Barrel cuts                                                                                                                                                                      
      //float sieie_cut_EB= 0.013;
      float hoe_EB = (eg_hcalHForHoverE[i])/(eg_energy[i]);
      float hoe_cut_EB = 0.30;
      ////end of barrel cuts                                                                                                                                                               

      ////endcap cuts                                                                                                                                                                       
      float hoe_EE = (eg_hgcalHForHoverE[i])/(eg_energy[i]);
      float hoe_cut_EE = 0.30;
      //float vv_cut_EE = 0.9*0.9;
      //float ww_cut_EE = 9*9;
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
      if (  fabs(eg_eta[i]) > 2.8 ) {
	passL1=true;
      }
      if ( (fabs(eg_eta[i])<2.8) &&  (eg_l1pho_et[i]>0) && (eg_l1pho_etThres[i]>51.) && (eg_l1pho_passQual[i]) ) {
	//if ( (eg_l1pho_et[i]>0) && (eg_l1pho_etThres[i]>51.) && (eg_l1pho_passQual[i]) ) {
        passL1=true;
      }

      ///endcap start    
      if ( (eg_gen_et[i]>10.0) && ((fabs(eg_eta[i]))>1.56) && ((fabs(eg_eta[i]))<2.70) ) {
	den_ele_genpt_EE->Fill(eg_gen_et[i]);

	//all cuts                                                                         
	if ( passL1 &&
	     (hoe_EE<hoe_cut_EE) 
	     //(eg_sigma2vv[i]<vv_cut_EE)
	     ) {
	  num_ele_genpt_all_EE->Fill(eg_gen_et[i]);
	}
      } // endcap end

      //barrel start
      if ( (eg_gen_et[i]>10.0) && ((fabs(eg_eta[i]))<1.44) ) {
	den_ele_genpt_EB->Fill(eg_gen_et[i]);

	//all cuts
	if ( passL1 && 
	     (hoe_EB<hoe_cut_EB)
	     //(eg_sigmaIEtaIEta[i]<sieie_cut_EB)
	     ) {
	  num_ele_genpt_all_EB->Fill(eg_gen_et[i]);
	}
      } //barrel end

      ///////////////////////////////////////////////
      /////// Barrel + Endcap //////////////////////
      //////////////////////////////////////////////
      if ( (eg_gen_et[i]>150.) &&  (eg_gen_et[i]<500.0) ) {
	den_ele_geneta->Fill(eg_gen_eta[i]);

	//all cuts
	if ( passL1 && 
	     (hoe<hoe_cut)  
	     // (eg_sigma2vv[i]<vv_cut_EE) &&  
	     //(eg_sigma2ww[i]<ww_cut_EE) && 
	     // (eg_sigmaIEtaIEta[i]<sieie_cut_EB) 
	     ) {
	  num_ele_geneta_all->Fill(eg_gen_eta[i]);
	}
      }



    } // egamma loop

  } //event loop

  den_ele_geneta->Write();
  num_ele_geneta_all->Write();
  num_ele_genpt_all_EE->Write();
  den_ele_genpt_EE->Write();
  den_ele_genpt_EB->Write();
  num_ele_genpt_all_EB->Write();

}
////////
// Single photon isolated, EB-only
////////
void CutOpt::SinglePhoIsoEBonly()
{
  if (fChain == 0) return;
  
  TFile* outputFile = new TFile("histSinglePhoIsoEBonly_Sig_Zprime.root","RECREATE");

  const Int_t NBINS_eta = 17;
  Double_t edges_eta[NBINS_eta + 1] = {-3.0,-2.7,-2.4,-2.0,-1.56,-1.44,-1.0,-0.6,-0.2,0.2,0.6,1.0,1.44,1.56,2.0,2.4,2.7,3.0};

  //  const Int_t NBINS_pt = 6;
  //  Double_t edges_pt[NBINS_pt + 1] = {200,250,300,400,600,800,1000};

  const Int_t NBINS_pt = 16;
  //  Double_t edges_pt[NBINS_pt + 1] = {0,20,40,60,80,100,120,140,160,180,200,250,300,350,400,450,500,600,700,800,900,1000,1200,1400,1600,1800,2000};
  Double_t edges_pt[NBINS_pt + 1] = {20,30,40,50,60,70,80,90,100,150,200,250,300,350,400,450,500};

  TH1D*  den_ele_geneta = new TH1D("den_ele_geneta", "den_ele_geneta", NBINS_eta, edges_eta); 
  TH1D*  num_ele_geneta_all = new TH1D("num_ele_geneta_all", "num_ele_geneta_all", NBINS_eta, edges_eta); 

  TH1D*  den_ele_genpt_EB = new TH1D("den_ele_genpt_EB", "den_ele_genpt_EB", NBINS_pt, edges_pt); 
  TH1D*  num_ele_genpt_all_EB = new TH1D("num_ele_genpt_all_EB", "num_ele_genpt_all_EB", NBINS_pt, edges_pt);

  Long64_t nentries = fChain->GetEntriesFast();
  
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
    for (int i=0; i<nrEgs; i++) {

      //// Barrel cuts                                                                                                                                                      
      float sieie_cut_EB= 0.01;
      float ecaliso_cut_EB= 3.0 + (0.02*eg_et[i]);
      float hcaliso_cut_EB= 5.3 + (0.02*eg_et[i]);
      float hoe_EB = (eg_hcalHForHoverE[i])/(eg_energy[i]);
      float hoe_cut_EB = 0.05;
      ////end of barrel cuts                                                                                                                                                

      ////define L1 pass/fail                                                                                                                                               

      bool passL1 = false;

      if (  fabs(eg_eta[i]) > 2.8 ) {
        passL1=true;
      }

      if ( (fabs(eg_eta[i])<2.8) &&  (eg_l1pho_et[i]>0) && (eg_l1pho_etThres[i]>51.) && (eg_l1pho_passQual[i]) ) {
        passL1=true;
      }

      if ( (fabs(eg_eta[i])<2.4) && (eg_l1pho_et[i]>0) && (eg_l1pho_etThresIso[i]>36.) && (eg_l1pho_passQual[i]) && (eg_l1pho_passIsol[i]) ) {
        passL1=true;
      }


      //barrel start
      if ( (eg_gen_et[i]>10.0) && ((fabs(eg_eta[i]))<1.44) ) {
	den_ele_genpt_EB->Fill(eg_gen_et[i]);

	//all cuts
	if ( 
	    passL1  && 
	    (hoe_EB<hoe_cut_EB)  &&
	    (eg_sigmaIEtaIEta[i]<sieie_cut_EB)  &&
	    (eg_ecaliso[i]<ecaliso_cut_EB ) &&
	    (eg_hcalPFIsol_default[i]<hcaliso_cut_EB) 
	     ) {
	  num_ele_genpt_all_EB->Fill(eg_gen_et[i]);
	}
      } //barrel end

      ///////////////////////////////////////////////
      /////// Barrel, for eta plot /////////////////
      //////////////////////////////////////////////
      if ( (eg_gen_et[i]>100.) &&  (eg_gen_et[i]<500.0) ) {
	den_ele_geneta->Fill(eg_gen_eta[i]);

	//all cuts
	if ( passL1 && 
	     ( fabs(eg_eta[i]) <1.479 ) &&
	     (hoe_EB<hoe_cut_EB) &&
	     (eg_ecaliso[i]<ecaliso_cut_EB ) &&
	     (eg_hcalPFIsol_default[i]<hcaliso_cut_EB) &&
	     (eg_sigmaIEtaIEta[i]<sieie_cut_EB)	
	     ) {
	  num_ele_geneta_all->Fill(eg_gen_eta[i]);
	}
      }



    } // egamma loop

  } //event loop

  den_ele_geneta->Write();
  num_ele_geneta_all->Write();
  den_ele_genpt_EB->Write();
  num_ele_genpt_all_EB->Write();

}



/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////// DoubleEle HLT
/////////////////////////////////////////////////////////////////////////////////////////////////////////
void CutOpt::DoubleEle()
{

  if (fChain == 0) return;
  
  TFile* outputFile  = new TFile("hist_Zprime_DoubleEle.root","RECREATE");

  const Int_t NBINS_eta = 17;
  Double_t edges_eta[NBINS_eta + 1] = {-3.0,-2.7,-2.4,-2.0,-1.56,-1.44,-1.0,-0.6,-0.2,0.2,0.6,1.0,1.44,1.56,2.0,2.4,2.7,3.0};

  const Int_t NBINS_pt = 16;
  //  Double_t edges_pt[NBINS_pt + 1] = {20,25,30,35,40,45,50,55,60,65,70,80,100,150};
  Double_t edges_pt[NBINS_pt + 1] = {20,30,40,50,60,70,80,90,100,150,200,250,300,350,400,450,500};

  TH1D*  den_ele_geneta = new TH1D("den_ele_geneta", "den_ele_geneta", NBINS_eta, edges_eta); 
  TH1D*  num_ele_geneta_all_leg1 = new TH1D("num_ele_geneta_all_leg1", "num_ele_geneta_all_leg1", NBINS_eta, edges_eta); 
  TH1D*  num_ele_geneta_all_leg2 = new TH1D("num_ele_geneta_all_leg2", "num_ele_geneta_all_leg2", NBINS_eta, edges_eta); 

  TH1D*  den_ele_genpt_EB = new TH1D("den_ele_genpt_EB", "den_ele_genpt_EB", NBINS_pt, edges_pt);
  TH1D*  den_ele_genpt_EE = new TH1D("den_ele_genpt_EE", "den_ele_genpt_EE", NBINS_pt, edges_pt);

  TH1D*  num_ele_genpt_EB_leg1 = new TH1D("num_ele_genpt_EB_leg1", "num_ele_genpt_EB_leg1", NBINS_pt, edges_pt);
  TH1D*  num_ele_genpt_EB_leg2 = new TH1D("num_ele_genpt_EB_leg2", "num_ele_genpt_EB_leg2", NBINS_pt, edges_pt);

  TH1D*  num_ele_genpt_EB_sieie = new TH1D("num_ele_genpt_EB_sieie", "num_ele_genpt_EB_sieie", NBINS_pt, edges_pt);

  TH1D*  num_ele_genpt_EB_L1_leg1 = new TH1D("num_ele_genpt_EB_L1_leg1", "num_ele_genpt_EB_L1_leg1", NBINS_pt, edges_pt);
  TH1D*  num_ele_genpt_EB_L1_leg2 = new TH1D("num_ele_genpt_EB_L1_leg2", "num_ele_genpt_EB_L1_leg2", NBINS_pt, edges_pt);

  TH1D*  num_ele_genpt_EE_L1_leg1 = new TH1D("num_ele_genpt_EE_L1_leg1", "num_ele_genpt_EE_L1_leg1", NBINS_pt, edges_pt);
  TH1D*  num_ele_genpt_EE_L1_leg2 = new TH1D("num_ele_genpt_EE_L1_leg2", "num_ele_genpt_EE_L1_leg2", NBINS_pt, edges_pt);

  TH1D*  num_ele_genpt_EE_leg1 = new TH1D("num_ele_genpt_EE_leg1", "num_ele_genpt_EE_leg1", NBINS_pt, edges_pt);
  TH1D*  num_ele_genpt_EE_leg2 = new TH1D("num_ele_genpt_EE_leg2", "num_ele_genpt_EE_leg2", NBINS_pt, edges_pt);

  Long64_t nentries = fChain->GetEntriesFast();
  
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
    for (int i=0; i<nrEgs; i++) {

      //// Barrel cuts // DoubleEle HLT
      float pms2_cut_EB= 55.0; 
      float sieie_cut_EB= 0.014; //0.013; 
      float hoe_EB = (eg_hcalHForHoverE[i])/(eg_energy[i]);
      float hoe_cut_EB = 0.19; 
      ////end of barrel cuts
      
      ////endcap cuts // DoubleEle HLT  
      float hoe_EE = (eg_hgcalHForHoverE[i])/(eg_energy[i]);
      float hoe_cut_EE = 0.19; //0.15 + (5.0/(eg_energy[i]));
      float vv_cut_EE = (0.8*0.8)+(0.0008*eg_et[i]); 
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
      if (  fabs(eg_eta[i]) > 2.8 ) {
	passL1_highpt=true;
      }
      if ( (fabs(eg_eta[i])<2.8) &&  (eg_l1pho_et[i]>0) && (eg_l1pho_etThres[i]>37.) && (eg_l1pho_passQual[i]) ) {
	passL1_highpt=true;
      }
      if ( (fabs(eg_eta[i])<2.4) && (eg_l1ele_et[i]>0) && (eg_l1ele_etThresNonIso[i]>25.) && (eg_l1ele_passQual[i]) ) {
	passL1_highpt=true;
      }
      /////////

      ////define L1 pass/fail, low pT leg
      bool passL1_lowpt = false;
      if (  fabs(eg_eta[i]) > 2.8 ) {
	passL1_lowpt=true;
      }
      if ( (fabs(eg_eta[i])<2.8) &&  (eg_l1pho_et[i]>0) && (eg_l1pho_etThres[i]>24.) && (eg_l1pho_passQual[i]) ) {
	passL1_lowpt=true;
      }
      if ( (fabs(eg_eta[i])<2.4) && (eg_l1ele_et[i]>0) && (eg_l1ele_etThresNonIso[i]>12.) && (eg_l1ele_passQual[i]) ) {
	passL1_lowpt=true;
      }

      ///endcap start    
      if ( (eg_gen_et[i]>10.0) && ((fabs(eg_eta[i]))>1.56) && ((fabs(eg_eta[i]))<2.70) ) {
	den_ele_genpt_EE->Fill(eg_gen_et[i]);

	if (passL1_highpt) {
	  num_ele_genpt_EE_L1_leg1->Fill(eg_gen_et[i]);
	}
	if (passL1_lowpt) {
	  num_ele_genpt_EE_L1_leg2->Fill(eg_gen_et[i]);
	}

	if (
            (passL1_highpt) &&
            (hoe<hoe_cut) &&
            (eg_sigma2vv[i]<vv_cut_EE) &&
            (eg_pms2[i]<pms2_cut) 
	    ) {
          num_ele_genpt_EE_leg1->Fill(eg_gen_et[i]);
        }

	if (
            (passL1_lowpt) &&
            (hoe<hoe_cut) &&
            (eg_sigma2vv[i]<vv_cut_EE) &&
            (eg_pms2[i]<pms2_cut) 
	    ) {
          num_ele_genpt_EE_leg2->Fill(eg_gen_et[i]);
        }
      } //endcap end

      //barrel start
      if ( (eg_gen_et[i]>10.0) && ((fabs(eg_eta[i]))<1.44) ) {
	den_ele_genpt_EB->Fill(eg_gen_et[i]);

	if ( eg_sigmaIEtaIEta[i]<sieie_cut_EB ) {
	  num_ele_genpt_EB_sieie->Fill(eg_gen_et[i]);
	}

	if (passL1_highpt) {
	  num_ele_genpt_EB_L1_leg1->Fill(eg_gen_et[i]);
	}

	if (passL1_lowpt) {
	  num_ele_genpt_EB_L1_leg2->Fill(eg_gen_et[i]);
	}
      
	if (
            (passL1_highpt) &&
            (hoe<hoe_cut) &&
            (eg_pms2[i]<pms2_cut) &&
            (eg_sigmaIEtaIEta[i]<sieie_cut_EB)
	    ) {
          num_ele_genpt_EB_leg1->Fill(eg_gen_et[i]);
        }

	if (
            (passL1_lowpt) &&
            (hoe<hoe_cut) &&
            (eg_pms2[i]<pms2_cut) &&
            (eg_sigmaIEtaIEta[i]<sieie_cut_EB)
	    ) {
          num_ele_genpt_EB_leg2->Fill(eg_gen_et[i]);
        }
      }

      ///////////////////////////////////////////////
      /////// Barrel + Endcap //////////////////////
      //////////////////////////////////////////////
      if ( (eg_gen_et[i]>40.) &&  (eg_gen_et[i]<500.0) ) {
	den_ele_geneta->Fill(eg_gen_eta[i]);
	//
	//all cuts leg1
	if ( 
	    (passL1_highpt) && 
	    (hoe<hoe_cut) &&  
	    (eg_sigma2vv[i]<vv_cut_EE) &&  
	    (eg_pms2[i]<pms2_cut) && 
	    (eg_sigmaIEtaIEta[i]<sieie_cut_EB) 
	       
	     ) {
	  num_ele_geneta_all_leg1->Fill(eg_gen_eta[i]);
	}

	if (
            (passL1_lowpt) &&
            (hoe<hoe_cut) &&
            (eg_sigma2vv[i]<vv_cut_EE) &&
            (eg_pms2[i]<pms2_cut) &&
            (eg_sigmaIEtaIEta[i]<sieie_cut_EB)
	    
	    ) {
          num_ele_geneta_all_leg2->Fill(eg_gen_eta[i]);
        }
	
      }//
    }//ele loop
    
  } //event loop

  den_ele_geneta->Write();
  num_ele_geneta_all_leg1->Write();
  num_ele_geneta_all_leg2->Write();

  den_ele_genpt_EB->Write();
  den_ele_genpt_EE->Write();

  num_ele_genpt_EB_leg1->Write();
  num_ele_genpt_EE_leg1->Write();

  num_ele_genpt_EB_leg2->Write();
  num_ele_genpt_EE_leg2->Write();

  num_ele_genpt_EB_sieie->Write();

  num_ele_genpt_EB_L1_leg1->Write();
  num_ele_genpt_EB_L1_leg2->Write();

  num_ele_genpt_EE_L1_leg1->Write();
  num_ele_genpt_EE_L1_leg2->Write();
}

/////// double photon HLT /////////

void CutOpt::DoublePhoton()
{

  if (fChain == 0) return;
  
  TFile* outputFile  = new TFile("hist_ZPrime_DoublePho.root","RECREATE");

  const Int_t NBINS_eta = 17;
  Double_t edges_eta[NBINS_eta + 1] = {-3.0,-2.7,-2.4,-2.0,-1.56,-1.44,-1.0,-0.6,-0.2,0.2,0.6,1.0,1.44,1.56,2.0,2.4,2.7,3.0};

  const Int_t NBINS_pt = 16;
  //  Double_t edges_pt[NBINS_pt + 1] = {20,25,30,35,40,45,50,55,60,65,70,80,100,150};
  //  Double_t edges_pt[NBINS_pt + 1] = {20,30,40,50,60,70,80,100,150,200,300,400,500,600,700,800,900,1000};
  Double_t edges_pt[NBINS_pt + 1] = {20,30,40,50,60,70,80,90,100,150,200,250,300,350,400,450,500};

  const Int_t NBINS_pu = 8;
  Double_t edges_pu[NBINS_pu + 1] = {150,170,180,190,200,210,220,230,250};

  TH1D*  den_ele_PU_EB = new TH1D("den_ele_PU_EB", "den_ele_PU_EB", NBINS_pu, edges_pu); 
  TH1D*  den_ele_PU_EE = new TH1D("den_ele_PU_EE", "den_ele_PU_EE", NBINS_pu, edges_pu); 

  TH1D*  num_ele_PU_EB_leg1_all = new TH1D("num_ele_PU_EB_leg1_all", "num_ele_PU_EB_leg1_all", NBINS_pu, edges_pu); 
  TH1D*  num_ele_PU_EE_leg1_all = new TH1D("num_ele_PU_EE_leg1_all", "num_ele_PU_EE_leg1_all", NBINS_pu, edges_pu); 

  TH1D*  num_ele_PU_EB_leg2_all = new TH1D("num_ele_PU_EB_leg2_all", "num_ele_PU_EB_leg2_all", NBINS_pu, edges_pu); 
  TH1D*  num_ele_PU_EE_leg2_all = new TH1D("num_ele_PU_EE_leg2_all", "num_ele_PU_EE_leg2_all", NBINS_pu, edges_pu); 

  TH1D*  den_ele_geneta = new TH1D("den_ele_geneta", "den_ele_geneta", NBINS_eta, edges_eta); 
  TH1D*  num_ele_geneta_all_leg1 = new TH1D("num_ele_geneta_all_leg1", "num_ele_geneta_all_leg1", NBINS_eta, edges_eta); 
  TH1D*  num_ele_geneta_all_leg2 = new TH1D("num_ele_geneta_all_leg2", "num_ele_geneta_all_leg2", NBINS_eta, edges_eta); 

  TH1D*  den_ele_genpt_EB = new TH1D("den_ele_genpt_EB", "den_ele_genpt_EB", NBINS_pt, edges_pt);
  TH1D*  den_ele_genpt_EE = new TH1D("den_ele_genpt_EE", "den_ele_genpt_EE", NBINS_pt, edges_pt);

  TH1D*  num_ele_genpt_EB_leg1_all = new TH1D("num_ele_genpt_EB_leg1_all", "num_ele_genpt_EB_leg1_all", NBINS_pt, edges_pt);
  TH1D*  num_ele_genpt_EB_leg2_all = new TH1D("num_ele_genpt_EB_leg2_all", "num_ele_genpt_EB_leg2_all", NBINS_pt, edges_pt);

  TH1D*  num_ele_genpt_EE_leg1_all = new TH1D("num_ele_genpt_EE_leg1_all", "num_ele_genpt_EE_leg1_all", NBINS_pt, edges_pt);
  TH1D*  num_ele_genpt_EE_leg2_all = new TH1D("num_ele_genpt_EE_leg2_all", "num_ele_genpt_EE_leg2_all", NBINS_pt, edges_pt);

  Long64_t nentries = fChain->GetEntriesFast();
  
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
    for (int i=0; i<nrEgs; i++) {

      //// Barrel cuts // DoublePhoton HLT
      float sieie_cut_EB= 0.0113;
      float hoe_EB = (eg_hcalHForHoverE[i])/(eg_energy[i]);
      float hoe_cut_EB = 0.18;
      float ecaliso_cut_EB= 7.5 + (0.02*eg_et[i]);
      float hcaliso_cut_EB=11 +(0.02*eg_et[i]);  //9.5;
      if ( ( fabs(eg_eta[i]) > 0.8 ) && ( fabs(eg_eta[i]) < 1.479 )  ) {
      	hcaliso_cut_EB = 12 +(0.02*eg_et[i]);
      }
      ////end of barrel cuts
      
      ////endcap cuts // DoublePhoton HLT  
      float hoe_EE = (eg_hgcalHForHoverE[i])/(eg_energy[i]);
      float hoe_cut_EE = 0.125; 
      float vv_cut_EE = (0.8*0.8) +(0.0008*eg_et[i]) ;
      float ww_cut_EE = (8*8) +(0.04*eg_et[i]) ;
      float hgcaliso_cut_EE = 140 +(0.05*eg_energy[i]);
      if ( ( fabs(eg_eta[i]) > 2.0 ) ) {
	hgcaliso_cut_EE = 370 +(0.05*eg_energy[i]); 
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
      
      ///
      //// Double StaEG 37,24 
      //// Double TkIsoPhoton 22, 12 
      ///
      ////define L1 pass/fail, high pT leg
      bool passL1_highpt = false;
      if (  fabs(eg_eta[i]) > 2.8 ) {
     	passL1_highpt=true;
      }
      if ( (fabs(eg_eta[i])<2.8) &&  (eg_l1pho_et[i]>0) && (eg_l1pho_etThres[i]>37.) && (eg_l1pho_passQual[i]) ) {
      //if ( (eg_l1pho_et[i]>0) && (eg_l1pho_etThres[i]>37.) && (eg_l1pho_passQual[i]) ) {
	passL1_highpt=true;
      }
      if ( (fabs(eg_eta[i])<2.4) && (eg_l1pho_et[i]>0) && (eg_l1pho_etThresIso[i]>22.) && (eg_l1pho_passQual[i]) && (eg_l1pho_passIsol[i]) ) {
	passL1_highpt=true;
      }
      /////////                                                                                                                                                                                               
      ////define L1 pass/fail, low pT leg                                                                                                                                                                     
      bool passL1_lowpt = false;
      if (  fabs(eg_eta[i]) > 2.8 ) {
	passL1_lowpt=true;
      }
      if ( (fabs(eg_eta[i])<2.8) &&  (eg_l1pho_et[i]>0) && (eg_l1pho_etThres[i]>24.) && (eg_l1pho_passQual[i]) ) {
	//if ( (eg_l1pho_et[i]>0) && (eg_l1pho_etThres[i]>24.) && (eg_l1pho_passQual[i]) ) {
	passL1_lowpt=true;
      }
      if ( (fabs(eg_eta[i])<2.4) && (eg_l1pho_et[i]>0) && (eg_l1pho_etThresIso[i]>12.) && (eg_l1pho_passQual[i]) && (eg_l1pho_passIsol[i]) ) {
	passL1_lowpt=true;
      }

      ///endcap start    
      if ( (eg_gen_et[i]>10.0) && ((fabs(eg_eta[i]))>1.56) && ((fabs(eg_eta[i]))<2.70) ) {
	den_ele_genpt_EE->Fill(eg_gen_et[i]);

	if (
            (passL1_highpt) &&
	    (hoe<hoe_cut) &&
	    (eg_sigma2vv[i]<vv_cut_EE) &&
	    (eg_sigma2ww[i]<ww_cut_EE) &&  
	    (eg_hgcaliso_layerclus[i]<hgcaliso_cut_EE) 
	    ) {
          num_ele_genpt_EE_leg1_all->Fill(eg_gen_et[i]);
        }

	if (
            (passL1_lowpt) &&
	    (hoe<hoe_cut) &&
	    (eg_sigma2vv[i]<vv_cut_EE) &&
	    (eg_sigma2ww[i]<ww_cut_EE) &&
	    (eg_hgcaliso_layerclus[i]<hgcaliso_cut_EE) 
	    ) {
          num_ele_genpt_EE_leg2_all->Fill(eg_gen_et[i]);
        }
      } //endcap end

      //////PU
      ///endcap start    
      if ( (eg_gen_et[i]>40.0) && ((fabs(eg_eta[i]))>1.56) && ((fabs(eg_eta[i]))<2.70) ) {
	den_ele_PU_EE->Fill(nrPtHats);

	if (
            (passL1_highpt) &&
            (hoe<hoe_cut) &&
            (eg_sigma2vv[i]<vv_cut_EE) &&
	    (eg_sigma2ww[i]<ww_cut_EE) &&  
	    (eg_hgcaliso_layerclus[i]<hgcaliso_cut_EE) 
	    ) {
          num_ele_PU_EE_leg1_all->Fill(nrPtHats);
        }

	if (
            (passL1_lowpt) &&
            (hoe<hoe_cut) &&
            (eg_sigma2vv[i]<vv_cut_EE) &&
	    (eg_sigma2ww[i]<ww_cut_EE)  &&
	    (eg_hgcaliso_layerclus[i]<hgcaliso_cut_EE) 
	    ) {
          num_ele_PU_EE_leg2_all->Fill(nrPtHats);
        }
      } //endcap end

      //barrel start
      if ( (eg_gen_et[i]>10.0) && ((fabs(eg_eta[i]))<1.44) ) {
	den_ele_genpt_EB->Fill(eg_gen_et[i]);

	if (
	    (passL1_highpt) &&
	    (hoe<hoe_cut) &&
	    (eg_sigmaIEtaIEta[i]<sieie_cut_EB) &&
	    (eg_ecaliso[i]<ecaliso_cut_EB ) &&
	    (eg_hcalPFIsol_default[i]<hcaliso_cut_EB) 
	    ) {
          num_ele_genpt_EB_leg1_all->Fill(eg_gen_et[i]);
        }

	if (
            (passL1_lowpt) &&
	    (hoe<hoe_cut) &&
	    (eg_sigmaIEtaIEta[i]<sieie_cut_EB) &&
	    (eg_ecaliso[i]<ecaliso_cut_EB ) &&
	    (eg_hcalPFIsol_default[i]<hcaliso_cut_EB) 
	    ) {
          num_ele_genpt_EB_leg2_all->Fill(eg_gen_et[i]);
        }
      }

      ////PU
      //barrel start
      if ( (eg_gen_et[i]>40.0) && ((fabs(eg_eta[i]))<1.44) ) {
	den_ele_PU_EB->Fill(nrPtHats);

	if (
            (passL1_highpt) &&
            (hoe<hoe_cut) &&
            (eg_sigmaIEtaIEta[i]<sieie_cut_EB) &&
	    (eg_ecaliso[i]<ecaliso_cut_EB ) &&
            (eg_hcalPFIsol_default[i]<hcaliso_cut_EB) 
	    ) {
          num_ele_PU_EB_leg1_all->Fill(nrPtHats);
        }

	if (
            (passL1_lowpt) &&
            (hoe<hoe_cut) &&
	    (eg_sigmaIEtaIEta[i]<sieie_cut_EB) &&
	    (eg_ecaliso[i]<ecaliso_cut_EB ) &&
            (eg_hcalPFIsol_default[i]<hcaliso_cut_EB) 
	    ) {
          num_ele_PU_EB_leg2_all->Fill(nrPtHats);
        }
      }


      ///////////////////////////////////////////////
      /////// Barrel + Endcap //////////////////////
      //////////////////////////////////////////////
      if ( (eg_gen_et[i]>40.) && (eg_gen_et[i]<500.0) ) {
	den_ele_geneta->Fill(eg_gen_eta[i]);
	//
	//all cuts leg1
	if ( 
	    (passL1_highpt) &&
	    // (eg_et[i]>40.0) &&
	    (hoe<hoe_cut) &&  
	    (eg_sigma2vv[i]<vv_cut_EE) &&  
	    (eg_sigma2ww[i]<ww_cut_EE) &&
            (eg_hgcaliso_layerclus[i]<hgcaliso_cut_EE) &&
            (eg_ecaliso[i]<ecaliso_cut_EB ) &&
            (eg_hcalPFIsol_default[i]<hcaliso_cut_EB) &&
	    (eg_sigmaIEtaIEta[i]<sieie_cut_EB) 
	       
	     ) {
	  num_ele_geneta_all_leg1->Fill(eg_gen_eta[i]);
	}

	if (
            (passL1_lowpt) &&
	    // (eg_et[i]>40) &&
            (hoe<hoe_cut) &&
            (eg_sigma2vv[i]<vv_cut_EE) &&
	    (eg_sigma2ww[i]<ww_cut_EE) &&
            (eg_hgcaliso_layerclus[i]<hgcaliso_cut_EE) &&
            (eg_ecaliso[i]<ecaliso_cut_EB ) &&
            (eg_hcalPFIsol_default[i]<hcaliso_cut_EB) &&
            (eg_sigmaIEtaIEta[i]<sieie_cut_EB)
	    
	    ) {
          num_ele_geneta_all_leg2->Fill(eg_gen_eta[i]);
        }
	
      }//
    }//ele loop
    
  } //event loop

  den_ele_geneta->Write();
  num_ele_geneta_all_leg1->Write();
  num_ele_geneta_all_leg2->Write();

  den_ele_genpt_EB->Write();
  den_ele_genpt_EE->Write();

  num_ele_genpt_EB_leg1_all->Write();
  num_ele_genpt_EE_leg1_all->Write();

  num_ele_genpt_EB_leg2_all->Write();
  num_ele_genpt_EE_leg2_all->Write();

  den_ele_PU_EB->Write();
  den_ele_PU_EE->Write();

  num_ele_PU_EB_leg1_all->Write();
  num_ele_PU_EE_leg1_all->Write();

  num_ele_PU_EB_leg2_all->Write();
  num_ele_PU_EE_leg2_all->Write();

}



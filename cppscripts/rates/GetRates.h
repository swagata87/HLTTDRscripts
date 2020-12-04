/////////////////
//file:   /eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/Upgrade/EGM_PhaseII/mc/11_1_4/ntup/QCD_Pt_170to300_TuneCP5_14TeV_pythia8__Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_Nov16thEM.root
//////////////////////////////////////////////////////////

#ifndef GetRates_h
#define GetRates_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>

// Header file for the classes stored in the TTree if any.

class GetRates {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          runnr;
   UInt_t          lumiSec;
   UInt_t          eventnr;
   Float_t         weightV1;
   Float_t         weightV2;
   Float_t         weightV2NoEM;
   UInt_t          nrEgs;
   Float_t         eg_hltisov72[8];   //[nrEgs]
   Float_t         eg_trkMissHits[8];   //[nrEgs]
   Float_t         eg_trkDEtaSeed[8];   //[nrEgs]
   Float_t         eg_ecaliso_validation[8];   //[nrEgs]
   Float_t         eg_hcalHForHoverE[8];   //[nrEgs]
   Float_t         eg_trkValidHits[8];   //[nrEgs]
   Float_t         eg_sigma2yz[8];   //[nrEgs]
   Float_t         eg_hcalH_dep1[8];   //[nrEgs]
   Float_t         eg_hSumForHoverE[8];   //[nrEgs]
   Float_t         eg_sigmaIEtaIEta[8];   //[nrEgs]
   Float_t         eg_trkIsolV6_default[8];   //[nrEgs]
   Float_t         eg_trkDEta[8];   //[nrEgs]
   Int_t           eg_nLayerIT[8];   //[nrEgs]
   Float_t         eg_sigma2xy[8];   //[nrEgs]
   Float_t         eg_hcalH_dep3[8];   //[nrEgs]
   Float_t         eg_sigma2xx[8];   //[nrEgs]
   Float_t         eg_invESeedInvP[8];   //[nrEgs]
   Float_t         eg_energy[8];   //[nrEgs]
   Float_t         eg_sigma2zz[8];   //[nrEgs]
   Float_t         eg_sigma2vv[8];   //[nrEgs]
   Float_t         eg_ecalPFIsol_default[8];   //[nrEgs]
   Float_t         eg_sigma2uu[8];   //[nrEgs]
   Float_t         eg_hgcaliso[8];   //[nrEgs]
   Float_t         eg_hltisov6[8];   //[nrEgs]
   Float_t         eg_hgcaliso_validation[8];   //[nrEgs]
   Float_t         eg_sigma2zx[8];   //[nrEgs]
   Float_t         eg_phi[8];   //[nrEgs]
   Float_t         eg_trkIsolV72_default[8];   //[nrEgs]
   Float_t         eg_hcalH_dep4[8];   //[nrEgs]
   Float_t         eg_invEInvP[8];   //[nrEgs]
   Int_t           eg_nLayerOT[8];   //[nrEgs]
   Float_t         eg_normChi2[8];   //[nrEgs]
   Float_t         eg_rVar[8];   //[nrEgs]
   Float_t         eg_trkChi2_default[8];   //[nrEgs]
   Float_t         eg_hgcalPFIsol_default[8];   //[nrEgs]
   Float_t         eg_hgcaliso_layerclus[8];   //[nrEgs]
   Float_t         eg_trkDPhi[8];   //[nrEgs]
   Float_t         eg_hcalH_dep2[8];   //[nrEgs]
   Float_t         eg_hgcalHForHoverE[8];   //[nrEgs]
   Float_t         eg_trkIsolV0[8];   //[nrEgs]
   Float_t         eg_hltisov6_validation[8];   //[nrEgs]
   Float_t         eg_hForHoverE[8];   //[nrEgs]
   Float_t         eg_ecaliso[8];   //[nrEgs]
   Float_t         eg_sigma2yy[8];   //[nrEgs]
   Int_t           eg_nGsf[8];   //[nrEgs]
   Float_t         eg_pms2_default[8];   //[nrEgs]
   Float_t         eg_hcaliso[8];   //[nrEgs]
   Float_t         eg_l1iso[8];   //[nrEgs]
   Float_t         eg_hcaliso_validation[8];   //[nrEgs]
   Float_t         eg_eta[8];   //[nrEgs]
   Float_t         eg_hltisov72_validation[8];   //[nrEgs]
   Float_t         eg_et[8];   //[nrEgs]
   Float_t         eg_hcalPFIsol_default[8];   //[nrEgs]
   Float_t         eg_pms2[8];   //[nrEgs]
   Float_t         eg_sigma2ww[8];   //[nrEgs]
   Float_t         eg_gen_phi[8];   //[nrEgs]
   Float_t         eg_gen_pt[8];   //[nrEgs]
   Float_t         eg_gen_energy[8];   //[nrEgs]
   Float_t         eg_gen_vz[8];   //[nrEgs]
   Float_t         eg_gen_et[8];   //[nrEgs]
   Float_t         eg_gen_eta[8];   //[nrEgs]
   Float_t         eg_l1pho_phi[8];   //[nrEgs]
   UChar_t         eg_l1pho_passIsol[8];   //[nrEgs]
   Float_t         eg_l1pho_etThresIso[8];   //[nrEgs]
   Float_t         eg_l1pho_trkIsol[8];   //[nrEgs]
   Float_t         eg_l1pho_etThresNonIso[8];   //[nrEgs]
   UChar_t         eg_l1pho_passQual[8];   //[nrEgs]
   Float_t         eg_l1pho_etThres[8];   //[nrEgs]
   Float_t         eg_l1pho_trkIsolPV[8];   //[nrEgs]
   Float_t         eg_l1pho_hwQual[8];   //[nrEgs]
   Float_t         eg_l1pho_et[8];   //[nrEgs]
   Float_t         eg_l1pho_eta[8];   //[nrEgs]
   Float_t         eg_l1ele_trkCurve[8];   //[nrEgs]
   Float_t         eg_l1ele_trkIsolPV[8];   //[nrEgs]
   Float_t         eg_l1ele_phi[8];   //[nrEgs]
   Float_t         eg_l1ele_etThresNonIso[8];   //[nrEgs]
   Float_t         eg_l1ele_etThresIso[8];   //[nrEgs]
   Float_t         eg_l1ele_hwQual[8];   //[nrEgs]
   Float_t         eg_l1ele_trkIsol[8];   //[nrEgs]
   UChar_t         eg_l1ele_passIsol[8];   //[nrEgs]
   Float_t         eg_l1ele_et[8];   //[nrEgs]
   Float_t         eg_l1ele_vz[8];   //[nrEgs]
   UChar_t         eg_l1ele_passQual[8];   //[nrEgs]
   Float_t         eg_l1ele_etThres[8];   //[nrEgs]
   Float_t         eg_l1ele_eta[8];   //[nrEgs]
   UChar_t         path_Gen_QCDMuGenFilter;
   UChar_t         path_Gen_QCDBCToEFilter;
   UChar_t         path_Gen_QCDEmEnrichingFilter;
   UChar_t         path_Gen_QCDEmEnrichingNoBCToEFilter;

   // List of branches
   TBranch        *b_runnr;   //!
   TBranch        *b_lumiSec;   //!
   TBranch        *b_eventnr;   //!
   TBranch        *b_weightV1;   //!
   TBranch        *b_weightV2;   //!
   TBranch        *b_weightV2NoEM;   //!
   TBranch        *b_nrEgs;   //!
   TBranch        *b_eg_hltisov72;   //!
   TBranch        *b_eg_trkMissHits;   //!
   TBranch        *b_eg_trkDEtaSeed;   //!
   TBranch        *b_eg_ecaliso_validation;   //!
   TBranch        *b_eg_hcalHForHoverE;   //!
   TBranch        *b_eg_trkValidHits;   //!
   TBranch        *b_eg_sigma2yz;   //!
   TBranch        *b_eg_hcalH_dep1;   //!
   TBranch        *b_eg_hSumForHoverE;   //!
   TBranch        *b_eg_sigmaIEtaIEta;   //!
   TBranch        *b_eg_trkIsolV6_default;   //!
   TBranch        *b_eg_trkDEta;   //!
   TBranch        *b_eg_nLayerIT;   //!
   TBranch        *b_eg_sigma2xy;   //!
   TBranch        *b_eg_hcalH_dep3;   //!
   TBranch        *b_eg_sigma2xx;   //!
   TBranch        *b_eg_invESeedInvP;   //!
   TBranch        *b_eg_energy;   //!
   TBranch        *b_eg_sigma2zz;   //!
   TBranch        *b_eg_sigma2vv;   //!
   TBranch        *b_eg_ecalPFIsol_default;   //!
   TBranch        *b_eg_sigma2uu;   //!
   TBranch        *b_eg_hgcaliso;   //!
   TBranch        *b_eg_hltisov6;   //!
   TBranch        *b_eg_hgcaliso_validation;   //!
   TBranch        *b_eg_sigma2zx;   //!
   TBranch        *b_eg_phi;   //!
   TBranch        *b_eg_trkIsolV72_default;   //!
   TBranch        *b_eg_hcalH_dep4;   //!
   TBranch        *b_eg_invEInvP;   //!
   TBranch        *b_eg_nLayerOT;   //!
   TBranch        *b_eg_normChi2;   //!
   TBranch        *b_eg_rVar;   //!
   TBranch        *b_eg_trkChi2_default;   //!
   TBranch        *b_eg_hgcalPFIsol_default;   //!
   TBranch        *b_eg_hgcaliso_layerclus;   //!
   TBranch        *b_eg_trkDPhi;   //!
   TBranch        *b_eg_hcalH_dep2;   //!
   TBranch        *b_eg_hgcalHForHoverE;   //!
   TBranch        *b_eg_trkIsolV0;   //!
   TBranch        *b_eg_hltisov6_validation;   //!
   TBranch        *b_eg_hForHoverE;   //!
   TBranch        *b_eg_ecaliso;   //!
   TBranch        *b_eg_sigma2yy;   //!
   TBranch        *b_eg_nGsf;   //!
   TBranch        *b_eg_pms2_default;   //!
   TBranch        *b_eg_hcaliso;   //!
   TBranch        *b_eg_l1iso;   //!
   TBranch        *b_eg_hcaliso_validation;   //!
   TBranch        *b_eg_eta;   //!
   TBranch        *b_eg_hltisov72_validation;   //!
   TBranch        *b_eg_et;   //!
   TBranch        *b_eg_hcalPFIsol_default;   //!
   TBranch        *b_eg_pms2;   //!
   TBranch        *b_eg_sigma2ww;   //!
   TBranch        *b_eg_gen_phi;   //!
   TBranch        *b_eg_gen_pt;   //!
   TBranch        *b_eg_gen_energy;   //!
   TBranch        *b_eg_gen_vz;   //!
   TBranch        *b_eg_gen_et;   //!
   TBranch        *b_eg_gen_eta;   //!
   TBranch        *b_eg_l1pho_phi;   //!
   TBranch        *b_eg_l1pho_passIsol;   //!
   TBranch        *b_eg_l1pho_etThresIso;   //!
   TBranch        *b_eg_l1pho_trkIsol;   //!
   TBranch        *b_eg_l1pho_etThresNonIso;   //!
   TBranch        *b_eg_l1pho_passQual;   //!
   TBranch        *b_eg_l1pho_etThres;   //!
   TBranch        *b_eg_l1pho_trkIsolPV;   //!
   TBranch        *b_eg_l1pho_hwQual;   //!
   TBranch        *b_eg_l1pho_et;   //!
   TBranch        *b_eg_l1pho_eta;   //!
   TBranch        *b_eg_l1ele_trkCurve;   //!
   TBranch        *b_eg_l1ele_trkIsolPV;   //!
   TBranch        *b_eg_l1ele_phi;   //!
   TBranch        *b_eg_l1ele_etThresNonIso;   //!
   TBranch        *b_eg_l1ele_etThresIso;   //!
   TBranch        *b_eg_l1ele_hwQual;   //!
   TBranch        *b_eg_l1ele_trkIsol;   //!
   TBranch        *b_eg_l1ele_passIsol;   //!
   TBranch        *b_eg_l1ele_et;   //!
   TBranch        *b_eg_l1ele_vz;   //!
   TBranch        *b_eg_l1ele_passQual;   //!
   TBranch        *b_eg_l1ele_etThres;   //!
   TBranch        *b_eg_l1ele_eta;   //!
   TBranch        *b_path_Gen_QCDMuGenFilter;   //!
   TBranch        *b_path_Gen_QCDBCToEFilter;   //!
   TBranch        *b_path_Gen_QCDEmEnrichingFilter;   //!
   TBranch        *b_path_Gen_QCDEmEnrichingNoBCToEFilter;   //!

   GetRates(TTree *tree=0);
   virtual ~GetRates();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef GetRates_cxx
GetRates::GetRates(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/Upgrade/EGM_PhaseII/mc/11_1_4/ntup/QCD_Pt_170to300_TuneCP5_14TeV_pythia8__Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_Nov16thEM.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/Upgrade/EGM_PhaseII/mc/11_1_4/ntup/QCD_Pt_170to300_TuneCP5_14TeV_pythia8__Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_Nov16thEM.root");
      }
      f->GetObject("egHLTTree",tree);
      std::cout << "/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/Upgrade/EGM_PhaseII/mc/11_1_4/ntup/QCD_Pt_170to300_TuneCP5_14TeV_pythia8__Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_Nov16thEM.root" << std::endl;

   }
   Init(tree);
}

GetRates::~GetRates()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t GetRates::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t GetRates::LoadTree(Long64_t entry)
{
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

void GetRates::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("runnr", &runnr, &b_runnr);
   fChain->SetBranchAddress("lumiSec", &lumiSec, &b_lumiSec);
   fChain->SetBranchAddress("eventnr", &eventnr, &b_eventnr);
   fChain->SetBranchAddress("weightV1", &weightV1, &b_weightV1);
   fChain->SetBranchAddress("weightV2", &weightV2, &b_weightV2);
   fChain->SetBranchAddress("weightV2NoEM", &weightV2NoEM, &b_weightV2NoEM);
   fChain->SetBranchAddress("nrEgs", &nrEgs, &b_nrEgs);
   fChain->SetBranchAddress("eg_hltisov72", eg_hltisov72, &b_eg_hltisov72);
   fChain->SetBranchAddress("eg_trkMissHits", eg_trkMissHits, &b_eg_trkMissHits);
   fChain->SetBranchAddress("eg_trkDEtaSeed", eg_trkDEtaSeed, &b_eg_trkDEtaSeed);
   fChain->SetBranchAddress("eg_ecaliso_validation", eg_ecaliso_validation, &b_eg_ecaliso_validation);
   fChain->SetBranchAddress("eg_hcalHForHoverE", eg_hcalHForHoverE, &b_eg_hcalHForHoverE);
   fChain->SetBranchAddress("eg_trkValidHits", eg_trkValidHits, &b_eg_trkValidHits);
   fChain->SetBranchAddress("eg_sigma2yz", eg_sigma2yz, &b_eg_sigma2yz);
   fChain->SetBranchAddress("eg_hcalH_dep1", eg_hcalH_dep1, &b_eg_hcalH_dep1);
   fChain->SetBranchAddress("eg_hSumForHoverE", eg_hSumForHoverE, &b_eg_hSumForHoverE);
   fChain->SetBranchAddress("eg_sigmaIEtaIEta", eg_sigmaIEtaIEta, &b_eg_sigmaIEtaIEta);
   fChain->SetBranchAddress("eg_trkIsolV6_default", eg_trkIsolV6_default, &b_eg_trkIsolV6_default);
   fChain->SetBranchAddress("eg_trkDEta", eg_trkDEta, &b_eg_trkDEta);
   fChain->SetBranchAddress("eg_nLayerIT", eg_nLayerIT, &b_eg_nLayerIT);
   fChain->SetBranchAddress("eg_sigma2xy", eg_sigma2xy, &b_eg_sigma2xy);
   fChain->SetBranchAddress("eg_hcalH_dep3", eg_hcalH_dep3, &b_eg_hcalH_dep3);
   fChain->SetBranchAddress("eg_sigma2xx", eg_sigma2xx, &b_eg_sigma2xx);
   fChain->SetBranchAddress("eg_invESeedInvP", eg_invESeedInvP, &b_eg_invESeedInvP);
   fChain->SetBranchAddress("eg_energy", eg_energy, &b_eg_energy);
   fChain->SetBranchAddress("eg_sigma2zz", eg_sigma2zz, &b_eg_sigma2zz);
   fChain->SetBranchAddress("eg_sigma2vv", eg_sigma2vv, &b_eg_sigma2vv);
   fChain->SetBranchAddress("eg_ecalPFIsol_default", eg_ecalPFIsol_default, &b_eg_ecalPFIsol_default);
   fChain->SetBranchAddress("eg_sigma2uu", eg_sigma2uu, &b_eg_sigma2uu);
   fChain->SetBranchAddress("eg_hgcaliso", eg_hgcaliso, &b_eg_hgcaliso);
   fChain->SetBranchAddress("eg_hltisov6", eg_hltisov6, &b_eg_hltisov6);
   fChain->SetBranchAddress("eg_hgcaliso_validation", eg_hgcaliso_validation, &b_eg_hgcaliso_validation);
   fChain->SetBranchAddress("eg_sigma2zx", eg_sigma2zx, &b_eg_sigma2zx);
   fChain->SetBranchAddress("eg_phi", eg_phi, &b_eg_phi);
   fChain->SetBranchAddress("eg_trkIsolV72_default", eg_trkIsolV72_default, &b_eg_trkIsolV72_default);
   fChain->SetBranchAddress("eg_hcalH_dep4", eg_hcalH_dep4, &b_eg_hcalH_dep4);
   fChain->SetBranchAddress("eg_invEInvP", eg_invEInvP, &b_eg_invEInvP);
   fChain->SetBranchAddress("eg_nLayerOT", eg_nLayerOT, &b_eg_nLayerOT);
   fChain->SetBranchAddress("eg_normChi2", eg_normChi2, &b_eg_normChi2);
   fChain->SetBranchAddress("eg_rVar", eg_rVar, &b_eg_rVar);
   fChain->SetBranchAddress("eg_trkChi2_default", eg_trkChi2_default, &b_eg_trkChi2_default);
   fChain->SetBranchAddress("eg_hgcalPFIsol_default", eg_hgcalPFIsol_default, &b_eg_hgcalPFIsol_default);
   fChain->SetBranchAddress("eg_hgcaliso_layerclus", eg_hgcaliso_layerclus, &b_eg_hgcaliso_layerclus);
   fChain->SetBranchAddress("eg_trkDPhi", eg_trkDPhi, &b_eg_trkDPhi);
   fChain->SetBranchAddress("eg_hcalH_dep2", eg_hcalH_dep2, &b_eg_hcalH_dep2);
   fChain->SetBranchAddress("eg_hgcalHForHoverE", eg_hgcalHForHoverE, &b_eg_hgcalHForHoverE);
   fChain->SetBranchAddress("eg_trkIsolV0", eg_trkIsolV0, &b_eg_trkIsolV0);
   fChain->SetBranchAddress("eg_hltisov6_validation", eg_hltisov6_validation, &b_eg_hltisov6_validation);
   fChain->SetBranchAddress("eg_hForHoverE", eg_hForHoverE, &b_eg_hForHoverE);
   fChain->SetBranchAddress("eg_ecaliso", eg_ecaliso, &b_eg_ecaliso);
   fChain->SetBranchAddress("eg_sigma2yy", eg_sigma2yy, &b_eg_sigma2yy);
   fChain->SetBranchAddress("eg_nGsf", eg_nGsf, &b_eg_nGsf);
   fChain->SetBranchAddress("eg_pms2_default", eg_pms2_default, &b_eg_pms2_default);
   fChain->SetBranchAddress("eg_hcaliso", eg_hcaliso, &b_eg_hcaliso);
   fChain->SetBranchAddress("eg_l1iso", eg_l1iso, &b_eg_l1iso);
   fChain->SetBranchAddress("eg_hcaliso_validation", eg_hcaliso_validation, &b_eg_hcaliso_validation);
   fChain->SetBranchAddress("eg_eta", eg_eta, &b_eg_eta);
   fChain->SetBranchAddress("eg_hltisov72_validation", eg_hltisov72_validation, &b_eg_hltisov72_validation);
   fChain->SetBranchAddress("eg_et", eg_et, &b_eg_et);
   fChain->SetBranchAddress("eg_hcalPFIsol_default", eg_hcalPFIsol_default, &b_eg_hcalPFIsol_default);
   fChain->SetBranchAddress("eg_pms2", eg_pms2, &b_eg_pms2);
   fChain->SetBranchAddress("eg_sigma2ww", eg_sigma2ww, &b_eg_sigma2ww);
   fChain->SetBranchAddress("eg_gen_phi", eg_gen_phi, &b_eg_gen_phi);
   fChain->SetBranchAddress("eg_gen_pt", eg_gen_pt, &b_eg_gen_pt);
   fChain->SetBranchAddress("eg_gen_energy", eg_gen_energy, &b_eg_gen_energy);
   fChain->SetBranchAddress("eg_gen_vz", eg_gen_vz, &b_eg_gen_vz);
   fChain->SetBranchAddress("eg_gen_et", eg_gen_et, &b_eg_gen_et);
   fChain->SetBranchAddress("eg_gen_eta", eg_gen_eta, &b_eg_gen_eta);
   fChain->SetBranchAddress("eg_l1pho_phi", eg_l1pho_phi, &b_eg_l1pho_phi);
   fChain->SetBranchAddress("eg_l1pho_passIsol", eg_l1pho_passIsol, &b_eg_l1pho_passIsol);
   fChain->SetBranchAddress("eg_l1pho_etThresIso", eg_l1pho_etThresIso, &b_eg_l1pho_etThresIso);
   fChain->SetBranchAddress("eg_l1pho_trkIsol", eg_l1pho_trkIsol, &b_eg_l1pho_trkIsol);
   fChain->SetBranchAddress("eg_l1pho_etThresNonIso", eg_l1pho_etThresNonIso, &b_eg_l1pho_etThresNonIso);
   fChain->SetBranchAddress("eg_l1pho_passQual", eg_l1pho_passQual, &b_eg_l1pho_passQual);
   fChain->SetBranchAddress("eg_l1pho_etThres", eg_l1pho_etThres, &b_eg_l1pho_etThres);
   fChain->SetBranchAddress("eg_l1pho_trkIsolPV", eg_l1pho_trkIsolPV, &b_eg_l1pho_trkIsolPV);
   fChain->SetBranchAddress("eg_l1pho_hwQual", eg_l1pho_hwQual, &b_eg_l1pho_hwQual);
   fChain->SetBranchAddress("eg_l1pho_et", eg_l1pho_et, &b_eg_l1pho_et);
   fChain->SetBranchAddress("eg_l1pho_eta", eg_l1pho_eta, &b_eg_l1pho_eta);
   fChain->SetBranchAddress("eg_l1ele_trkCurve", eg_l1ele_trkCurve, &b_eg_l1ele_trkCurve);
   fChain->SetBranchAddress("eg_l1ele_trkIsolPV", eg_l1ele_trkIsolPV, &b_eg_l1ele_trkIsolPV);
   fChain->SetBranchAddress("eg_l1ele_phi", eg_l1ele_phi, &b_eg_l1ele_phi);
   fChain->SetBranchAddress("eg_l1ele_etThresNonIso", eg_l1ele_etThresNonIso, &b_eg_l1ele_etThresNonIso);
   fChain->SetBranchAddress("eg_l1ele_etThresIso", eg_l1ele_etThresIso, &b_eg_l1ele_etThresIso);
   fChain->SetBranchAddress("eg_l1ele_hwQual", eg_l1ele_hwQual, &b_eg_l1ele_hwQual);
   fChain->SetBranchAddress("eg_l1ele_trkIsol", eg_l1ele_trkIsol, &b_eg_l1ele_trkIsol);
   fChain->SetBranchAddress("eg_l1ele_passIsol", eg_l1ele_passIsol, &b_eg_l1ele_passIsol);
   fChain->SetBranchAddress("eg_l1ele_et", eg_l1ele_et, &b_eg_l1ele_et);
   fChain->SetBranchAddress("eg_l1ele_vz", eg_l1ele_vz, &b_eg_l1ele_vz);
   fChain->SetBranchAddress("eg_l1ele_passQual", eg_l1ele_passQual, &b_eg_l1ele_passQual);
   fChain->SetBranchAddress("eg_l1ele_etThres", eg_l1ele_etThres, &b_eg_l1ele_etThres);
   fChain->SetBranchAddress("eg_l1ele_eta", eg_l1ele_eta, &b_eg_l1ele_eta);
   fChain->SetBranchAddress("path_Gen_QCDMuGenFilter", &path_Gen_QCDMuGenFilter, &b_path_Gen_QCDMuGenFilter);
   fChain->SetBranchAddress("path_Gen_QCDBCToEFilter", &path_Gen_QCDBCToEFilter, &b_path_Gen_QCDBCToEFilter);
   fChain->SetBranchAddress("path_Gen_QCDEmEnrichingFilter", &path_Gen_QCDEmEnrichingFilter, &b_path_Gen_QCDEmEnrichingFilter);
   fChain->SetBranchAddress("path_Gen_QCDEmEnrichingNoBCToEFilter", &path_Gen_QCDEmEnrichingNoBCToEFilter, &b_path_Gen_QCDEmEnrichingNoBCToEFilter);
   Notify();
}

Bool_t GetRates::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void GetRates::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t GetRates::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef GetRates_cxx

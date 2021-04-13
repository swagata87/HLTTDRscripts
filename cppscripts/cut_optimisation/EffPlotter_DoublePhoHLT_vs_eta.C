#include <iostream>
#include <algorithm>
#include "TLatex.h"
#include <iomanip>
#include <vector>
#include "TPad.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TSystem.h"
#include "TImage.h"
#include "TKey.h"
#include "TH1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TObjArray.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TPostScript.h"
#include <TPaveStats.h>
#include "TLegend.h"
#include <TProfile.h>
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "THStack.h"
#include "TEfficiency.h"
#include "TGraphErrors.h"

void EffPlotter_DoublePhoHLT_vs_eta() {

  TFile *file_sig = new TFile("hist_ZPrime_DoublePho_L1seeded.root");
  
  //--Plotting Styles//
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadBottomMargin(0.12);  
  gStyle->SetPadTopMargin(0.05);   
  gStyle->SetTitleXSize(0.05);
  gStyle->SetTitleXOffset(1.05);
  gStyle->SetTitleYSize(0.05);
  gStyle->SetTitleYOffset(1.05);
  gStyle->SetOptStat();

  TH1D* den_sig = (TH1D*)file_sig->Get("den_ele_geneta");
  TH1D* num_sig_leg1 = (TH1D*)file_sig->Get("num_ele_geneta_all_leg1");
  TH1D* num_sig_leg2 = (TH1D*)file_sig->Get("num_ele_geneta_all_leg2");
  
  TEfficiency* pEff_sig_leg1 = 0;
  TEfficiency* pEff_sig_leg2 = 0;

  TCanvas* my_canvas1 = new TCanvas("canvas","canvas",800,700);
  my_canvas1->cd();

  if(TEfficiency::CheckConsistency(*num_sig_leg1,*den_sig)) {
    pEff_sig_leg1 = new TEfficiency(*num_sig_leg1,*den_sig);
    pEff_sig_leg1->SetLineColor(kAzure+3);
    pEff_sig_leg1->SetMarkerColor(kAzure+3);
    pEff_sig_leg1->SetMarkerStyle(21);
    pEff_sig_leg1->SetLineWidth(2);
    pEff_sig_leg1->SetTitle(";gen #eta ;per electron efficiency");
    pEff_sig_leg1->Draw();
    gPad->Update();
    auto graph = pEff_sig_leg1->GetPaintedGraph(); 
    graph->SetMinimum(0.00);
    graph->SetMaximum(1.01); 
    //    graph->GetXaxis()->SetLimits(20.,110.);
  }
  //

  if(TEfficiency::CheckConsistency(*num_sig_leg2,*den_sig)) {
    pEff_sig_leg2 = new TEfficiency(*num_sig_leg2,*den_sig);
    pEff_sig_leg2->SetLineColor(kAzure+1);
    pEff_sig_leg2->SetMarkerColor(kAzure+1);
    pEff_sig_leg2->SetMarkerStyle(22);
    pEff_sig_leg2->SetLineWidth(1);
    pEff_sig_leg2->Draw("same");
    gPad->Update();
  }


  TLegend *leg_example2 = new TLegend(0.25,0.2,0.85,0.55);
  leg_example2->SetHeader("L1+HLT Efficiency","C"); // option "C" allows to center the header
  leg_example2->AddEntry((TObject*)0, "HLT_Diphoton30_23_IsoCaloId", "");
  leg_example2->AddEntry((TObject*)0, "40<pT<500 GeV", "");
  leg_example2->AddEntry((TObject*)0, "Z'(ee) M=6 TeV, PU=200", "");
  // leg_example2->AddEntry((TObject*)0, "PU=200", "");
  leg_example2->SetFillColor(0);
  leg_example2->SetTextFont(42);
  leg_example2->SetBorderSize(0);
  leg_example2->AddEntry(pEff_sig_leg1, "leg1", "lp");
  leg_example2->AddEntry(pEff_sig_leg2, "leg2", "lp");
  leg_example2->Draw("same");

  TPaveText *t = new TPaveText(-3.8,1.015,1,1.05);
  t->AddText("#scale[1.2]{CMS} #bf{#it{Simulation Preliminary}}"); ((TText*)t->GetListOfLines()->Last())->SetTextColor(kBlack);
  t->SetTextAlign(11);
  t->SetFillColor(0);
  t->SetLineColor(0);
  t->SetBorderSize(1);
  t->Draw("same");
  
  my_canvas1->SetGrid();

  my_canvas1->Update();
  my_canvas1->SaveAs("DoublePho_fullEff_vs_geneta_leg_1_2.png");
}

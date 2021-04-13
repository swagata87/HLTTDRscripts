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

void EffPlotter_DoublePhoHLT_vs_pt() {

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

  TH1D* den_sig = (TH1D*)file_sig->Get("den_ele_genpt_EB");
  TH1D* num_sig_leg1 = (TH1D*)file_sig->Get("num_ele_genpt_EB_leg1_all");
  TH1D* num_sig_leg2 = (TH1D*)file_sig->Get("num_ele_genpt_EB_leg2_all");
  
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
    pEff_sig_leg1->SetTitle(";gen p_{T} [GeV];per electron efficiency");
    pEff_sig_leg1->Draw();
    gPad->Update();
    auto graph = pEff_sig_leg1->GetPaintedGraph(); 
    graph->SetMinimum(0.0);
    //    graph->SetMaximum(1.05); 
    graph->SetMaximum(1.01); 
    //    graph->GetXaxis()->SetLimits(20.,110.);
  }

  if(TEfficiency::CheckConsistency(*num_sig_leg2,*den_sig)) {
    pEff_sig_leg2 = new TEfficiency(*num_sig_leg2,*den_sig);
    pEff_sig_leg2->SetLineColor(kAzure+1);
    pEff_sig_leg2->SetMarkerColor(kAzure+1);
    pEff_sig_leg2->SetMarkerStyle(22);
    pEff_sig_leg2->SetLineWidth(1);
    //    pEff_sig_leg2->SetTitle("Signal efficiency vs gen pT;gen pT [GeV];efficiency");
    pEff_sig_leg2->Draw("same");
    gPad->Update();
  }

  TLegend *leg_example2 = new TLegend(0.4,0.2,0.9,0.5);
  leg_example2->SetHeader("L1+HLT efficiency","C"); // option "C" allows to center the header
  leg_example2->AddEntry((TObject*)0, "HLT_Diphoton30_23_IsoCaloId", "");
  //leg_example2->AddEntry((TObject*)0, "Endcap, 1.56<|#eta|<2.70", "");
  leg_example2->AddEntry((TObject*)0, "Barrel, |#eta|<1.44", "");
  leg_example2->AddEntry((TObject*)0, "Z'(ee) M=6 TeV, PU=200", "");
  // leg_example2->AddEntry((TObject*)0, "PU=200", "");
  leg_example2->SetFillColor(0);
  leg_example2->SetTextFont(42);
  leg_example2->SetBorderSize(0);
  leg_example2->AddEntry(pEff_sig_leg1, "leg1", "lp");
  leg_example2->AddEntry(pEff_sig_leg2, "leg2", "lp");
  leg_example2->Draw("same");
  
  TPaveText *t1 = new TPaveText(-6,1.015,500,1.05);
  t1->AddText("#scale[1.2]{CMS} #bf{#it{Simulation Preliminary}}"); ((TText*)t1->GetListOfLines()->Last())->SetTextColor(kBlack);
  t1->SetTextAlign(11);
  t1->SetFillColor(0);
  t1->SetLineColor(0);
  t1->SetBorderSize(1);
  t1->Draw("same");

  my_canvas1->SetGrid(1,1);

  my_canvas1->Update();

  ///////
  /////endcap
  ///////

  TH1D* den_sig_EE = (TH1D*)file_sig->Get("den_ele_genpt_EE");
  TH1D* num_sig_leg1_EE = (TH1D*)file_sig->Get("num_ele_genpt_EE_leg1_all");
  TH1D* num_sig_leg2_EE = (TH1D*)file_sig->Get("num_ele_genpt_EE_leg2_all");
  
  TEfficiency* pEff_sig_leg1_EE = 0;
  TEfficiency* pEff_sig_leg2_EE = 0;

  TCanvas* my_canvas2 = new TCanvas("canvas2","canvas2",800,700);
  my_canvas2->cd();

  if(TEfficiency::CheckConsistency(*num_sig_leg1_EE,*den_sig_EE)) {
    pEff_sig_leg1_EE = new TEfficiency(*num_sig_leg1_EE,*den_sig_EE);
    pEff_sig_leg1_EE->SetLineColor(kAzure+3);
    pEff_sig_leg1_EE->SetMarkerColor(kAzure+3);
    pEff_sig_leg1_EE->SetMarkerStyle(21);
    pEff_sig_leg1_EE->SetLineWidth(2);
    pEff_sig_leg1_EE->SetTitle(";gen p_{T} [GeV];per electron efficiency");
    pEff_sig_leg1_EE->Draw();
    gPad->Update();
    auto graph_EE = pEff_sig_leg1_EE->GetPaintedGraph(); 
    graph_EE->SetMinimum(0.0);
    //    graph->SetMaximum(1.05); 
    graph_EE->SetMaximum(1.01); 
    //    graph->GetXaxis()->SetLimits(20.,110.);
  }

  if(TEfficiency::CheckConsistency(*num_sig_leg2_EE,*den_sig_EE)) {
    pEff_sig_leg2_EE = new TEfficiency(*num_sig_leg2_EE,*den_sig_EE);
    pEff_sig_leg2_EE->SetLineColor(kAzure+1);
    pEff_sig_leg2_EE->SetMarkerColor(kAzure+1);
    pEff_sig_leg2_EE->SetMarkerStyle(22);
    pEff_sig_leg2_EE->SetLineWidth(1);
    //    pEff_sig_leg2->SetTitle("Signal efficiency vs gen pT;gen pT [GeV];efficiency");
    pEff_sig_leg2_EE->Draw("same");
    gPad->Update();
  }

  TLegend *leg_example22 = new TLegend(0.4,0.2,0.9,0.5);
  leg_example22->SetHeader("L1+HLT efficiency","C"); // option "C" allows to center the header
  leg_example22->AddEntry((TObject*)0, "HLT_Diphoton30_23_IsoCaloId", "");
  leg_example22->AddEntry((TObject*)0, "Endcap, 1.56<|#eta|<2.40", "");
  //leg_example2->AddEntry((TObject*)0, "Barrel, |#eta|<1.44", "");
  leg_example22->AddEntry((TObject*)0, "Z'(ee) M=6 TeV, PU=200", "");
  // leg_example2->AddEntry((TObject*)0, "PU=200", "");
  leg_example22->SetFillColor(0);
  leg_example22->SetTextFont(42);
  leg_example22->SetBorderSize(0);
  leg_example22->AddEntry(pEff_sig_leg1_EE, "leg1", "lp");
  leg_example22->AddEntry(pEff_sig_leg2_EE, "leg2", "lp");
  leg_example22->Draw("same");
  
  TPaveText *t = new TPaveText(-6,1.015,500,1.05);
  t->AddText("#scale[1.2]{CMS} #bf{#it{Simulation Preliminary}}"); ((TText*)t->GetListOfLines()->Last())->SetTextColor(kBlack);
  t->SetTextAlign(11);
  t->SetFillColor(0);
  t->SetLineColor(0);
  t->SetBorderSize(1);
  t->Draw("same");  

  my_canvas2->SetGrid(1,1);

  my_canvas2->Update();

  my_canvas1->SaveAs("DoublePho_fullEff_vs_genpt_EB_leg_1_2.png");
  my_canvas2->SaveAs("DoublePho_fullEff_vs_genpt_EE_leg_1_2.png");

}

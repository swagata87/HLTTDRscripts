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

void EffPlotter_var_vs_pt_tot_70_80() {

  TFile *file_sig80 = new TFile("hist_Zprime_WP80_L1seeded.root");
  TFile *file_sig70 = new TFile("hist_Zprime_WP70_L1seeded.root");
  
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

  TH1D* den_sig_EE80 = (TH1D*)file_sig80->Get("den_ele_genpt_EE");
  TH1D* num_sig_EE80 = (TH1D*)file_sig80->Get("num_ele_genpt_all_EE");
  TH1D* den_sig_EB80 = (TH1D*)file_sig80->Get("den_ele_genpt_EB");
  TH1D* num_sig_EB80 = (TH1D*)file_sig80->Get("num_ele_genpt_all_EB");

  TH1D* den_sig_EE70 = (TH1D*)file_sig70->Get("den_ele_genpt_EE");
  TH1D* num_sig_EE70 = (TH1D*)file_sig70->Get("num_ele_genpt_all_EE");
  TH1D* den_sig_EB70 = (TH1D*)file_sig70->Get("den_ele_genpt_EB");
  TH1D* num_sig_EB70 = (TH1D*)file_sig70->Get("num_ele_genpt_all_EB");
  
  TEfficiency* pEff_sig_EE80 = 0;
  TEfficiency* pEff_sig_EB80 = 0;

  TEfficiency* pEff_sig_EE70 = 0;
  TEfficiency* pEff_sig_EB70 = 0;

  TCanvas* my_canvas1 = new TCanvas("canvas","canvas",800,700);
  my_canvas1->cd();

  if(TEfficiency::CheckConsistency(*num_sig_EE80,*den_sig_EE80)) {
    pEff_sig_EE80 = new TEfficiency(*num_sig_EE80,*den_sig_EE80);
    pEff_sig_EE80->SetLineColor(kGreen+2);
    pEff_sig_EE80->SetMarkerColor(kGreen+2);
    pEff_sig_EE80->SetMarkerStyle(21);
    pEff_sig_EE80->SetLineWidth(2);
    pEff_sig_EE80->SetTitle(";gen p_{T} [GeV];per electron efficiency");
    pEff_sig_EE80->Draw();
    gPad->Update();
    auto graph = pEff_sig_EE80->GetPaintedGraph(); 
    graph->SetMinimum(0.0);
    graph->SetMaximum(1.01); 
    //    graph->GetXaxis()->SetLimits(20.,110.);
  }

  if(TEfficiency::CheckConsistency(*num_sig_EE70,*den_sig_EE70)) {
    pEff_sig_EE70 = new TEfficiency(*num_sig_EE70,*den_sig_EE70);
    pEff_sig_EE70->SetLineColor(kGreen);
    pEff_sig_EE70->SetMarkerColor(kGreen);
    pEff_sig_EE70->SetMarkerStyle(21);
    pEff_sig_EE70->SetLineWidth(2);
    //    pEff_sig_EE70->SetTitle("Signal efficiency vs gen pT;gen pT [GeV];per electron efficiency");
    pEff_sig_EE70->Draw("same");
    gPad->Update();
    //    auto graph = pEff_sig_EE70->GetPaintedGraph(); 
    //graph->SetMinimum(0.0);
    //    graph->SetMaximum(1.05); 
    // graph->SetMaximum(1.0); 
    //    graph->GetXaxis()->SetLimits(20.,110.);
  }

  TLegend *leg_example2 = new TLegend(0.5,0.25,0.9,0.5);
  leg_example2->SetHeader("L1+HLT efficiency","C"); // option "C" allows to center the header
  //leg_example2->AddEntry((TObject*)0, "SingleEle_Isolated", "");
  leg_example2->AddEntry((TObject*)0, "Endcap, 1.56<|#eta|<2.40", "");
  leg_example2->AddEntry((TObject*)0, "Z'(ee) M=6 TeV, PU=200", "");
  // leg_example2->AddEntry((TObject*)0, "PU=200", "");
  leg_example2->SetFillColor(0);
  leg_example2->SetTextFont(42);
  leg_example2->SetBorderSize(0);
  leg_example2->AddEntry(pEff_sig_EE80, "HLT_Ele32_WP80", "lp");
  leg_example2->AddEntry(pEff_sig_EE70, "HLT_Ele26_WP75", "lp");
  leg_example2->Draw("same");

  TPaveText *t1 = new TPaveText(-6,1.015,250,1.05);
  t1->AddText("#scale[1.2]{CMS} #bf{#it{Simulation Preliminary}}"); ((TText*)t1->GetListOfLines()->Last())->SetTextColor(kBlack);
  t1->SetTextAlign(11);
  t1->SetFillColor(0);
  t1->SetLineColor(0);
  t1->SetBorderSize(1);
  t1->Draw("same");


  TCanvas* my_canvas2 = new TCanvas("canvas2","canvas2",800,700);
  my_canvas2->cd();

  if(TEfficiency::CheckConsistency(*num_sig_EB80,*den_sig_EB80)) {
    pEff_sig_EB80 = new TEfficiency(*num_sig_EB80,*den_sig_EB80);
    pEff_sig_EB80->SetLineColor(kGreen+2);
    pEff_sig_EB80->SetMarkerColor(kGreen+2);
    pEff_sig_EB80->SetMarkerStyle(21);
    pEff_sig_EB80->SetLineWidth(2);
    pEff_sig_EB80->SetTitle(";gen p_{T} [GeV];per electron efficiency");
    pEff_sig_EB80->Draw();
    gPad->Update();
    auto graph2 = pEff_sig_EB80->GetPaintedGraph();
    graph2->SetMinimum(0.0);
    //    graph->GetXaxis()->SetLimits(20.,110.);   
    graph2->SetMaximum(1.01);
    gPad->Update();
  }
  //////


  if(TEfficiency::CheckConsistency(*num_sig_EB70,*den_sig_EB70)) {
    pEff_sig_EB70 = new TEfficiency(*num_sig_EB70,*den_sig_EB70);
    pEff_sig_EB70->SetLineColor(kGreen);
    pEff_sig_EB70->SetMarkerColor(kGreen);
    pEff_sig_EB70->SetMarkerStyle(21);
    pEff_sig_EB70->SetLineWidth(2);
    pEff_sig_EB70->SetTitle(";gen p_{T} [GeV];per electron efficiency");
    pEff_sig_EB70->Draw("same");
    gPad->Update();
  }


  TLegend *leg_example22 = new TLegend(0.5,0.25,0.9,0.5);
  leg_example22->SetHeader("L1+HLT efficiency","C"); // option "C" allows to center the header
  //leg_example22->AddEntry((TObject*)0, "SingleEle_Isolated", "");
  leg_example22->AddEntry((TObject*)0, "Barrel, |#eta|<1.44", "");
  leg_example22->AddEntry((TObject*)0, "Z'(ee) M=6 TeV, PU=200", "");
  // leg_example2->AddEntry((TObject*)0, "PU=200", "");
  leg_example22->SetFillColor(0);
  leg_example22->SetTextFont(42);
  leg_example22->SetBorderSize(0);
  leg_example22->AddEntry(pEff_sig_EB80, "HLT_Ele32_WP80", "lp");
  leg_example22->AddEntry(pEff_sig_EB70, "HLT_Ele26_WP75", "lp");
  leg_example22->Draw("same");
  
  TPaveText *t = new TPaveText(-6,1.015,250,1.05); 
  t->AddText("#scale[1.2]{CMS} #bf{#it{Simulation Preliminary}}"); ((TText*)t->GetListOfLines()->Last())->SetTextColor(kBlack); 
  t->SetTextAlign(11);
  t->SetFillColor(0); 
  t->SetLineColor(0); 
  t->SetBorderSize(1); 
  t->Draw("same");

  my_canvas1->SetGrid(1,1);
  my_canvas2->SetGrid(1,1);

  my_canvas1->Update();
  my_canvas1->SaveAs("SingleEle_fullEff_vs_genpt_EE_wp_70_80_zprime.png");

  my_canvas2->Update();
  my_canvas2->SaveAs("SingleEle_fullEff_vs_genpt_EB_wp_70_80_zprime.png");
}

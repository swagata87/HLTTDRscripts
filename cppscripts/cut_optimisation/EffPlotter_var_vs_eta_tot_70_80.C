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

void EffPlotter_var_vs_eta_tot_70_80() {

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

  TH1D* den_sig80 = (TH1D*)file_sig80->Get("den_ele_geneta");
  TH1D* num_sig80 = (TH1D*)file_sig80->Get("num_ele_geneta_all");

  TH1D* den_sig70 = (TH1D*)file_sig70->Get("den_ele_geneta");
  TH1D* num_sig70 = (TH1D*)file_sig70->Get("num_ele_geneta_all");
  
  TEfficiency* pEff_sig80 = 0;
  TEfficiency* pEff_sig70 = 0;

  TCanvas* my_canvas1 = new TCanvas("canvas","canvas",800,700);
  my_canvas1->cd();

  if(TEfficiency::CheckConsistency(*num_sig80,*den_sig80)) {
    pEff_sig80 = new TEfficiency(*num_sig80,*den_sig80);
    pEff_sig80->SetLineColor(kGreen+2);
    pEff_sig80->SetMarkerColor(kGreen+2);
    pEff_sig80->SetMarkerStyle(21);
    pEff_sig80->SetLineWidth(3);
    pEff_sig80->SetTitle(";gen #eta ;per electron efficiency");
    pEff_sig80->Draw();
    gPad->Update();
    auto graph = pEff_sig80->GetPaintedGraph(); 
    graph->SetMinimum(0.0);
    //    graph->SetMaximum(1.05); 
    graph->SetMaximum(1.01); 
    //    graph->GetXaxis()->SetLimits(20.,110.);
  }


  if(TEfficiency::CheckConsistency(*num_sig70,*den_sig70)) {
    pEff_sig70 = new TEfficiency(*num_sig70,*den_sig70);
    pEff_sig70->SetLineColor(kGreen);
    pEff_sig70->SetMarkerColor(kGreen);
    pEff_sig70->SetMarkerStyle(20);
    pEff_sig70->SetLineWidth(2);
    //pEff_sig70->SetTitle("Signal efficiency vs #eta;gen #eta ;efficiency");
    pEff_sig70->Draw("same");
    gPad->Update();
    //    auto graph = pEff_sig70->GetPaintedGraph(); 
    //graph->SetMinimum(0.0);
    //graph->SetMaximum(1.0); 
  }


  TLegend *leg_example2 = new TLegend(0.3,0.25,0.7,0.5);
  leg_example2->SetHeader("L1+HLT efficiency","C"); // option "C" allows to center the header
  //leg_example2->AddEntry((TObject*)0, "SingleEle_Isolated", "");
  leg_example2->AddEntry((TObject*)0, "40<pT<500 GeV", "");
  leg_example2->AddEntry((TObject*)0, "Z'(ee) M=6 TeV, PU=200", "");
  // leg_example2->AddEntry((TObject*)0, "PU=200", "");
  leg_example2->SetFillColor(0);
  leg_example2->SetTextFont(42);
  leg_example2->SetBorderSize(0);
  leg_example2->AddEntry(pEff_sig80, "HLT_Ele32_WP80", "lp");
  leg_example2->AddEntry(pEff_sig70, "HLT_Ele32_WP75", "lp");
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
  my_canvas1->SaveAs("SingleEle_fullEff_vs_geneta_wp_70_80_Zprime.png");
}

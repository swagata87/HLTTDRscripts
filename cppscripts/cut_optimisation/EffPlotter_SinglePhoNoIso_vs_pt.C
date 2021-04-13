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

void EffPlotter_SinglePhoNoIso_vs_pt() {

  TFile *file_sig = new TFile("histSinglePhoNoIso_Sig_Zprime_L1seeded.root");
  
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
  TH1D* num_sig = (TH1D*)file_sig->Get("num_ele_genpt_all_EB");
    
  TEfficiency* pEff_sig = 0;

  TCanvas* my_canvas1 = new TCanvas("canvas","canvas",800,700);
  my_canvas1->cd();

  if(TEfficiency::CheckConsistency(*num_sig,*den_sig)) {
    pEff_sig = new TEfficiency(*num_sig,*den_sig);
    pEff_sig->SetLineColor(kAzure+3);
    pEff_sig->SetMarkerColor(kAzure+3);
    pEff_sig->SetMarkerStyle(21);
    pEff_sig->SetLineWidth(3);
    pEff_sig->SetTitle(";gen p_{T} [GeV];per electron efficiency");
    pEff_sig->Draw();
    gPad->Update();
    auto graph = pEff_sig->GetPaintedGraph(); 
    graph->SetMinimum(0.0);
    //    graph->SetMaximum(1.05); 
    graph->SetMaximum(1.01); 
    //    graph->GetXaxis()->SetLimits(0.,1000);
  }

  //  TLegend *leg_example22 = new TLegend(0.4,0.25,0.9,0.5);
  //  TLegend *leg_example2 = new TLegend(0.68,0.8,0.94,0.94);
  TLegend *leg_example2 = new TLegend(0.5,0.30,0.9,0.6);
  leg_example2->SetHeader("L1+HLT efficiency","C"); // option "C" allows to center the header
  //  leg_example2->SetHeader("L1 Photon 51 GeV + hwQual","C"); // option "C" allows to center the header
  leg_example2->AddEntry((TObject*)0, "HLT_Photon187", "");
  //leg_example2->AddEntry((TObject*)0, "Endcap, 1.56<|#eta|<2.70", "");
  leg_example2->AddEntry((TObject*)0, "Barrel, |#eta|<1.44", "");
  leg_example2->AddEntry((TObject*)0, "Z'(ee) M=6 TeV, PU=200", "");
  // leg_example2->AddEntry((TObject*)0, "PU=200", "");
  leg_example2->SetFillColor(0);
  leg_example2->SetTextFont(42);
  leg_example2->SetBorderSize(0);
  //leg_example2->AddEntry(pEff_sig, "efficiency", "lp");
  //  leg_example2->AddEntry(pEff_sig_leg2, "leg2", "lp");
  leg_example2->Draw("same");

  TPaveText *t1 = new TPaveText(-6,1.015,250,1.05);
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
  TH1D* num_sig_EE = (TH1D*)file_sig->Get("num_ele_genpt_all_EE");
  //  TH1D* num_sig_leg2_EE = (TH1D*)file_sig->Get("num_ele_genpt_EE_leg2_all");
  
  TEfficiency* pEff_sig_EE = 0;
  //  TEfficiency* pEff_sig_leg2_EE = 0;

  TCanvas* my_canvas2 = new TCanvas("canvas2","canvas2",800,700);
  my_canvas2->cd();

  if(TEfficiency::CheckConsistency(*num_sig_EE,*den_sig_EE)) {
    pEff_sig_EE = new TEfficiency(*num_sig_EE,*den_sig_EE);
    pEff_sig_EE->SetLineColor(kAzure+3);
    pEff_sig_EE->SetMarkerColor(kAzure+3);
    pEff_sig_EE->SetMarkerStyle(21);
    pEff_sig_EE->SetLineWidth(3);
    pEff_sig_EE->SetTitle(";gen p_{T} [GeV];per electron efficiency");
    pEff_sig_EE->Draw();
    gPad->Update();
    auto graph_EE = pEff_sig_EE->GetPaintedGraph(); 
    graph_EE->SetMinimum(0.0);
    graph_EE->SetMaximum(1.01); 
  }
  TLegend *leg_example22 = new TLegend(0.5,0.30,0.9,0.6);
  // TLegend *leg_example22 = new TLegend(0.4,0.25,0.9,0.5);
  leg_example22->SetHeader("L1+HLT efficiency","C"); // option "C" allows to center the header
  leg_example22->AddEntry((TObject*)0, "HLT_Photon187", "");
  leg_example22->AddEntry((TObject*)0, "Endcap, 1.56<|#eta|<2.40", "");
  leg_example22->AddEntry((TObject*)0, "Z'(ee) M=6 TeV, PU=200", "");
  leg_example22->SetFillColor(0);
  leg_example22->SetTextFont(42);
  leg_example22->SetBorderSize(0);
  //  leg_example22->AddEntry(pEff_sig_EE, "efficiency", "lp");
  //  leg_example22->AddEntry(pEff_sig_leg2_EE, "leg2", "lp");
  leg_example22->Draw("same");


  TPaveText *t = new TPaveText(-6,1.015,250,1.05);
  t->AddText("#scale[1.2]{CMS} #bf{#it{Simulation Preliminary}}"); ((TText*)t->GetListOfLines()->Last())->SetTextColor(kBlack);
  t->SetTextAlign(11);
  t->SetFillColor(0);
  t->SetLineColor(0);
  t->SetBorderSize(1);
  t->Draw("same");
  
  my_canvas2->SetGrid(1,1);

  my_canvas2->Update();

  //  my_canvas1->SaveAs("SinglePhoNoIso_fullEff_vs_genpt_EB.png");
  // my_canvas2->SaveAs("SinglePhoNoIso_fullEff_vs_genpt_EE.png");

  my_canvas1->SaveAs("SinglePho_vs_genpt_EB.png");
  my_canvas2->SaveAs("SinglePho_vs_genpt_EE.png");

}

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

void compare_sigww_EE() {

  std::cout << "Get the root files " << std::endl;
  TFile *file_sig = new TFile("redef_sig_Dec.root");
  TFile *file_bkg = new TFile("redef_bkg_Dec.root");

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

  TH1D* h1_sig = (TH1D*)file_sig->Get("sigmaww_sig_EE");
  h1_sig->Rebin(2);
  h1_sig->SetLineColor(kGreen+2);
  h1_sig->SetFillColor(kGreen+2);
  h1_sig->SetFillStyle(3144);
  h1_sig->SetLineWidth(3);
  h1_sig->Scale(1/h1_sig->Integral());
  h1_sig->GetXaxis()->SetTitle("#sigma_{ww}");
  h1_sig->GetYaxis()->SetTitle("Arbitrary unit");
  h1_sig->SetTitle("");
  h1_sig->SetStats(0);
  h1_sig->GetXaxis()->SetRangeUser(1, 19);

  TH1D* h1_bkg = (TH1D*)file_bkg->Get("sigmaww_bkg_EE");
  h1_bkg->Rebin(2);
  h1_bkg->SetLineColor(kRed+1);
  h1_bkg->SetFillColor(kRed+1);
  h1_bkg->SetFillStyle(3350);
  h1_bkg->SetLineWidth(3);
  h1_bkg->Scale(1/h1_bkg->Integral());
  //  h1_bkg->GetXaxis()->SetTitle("#sigmaww");
  //h1_bkg->SetTitle("");
  //h1_bkg->SetStats(0);
  //h1_bkg->GetXaxis()->SetRangeUser(0.4, 2.0);


  TCanvas* my_canvas1 = new TCanvas("canvas1","canvas1",800,600);
  my_canvas1->cd();
  gPad->SetLogy();
  h1_sig->Draw("hist");
  h1_bkg->Draw("hist same");

  TLegend *leg_example = new TLegend(0.56,0.75,0.94,0.94);
  leg_example->SetHeader("Electron #sigma_{ww}, HGCAL","C"); // option "C" allows to center the header
  leg_example->AddEntry((TObject*)0, "pT>30 GeV, PU=200", "");
  leg_example->SetFillColor(0);
  leg_example->SetTextFont(42);
  leg_example->SetBorderSize(0);
  leg_example->AddEntry(h1_sig, "Signal: Drell-Yan, W+Jets", "f");
  leg_example->AddEntry(h1_bkg, "Background: QCD", "f");
  //leg_example->AddEntry(h1_dphiin_EE_ticl_v3_trkv6, "TICL V3 HLT tracks", "lp");
  leg_example->Draw("same");
  my_canvas1->SetGrid();

  TPaveText *t = new TPaveText(0.8,0.26,8,0.4); 
  t->AddText("#scale[1.2]{CMS} #bf{#it{Simulation Preliminary}}"); ((TText*)t->GetListOfLines()->Last())->SetTextColor(kBlack); 
  t->SetTextAlign(11); 
  t->SetFillColor(0); 
  t->SetLineColor(0); 
  t->SetBorderSize(1); 
  t->Draw("same");

  my_canvas1->SaveAs("sigmaww_EE_sig_bkg.png");

}

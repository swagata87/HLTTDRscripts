#include <algorithm>
#include "TLatex.h"
#include "TROOT.h"
#include <iomanip>
#include "TObjArray.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>
#include <vector>
#include "TGraph.h"
#include "TFile.h"
#include "TLegend.h"
#include "TObject.h"
#include "TAttLine.h"
#include "TAttBBox2D.h"
#include "TPoint.h"
#include "TPad.h"
#include "TColor.h"
#include "TStyle.h"
#include "TPostScript.h"
#include "TAttBBox2D.h"
#include "TArrow.h"
#include <stdio.h>
#include "TAttText.h"
#include "TPaveText.h"

int pms2_ROC_EE() {

  TFile *file_bkg = new TFile("redef_bkg_Dec.root");
  TFile *file_sig = new TFile("redef_sig_Dec.root");
  /////

  TH1D* h1_sig_old = (TH1D*)file_sig->Get("pms2_old_sig_EE");
  TH1D* h1_bkg_old = (TH1D*)file_bkg->Get("pms2_old_bkg_EE");

  TH1D* h1_sig_new = (TH1D*)file_sig->Get("pms2_new_sig_EE");
  TH1D* h1_bkg_new = (TH1D*)file_bkg->Get("pms2_new_bkg_EE");

  int nbins = h1_sig_old->GetNbinsX();

  std::cout << "number of bins" << std::endl;
  std::cout << h1_sig_old->GetNbinsX() << std::endl ;
  std::cout << h1_bkg_old->GetNbinsX() << std::endl ;
  std::cout << h1_sig_new->GetNbinsX() << std::endl ;
  std::cout << h1_bkg_new->GetNbinsX() << std::endl ;

  /*
  std::cout << "max value" << std::endl;
  std::cout << h1_sig_old->GetMaximum() << std::endl ;
  std::cout << h1_bkg_old->GetMaximum() << std::endl ;
  std::cout << h1_sig_new->GetMaximum() << std::endl ;
  std::cout << h1_bkg_new->GetMaximum() << std::endl ;
  */
  // get the total integrals for each histogram
  float sig_integral_old = h1_sig_old->Integral(1,nbins);
  float bkg_integral_old = h1_bkg_old->Integral(1,nbins);

  float sig_integral_new = h1_sig_new->Integral(1,nbins);
  float bkg_integral_new = h1_bkg_new->Integral(1,nbins);

  // create containers sig = x points, bkg = y points

  std::vector<float> sigPoints_old(nbins);
  std::vector<float> bkgPoints_old(nbins);

  std::vector<float> sigPoints_new(nbins);
  std::vector<float> bkgPoints_new(nbins);

  // in the loop, fill the containers with the x and y points
  // each x point is this loop's slice integral over total integral (signal)
  // each y point is this loop's slice integral over total integral (background)
  //  std::cout << "print only old_pu0 \n " ;
  //  std::cout << "files " << file_sig_old_pu0 << " " << file_bkg_old_pu0 << " " << file_sig_new << " " << file_bkg_new << std::endl;
  //myfile << "cut / sig eff / bkg rej \n";

  for ( int i = 1; i < nbins+1; ++i ) {
    // notice the slice integral is dependent on i!
    // on each iteration we take a larger and larger slice of the histograms
    // eventually the slice will be the total integral (from bin 1 to bin nbins)
    // that point is (1,1) on the ROC curve.
    //    float sig_slice_integral_old_pu0 = h1_sig_old_pu0->Integral(nbins-i,nbins);
    // float bkg_slice_integral_old_pu0 = h1_bkg_old_pu0->Integral(nbins-i,nbins);

    //float sig_slice_integral_new = h1_sig_new->Integral(nbins-i,nbins);
    //float bkg_slice_integral_new = h1_bkg_new->Integral(nbins-i,nbins);

    float sig_slice_integral_old = h1_sig_old->Integral(1,i);
    float bkg_slice_integral_old = h1_bkg_old->Integral(1,i);

    float sig_slice_integral_new = h1_sig_new->Integral(1,i);
    float bkg_slice_integral_new = h1_bkg_new->Integral(1,i);

    //    if ( (h1_sig_old_pu0->GetBinCenter(i) > 0.0110) && (h1_sig_old_pu0->GetBinCenter(i) < 0.0112) ) {
    // std::cout << "HoE cut = " << h1_sig_old_pu0->GetBinCenter(i) << " = " << h1_bkg_old_pu0->GetBinCenter(i) << std::endl;
    // std::cout << "sig eff old_pu0 " << sig_slice_integral_old_pu0/sig_old_pu0_integral << " bkg rej old_pu0 " << 1-(bkg_slice_integral_old_pu0/bkg_old_pu0_integral) << std::endl;

    //}

    sigPoints_old.push_back(sig_slice_integral_old/sig_integral_old);
    bkgPoints_old.push_back(1-(bkg_slice_integral_old/bkg_integral_old));

    sigPoints_new.push_back(sig_slice_integral_new/sig_integral_new);
    bkgPoints_new.push_back(1-(bkg_slice_integral_new/bkg_integral_new));


    //myfile << h1_sig_pu200_default->GetBinCenter(i) << " / " << sig_slice_integral_pu200_default/sig_pu200_integral_default << " / " << 1-(bkg_slice_integral_pu200_default/bkg_pu200_integral_default) << "\n"  ;

  }

  TCanvas *c1 = new TCanvas("c1","eff", 200,10,600,400);

  // create a TGraph from the containers
  // this graph will have N (=nbins) number of points forming the curve.

  TGraph *g_old = new TGraph(sigPoints_old.size(),&sigPoints_old[0],&bkgPoints_old[0]);
  g_old->SetMarkerColor(kGreen+1);
  g_old->SetLineColor(kGreen+1);
  g_old->SetLineWidth(3);
  g_old->SetMarkerStyle(20);
  g_old->Draw("AP");
  g_old->SetTitle("Pixel matching S^{2} (Endcap)");
  g_old->GetYaxis()->SetTitle("Background Rejection");
  g_old->GetXaxis()->SetTitle("Signal Efficiency");
  g_old->GetHistogram()->SetMaximum(0.85);   // along                                                                      
  g_old->GetHistogram()->SetMinimum(0.0);  //   Y                                                                                
  g_old->GetXaxis()->SetLimits(0.88,1.002);


  TGraph *g_new = new TGraph(sigPoints_new.size(),&sigPoints_new[0],&bkgPoints_new[0]);
  g_new->SetMarkerColor(kBlue+1);
  g_new->SetLineColor(kBlue+1);
  g_new->SetLineWidth(3);
  g_new->SetMarkerStyle(21);
  g_new->Draw("P: same");


  /*
  int sig200_size =  sigPoints_pu200.size();
  float q1,q2,p1,p2;
  q1 = 1-bkgPoints_pu200[0];
  q2 = sigPoints_pu200[0];
  float area = 0.0;
  for(int i=1;i < sig200_size ;++i){
    p1 = 1-bkgPoints_pu200[i];
    p2 = sigPoints_pu200[i];
    area += sqrt(pow( ((1-q1)+(1-p1))/2 * (q2-p2),2));
    q1=p1;
    q2=p2;   
  }
  std::cout << "auc pu200 roc " << area << std::endl;
  */
  //  char myarea[40];

  // std::sprintf(myarea,"AUC = %g",area);
   
  //using normal text in NDC
  //  TText *t2 = new TText(0.13,0.4,myarea);
  // t2->SetNDC();
  //  t2->SetTextColor(kBlue);
  //  t2->SetTextSize(0.04);
  // t2->Draw("same");

  TLegend *leg_example = new TLegend(0.17,0.17,0.65,0.55);
  //  TLegend *leg_example = new TLegend(0.14,0.14,0.55,0.4);
  leg_example->SetHeader("p_{T}>30 GeV, PU=200","C"); // option "C" allows to center the header              
  //leg_example->AddEntry((TObject*)0, "H/E<0.10+5.0/pT", "");
  leg_example->AddEntry((TObject*)0, "Signal: Drell-Yan, W+Jets", "");
  leg_example->AddEntry((TObject*)0, "Background: QCD", "");
  //leg_example->AddEntry((TObject*)0, "Using 3D TICL clusters", "");
  leg_example->SetFillColor(0);
  leg_example->SetTextFont(42);
  leg_example->SetBorderSize(0);
  leg_example->AddEntry(g_old, "Run-2 definition","pl");
  leg_example->AddEntry(g_new, "Phase-2 definition","pl");
  //  leg_example->AddEntry(g_pu200_layerclus, "2D cluster","pl");
  leg_example->Draw("same");

  TPaveText *t = new TPaveText(0.95,0.75,1.00,0.82);
  t->AddText("#scale[1.2]{CMS} #bf{#it{Simulation Preliminary}}"); ((TText*)t->GetListOfLines()->Last())->SetTextColor(kBlack);
  t->SetTextAlign(11);
  t->SetFillColor(0);
  t->SetLineColor(0);
  t->SetBorderSize(1);
  t->Draw("same");

  c1->SetGrid();
  c1->SaveAs("ROC_pms2_new_vs_old_Endcap.png");

  //myfile.close(); 

  return 0;
}

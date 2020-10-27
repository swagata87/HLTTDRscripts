#!/usr/bin/env python

#RUN it like this
#python pixelmatch_S2.py  inputFile.root    -o=pms2.root  -g=0/1   -weights=path/to/weight_json 

from array import array
import re
import argparse
import sys
import math
from DataFormats.FWLite import Events, Handle
import ROOT
from ROOT import TCanvas, TFile, TProfile, TNtuple, TH1F, TH2F, TGraph
import json
import numpy

#These 'a' parameters fit values are taken from here: 
#https://gitlab.cern.ch/sharper/EgHLTPhase2/-/blob/master/EDProducers/hltEgammaPixelMatchVarsUnseeded_cfi.py
def getValFromFit_a_dphi1(abseta,nclus):
   fx=999
   x=abseta
   if (nclus==1):
      if (x>=0 and x<1.5): fx = 0.00112 + (0.000752 * x) - (0.00122 *x*x) + (0.00109*x*x*x)
      if (x>=1.5 and x<2): fx = 0.00823 - (0.0029 * x )
      if (x>=2 and x<=3): fx = 0.00282
   if (nclus==2):
      if (x>=0 and x<1.5): fx = 0.00222 + (0.000196 * x) - (0.000203 * x*x) + (0.000447*x*x*x)
      if (x>=1.5 and x<2): fx = 0.010838  - (0.00345*x)
      if (x>=2 and x<=3): fx = 0.0043 
   if (nclus>2):
      if (x>=0 and x<1.5): fx = 0.00236 + (0.000691*x) + (0.000199*x*x) + (0.000416*x*x*x) 
      if (x>=1.5 and x<=3): fx = 0.0208 - 0.0125*x + 0.00231*x*x 

   return fx

####
#In EDProducers/hltEgammaPixelMatchVarsUnseeded_cfi.py, there is a overlap of eta bin for a_dphi2, 0-1.6 and 1.5-1.9.
#I believe this is not intended. So I correct it here: 0-1.5 and 1.5-1.9 
def getValFromFit_a_dphi2(abseta,nclus):
   fx=999
   x=abseta
   if (nclus>0):
      if (x>=0 and x<1.5): fx = 0.00013
      if (x>=1.5 and x<1.9): fx = 0.00045 - 0.000199*x
      if (x>=1.9 and x<=3): fx = 7.94e-05

   return fx

####
def getValFromFit_a_drz2(abseta,nclus):
   fx=999
   x=abseta
   if (nclus>0):
      if (x>=0 and x<1.5): fx = 0.00299 + 0.000299*x - 4.13e-06*x*x + 0.00191*x*x*x
      if (x>=1.5 and x<=3): fx = 0.248 - 0.329*x + 0.148*x*x - 0.0222*x*x*x

   return fx

## class for event weight from Sam ##
class EvtWeights:
    def __init__(self,input_filename,lumi=0.075):
        if input_filename:
            with open(input_filename,'r') as f:
                self.data = json.load(f)
        else:
            self.data = {}
        self.warned = []
        self.lumi = lumi #luminosity to weight to in pb                                         

    def weight_from_name(self,dataset_name):
        if dataset_name in self.data:
            val = self.data[dataset_name]
            return val['xsec']/val['nrtot']*self.lumi
        else:
            if dataset_name not in self.warned:
                self.warned.append(dataset_name)
                print "{} not in weights file, returning weight 1".format(dataset_name)
            return 1.

    def weight_from_evt(self,event):
        filename = event.getTFile().GetName().split("/")[-1]
        dataset_name = re.search(r'(.+)(_\d+_EDM.root)',filename).groups()[0]
        #print dataset_name
        return self.weight_from_name(dataset_name)

    def dataset_from_evt(self,event):
        filename = event.getTFile().GetName().split("/")[-1]
        dataset_name = re.search(r'(.+)(_\d+_EDM.root)',filename).groups()[0]
        return dataset_name



if __name__ == "__main__":

    oldargv = sys.argv[:]
    sys.argv = [ '-b-' ]
    sys.argv = oldargv
    ROOT.gSystem.Load("libFWCoreFWLite.so");
    ROOT.gSystem.Load("libDataFormatsFWLite.so");
    ROOT.FWLiteEnabler.enable()

    parser = argparse.ArgumentParser(description='example e/gamma HLT analyser')
    parser.add_argument('in_filename',nargs="+",help='input filename')
    parser.add_argument("-o", "--output", help="Directs the output to a name of your choice")
    parser.add_argument("-g", "--genmatched", help="1 for signal, 0 for bkg")
    parser.add_argument("-weights", "--weights", help="json filename having weights / 0 if weights not needed")

    args = parser.parse_args()
   
    print 'flag = ', args.genmatched
    weights = EvtWeights(args.weights)
    print weights

    trigBitsHandle, trigBitsLabel = Handle("edm::TriggerResults"), "TriggerResults"
    trigobjs_handle, trigobjs_label = Handle("std::vector<reco::EgTrigSumObj>"), "hltEgammaHLTExtra"
    if (args.genmatched=='1'): gen_handle, gen_label = Handle("std::vector<reco::GenParticle>"), "genParticles"

    to_remove=[]
    for file in args.in_filename:
        file_temp = ROOT.TFile(file)
        if ( file_temp.IsZombie() ) :
            print "******** remove file **********"
            to_remove.append(file)

    new_list = [x for x in args.in_filename if x not in to_remove]
    #print "AFTER REMOVING ALL ZOMBIES, this is the new list of files "
    #print new_list

    events = Events(new_list)    

    #saved pms2 and my pms2 should match
    saved_PMS2 = ROOT.TH1D("saved_PMS2",";PMS2;#entries",100,0,50)
    my_PMS2 = ROOT.TH1D("my_PMS2",";PMS2;#entries",100,0,50)

    #if they do not match perfectly, then check if the mismatch is concentrated in eta
    saved_PMS2_vs_abseta = ROOT.TProfile("saved_PMS2_vs_abseta",";|#eta|;saved PMS2",30,0,3.0,0,4000)
    my_PMS2_vs_abseta = ROOT.TProfile("my_PMS2_vs_abseta",";|#eta|;my PMS2",30,0,3.0,0,4000)

    #slide 14 of https://indico.cern.ch/event/482680/contributions/2296717/attachments/1332989/2004104/hltUpdates_070916_V2.pdf
    dphi1_vs_eta = ROOT.TH2D("dphi1_vs_eta",";#eta;d#phi1",640,-3.2,3.2,100,-0.05,0.05)
    dphi2_vs_eta = ROOT.TH2D("dphi2_vs_eta",";#eta;d#phi2",640,-3.2,3.2,60,-0.003,0.003)
    drz2_vs_eta = ROOT.TH2D("drz2_vs_eta",";#eta;drz2",640,-3.2,3.2,180,-0.09,0.09)

    #    
    prof_dphi1_vs_eta_nclus1 = ROOT.TProfile("prof_dphi1_vs_eta_nclus1",";#eta;mean d#phi1",60,-3.0,3.0,-0.09,0.09,"s")
    prof_dphi1_vs_eta_nclus2 = ROOT.TProfile("prof_dphi1_vs_eta_nclus2",";#eta;mean d#phi1",60,-3.0,3.0,-0.09,0.09,"s")
    prof_dphi1_vs_eta_nclus3 = ROOT.TProfile("prof_dphi1_vs_eta_nclus3",";#eta;mean d#phi1",60,-3.0,3.0,-0.09,0.09,"s")

    #slide 5-6 of https://indico.cern.ch/event/613828/contributions/2561829/attachments/1447465/2230535/pixel_matching_hlt_21apr2017.pdf
    dphi1_vs_nclus = ROOT.TH2D("dphi1_vs_nclus",";nclus;d#phi1",10,0,10,100,-0.05,0.05)
    dphi2_vs_nclus = ROOT.TH2D("dphi2_vs_nclus",";nclus;d#phi2",10,0,10,120,-0.003,0.003)
    drz2_vs_nclus = ROOT.TH2D("drz2_vs_nclus",";nclus;drz2",10,0,10,180,-0.09,0.09)

    #slide 6-8 of https://indico.cern.ch/event/379088/contributions/1803928/attachments/756169/1037305/aidan.pdf
    dphi1_hist = ROOT.TH1D("dphi1_hist",";d#phi1;nEntries",100,-0.05,0.05)
    dphi2_hist = ROOT.TH1D("dphi2_hist",";d#phi2;nEntries",60,-0.003,0.003)
    drz2_hist = ROOT.TH1D("drz2_hist",";drz2;nEntries",180,-0.09,0.09)

    for event_nr,event in enumerate(events):

        event.getByLabel(trigobjs_label,trigobjs_handle)
        if (args.genmatched=='1'):  event.getByLabel(gen_label,gen_handle)

        trigobjs = trigobjs_handle.product()
        if (args.genmatched=='1'): genobjs=gen_handle.product()

        if (args.genmatched=='1'):
            weight = 1

        if (args.genmatched=='0'):
            weight = weights.weight_from_evt(event.object())
       
        for eg in trigobjs:
          min_dr=999.9
          if (args.genmatched=='1'):
              for gen in genobjs:
                  if math.fabs(gen.pdgId())==11:
                      dr2 = ROOT.reco.deltaR2(eg.eta(),eg.phi(),gen.eta(),gen.phi())
                      dr=math.sqrt(dr2)
                      if dr<min_dr: min_dr=dr
          elif (args.genmatched=='0'): min_dr=-1

          if (min_dr<=0.1) : 
             if eg.et()>=30.0 and abs(eg.superCluster().eta())<3.0:

                saved_PMS2.Fill(eg.var("hltEgammaPixelMatchVarsUnseeded_s2"))
                saved_PMS2_vs_abseta.Fill(abs(eg.superCluster().eta()),eg.var("hltEgammaPixelMatchVarsUnseeded_s2"))

                abseta=abs(eg.superCluster().eta())
                nclus=eg.superCluster().clustersSize()

                dPhi1Const =  getValFromFit_a_dphi1(abseta,nclus)
                dPhi2Const =  getValFromFit_a_dphi2(abseta,nclus)
                dRZ2Const =  getValFromFit_a_drz2(abseta,nclus)

                if (dPhi1Const==999 or dPhi2Const==999 or dRZ2Const==999):
                   print '\nWARNING !!!!  a param value 999 ... something is wrong ... \n'
                
                min_my_s2=999999
                best_dphi1=999999
                best_dphi2=999999
                best_drz2=999999

                for seed in eg.seeds():
                   #taking help from https://github.com/cms-sw/cmssw/blob/master/RecoEgamma/EgammaHLTProducers/plugins/EgammaHLTPixelMatchVarProducer.cc
                   charge = seed.getCharge()

                   #dPhiPos(hitNr) etc are defined here: https://github.com/cms-sw/cmssw/blob/master/DataFormats/EgammaReco/interface/ElectronSeed.h
                   #dPhi1 = difference in phi of first hit from SC traj, similar for dPhi2, dRZ2
                   if (charge < 0): 
                      dPhi1=seed.dPhiNeg(0) 
                      dPhi2=seed.dPhiNeg(1) 
                      dRZ2=seed.dRZNeg(1)

                   elif (charge > 0): 
                      dPhi1=seed.dPhiPos(0)
                      dPhi2=seed.dPhiPos(1) 
                      dRZ2=seed.dRZPos(1) 
           
                   dPhi1_term = dPhi1/dPhi1Const
                   dPhi2_term = dPhi2/dPhi2Const
                   dRZ2_term = dRZ2/dRZ2Const
                
                   #there are many pixel seeds for each supercluster. Take the one which has lowest PM_S2. 
                   my_s2 = dPhi1_term*dPhi1_term + dPhi2_term*dPhi2_term + dRZ2_term*dRZ2_term
                   if (my_s2 < min_my_s2): 
                      min_my_s2=my_s2
                      best_dphi1=dPhi1
                      best_dphi2=dPhi2
                      best_drz2=dRZ2

                #histograms are filled w/o weight because for this study I expect to run only on signal (ele gun or DY).
                #To get 'a' params, we do not need to run on QCD, just run on real electrons
                #But later, if this code is to be used for ROC etc, then weights will be needed
                my_PMS2.Fill(min_my_s2)
                my_PMS2_vs_abseta.Fill(abs(eg.superCluster().eta()),min_my_s2)

                #
                dphi1_vs_eta.Fill(eg.superCluster().eta(),best_dphi1)
                dphi2_vs_eta.Fill(eg.superCluster().eta(),best_dphi2)
                drz2_vs_eta.Fill(eg.superCluster().eta(),best_drz2) 
                #
                # In the past it was found that some delta params width have stron ncluster dependence
                dphi1_vs_nclus.Fill(nclus,best_dphi1)
                dphi2_vs_nclus.Fill(nclus,best_dphi2)
                drz2_vs_nclus.Fill(nclus,best_drz2) 
                #
                dphi1_hist.Fill(best_dphi1)
                dphi2_hist.Fill(best_dphi2)
                drz2_hist.Fill(best_drz2) 

                if (nclus==1): prof_dphi1_vs_eta_nclus1.Fill(eg.superCluster().eta(),best_dphi1)
                if (nclus==2): prof_dphi1_vs_eta_nclus2.Fill(eg.superCluster().eta(),best_dphi1)
                if (nclus>2): prof_dphi1_vs_eta_nclus3.Fill(eg.superCluster().eta(),best_dphi1)

    output_file = TFile( args.output, 'recreate' )

    saved_PMS2.Write()
    my_PMS2.Write()

    dphi1_vs_eta.Write()
    dphi2_vs_eta.Write()
    drz2_vs_eta.Write()

    dphi1_hist.Write()
    dphi2_hist.Write()
    drz2_hist.Write()

    dphi1_vs_nclus.Write()
    dphi2_vs_nclus.Write()
    drz2_vs_nclus.Write()

    prof_dphi1_vs_eta_nclus1.Write()
    prof_dphi1_vs_eta_nclus2.Write()
    prof_dphi1_vs_eta_nclus3.Write()

    my_PMS2_vs_abseta.Write()
    saved_PMS2_vs_abseta.Write()

    output_file.Close()

    #whats next?
    #plot spread of dphi1/2, dz2 as a function of eta (for nclus = 1, 2, ...) like this:
    #slide 15 of https://indico.cern.ch/event/482680/contributions/2296717/attachments/1332989/2004104/hltUpdates_070916_V2.pdf
    #check the trends and fit accordingly

    #Further references:
    #https://indico.cern.ch/event/626270/contributions/2529268/attachments/1434975/2206722/egPhase1XPOG.pdf
    #https://indico.cern.ch/event/604937/contributions/2596214/attachments/1461173/2257099/egHLTPMTSG170517.pdf

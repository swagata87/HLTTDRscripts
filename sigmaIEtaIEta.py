#!/usr/bin/env python

from array import array
import re
import argparse
import sys
import math
from DataFormats.FWLite import Events, Handle
import ROOT
from ROOT import TCanvas, TFile, TProfile, TNtuple, TH1F, TH2F
import json
import FWCore.ParameterSet.Config as cms
from CondCore.CondDB.CondDB_cfi import *
from PhysicsTools.PythonAnalysis import *

def get_hitEnergy(hit):
    En=hit.energy()
    #if ( hit.checkFlag(hit.kTowerRecovered) or hit.checkFlag(hit.kWeird) or hit.checkFlag(hit.kDiWeird) ):
        #En=0.0
    return En

def get_E5x5(ele,ebhits):

    seed_id = eg.superCluster().seed().seed()
    if seed_id.subdetId()==1: #not ecal barrel, return
        seed_ebid = ROOT.EBDetId(seed_id)
        seed_crys_ieta = seed_ebid.ieta()
        seed_crys_iphi = seed_ebid.iphi()
        E5x5=0
        my_ieta=[0,1,-1,2,-2]
        my_iphi=[0,1,-1,2,-2]
        for hit in ebhits:
            ecalhit=ROOT.EBDetId(hit.id())
            for i in my_ieta:
                for j in my_iphi:
                    if ecalhit.ieta()==seed_crys_ieta+i and ecalhit.iphi()==seed_crys_iphi+j:
                        #print 'ieta/ iphi ', i, j 
                        E5x5=E5x5+get_hitEnergy(hit)
        
    return E5x5

#
def get_wi(Ei,E5x5):
    if Ei>0 and E5x5>0:
        wi=max(0,4.7+math.log(Ei/E5x5))
    return wi

#
def getNrCrysDiffInEta(crysIEta,orginIEta):
    nrCrysDiff = crysIEta-orginIEta
    if (crysIEta*orginIEta<0):
        if (crysIEta>0): nrCrysDiff=nrCrysDiff-1
        else: nrCrysDiff=nrCrysDiff+1
    return nrCrysDiff
#
def EtaBar5x5(eg,E5x5,ebhits):
    seed_id = eg.superCluster().seed().seed()
    if seed_id.subdetId()==1: #not ecal barrel, return
        seed_ebid = ROOT.EBDetId(seed_id)
        seed_crys_ieta = seed_ebid.ieta()
        seed_crys_iphi = seed_ebid.iphi()
        meanDEta=0
        for hit in ebhits:
            ecalhit=ROOT.EBDetId(hit.id())
            my_ieta=[0,1,-1,2,-2]
            my_iphi=[0,1,-1,2,-2]
            for i in my_ieta:
                for j in my_iphi:
                    if ecalhit.ieta()==seed_crys_ieta+i and ecalhit.iphi()==seed_crys_iphi+j:
                        Ei=get_hitEnergy(hit)
                        crysIEta = ecalhit.ieta()
                        orginIEta=seed_crys_ieta
                        if (Ei>0):
                            meanDEta = meanDEta + (Ei * getNrCrysDiffInEta(crysIEta,orginIEta) )
        meanDEta = meanDEta/E5x5
    return meanDEta
                      
#

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
   
    ebhit_handle, ebhit_label = Handle("edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >"), "hltEgammaHLTExtra:EcalRecHitsEB"
    ele_handle, ele_label = Handle("std::vector<reco::EgTrigSumObj>"), "hltEgammaHLTExtra"
    if (args.genmatched=='1'): gen_handle, gen_label = Handle("std::vector<reco::GenParticle>"), "genParticles"
    weights = EvtWeights(args.weights)

    to_remove=[]
    for file in args.in_filename:
        file_temp = ROOT.TFile(file)
        if ( file_temp.IsZombie() ) :
            print "******** remove file **********"
            to_remove.append(file)

    new_list = [x for x in args.in_filename if x not in to_remove]
    print "AFTER REMOVING ALL ZOMBIES, this is the new list of files "
    print new_list

    events = Events(new_list)    
    #events = Events(args.in_filename)

    #ele_seed_time_EB = ROOT.TH1D("ele_seed_time_EB",";ele seed time;#entries",600,-6,6)
    ele_my_sieie_default_EB = ROOT.TH1D("ele_my_sieie_default_EB",";my_sieie_default;#entries",2000,0,2) 
    ele_my_sieie_NC_EB = ROOT.TH1D("ele_my_sieie_NC_EB",";my_sieie_NC;#entries",2000,0,2)
    ele_saved_sieie_default_EB = ROOT.TH1D("ele_saved_sieie_default_EB",";saved_sieie_default;#entries",2000,0,2) 
    ele_E5x5 = ROOT.TH1D("ele_E5x5","E5x5;#entries",100,0,100)

    noise_vs_ieta = ROOT.TH2D("noise_vs_ieta",";ieta;PF recHit threshold",200,-100,100,500,0,5)


    crysSize=0.01745
    
    #now we'll loop over the events rather than one by one
    for event_nr,event in enumerate(events):

        #recosf_file  = ROOT.TFile(args.recosf_file,"READ") 
        cond_file = ROOT.TFile("egamma__conditions__111X_mcRun4_realistic_T15_v1.root","READ")  
        rec_hit_thres = ROOT.EcalCondObjectContainer("float")()
        cond_file.EcalPFRecHitThresholdsRcd.SetBranchAddress("EcalCondObjectContainer_9float_0__",rec_hit_thres)
        cond_file.EcalPFRecHitThresholdsRcd.GetEntry(0)
        
        event.getByLabel(ebhit_label,ebhit_handle)
        event.getByLabel(ele_label,ele_handle)
        if (args.genmatched=='1'):  event.getByLabel(gen_label,gen_handle)

        ebhits = ebhit_handle.product()
        eles = ele_handle.product()
        if (args.genmatched=='1'): genobjs=gen_handle.product()

        if (args.genmatched=='1'):
            weight = 1

        if (args.genmatched=='0'):
            weight = weights.weight_from_evt(event.object())

        for eg in eles:
          min_dr=999.9

          if (args.genmatched=='1'):
              for gen in genobjs:
                  if math.fabs(gen.pdgId())==11:
                      dr2 = ROOT.reco.deltaR2(eg.eta(),eg.phi(),gen.eta(),gen.phi())
                      dr=math.sqrt(dr2)
                      if dr<min_dr: min_dr=dr
          elif (args.genmatched=='0'): min_dr=-1

          if (min_dr<=0.1) : 
              if eg.et()>=30.0 and eg.et()<50.0 and abs(eg.superCluster().eta())<1.44:
                  #print '\n\n\n***** new ele'
                  seed_id = eg.superCluster().seed().seed()
                  if seed_id.subdetId()==1: #not ecal barrel, return
                      seed_ebid = ROOT.EBDetId(seed_id)
                      seed_crys_ieta = seed_ebid.ieta()
                      seed_crys_iphi = seed_ebid.iphi()
                      E5x5=get_E5x5(eg,ebhits)
                      ele_E5x5.Fill(E5x5)
                      #print 'E5x5=', E5x5
                      #print '0.9% of E5x5 = ', E5x5*(0.9/100.0)
                      if E5x5>0: meanDEta=EtaBar5x5(eg,E5x5,ebhits)
                      numeratorEtaEta = 0
                      den=0
                      numeratorEtaEtaNC = 0
                      denNC=0
                      if (E5x5>0):
                        for hit in ebhits:
                          ecalhit=ROOT.EBDetId(hit.id())
                          
                          my_ieta=[0,1,-1,2,-2]
                          my_iphi=[0,1,-1,2,-2]
                          for i in my_ieta:
                              for j in my_iphi:
                                  if ecalhit.ieta()==seed_crys_ieta+i and ecalhit.iphi()==seed_crys_iphi+j:
                                      if ecalhit.subdetId()!=1: continue 
                                      Ei=get_hitEnergy(hit)

                                      crysIEta = ecalhit.ieta()
                                      orginIEta=seed_crys_ieta
                                      if (Ei>0): #default
                         
                                           dEta = getNrCrysDiffInEta(crysIEta,orginIEta) - meanDEta
                                           w=get_wi(Ei,E5x5) 
                                           
                                           numeratorEtaEta = numeratorEtaEta + (w * dEta * dEta)
                                           
                                           den=den+w
                                           
                                      noise=rec_hit_thres[ecalhit.rawId()]
                                      noise_vs_ieta.Fill(ecalhit.ieta(),noise)
                                     
                              
                                      if (Ei>0.0 and Ei>noise): #noise cleaned
                                  
                                           dEtaNC = getNrCrysDiffInEta(crysIEta,orginIEta) - meanDEta
                                  
                                           wNC=get_wi(Ei,E5x5) 
                                  
                                           numeratorEtaEtaNC = numeratorEtaEtaNC + (wNC * dEtaNC * dEtaNC)
                                  
                                           denNC=denNC+wNC
                                  
                      
        
        
                  if (den!=0):
                    covEtaEta2 = (crysSize * crysSize * numeratorEtaEta) / den
                    covEtaEta = math.sqrt(covEtaEta2)
                  else:
                    covEtaEta = 1


                  if (denNC!=0):
                    covEtaEtaNC2 = (crysSize * crysSize * numeratorEtaEtaNC) / denNC
                    covEtaEtaNC = math.sqrt(covEtaEtaNC2)
                  else:
                    covEtaEtaNC = 1

        
                  ele_my_sieie_default_EB.Fill(covEtaEta,weight)
                  ele_my_sieie_NC_EB.Fill(covEtaEtaNC,weight)
                  ele_saved_sieie_default_EB.Fill(eg.var("hltEgammaClusterShapeUnseeded_sigmaIEtaIEta5x5"),weight)
                  
    output_file = TFile( args.output, 'recreate' )
    #ele_seed_time_EB.Write()
    ele_my_sieie_default_EB.Write()
    ele_saved_sieie_default_EB.Write()
    ele_my_sieie_NC_EB.Write()
    noise_vs_ieta.Write()
    ele_E5x5.Write()
    output_file.Close()

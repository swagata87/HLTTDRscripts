#!/usr/bin/env python

#RUN it like this
#python CheckFilter.py  input.root  

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
import time

def getListFilterPassedObj(filterName,hltsevt):
        eg_trig_objs = []
        #get a list of trigger objects that passed a given filter
        filterIndex = getFilterIndex(hltsevt,filterName)
        if (filterIndex < hltsevt.sizeFilters() ):
                for filterKey in hltsevt.filterKeys(filterIndex):
                        obj = hltsevt.getObjects()[filterKey]
                        eg_trig_objs.append(obj)
        return eg_trig_objs


#from https://github.com/cms-egamma/EgammaDAS2020/blob/solutions/test/egFWLiteExercise2a.py
def match_trig_objs(eta,phi,trig_objs,max_dr=0.1):    
    max_dr2 = max_dr*max_dr
    matched_objs = [obj for obj in trig_objs if ROOT.reco.deltaR2(eta,phi,obj.eta(),obj.phi()) < max_dr2]
    return matched_objs

#from https://github.com/Sam-Harper/usercode/blob/100XNtup/SHNtupliser/test/checkTrigsAOD.py
def getFilterIndex(trigEvt,filterName):
    for index in range(0,trigEvt.sizeFilters()):
        print trigEvt.filterLabel(index)
        if filterName==trigEvt.filterLabel(index):
                return index
    return trigEvt.sizeFilters()

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

    args = parser.parse_args()

    ele_handle, ele_label = Handle("std::vector<reco::EgTrigSumObj>"), "hltEgammaHLTExtra"
    hlt_handle, hlt_label = Handle("edm::TriggerResults"), "TriggerResults::HLTX"
    hltevt_handle, hltevt_label = Handle("trigger::TriggerEvent"), "hltTriggerSummaryAOD"

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

    for event_nr,event in enumerate(events):

        event.getByLabel(ele_label,ele_handle)
        event.getByLabel(hlt_label,hlt_handle)
        event.getByLabel(hltevt_label,hltevt_handle)

        eles = ele_handle.product()
        hlts = hlt_handle.product()
        hltsevt = hltevt_handle.product()

        trigdict=event.object().triggerNames(hlts).triggerNames()

        #print "\nsieie"
        eg_trig_objs_sieieFilter = getListFilterPassedObj("hltEle5WPTightClusterShapeUnseededFilter",hltsevt)

        #print "\nH/E"
        eg_trig_objs_hoeFilter = getListFilterPassedObj("hltEle5WPTightHEUnseededFilter",hltsevt)

        #print "\necal iso"
        eg_trig_objs_ecalisoFilter = getListFilterPassedObj("hltEle5WPTightEcalIsoUnseededFilter",hltsevt)

        #print "\nhcal iso"
        eg_trig_objs_hcalisoFilter = getListFilterPassedObj("hltEle5WPTightHcalIsoUnseededFilter",hltsevt)

        #print "\nPM"
        eg_trig_objs_pmFilter = getListFilterPassedObj("hltEle5WPTightPixelMatchUnseededFilter",hltsevt)

        #print "\nPM S2"
        eg_trig_objs_pms2Filter = getListFilterPassedObj("hltEle5WPTightPMS2UnseededFilter",hltsevt)

        #print "\n1/E-1/p"
        eg_trig_objs_ooemoopFilter = getListFilterPassedObj("hltEle5WPTightGsfOneOEMinusOneOPUnseededFilter",hltsevt)

        #print "\nmissing hit"
        eg_trig_objs_mhitFilter = getListFilterPassedObj("hltEle5WPTightGsfMissingHitsUnseededFilter",hltsevt)

        #print "\ndEta"
        eg_trig_objs_detaFilter = getListFilterPassedObj("hltEle5WPTightGsfDetaUnseededFilter",hltsevt)

        #print "\ndPhi"
        eg_trig_objs_dphiFilter = getListFilterPassedObj("hltEle5WPTightGsfDphiUnseededFilter",hltsevt)

        #print "\n nLayerIT"
        eg_trig_objs_npixFilter = getListFilterPassedObj("hltEle5WPTightBestGsfNLayerITUnseededFilter",hltsevt)

        #print "\nchi2"
        eg_trig_objs_chi2Filter = getListFilterPassedObj("hltEle5WPTightBestGsfChi2UnseededFilter",hltsevt)

        #print "\ntrack iso"
        eg_trig_objs_trkisoFilter = getListFilterPassedObj("hltEle5WPTightGsfTrackIsoUnseededFilter",hltsevt)

        for eg in eles:
          
          print '\n\nnew ele'

          #print eg.varNames()

          #check if the electron matches with any trigger object that passed a given filter 
          matched_objs_sieieFilter = match_trig_objs(eg.superCluster().eta(),eg.superCluster().phi(),eg_trig_objs_sieieFilter)
          nmatch_sieieFilter = len(matched_objs_sieieFilter)

          matched_objs_hoeFilter = match_trig_objs(eg.superCluster().eta(),eg.superCluster().phi(),eg_trig_objs_hoeFilter)
          nmatch_hoeFilter = len(matched_objs_hoeFilter)

          matched_objs_ecalisoFilter = match_trig_objs(eg.superCluster().eta(),eg.superCluster().phi(),eg_trig_objs_ecalisoFilter)
          nmatch_ecalisoFilter = len(matched_objs_ecalisoFilter)

          matched_objs_hcalisoFilter = match_trig_objs(eg.superCluster().eta(),eg.superCluster().phi(),eg_trig_objs_hcalisoFilter)
          nmatch_hcalisoFilter = len(matched_objs_hcalisoFilter)

          matched_objs_pmFilter = match_trig_objs(eg.superCluster().eta(),eg.superCluster().phi(),eg_trig_objs_pmFilter)
          nmatch_pmFilter = len(matched_objs_pmFilter)

          matched_objs_pms2Filter = match_trig_objs(eg.superCluster().eta(),eg.superCluster().phi(),eg_trig_objs_pms2Filter)
          nmatch_pms2Filter = len(matched_objs_pms2Filter)

          matched_objs_ooemoopFilter = match_trig_objs(eg.superCluster().eta(),eg.superCluster().phi(),eg_trig_objs_ooemoopFilter)
          nmatch_ooemoopFilter = len(matched_objs_ooemoopFilter)

          matched_objs_mhitFilter = match_trig_objs(eg.superCluster().eta(),eg.superCluster().phi(),eg_trig_objs_mhitFilter)
          nmatch_mhitFilter = len(matched_objs_mhitFilter)

          matched_objs_detaFilter = match_trig_objs(eg.superCluster().eta(),eg.superCluster().phi(),eg_trig_objs_detaFilter)
          nmatch_detaFilter = len(matched_objs_detaFilter)

          matched_objs_dphiFilter = match_trig_objs(eg.superCluster().eta(),eg.superCluster().phi(),eg_trig_objs_dphiFilter)
          nmatch_dphiFilter = len(matched_objs_dphiFilter)

          matched_objs_npixFilter = match_trig_objs(eg.superCluster().eta(),eg.superCluster().phi(),eg_trig_objs_npixFilter)
          nmatch_npixFilter = len(matched_objs_npixFilter)

          matched_objs_chi2Filter = match_trig_objs(eg.superCluster().eta(),eg.superCluster().phi(),eg_trig_objs_chi2Filter)
          nmatch_chi2Filter = len(matched_objs_chi2Filter)

          matched_objs_trkisoFilter = match_trig_objs(eg.superCluster().eta(),eg.superCluster().phi(),eg_trig_objs_trkisoFilter)
          nmatch_trkisoFilter = len(matched_objs_trkisoFilter)

          
                      
          print 'nmatch_trkisoFilter ', nmatch_trkisoFilter
          print 'nmatch_chi2Filter ', nmatch_chi2Filter
          print 'nmatch_npixFilter ', nmatch_npixFilter
          print 'nmatch_dphiFilter ', nmatch_dphiFilter
          print 'nmatch_detaFilter ', nmatch_detaFilter
          print 'nmatch_sieieFilter ', nmatch_sieieFilter

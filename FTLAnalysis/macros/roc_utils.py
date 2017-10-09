#!/bin/python

import sys
import re
import time
import math
import argparse
import os
import subprocess
import bisect
import ROOT
from string import maketrans

from ExternalTools.FuriousPlotter.fp_utils import *
from array import array
import numpy 

def roc(args, srcs):
    sig_tree = srcs[args[0]]
    bkg_tree = srcs[args[1]]

    table = maketrans('[]', '()')
    den_cut_sig = args[2].translate(table, '"')
    #num_cut_sig = den_cut_sig+" && "+args[4].translate(translate_table, '"')
    den_cut_bkg = args[3].translate(table, '"')
    #num_cut_bkg = den_cut_bkg+" && "+args[4].translate(translate_table, '"')

    ###---define arrays for sig & bkg values and errors
    N_num_sig = array('d', (0,)*int(args[5]))
    N_num_bkg = array('d', (0,)*int(args[5]))    
    eff_sig = array('d', (0,)*int(args[5]))
    eff_sig_erhi = array('d', (0,)*int(args[5]))
    eff_sig_erlo = array('d', (0,)*int(args[5]))
    eff_bkg = array('d', (0,)*int(args[5]))
    eff_bkg_erhi = array('d', (0,)*int(args[5]))
    eff_bkg_erlo = array('d', (0,)*int(args[5]))
    
    ###---get tot events
    elistsig = ROOT.TEventList("elistsig")
    elistbkg = ROOT.TEventList("elistbkg")    
    N_den_sig = float(sig_tree.Draw(">>elistsig", den_cut_sig))
    N_den_bkg = float(bkg_tree.Draw(">>elistbkg", den_cut_bkg))

    elistsig.Print()
    print(sig_tree.GetEntriesFast(), elistsig.GetN())
    print(bkg_tree.GetEntriesFast(), elistbkg.GetN())    
    
    ###---selections formulas
    test_values = numpy.linspace(float(args[6]), float(args[7]), int(args[5]), False)

    selection_histo = ROOT.TH1F("selection_histo", "", int(args[5]), float(args[6]), float(args[7]))

    for iEntry in range(elistsig.GetN()):
        ###---read the tre
        tEntry = elistsig.GetEntry(iEntry)
        if int(iEntry)%int(50000) == 0:
            print(iEntry, tEntry)            
        
        ###---get number of selected events
        sig_tree.Draw(args[4].translate(None, '"')+">>selection_histo", "", "goff", 1, tEntry)
        ibin = selection_histo.FindFirstBinAbove(0.5)
        if ibin > -1:
            for i in range(ibin, int(args[5])):        
                N_num_sig[i] = N_num_sig[i]+1
        selection_histo.Reset()

    for iEntry in range(elistbkg.GetN()):
        ###---read the tree
        tEntry = elistbkg.GetEntry(iEntry)
        if int(iEntry)%int(50000) == 0:
            print(iEntry, tEntry)            
        
        ###---get number of selected events
        bkg_tree.Draw(args[4].translate(None, '"')+">>selection_histo", "", "goff", 1, tEntry)
        ibin = selection_histo.FindFirstBinAbove(0.5)
        if ibin > -1:
            for i in range(ibin, int(args[5])):
                N_num_bkg[i] = N_num_bkg[i]+1
        selection_histo.Reset()
        
    for i, value in enumerate(test_values):                
        ###---compute eff and errors for sig
        eff_sig[i] = N_num_sig[i]/N_den_sig
        eff_sig_erhi[i] = ROOT.TEfficiency.ClopperPearson(N_den_sig, N_num_sig[i], 0.683, True) - eff_sig[i]
        eff_sig_erlo[i] = eff_sig[i] - ROOT.TEfficiency.ClopperPearson(N_den_sig, N_num_sig[i], 0.683, False)

        ###---compute eff and errors for bkg
        eff_bkg[i] = N_num_bkg[i]/N_den_bkg
        eff_bkg_erhi[i] = ROOT.TEfficiency.ClopperPearson(N_den_bkg, N_num_bkg[i], 0.683, True) - eff_bkg[i]
        eff_bkg_erlo[i] = eff_bkg[i] - ROOT.TEfficiency.ClopperPearson(N_den_bkg, N_num_bkg[i], 0.683, False)

    ###---create curve as TGraphAsymmErrors
    tmp = ROOT.TGraphAsymmErrors(int(args[5]), eff_sig, eff_bkg, eff_sig_erlo, eff_sig_erhi, eff_bkg_erlo, eff_bkg_erhi)

    return tmp

def abs_rate_vs_eff(args, srcs):
    sig_tree = srcs[args[0]]
    bkg_tree = srcs[args[1]]

    table = maketrans('[]', '()')
    den_cut_sig = args[2].translate(table, '"')
    #num_cut_sig = den_cut_sig+" && "+args[4].translate(translate_table, '"')
    den_cut_bkg = args[3].translate(table, '"')
    #num_cut_bkg = den_cut_bkg+" && "+args[4].translate(translate_table, '"')

    ###---define arrays for sig & bkg values and errors
    N_num_sig = array('d', (0,)*int(args[5]))
    N_num_bkg = array('d', (0,)*int(args[5]))    
    eff_sig = array('d', (0,)*int(args[5]))
    eff_sig_erhi = array('d', (0,)*int(args[5]))
    eff_sig_erlo = array('d', (0,)*int(args[5]))
    eff_bkg = array('d', (0,)*int(args[5]))
    eff_bkg_erhi = array('d', (0,)*int(args[5]))
    eff_bkg_erlo = array('d', (0,)*int(args[5]))
    
    ###---get tot events
    elistsig = ROOT.TEventList("elistsig")
    elistbkg = ROOT.TEventList("elistbkg")    
    N_den_sig = float(sig_tree.Draw(">>elistsig", den_cut_sig))
    bkg_tree.Draw(">>elistbkg", den_cut_bkg)
    N_den_bkg = 0

    elistsig.Print()
    print(sig_tree.GetEntriesFast(), elistsig.GetN())
    print(bkg_tree.GetEntriesFast(), elistbkg.GetN())    
    
    ###---selections formulas
    test_values = numpy.linspace(float(args[6]), float(args[7]), int(args[5]), False)

    selection_histo = ROOT.TH1F("selection_histo", "", int(args[5]), float(args[6]), float(args[7]))

    for iEntry in range(elistsig.GetN()):
        ###---read the tre
        tEntry = elistsig.GetEntry(iEntry)
        if int(iEntry)%int(50000) == 0:
            print(iEntry, tEntry)            
        
        ###---get number of selected events
        sig_tree.Draw(args[4].translate(None, '"')+">>selection_histo", "", "goff", 1, tEntry)
        ibin = selection_histo.FindFirstBinAbove(0.5)
        if ibin > -1:
            for i in range(ibin, int(args[5])):        
                N_num_sig[i] = N_num_sig[i]+1
        selection_histo.Reset()

    ###---read the tree
    crun = clumi = cevent = -1
    for iEntry in range(min(bkg_tree.GetEntriesFast(), 1000)):
        if int(iEntry)%int(50000) == 0:
            print(iEntry)
        bkg_tree.Draw("event:lumi:run:"+args[4].translate(None, '"'), "", "", 1, iEntry)
        events = bkg_tree.GetV1()
        lumis = bkg_tree.GetV2()
        runs = bkg_tree.GetV3()    
        if events[0] != cevent or lumis[0] != clumi or runs[0] != crun:
            cevent = events[0]
            clumi = lumis[0]
            crun = runs[0]
            N_den_bkg = N_den_bkg+1
        
        ###---get number of selected events
        if elistbkg.Contains(iEntry):
            selection_histo.Fill(bkg_tree.GetV4()[0])
            ibin = selection_histo.FindFirstBinAbove(0.5)
            if ibin > -1:
                for i in range(ibin, int(args[5])):
                    N_num_bkg[i] = N_num_bkg[i]+1
            selection_histo.Reset()

    for i, value in enumerate(test_values):                
        ###---compute eff and errors for sig
        eff_sig[i] = N_num_sig[i]/N_den_sig
        eff_sig_erhi[i] = ROOT.TEfficiency.ClopperPearson(N_den_sig, N_num_sig[i], 0.683, True) - eff_sig[i]
        eff_sig_erlo[i] = eff_sig[i] - ROOT.TEfficiency.ClopperPearson(N_den_sig, N_num_sig[i], 0.683, False)

        ###---compute eff and errors for bkg
        eff_bkg[i] = N_num_bkg[i]/N_den_bkg
        eff_bkg_erhi[i] = eff_bkg_erlo[i] = math.sqrt(N_num_bkg[i])/N_den_bkg

    ###---create curve as TGraphAsymmErrors
    tmp = ROOT.TGraphAsymmErrors(int(args[5]), eff_sig, eff_bkg, eff_sig_erlo, eff_sig_erhi, eff_bkg_erlo, eff_bkg_erhi)

    return tmp

def eff_scan(args, srcs):
    """Scan efficiency as a function of x-axis variable:
    - takes 1 argument which is a TH1
    - compute cumulative distribution and scales
    """
    
    tmp = srcs[args[0]].GetCumulative()
    tmp.Scale(1./srcs[args[0]].Integral())

    print srcs[args[0]].GetName(), tmp.GetBinCenter(tmp.FindFirstBinAbove(0.95))
    
    return tmp
    
dictionary = dict(ROC=roc, AbsRateVsEff=abs_rate_vs_eff, EffScan=eff_scan)

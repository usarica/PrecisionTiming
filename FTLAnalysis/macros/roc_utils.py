#!/bin/python

import sys
import re
import time
import argparse
import os
import subprocess
import ROOT

from ExternalTools.FuriousPlotter.fp_utils import *
from array import array
import numpy 

def roc(args, srcs):
    sig_tree = srcs[args[0]]
    bkg_tree = srcs[args[1]]
    
    den_cut_sig = args[2].translate(None, '"')
    num_cut_sig = den_cut_sig+" && "+args[4].translate(None, '"')
    den_cut_bkg = args[3].translate(None, '"')
    num_cut_bkg = den_cut_bkg+" && "+args[4].translate(None, '"')

    ###---define arrays for sig & bkg values and errors
    eff_sig = array('d', (0,)*int(args[5]))
    eff_sig_erhi = array('d', (0,)*int(args[5]))
    eff_sig_erlo = array('d', (0,)*int(args[5]))
    eff_bkg = array('d', (0,)*int(args[5]))
    eff_bkg_erhi = array('d', (0,)*int(args[5]))
    eff_bkg_erlo = array('d', (0,)*int(args[5]))
    
    ###---get tot events
    N_den_sig = float(sig_tree.Draw("Entry$", den_cut_sig, "goff"))
    N_den_bkg = float(bkg_tree.Draw("Entry$", den_cut_bkg, "goff"))

    for i, value in enumerate(numpy.linspace(float(args[6]), float(args[7]), int(args[5]), False)):
        cut_sig = num_cut_sig+str(value)
        cut_bkg = num_cut_bkg+str(value)

        ###---get selected events
        N_num_sig = float(sig_tree.Draw("Entry$", cut_sig, "goff"))
        N_num_bkg = float(bkg_tree.Draw("Entry$", cut_bkg, "goff"))

        ###---compute eff and errors for sig
        eff_sig[i] = N_num_sig/N_den_sig
        eff_sig_erhi[i] = ROOT.TEfficiency.ClopperPearson(N_den_sig, N_num_sig, 0.683, True) - eff_sig[i]
        eff_sig_erlo[i] = eff_sig[i] - ROOT.TEfficiency.ClopperPearson(N_den_sig, N_num_sig, 0.683, False)

        ###---compute eff and errors for bkg
        eff_bkg[i] = N_num_bkg/N_den_bkg
        eff_bkg_erhi[i] = ROOT.TEfficiency.ClopperPearson(N_den_bkg, N_num_bkg, 0.683, True) - eff_bkg[i]
        eff_bkg_erlo[i] = eff_bkg[i] - ROOT.TEfficiency.ClopperPearson(N_den_bkg, N_num_bkg, 0.683, False)

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
    
dictionary = dict(ROC=roc, EffScan=eff_scan)

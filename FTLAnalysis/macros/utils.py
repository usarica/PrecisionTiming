#!/bin/python

import sys
import re
import time
import argparse
import os
import subprocess
import ROOT

###---Effective sigma: args[0] = fraction of events to be contained inside the effective sigma interval
###   the distribution is assumed to be peaked around its mean
def effective_sigma(args, srcs):
    """Compute effective sigma of TH1"""

    bins = srcs.GetNbinsX()-1
    quantile = float(args[0])*srcs.Integral(1, bins) 
    left = 1    
    right = bins

    for lbin in range(1, bins):
        for rbin in range(bins, lbin, -1):
            if srcs.Integral(lbin, rbin) >= quantile:
                if right-left < rbin-lbin:
                    left = lbin
                    right = rbin
            else:
                break

    return (srcs.GetBinCenter(right)-srcs.GetBinCenter(left))/2
    
def effective_sigma_2d(args, srcs):
    """Compute effective sigma of TH2"""

    quantile = args[1] if len(args)>1 else "0.68"
    tmp = srcs[args[0]].ProjectionX("tmp_x", 1, 2, "o")
    tmp.Reset()
    for ibin in range(1, srcs[args[0]].GetNbinsX()-1):
        tmp_y = srcs[args[0]].ProjectionY("tmp_y", ibin, ibin, "eo")
        if tmp_y.GetEntries() > 0:
            tmp.SetBinContent(ibin, effective_sigma(quantile, tmp_y))
        tmp_y.Delete()        

    return tmp                               

def rms_projection(args, srcs):
    """Compute RMS of Y distribution as a function of X"""

    tmp = srcs[args[0]].ProjectionX("tmp_x", 1, 2, "o")
    tmp.Reset()
    for ibin in range(1, srcs[args[0]].GetNbinsX()-1):
        tmp_y = srcs[args[0]].ProjectionY("tmp_y", ibin, ibin, "eo")
        if tmp_y.GetEntries() > 0:
            tmp.SetBinContent(ibin, tmp_y.GetRMS())
            tmp.SetBinError(ibin, tmp_y.GetRMSError())
        tmp_y.Delete()        

    return tmp                               

def rms_global_projection(args, srcs):
    """Compute global RMS of crystal history: 
    this is done using a standard TProfile2D and swapping the bin error with the bin content"""

    h_tmp = srcs[args[0]]
    for xbin in range(1, srcs[args[0]].GetNbinsX()):
        for ybin in range(1, srcs[args[0]].GetNbinsY()):
            h_tmp.SetBinContent(xbin, ybin, srcs[args[0]].GetBinError(xbin, ybin))

    return h_tmp

def make_graph_with_errors(args, srcs):
    """Combine two histograms: the first one set the point value, the second the point error"""

    h_tmp = ROOT.TGraphErrors()
    for xbin in range(1, srcs[args[0]].GetNbinsX()):
        h_tmp.SetPoint(xbin-1, srcs[args[0]].GetBinCenter(xbin), srcs[args[0]].GetBinContent(xbin))
        h_tmp.SetPointError(xbin-1, 0, srcs[args[1]].GetBinContent(xbin))

    return h_tmp
                          
dictionary = dict(EffSigma=effective_sigma, EffSigma2D=effective_sigma_2d,
                  RMSProj=rms_projection, RMSMap=rms_global_projection,
                  MakeHistoErrors=make_graph_with_errors)


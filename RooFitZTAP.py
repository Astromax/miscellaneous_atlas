#!/usr/bin/env python
"""
    NAME: RooFitZTAP.py: this script is intended to look over data which has run through the ZTAP selector

    OPTIONS: fill in later
"""
#--------------------------------------------------------------------------------------
### Imports
#--------------------------------------------------------------------------------------
## std                                                                                                                                                                                                           
import os, argparse, sys, glob, time

#special                                                                                                                                                                                                              
import ROOT as root
import math
from array import array
import ast

root.gROOT.SetBatch(root.kTRUE)

## Globals                                                                                                                                                                                                    
v = 'v13'
time = time.strftime('%m-%d')
path = 'RFZTAP_storage_%s_%s' % (time, v)
os.mkdir(path)
RooFitZTAP_log = open('%s/RFZTAP_log_data_%s_%s.txt' % (path,time,v), 'w')
#--------------------------------------------------------------------------------------
### Options
#--------------------------------------------------------------------------------------
def options():
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--target', default='ZTAP_Full.txt',
            help='The list of target files')
    parser.add_argument('-sl', '--selections', default='ZTAP_electron_selections.txt',
            help='The list of numerator, denominator, and common selections, all independent of the bin-slicing')
    parser.add_argument('-s', '--slicings', default='ZTAP_slicings.txt',
            help='The list of variable slicings to use')    
    return parser.parse_args()

#--------------------------------------------------------------------------------------
### Main
#--------------------------------------------------------------------------------------
def main():
    ops = options()

    #Particulars 
    #Note to self: modify this kw0/mode thing to make it smoother
    fittype = 'BW+L'
    tgt = ops.target
    tgtwords = tgt.split("_")
    if tgtwords[1] == 'Test.txt':
        kw0 = 'tree_NoSys'
        mode = 'Test'
    elif tgtwords[1] == 'Full.txt':
        kw0 = 'data'
        mode = 'Full'
    elif tgtwords[-1] == 'MC.txt':
        kw0 = '_NoSys'
        mode = 'MC'
        
    kw1 = 'zmass'

    RooFitZTAP_log.write('Running with fittype %s on %s' % (fittype, mode) + ' \n')

    #First: open the list of files & save them into a list
    samplelist = list()
    samples = open(ops.target, 'r')
    for s in samples:
        samplelist.append(s)

    for sl in samplelist:
        print 'We have this sample in samplelist'
        print(sl)

    fpairs = NameTester(samplelist)
    fsams = [root.TFile(a) for a,b in fpairs]
    fnames = [b for a,b in fpairs]

    for f in fpairs:
        print 'We have this fpair'
        print(f)

    #Second: define the selections & the bin slicings
    selwords = ops.selections.split('_')
    lepmode = selwords[1]
    RooFitZTAP_log.write('Running in lepmode %s' % lepmode + ' \n')
    selections = list()
    selectionset = open(ops.selections, 'r')
    for line in selectionset:
        selection = ast.literal_eval(line)
        if selection['Category'] == 'Numerator':
            numselbase = selection['Criteria']
        elif selection['Category'] == 'Denominator':
            denomselbase = selection['Criteria']
        elif selection['Category'] == 'Common':
            selections.append(selection['Criteria'])

    slicings = list()
    slics = open(ops.slicings, 'r')
    for line in slics:
        slicing = ast.literal_eval(line)
        slicings.append(slicing)

    #Third: loop over the samples, apply the fit, save to pngs and root file
    ZTAP_Efficiencies = root.TFile('%s/ZTAP_%s_efficiencies_%s_%s_%s.root' % (path, lepmode, fittype, time, v) , 'RECREATE')
    for i,f in enumerate(fsams):
        keys = f.GetListOfKeys()
        for key in keys:
            if kw0 in key.GetName():
                t = key.GetName()
        mytree = f.Get(t)
        fname = fnames[i]
        for sel in selections:
            numselbase +='&&'+sel
            denomselbase +='&&'+sel
            for slic in slicings:
                if slic['Dimensionality'] == 1:
                    dim = 1
                    var = slic['Branch']
                    bin_edges = slic['bin_edges']
                    bins = len(bin_edges) - 1
                    Tracklets    = root.TH1D("Tracklets", "Tracklets vs %s" % var, bins, min(bin_edges), max(bin_edges))
                    Leptons      = root.TH1D("Leptons", "Leptons vs %s" % var, bins, min(bin_edges), max(bin_edges))
                    Efficiencies = root.TH1D("Efficiencies", "Efficiencies vs %s" % var, bins, min(bin_edges), max(bin_edges))
                    for j, be in enumerate(bin_edges):
                        if be == max(bin_edges):
                            continue
                        low = be
                        high = bin_edges[j+1]
                        limits = [low, high]
                        numsel = numselbase+'&& %s>%f && %s<%f' % (var, low, var, high)
                        denomsel = denomselbase+'&& %s>%f && %s<%f' % (var, low, var, high)
                        sels = [numsel, denomsel]
                        output  = General_Fitter(mode, lepmode, fname, mytree, kw1, dim, var, sels, limits, time, fittype)
                        numPair = output['Tracklets']
                        denomPair = output['Leptons']
                        effPair = output['Efficiency']
                        Tracklets.SetBinContent(j+1, numPair[0])
                        Tracklets.SetBinError(j+1, numPair[1])
                        Leptons.SetBinContent(j+1, denomPair[0])
                        Leptons.SetBinError(j+1, denomPair[1])
                        Efficiencies.SetBinContent(j+1, effPair[0])
                        Efficiencies.SetBinError(j+1, effPair[1])
                    effC = root.TCanvas("effHist", "Efficiency vs %s" % var, 800, 400)
                    Efficiencies.SetStats(0)
                    Efficiencies.GetXaxis().SetTitle('%s' % var)
                    effC.SetLogy()
                    Efficiencies.Draw("E1")
                    Efficiencies.SetName('%s_efficiency' % var)
                    Efficiencies.SetTitle('%s_efficiency' % var)
                    effC.SaveAs("%s/%s_%s_%s_Efficiency_vs_%s_%s.png" % (path, mode, lepmode, fname, var, time))
                elif slic['Dimensionality'] == 2:
                    dim = 2
                    varx = slic['xBranch']
                    vary = slic['yBranch']
                    varpair = [varx, vary]
                    xbin_edges = slic['xbin_edges']
                    ybin_edges = slic['ybin_edges']
                    xbins = len(xbin_edges) - 1
                    ybins = len(ybin_edges) - 1
                    Tracklets    = root.TH2D("Tracklets", "Tracklet Transfer Map vs %s and %s" % (varx, vary), xbins, min(xbin_edges), max(xbin_edges), ybins, min(ybin_edges), max(ybin_edges))
                    Leptons      = root.TH2D("Leptons", "Lepton Transfer Map vs %s and %s" % (varx, vary), xbins, min(xbin_edges), max(xbin_edges), ybins, min(ybin_edges), max(ybin_edges))
                    Efficiencies = root.TH2D("Efficiencies", "Efficiency Transfer Map vs %s and %s" % (varx, vary), xbins, min(xbin_edges), max(xbin_edges), ybins, min(ybin_edges), max(ybin_edges))
                    for j, x in enumerate(xbin_edges):
                        if x == max(xbin_edges):
                            continue
                        xhigh = xbin_edges[j+1]
                        xlimits = [x, xhigh]
                        for k, y in enumerate(ybin_edges):
                            if y == max(ybin_edges):
                                continue
                            yhigh = ybin_edges[k+1]
                            ylimits = [y, yhigh]
                            limits = [xlimits, ylimits]
                            numsel = numselbase+'&& %s>%f && %s<%f && %s>%f && %s<%f' % (varx, x, varx, xhigh, vary, y, vary, yhigh)
                            denomsel = denomselbase+'&& %s>%f && %s<%f && %s>%f && %s<%f' % (varx, x, varx, xhigh, vary, y, vary, yhigh)
                            sels = [numsel, denomsel]
                            output  = General_Fitter(mode, lepmode, fname, mytree, kw1, dim, varpair, sels, limits, time, fittype)
                            numPair = output['Tracklets']
                            denomPair = output['Leptons']
                            effPair = output['Efficiency']
                            Tracklets.SetBinContent(j+1, k+1, numPair[0])
                            Tracklets.SetBinError(j+1, k+1, numPair[1])
                            Leptons.SetBinContent(j+1, k+1, denomPair[0])
                            Leptons.SetBinError(j+1, k+1, denomPair[1])
                            Efficiencies.SetBinContent(j+1, k+1, effPair[0])
                            Efficiencies.SetBinError(j+1, k+1, effPair[1])
                    effC = root.TCanvas("effHist", 'Efficiency Transfer Map vs %s and %s' % (varx, vary), 800, 400)
                    Efficiencies.SetStats(0)
                    Efficiencies.GetXaxis().SetTitle('%s' % varx)
                    Efficiencies.GetYaxis().SetTitle('%s' % vary)
                    Efficiencies.Draw("COLZ")
                    Efficiencies.SetName('Efficiency_TransferMap_vs_%s_and_%s' % (varx, vary))
                    Efficiencies.SetTitle('Efficiency_TransferMap_vs_%s_and_%s' % (varx, vary))
                    effC.SaveAs('%s/%s_%s_%s_Efficiency_TransferMap_vs_%s_and_%s_%s.png' % (path, mode, lepmode, fname, varx, vary, time))
                elif slic['Dimensionality'] == 3:
                    dim = 3
                    varx = slic['xBranch']
                    vary = slic['yBranch']
                    varz = slic['zBranch']
                    varset = [varx,vary,varz]
                    xbin_edges = slic['xbin_edges']
                    ybin_edges = slic['ybin_edges']
                    zbin_edges = slic['zbin_edges']
                    xbins = len(xbin_edges) - 1
                    ybins = len(ybin_edges) - 1
                    zbins = len(zbin_edges) - 1
                    Tracklets = root.TH3D("Tracklets", "Tracklet Cube vs %s and %s and %s" % (varx, vary, varz), xbins, min(xbin_edges), max(xbin_edges), ybins, min(ybin_edges), max(ybin_edges), zbins, min(zbin_edges), max(zbin_edges))
                    Leptons  = root.TH3D("Leptons", "Lepton Cube vs %s and %s and %s" % (varx, vary, varz), xbins, min(xbin_edges), max(xbin_edges), ybins, min(ybin_edges), max(ybin_edges), zbins, min(zbin_edges), max(zbin_edges))
                    Efficiencies = root.TH3D("Efficiencies", "Efficiency Cube vs %s and %s and %s" % (varx, vary, varz), xbins, min(xbin_edges), max(xbin_edges), ybins, min(ybin_edges), max(ybin_edges), zbins, min(zbin_edges), max(zbin_edges))
                    for j,x in enumerate(xbin_edges):
                        if x == max(xbin_edges):
                            continue
                        xhigh = xbin_edges[j+1]
                        xlimits = [x, xhigh]
                        for k,y in enumerate(ybin_edges):
                            if y == max(ybin_edges):
                                continue
                            yhigh = ybin_edges[k+1]
                            ylimits = [y, yhigh]
                            for l,z in enumerate(zbin_edges):
                                if z == max(zbin_edges):
                                    continue
                                zhigh = zbin_edges[l+1]
                                zlimits = [z, zhigh]
                                limits = [xlimits, ylimits, zlimits]
                                numsel = numselbase+'&& %s>%f && %s<%f && %s>%f && %s<%f && %s>%f && %s<%f' % (varx, x, varx, xhigh, vary, y, vary, yhigh, varz, z, varz, zhigh)
                                denomsel = denomselbase+'&& %s>%f && %s<%f && %s>%f && %s<%f && %s>%f && %s<%f' % (varx, x, varx, xhigh, vary, y, vary, yhigh, varz, z, varz, zhigh)
                                sels = [numsel, denomsel]
                                output  = General_Fitter(mode, lepmode, fname, mytree, kw1, dim, varset, sels, limits, time, fittype)
                                numPair = output['Tracklets']
                                denomPair = output['Leptons']
                                effPair = output['Efficiency']
                                Tracklets.SetBinContent(j+1, k+1, l+1, numPair[0])
                                Tracklets.SetBinError(j+1, k+1, l+1, numPair[1])
                                Leptons.SetBinContent(j+1, k+1, l+1, denomPair[0])
                                Leptons.SetBinError(j+1, k+1, l+1, denomPair[1])
                                Efficiencies.SetBinContent(j+1, k+1, l+1, effPair[0])
                                Efficiencies.SetBinError(j+1, k+1, l+1, effPair[1])
                ZTAP_Efficiencies.Write()
                del Efficiencies
                if dim==1 or dim==2:
                    del effC

    return
#--------------------------------------------------------------------------------------
### General Fitter
#--------------------------------------------------------------------------------------
def General_Fitter(mode, lepmode, fname, mytree, kw1, dim, var, sels, limits, time, fittype):
    efflist = list()

    numsel = sels[0]
    denomsel = sels[1]

    RooFitZTAP_log.write('Tracklet Selection: %s' % numsel+ ' \n')
    RooFitZTAP_log.write('Lepton Selection: %s' % denomsel+ ' \n')

    # Eliminate decimal points from ranges (decimal points sometimes cause confusion in output files, ex: Fit_eta_1.0_2.0.png may cause an error but Fit_eta_1p0_2p0.png will not)
    if dim==1:
        outlims = list()
        for l in limits:
            if '.' in str(l):
                outlwrds = str(l).split('.')
                outl = outlwrds[0]+'p'+outlwrds[1]
            else:
                outl = l
            outlims.append(outl)
    elif dim==2:
        outxlims = list()
        outylims = list()
        for x in limits[0]:
            if '.' in str(x):
                outxwrds = str(x).split('.')
                outx = outxwrds[0]+'p'+outxwrds[1]
            else:
                outx = x
            outxlims.append(outx)
        for y in limits[1]:
            if '.' in str(y):
                outywrds = str(y).split('.')
                outy = outywrds[0]+'p'+outywrds[1]
            else:
                outy = y
            outylims.append(outy)
    elif dim==3:
        outxlims = list()
        outylims = list()
        outzlims = list()
        for x in limits[0]:
            if '.' in str(x):
                outxwrds = str(x).split('.')
                outx = outxwrds[0]+'p'+outxwrds[1]
            else:
                outx = x
            outxlims.append(outx)
        for y in limits[1]:
            if '.' in str(y):
                outywrds = str(y).split('.')
                outy = outywrds[0]+'p'+outywrds[1]
            else:
                outy = y
            outylims.append(outy)
        for z in limits[2]:
            if '.' in str(z):
                outzwrds = str(z).split('.')
                outz = outzwrds[0]+'p'+outzwrds[1]
            else:
                outz = z
            outzlims.append(outz)


    # Note to self: generalize by having it read in a dictionary with the relevant parameters, like selections.txt for the scanner
    temptracklethist = root.TH1D("temptracklethist", "temptracklethist", 30, 60, 120)
    templeptonhist = root.TH1D("templeptonhist", "templeptonhist", 30, 60, 120)

    mytree.Draw('%s>>temptracklethist' % kw1, numsel)
    mytree.Draw('%s>>templeptonhist' % kw1, denomsel)

    # Establish the basic variables needed to do this Fit (these would also need to be in the dictionary)
    # Pattern: trk_means = [trk_mean_guess, trk_mean_low, trk_mean_high], same order as RooRealVar
    trk_means = [91.19, 70, 112]
    trk_sigmas = [5, 2.5, 35]
    trk_lambdas = [0, -5, 5]
    trk_bkg_fracs = [0.1, 0, 1]
    trk_sig_fracs = [0.9, 0, 1]

    lep_means = [91.19, 70, 112]
    lep_sigmas = [3, 2.5, 30]
    lep_lambdas = [0, -5, 5]
    lep_bkg_fracs = [0.1, 0, 1]
    lep_sig_fracs = [0.9, 0, 1]

    m = root.RooRealVar("m", "Invariant Mass", 60, 120)
    nummean = root.RooRealVar("nummean", "numerator mean mass", trk_means[0], trk_means[1], trk_means[2])
    numsigma = root.RooRealVar("numsigma", "numerator sigma", trk_sigmas[0], trk_sigmas[1], trk_sigmas[2])
    
    denommean = root.RooRealVar("denommean", "denominator mean mass",  lep_means[0],  lep_means[1],  lep_means[2])
    denomsigma = root.RooRealVar("denomsigma", "denominator sigma",  lep_sigmas[0],  lep_sigmas[1],  lep_sigmas[2])
        
    tlam = root.RooRealVar("tlam", "tracklet lambda", trk_lambdas[0], trk_lambdas[1], trk_lambdas[2])
    llam = root.RooRealVar("llam", "lepton lambda", lep_lambdas[0], lep_lambdas[1], lep_lambdas[2])

    nt = temptracklethist.Integral()
    nl = templeptonhist.Integral()
    tbkg = root.RooExponential("tbkg","Background for Tracklets", m, tlam)
    ntbfrac = root.RooRealVar("ntbfrac", "Tracklet Background Fraction", trk_bkg_fracs[0], trk_bkg_fracs[1], trk_bkg_fracs[2])
    ntsfrac = root.RooRealVar("ntsfrac", "Tracklet Signal Fraction", trk_sig_fracs[0], trk_sig_fracs[1], trk_sig_fracs[2])

    lbkg = root.RooExponential("lbkg", "Background for Leptons", m, llam)
    nlbfrac = root.RooRealVar("nlbfrac", "Lepton Background Fraction", lep_bkg_fracs[0], lep_bkg_fracs[1], lep_bkg_fracs[2])
    nlsfrac = root.RooRealVar("nlsfrac", "Lepton Signal Fraction", lep_sig_fracs[0], lep_sig_fracs[1], lep_sig_fracs[2])

    RooFitZTAP_log.write('Below Parameters are arranged as in RooRealVar: (Initial Guess, Low End of Range, High End of Range) \n')
    RooFitZTAP_log.write('Tracklet Mean Parameters: (%f, %f, %f)' % (trk_means[0], trk_means[1], trk_means[2])+ ' \n')
    RooFitZTAP_log.write('Tracklet Sigma Parameters: (%f, %f, %f)' % (trk_sigmas[0], trk_sigmas[1], trk_sigmas[2])+ ' \n')
    RooFitZTAP_log.write('Tracklet Lambda Parameters: (%f, %f, %f)' % (trk_lambdas[0], trk_lambdas[1], trk_lambdas[2])+' \n')
    RooFitZTAP_log.write('Tracklet Background Fraction Parameters: (%f, %f, %f)' % (trk_bkg_fracs[0], trk_bkg_fracs[1], trk_bkg_fracs[2])+' \n')
    RooFitZTAP_log.write('Tracklet Signal Fraction Parameters: (%f, %f, %f)' % (trk_sig_fracs[0], trk_sig_fracs[1], trk_sig_fracs[2])+' \n')

    RooFitZTAP_log.write('Lepton Mean Parameters: (%f, %f, %f)' % (lep_means[0], lep_means[1], lep_means[2])+ ' \n')
    RooFitZTAP_log.write('Lepton Sigma Parameters: (%f, %f, %f)' % (lep_sigmas[0], lep_sigmas[1], lep_sigmas[2])+ ' \n')
    RooFitZTAP_log.write('Lepton Lambda Parameters: (%f, %f, %f)' % (lep_lambdas[0], lep_lambdas[1], lep_lambdas[2])+' \n')
    RooFitZTAP_log.write('Lepton Background Fraction Parameters: (%f, %f, %f)' % (lep_bkg_fracs[0], lep_bkg_fracs[1], lep_bkg_fracs[2])+' \n')
    RooFitZTAP_log.write('Lepton Signal Fraction Parameters: (%f, %f, %f)' % (lep_sig_fracs[0], lep_sig_fracs[1], lep_sig_fracs[2])+' \n')

    if dim==1:
        RooFitZTAP_log.write('The tracklet histogram integral for %s from %s to %s is %f' % (var, outlims[0], outlims[1], nt)+' \n')
        RooFitZTAP_log.write('The lepton histogram integral for %s from %s to %s is %f' % (var, outlims[0], outlims[1], nl)+' \n')
    elif dim==2:
        RooFitZTAP_log.write('The tracklet histogram integral for %s from %s to %s and %s from %s to %s is %f' % (var[0], outxlims[0], outxlims[1], var[1], outylims[0], outylims[1], nt)+' \n')
        RooFitZTAP_log.write('The lepton histogram integral for %s from %s to %s and %s from %s to %s is %f' % (var[0], outxlims[0], outxlims[1], var[1], outylims[0], outylims[1], nl)+' \n')

    if fittype == 'Breit-Wigner':
        trackletsig = root.RooBreitWigner("trkbreitwigner", "trkbreitwigner", m, nummean, numsigma)
        leptonsig   = root.RooBreitWigner("lepbreitwigner", "lepbreitwigner", m, denommean, denomsigma)
    elif fittype == 'Bukin':
        #trackletsig = root.RooBukin("trkbukin", "trkbukin", m, nummean, ...)
        #leptonsig   = root.RooBukin("lepbukin", "lepbukin", m, denommean, ...)
        print 'Bukin mode unavailable, please try again later'
        return (0,0)
    elif fittype == 'Landau':
        trackletsig = root.RooLandau("trklandau", "trklandau", m, nummean, numsigma)
        leptonsig   = root.RooLandau("leplandau", "leplandau", m, denommean, denomsigma)
    elif fittype == 'BW+L':
        trk_landau_sigmas = [3, 3.0, 35]
        trk_landau_fracs = [0.2, 0, 1]
        trk_BW_fracs = [0.8, 0, 1]
        lep_landau_sigmas = [3, 3.0, 35]
        lep_landau_fracs = [0.2, 0, 1]
        lep_BW_fracs = [0.8, 0, 1]

        tlausigma = root.RooRealVar("tlausigma", "tracklet landau sigma", trk_landau_sigmas[0], trk_landau_sigmas[1], trk_landau_sigmas[2])
        llausigma = root.RooRealVar("llausigma", "lepton landau sigma", lep_landau_sigmas[0], lep_landau_sigmas[1], lep_landau_sigmas[2])
        tlaufrac = root.RooRealVar("tlaufrac", "Tracklet Signal Landau component", trk_landau_fracs[0], trk_landau_fracs[1], trk_landau_fracs[2])
        tbrwfrac = root.RooRealVar("tbwfrac", "Tracklet Signal BW component", trk_BW_fracs[0], trk_BW_fracs[1], trk_BW_fracs[2])
        llaufrac = root.RooRealVar("llaufrac", "Lepton Signal Landau component", lep_landau_fracs[0], lep_landau_fracs[1], lep_landau_fracs[2])
        lbrwfrac = root.RooRealVar("lbwfrac", "Lepton Signal BW component", lep_BW_fracs[0], lep_BW_fracs[1], lep_BW_fracs[2])
        trackletsig1 = root.RooBreitWigner("trkbreitwigner", "trkbreitwigner", m, nummean, numsigma)
        trackletsig2 = root.RooLandau("trklandau", "trklandau", m, nummean, tlausigma)
        trackletsig = root.RooAddPdf("trackletsig", "Breit-Wigner + Landau Signal", root.RooArgList(trackletsig1, trackletsig2), root.RooArgList(tbrwfrac))

        leptonsig1 = root.RooBreitWigner("lepbreitwigner", "lepbreitwigner", m, denommean, denomsigma)
        leptonsig2 = root.RooLandau("leplandau", "leplandau", m, denommean, llausigma)
        leptonsig = root.RooAddPdf("leptonsig", "Breit-Wigner + Landau Signal", root.RooArgList(leptonsig1, leptonsig2), root.RooArgList(lbrwfrac))

        RooFitZTAP_log.write('Tracklet Landau Parameters: (%f, %f, %f)' % (trk_landau_sigmas[0], trk_landau_sigmas[1], trk_landau_sigmas[2])+ ' \n')
        RooFitZTAP_log.write('Tracklet Landau Fraction: (%f, %f, %f)' % (trk_landau_fracs[0], trk_landau_fracs[1], trk_landau_fracs[2])+ ' \n')
        RooFitZTAP_log.write('Tracklet Breit-Wigner Fraction: (%f, %f, %f)' % (trk_BW_fracs[0], trk_BW_fracs[1], trk_BW_fracs[2])+ ' \n')
        RooFitZTAP_log.write('Lepton Landau Parameters: (%f, %f, %f)' % (lep_landau_sigmas[0], lep_landau_sigmas[1], lep_landau_sigmas[2])+ ' \n')
        RooFitZTAP_log.write('Lepton Landau Fraction: (%f, %f, %f)' % (lep_landau_fracs[0], lep_landau_fracs[1], lep_landau_fracs[2])+ ' \n')
        RooFitZTAP_log.write('Lepton Breit-Wigner Fraction: (%f, %f, %f)' % (lep_BW_fracs[0], lep_BW_fracs[1], lep_BW_fracs[2])+ ' \n')


    trackletmodel = root.RooAddPdf("trackletmodel", "%s Signal + Falling Exponential Background, %s" % (fittype ,numsel), root.RooArgList(trackletsig, tbkg), root.RooArgList(ntsfrac))
    leptonmodel = root.RooAddPdf("leptonmodel", "%s Peak + Falling Exponential Background, %s" % (fittype, denomsel), root.RooArgList(leptonsig, lbkg), root.RooArgList(nlsfrac))

    # Now fit the two models to their respective histograms
    trkdh = root.RooDataHist("trkdh", "trkdh", root.RooArgList(m), root.RooFit.Import(temptracklethist))
    lepdh = root.RooDataHist("lepdh", "lepdh", root.RooArgList(m), root.RooFit.Import(templeptonhist))

    trackletmodel.fitTo(trkdh)
    leptonmodel.fitTo(lepdh)

    # Plot the fits on frames
    trkframe = m.frame(root.RooFit.Title("tracklets, %s" % numsel))
    lepframe = m.frame(root.RooFit.Title("leptons, %s" % denomsel))

    trkdh.plotOn(trkframe)
    lepdh.plotOn(lepframe)
    trackletmodel.plotOn(trkframe)
    leptonmodel.plotOn(lepframe)

    trkC = root.TCanvas("trkBWFE", "Tracklet Breit-Wigner + Falling Exponential", 800, 400)
    trkframe.Draw()

    lepC = root.TCanvas("lepBWFE", "Tracklet Breit-Wigner + Falling Exponential", 800, 400)
    lepframe.Draw()

    if dim==1:
        trkC.SaveAs('%s/%s_%s_%s_%s_Tracklet_Fit_%s_%s_%s_%s_%s.png' % (path, mode, lepmode, fittype, fname, var, outlims[0], outlims[1], time, v))
        lepC.SaveAs('%s/%s_%s_%s_%s_Lepton_Fit_%s_%s_%s_%s_%s.png' % (path, mode, lepmode, fittype, fname, var, outlims[0], outlims[1], time, v))
    elif dim==2:
        varx = var[0]
        vary = var[1]
        trkC.SaveAs('%s/%s_%s_%s_%s_Tracklet_Fit_%s_%s_%s_%s_%s_%s_%s_%s.png' % (path, mode, lepmode, fittype, fname, varx, outxlims[0], outxlims[1], vary, outylims[0], outylims[1], time, v))
        lepC.SaveAs('%s/%s_%s_%s_%s_Lepton_Fit_%s_%s_%s_%s_%s_%s_%s_%s.png' % (path, mode, lepmode, fittype, fname, varx, outxlims[0], outxlims[1], vary, outylims[0], outylims[1], time, v))
    elif dim==3:
        varx = var[0]
        vary = var[1]
        varz = var[2]
        trkC.SaveAs('%s/%s_%s_%s_Tracklet_Fit_%s_%s_%s_%s_%s_%s_%s_%s_%s_%s_%s.png' % (path, mode, lepmode, fname, varx, outxlims[0], outxlims[1], vary, outylims[0], outylims[1], varz, outzlims[0], outzlims[1], time, v))
        lepC.SaveAs('%s/%s_%s_%s_Lepton_Fit_%s_%s_%s_%s_%s_%s_%s_%s_%s_%s_%s.png' % (path, mode, lepmode, fname, varx, outxlims[0], outxlims[1], vary, outylims[0], outylims[1], varz, outzlims[0], outzlims[1], time, v))

    numfrac = ntsfrac.getVal()
    num = numfrac * nt
    numerr = ntsfrac.getError() * nt

    denomfrac = nlsfrac.getVal()
    denom = denomfrac * nl
    denomerr = nlsfrac.getError() * nl

    cut_eff = 0
    if denom == 0 or num == 0:
        simple_eff = 0
        simple_efferr = 0
        bruce_efferr = 0
    else:
        simple_eff = (num / denom)
        cut_eff = (nt/nl)
        simple_efferr = Error(num, denom, numerr, denomerr)
        bruce_efferr = simple_eff * math.sqrt((1 / num) + (1 / denom))

    RooFitZTAP_log.write('Raw Tracklet Integral: %f' % nt+' \n')
    RooFitZTAP_log.write('Raw Lepton Integral: %f' % nl+' \n')
    RooFitZTAP_log.write('Numerator: %f' % num+' \n')
    RooFitZTAP_log.write('Numerator Uncertainty: %f' % numerr+' \n')
    RooFitZTAP_log.write('Denominator: %f' % denom+' \n')
    RooFitZTAP_log.write('Denominator Uncertainty: %f' % denomerr+' \n')
    RooFitZTAP_log.write('Efficiency: %f' % simple_eff+' \n')
    RooFitZTAP_log.write('Efficiency Uncertainty: %f' % simple_efferr+' \n')
    RooFitZTAP_log.write('Efficiency Error (Bruce): %f' % bruce_efferr+' \n')

    if mode == 'Test':
        print 'Raw tracklet integral: %f' % nt
        print 'Raw lepton integral: %f' % nl
        print 'Numerator: %f' % num
        print 'Numerator Error: %f' % numerr
        print 'Denominator: %f' % denom
        print 'Denominator Error: %f' % denomerr
        print 'Simple Efficiency: %f' % simple_eff
        print 'Simple Efficiency Error: %f' % simple_efferr
        print 'Cut based eff: %f' % cut_eff

    efflist.append(simple_eff)
    efflist.append(simple_efferr)

    output = {'Tracklets': [num, numerr], 'Leptons': [denom, denomerr], 'Efficiency': [simple_eff, bruce_efferr]}

    del temptracklethist
    del templeptonhist
    del trkdh
    del lepdh
    del trackletmodel
    del leptonmodel
    del trkframe
    del lepframe
    del trkC
    del lepC

    return output
#-------------------------------------------------------------------------------------
# Name Tester: takes in the list of files & outputs appropriate fnames for them (my attempt at idiot-proofing)
#-------------------------------------------------------------------------------------
def NameTester(samples):
    uniquefiles = list()
    filepairs = list()
    truepairs = list()
    truepaths = list()

    for i, file in enumerate(samples):
        file = file.rstrip()
        if file not in uniquefiles and os.path.isfile(file):
            uniquefiles.append(file)
            bname = file.split('/')[-1]
            bname = bname[:-5]         #drop the .root portion
            filepair = (file, bname)
            filepairs.append(filepair)

    for fp in filepairs:
        if fp[0] in truepaths:
            continue
        path = fp[0]
        bn = fp[1]
        dups = [f for f in filepairs if f[1] == bn and f[0] != path]
        if len(dups) > 0:
            initpathwords = path.split('/')
            for i,d in enumerate(dups):
                pathwords = d[0].split('/')
                for j,pw in enumerate(pathwords):
                    if pw == initpathwords[j]:
                        continue
                    else:
                        distinct = pw
                        distinct_index = j
                        continue
                new_bname = distinct+'_'+d[1]
                newpair = (d[0], new_bname)
                truepairs.append(newpair)
            distinct = initpathwords[distinct_index]
            new_fname = distinct+'_'+bn
            new_pair = (path, new_fname)
            truepairs.append(new_pair)
        else:
            truepairs.append(fp)
        truepaths = [t[0] for t in truepairs]

    return truepairs
##------------------------------------------------------------------------------------                                                                                                                      
# Error Calculation                                                                                                                                                                                                                   
##------------------------------------------------------------------------------------                                                                                                                                        
def Error(v1, v2, e1, e2):

    Error = math.sqrt( (v1*e2)**2 + (v2*e1)**2 )/(v2**2)

    return Error

#--------------------------------------------------------------------------------------
if __name__ == '__main__': main()

# EOF
    

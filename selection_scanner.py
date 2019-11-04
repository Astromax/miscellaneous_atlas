#!/usr/bin/env python 
"""
NAME
    selection_scanner.py: reads in root files & tries out assorted selection criteria
OPTIONS
    
DESCRIPTION
    Reads in signal & background root files, tries out assorted selection criteria, and produces a csv file
    containing, along the top row, the selection criteria (independent except for the base), with each row showing the 
    number of entries for that particular signal/background given the criterion above, with the bottom rows showing a proxy 
    Figure of Merit for each signal, calculated by taking the signal value & dividing it by the square root of the sum of the 
    background numbers within it's column.  This proxy is NOT the actual FoM expected from the experiment, but it is *A* figure of merit

"""
#------------------------------------------------------------------------------                                                                                                                                  
### Imports                                                                                                                                                                                          
#------------------------------------------------------------------------------                                                                                                                        
## std                                                                                                                                                                                                  
import os, argparse, sys, glob, time

#special                                                                                                                                                                                                        
import ROOT as root
import math
import csv
import ast

root.gROOT.SetBatch(root.kFALSE)

## Globals
time = time.strftime('%m-%d')
#--------------------------------------------------------------------------------                                                                                                            
### Options                                                                                                                                                                                  
#--------------------------------------------------------------------------------             
def options():
    parser = argparse.ArgumentParser()
    parser.add_argument('-sg', '--signals', default='signals.txt',
            help='List of signal files')
    parser.add_argument('-bkg', '--backgrounds', default='backgrounds.txt',
            help='List of background files')
    parser.add_argument('-sel', '--selections', default='selections.txt',
            help='List of selection criteria')
    return parser.parse_args()

#------------------------------------------------------------------------------
### Main
#------------------------------------------------------------------------------
def main():
    ops = options()

    v = 'v0'
    mode = 'GO'
    lumi = 140000

    ##Particulars (things which may be different in your ROOT files than in mine)
    kw0 = '_NoSys'
    kw1 = 'nVtx'
    
    ##First section: acquire all of the relevant files 
    output = list()
    fsigs = list()
    fbkgs = list()

    signals = open(ops.signals, 'r')
    for file in signals:
        file = file.rstrip()
        sig = root.TFile(file)
        fsigs.append(sig)

    backgrounds = open(ops.backgrounds, 'r')
    for file in backgrounds:
        file = file.rstrip()
        bkg = root.TFile(file)
        fbkgs.append(bkg)

    ##Second section: import selection options, set the weight & base selection and the top row of the output csv file
    weight = 'genWeight*eventWeight*leptonWeight*jvtWeight*bTagWeight' 
    basesel = '((leadjet_MET_dPhi!=999 && abs(leadjet_MET_dPhi)>1.0) || leadjet_MET_dPhi == 999) && ((jet2_MET_dPhi!=999 && abs(jet2_MET_dPhi)>1.0) || jet2_MET_dPhi==999) && ((jet3_MET_dPhi!=999 && abs(jet3_MET_dPhi)>1.0) || jet3_MET_dPhi==999) && ((jet4_MET_dPhi!=999 && abs(jet4_MET_dPhi)>1.0) || jet4_MET_dPhi==999)'
    usedbases = list()
    seltypes = list()
    selections = open(ops.selections, 'r')
    for line in selections:
        sel = ast.literal_eval(line)
        seltypes.append(sel)

    if mode == 'DEBUG':
        print 'This is the full list of seltypes'
        for st in seltypes:
            print(st)

    for st in seltypes:
        print 'Now checking seltype %s, which is a %s' % (st['Name'], st['type'])
        if st in usedbases:
            continue
        if basesel == '':
            if st['type'] == 'baseequal':
                basesel = '%s==%f' % (st['Name'], st['lower_bound'])
                usedbases.append(st)
            elif st['type'] == 'basefloor':
                basesel = '%s>%f' % (st['Name'], st['lower_bound'])
                usedbases.append(st)
            elif st['type'] == 'baseceiling':
                basesel = '%s<%f' % (st['Name'], st['upper_bound'])
                usedbases.append(st)
            elif st['type'] == 'baseabsfloor':
                basesel = 'abs(%s)>%f' % (st['Name'], st['lower_bound'])
                usedbases.append(st)
            elif st['type'] == 'baseabsceiling':
                basesel = 'abs(%s)<%f' % (st['Name'], st['upper_bound'])
                usedbases.append(st)
        else:
            if st['type'] == 'baseequal':
                basesel += ' && %s==%f' % (st['Name'], st['lower_bound'])
                usedbases.append(st)
            elif st['type'] == 'basefloor':
                basesel += ' && %s>%f' % (st['Name'], st['lower_bound'])
                usedbases.append(st)
            elif st['type'] == 'baseceiling':
                basesel += ' && %s<%f' % (st['Name'], st['upper_bound'])
                usedbases.append(st)
            elif st['type'] == 'baseabsfloor':
                basesel += ' && abs(%s)>%f' % (st['Name'], st['lower_bound'])
                usedbases.append(st)
            elif st['type'] == 'baseabsceiling':
                basesel += ' && abs(%s)<%f' % (st['Name'], st['upper_bound'])
                usedbases.append(st)


    scantypes = [s for s in seltypes if s not in usedbases]

    if mode == 'DEBUG':
        print 'This is the list of scantypes'
        for st in scantypes:
            print(st)

    toprow = ['Sample', 'Raw Histogram Integral', basesel]
    for st in scantypes:
        histmin = st['lower_bound']
        histmax = st['upper_bound']
        ntest = st['ntest']
        gap = histmax - histmin
        interval = float(gap)/ntest
        for i in range(ntest+1):
            if st['type'] == 'floor':
                minimum = histmin + i*interval
                addsel = '%s>%f' % (st['Name'], minimum)
            elif st['type'] == 'ceiling':
                maximum = histmax - i*interval
                addsel = '%s<%f' % (st['Name'], maximum)
            elif st['type'] == 'absfloor':
                minimum = histmin + i*interval
                addsel = 'abs(%s)>%f' % (st['Name'], minimum)
            elif st['type'] == 'absceiling':
                maximum = histmax - i*interval
                addsel = 'abs(%s)<%f' % (st['Name'], maximum)
            elif st['type'] == 'equal':
                equality = histmin + i*interval
                addsel = '%s==%f' % (st['Name'], equality)
            toprow.append(addsel)
    output.append(toprow)
    if mode == 'DEBUG':
        print 'This is the top row'
        print(toprow)

    ## Third section: loop over the files & get the numbers for output 
    for f in fsigs:
        keys = f.GetListOfKeys()
        for key in keys:
            if kw0 in key.GetName():
                t = key.GetName()
        mytree = f.Get(t)
        sigrow = segmenter(t, mytree, scantypes, basesel, kw1, lumi, weight, mode)
        output.append(sigrow)

    for f in fbkgs:
        keys = f.GetListOfKeys()
        for key in keys:
            if kw0 in key.GetName():
                t = key.GetName()
        mytree = f.Get(t)
        bkgrow = segmenter(t, mytree, scantypes, basesel, kw1, lumi, weight, mode)
        output.append(bkgrow)

    for i, s in enumerate(fsigs):
        keys = s.GetListOfKeys()
        for key in keys:
            if kw0 in key.GetName():
                t = key.GetName()
        tname = t[:-6]
        proxy_row = ['Proxy FoM %s' % tname]
        for j in range(len(toprow)):
            if j==0:
                continue
            signal = 0
            totalbkg = 0            
            for k, r in enumerate(output):
                if k<i+1:
                    continue
                elif k==i+1:
                    signal += float(r[j])
                elif k>i+1 and k<=len(fsigs):
                    continue
                elif k<=(len(fsigs) + len(fbkgs)):
                    if r[j] < 0.0:
                        print 'WARNING: Background value below zero, you are likely in a low-stats region, resetting this term to zero'
                        r[j] = 0
                    totalbkg += r[j]
                else:
                    continue
            FoM_proxy = significance(signal, totalbkg, 'complete')
            proxy_row.append(FoM_proxy)
        output.append(proxy_row)

    ## Fourth section: dump everything into a csv file
    lumistring = str(lumi)
    if '.' in lumistring:
        lumiwords = lumistring.split('.')
        outlumi = lumiwords[0]+'p'+lumiwords[1]
    else:
        outlumi = lumistring

    with open('test_selection_scan_lumi_%s_ipb_%s_%s.csv' % (outlumi, time, v), 'wb') as f:
        writer = csv.writer(f)
        writer.writerows(output)

    return
#--------------------------------------------------------------------------------
### engine
#--------------------------------------------------------------------------------
def segmenter(t, mytree, scantypes, basesel, kw1, lumi, weight, mode):

    branches = mytree.GetListOfBranches()
    if kw1 not in branches:
        print 'WARNING: missing the dummy variable branch, exiting segmenter'
        return
    
    raw_hist = root.TH1D('raw_hist', 'raw_hist', 35, 0, 35)
    mytree.Draw('%s>>raw_hist' % kw1, '%f*%s' % (lumi, weight))

    base_hist = root.TH1D('base_hist', 'base_hist', 35, 0, 35)
    mytree.Draw('%s>>base_hist' % kw1, '('+basesel+')'+'*%f*%s' % (lumi, weight))
    tname = t[:-6]

    row = list()
    row.append(tname)

    raw_integral = raw_hist.Integral()
    row.append(raw_integral)

    base_integral = base_hist.Integral()
    row.append(base_integral)

    if mode == 'DEBUG':
        print 'The tname is %s' % tname
        print 'The integral of the raw histogram with luminosity %f ipb is %0.5f' % (lumi, raw_integral)
        print 'The integral of the base selection histogram with luminosity %f ipb is %0.5f' % (lumi, base_integral)

    for st in scantypes:
        histmin = st['lower_bound']
        histmax = st['upper_bound']
        ntest = st['ntest']
        gap = histmax - histmin
        interval = float(gap)/ntest
        for i in range(ntest+1):
            if st['type'] == 'floor':
                minimum = histmin + i*interval
                testsel = basesel+' && %s>%f' % (st['Name'], minimum)
            elif st['type'] == 'ceiling':
                maximum = histmax - i*interval
                testsel = basesel+' && %s<%f' % (st['Name'], maximum)
            elif st['type'] == 'absfloor':
                minimum = histmin + i*interval
                testsel = basesel+' && abs(%s)>%f' % (st['Name'], minimum)
            elif st['type'] == 'absceiling':
                maximum = histmax - i*interval
                testsel = basesel+' && abs(%s)<%f' % (st['Name'], maximum)
            elif st['type'] == 'equal':
                equality = histmin + i*interval
                testsel = basesel+' && %s==%f' % (st['Name'], equality)
            test_hist = root.TH1D('test_hist', 'test_hist', 35, 0, 35)
            mytree.Draw('%s>>test_hist' % kw1, '('+testsel+')'+'*%f*%s' % (lumi, weight))
            row.append(test_hist.Integral())
            if mode == 'DEBUG':
                print 'The integral of testsel %s is %f' % (testsel, test_hist.Integral())
            del test_hist

    del raw_hist
    del base_hist

    return row
#--------------------------------------------------------------------------------
## Figure of Merit Calculation
#--------------------------------------------------------------------------------
def significance(signal, background, mode):

    if signal<=0 or background <=0:
        return 0

    if mode == 'simple':
        significance = signal/math.sqrt(background)
    elif mode == 'complex':
        significance = math.sqrt(2 * ((signal + background) * math.log(1 + (signal/background)) - signal))
    elif mode == 'complete':
        n = signal + background
        sigma = 0.3 * background
        #sigma = math.sqrt(background)
        term1 = n * (background + sigma**2)/(background**2 + n * sigma**2)
        term2 = (signal * sigma**2)/(background*(background + sigma**2))
        significance = math.sqrt(2 * (n * math.log(term1) - math.pow((background/sigma), 2) * math.log(1 + term2)))

    return significance
#--------------------------------------------------------------------------------                                                                                                                                           
if __name__ == '__main__': main()

# EOF

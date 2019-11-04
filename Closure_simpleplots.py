#!/usr/bin/env python 


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
time = time.strftime('%m-%d')
#--------------------------------------------------------------------------------------                                                                                                              
### Options                                                                                                                                                                                              
#--------------------------------------------------------------------------------------                                                                                                               
def options():
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--target', default='Closure_sample.txt',
            help='The list of target files')
    parser.add_argument('-sl', '--selections', default='Closure_electron_selections.txt',
            help='The list of numerator, denominator, and common selections, all independent of the bin-slicing')
    parser.add_argument('-s', '--slicings', default='Closure_slicings.txt',
            help='The list of variable slicings to use')
    return parser.parse_args()

#--------------------------------------------------------------------------------------                                                                                                              
### Main                                                                               
#--------------------------------------------------------------------------------------
def main():
    ops = options()

    v = 'v0'
    kw0 = '_NoSys'

    # First section: get the samples 
    fsams = list()
    for file in open(ops.target, 'r'):
        file = file.rstrip()
        fsam = root.TFile(file)
        fsams.append(fsam)

    # Second section: get the selections & slicings
    selwords = ops.selections.split('_')
    lepmode = selwords[1]
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

    # Third section: scroll through the sample, produce the efficiency & dump to root file
    Closure_Efficiencies = root.TFile('Closure_Efficiencies_%s_%s_%s.root' % (lepmode, time, v) , 'RECREATE')
    for i,f in enumerate(fsams):
        keys = f.GetListOfKeys()
        for key in keys:
            if kw0 in key.GetName():
                t = key.GetName()
                break
        mytree = f.Get(t)
        for sel in selections:
            numselbase +='&&'+sel
            denomselbase +='&&'+sel
            for slic in slicings:
                if slic['Dimensionality'] == 1:
                    var = slic['Branch']
                    bin_edges = slic['bin_edges']
                    bins = len(bin_edges) - 1
                    tempnumhist = root.TH1D('tempnumhist', '%s_numerator' % var, bins, min(bin_edges), max(bin_edges))
                    tempdenomhist = root.TH1D('tempdenomhist', '%s_denominator' % var, bins, min(bin_edges), max(bin_edges))
                    tempeffhist = root.TH1D('tempeffhist', '%s_efficiency' % var, bins, min(bin_edges), max(bin_edges))

                    numC = root.TCanvas('numerator', 'Numerator', 800, 400)
                    tempnumhist.SetTitle('Numerator: %s' % numselbase)
                    tempnumhist.SetXTitle(var)
                    tempnumhist.SetYTitle('Raw Number of Events')
                    tempnumhist.SetStats(0)
                    mytree.Draw('%s>>tempnumhist' % var, numselbase)
                    numC.SaveAs('Closure_W_%s_numerator_vs_%s_%s_%s.png' % (lepmode, var, time, v))

                    denomC = root.TCanvas('denominator', 'Denominator', 800, 400)
                    tempdenomhist.SetTitle('Denominator: %s' % denomselbase)
                    tempdenomhist.SetXTitle(var)
                    tempdenomhist.SetYTitle('Raw Number of Events')
                    tempdenomhist.SetStats(0)
                    mytree.Draw('%s>>tempdenomhist' % var, denomselbase)
                    denomC.SaveAs('Closure_W_%s_denominator_vs_%s_%s_%s.png' % (lepmode, var, time, v))

                    effC = root.TCanvas('efficiency', 'Efficiency', 800, 400)
                    tempeffhist.SetTitle('Efficiency vs %s' % var)
                    tempeffhist.SetXTitle(var)
                    tempeffhist.SetYTitle('Efficiency')
                    tempeffhist.SetStats(0)
                    effC.SetLogy()
                    tempeffhist.Divide(tempnumhist, tempdenomhist, 1, 1, 'B')
                    for bin in range(bins):
                        num = tempnumhist.GetBinContent(bin)
                        numerr = tempnumhist.GetBinError(bin)
                        denom = tempdenomhist.GetBinContent(bin)
                        denomerr = tempdenomhist.GetBinError(bin)
                        if denom == 0:
                            efferr = 0
                        else:
                            efferr = Error(num, denom, numerr, denomerr)
                        tempeffhist.SetBinError(bin, efferr)

                    tempeffhist.Draw('E1')
                    tempnumhist.SetName('Closure_%s_numerator' % var)
                    tempdenomhist.SetName('Closure_%s_denominator' % var)
                    tempeffhist.SetName('Closure_%s_efficiency' % var)
                    Closure_Efficiencies.Write()
                    effC.SaveAs('Closure_W_%s_Efficiency_vs_%s_%s_%s.png' % (lepmode, var, time, v))

                    del tempnumhist
                    del tempdenomhist
                    del tempeffhist
                    del numC
                    del denomC
                    del effC
                else:
                    print '2 Dimensional functionality not yet available'
                    varx = slic['xBranch']
                    vary = slic['yBranch']
                    xbin_edges = slic['xbin_edges']
                    ybin_edges = slic['ybin_edges']
                    xbins = len(xbin_edges) - 1
                    ybins = len(ybin_edges) - 1
                    tempnumhist = root.TH2D('tempnumhist', '%s_vs_%s_numerator' % (varx, vary), xbins, min(xbin_edges), max(xbin_edges), ybins, min(ybin_edges), max(ybin_edges))
                    tempdenomhist = root.TH2D('tempdenomhist', '%s_vs_%s_denominator' % (varx, vary), xbins, min(xbin_edges), max(xbin_edges), ybins, min(ybin_edges), max(ybin_edges))
                    tempeffhist = root.TH2D('tempeffhist', '%s_vs_%s_efficiency' % (varx, vary), xbins, min(xbin_edges), max(xbin_edges), ybins, min(ybin_edges), max(ybin_edges))

                    numC = root.TCanvas('numerator', 'Numerator', 800, 400)
                    tempnumhist.SetTitle('Numerator: %s' % numselbase)
                    tempnumhist.SetXTitle(varx)
                    tempnumhist.SetYTitle(vary)
                    tempnumhist.SetStats(0)
                    mytree.Draw('%s:%s>>tempnumhist' % (vary, varx), numselbase)
                    numC.SaveAs('Closure_W_%s_numerator_vs_%s_and_%s_%s_%s.png' % (lepmode, varx, vary, time, v))

                    denomC = root.TCanvas('denominator', 'Denominator', 800, 400)
                    tempdenomhist.SetTitle('Denominator: %s' % denomselbase)
                    tempdenomhist.SetXTitle(varx)
                    tempdenomhist.SetYTitle(vary)
                    tempdenomhist.SetStats(0)
                    mytree.Draw('%s:%s>>tempdenomhist' % (vary, varx), denomselbase)
                    denomC.SaveAs('Closure_W_%s_denominator_vs_%s_and_%s_%s_%s.png' % (lepmode, varx, vary, time, v))

                    effC = root.TCanvas('efficiency', 'Efficiency', 800, 400)
                    tempeffhist.SetTitle('Efficiency vs %s and %s' % (varx,vary))
                    tempeffhist.SetXTitle(varx)
                    tempeffhist.SetYTitle(vary)
                    tempeffhist.SetStats(0)
                    effC.SetLogy()
                    tempeffhist.Divide(tempnumhist, tempdenomhist, 1, 1, 'B')
                    for xbin in range(xbins):
                        for ybin in range(ybins):
                            num = tempnumhist.GetBinContent(bin)
                            numerr = tempnumhist.GetBinError(bin)
                            denom = tempdenomhist.GetBinContent(bin)
                            denomerr = tempdenomhist.GetBinError(bin)
                            if denom == 0:
                                efferr = 0
                            else:
                                efferr = Error(num, denom, numerr, denomerr)
                            tempeffhist.SetBinError(bin, efferr)

                    tempnumhist.SetName('Closure_%s_%s_numerator' % (varx,vary))
                    tempdenomhist.SetName('Closure_%s_%s_denominator' % (varx,vary))
                    tempeffhist.SetName('Closure_%s_%s_efficiency' % (varx,vary))
                    Closure_Efficiencies.Write()
                    effC.SaveAs('Closure_W_%s_Efficiency_vs_%s_and_%s_%s_%s.png' % (lepmode, varx, vary, time, v))

                    del tempnumhist
                    del tempdenomhist
                    del tempeffhist
                    del numC
                    del denomC
                    del effC
            print 'All slicings for selection %s should be done' % sel
        print 'All selections should be done'
    print 'Done!'

    return
##------------------------------------------------------------------------------------                                                                                                                                                            
# Error Calculation                                                                                                                                                                                                                                 
##------------------------------------------------------------------------------------                                                                                                                                                     
def Error(v1, v2, e1, e2):
    Error = math.sqrt( (v1*e2)**2 + (v2*e1)**2 )/(v2**2)
    return Error
#--------------------------------------------------------------------------------------                                                                                                                        
if __name__ == '__main__': main()

# EOF

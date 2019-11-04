#!/usr/bin/env python
"""
NAME
      heatmap_postprocess.py: the purpose of this script is to fill in the holes in the fake rate heatmaps
OPTIONS
     -h, --help: Print manual & exits

2019-07-06

"""
#------------------------------------------------------------------------------
### Imports
#------------------------------------------------------------------------------                                                    

## std
import os, argparse, sys, glob, time

#special
import ROOT as root
import math
import ast

root.gROOT.SetBatch(root.kFALSE)

## Globals
v = 'v5'
time = time.strftime('%m-%d')
heatmap_post_process_log = open('HMPP_log_%s_%s.txt' % (time,v), 'w')
#-------------------------------------------------------------------------------
### Options
#-------------------------------------------------------------------------------
def options():
    parser = argparse.ArgumentParser()
    parser.add_argument('-rs', '--ratiosource', default='Closure_vs_ZTAP_muon_comparison_NoTM_11Oct.root',
                        help='The root file with the Ratio Map inside')
    parser.add_argument('-r', '--rootfile', default='RFZTAP_storage_10-11_v13/ZTAP_electron_efficiencies_BW+L_10-11_v13.root',
                        help='The rootfile with the relevant heatmap')
    parser.add_argument('-hm', '--heatmaps', default='muon_heatmaps.txt',
                        help='The list of relevant heatmaps')
    return parser.parse_args()
#-------------------------------------------------------------------------------
### Main
#-------------------------------------------------------------------------------
def main():
    ops = options()

#    rmaps = list()
#    heatmap_post_process_log.write('Now opening the ratio source file... \n')
#    R = root.TFile(ops.ratiosource, 'update')
#    if R.GetListOfKeys().Contains('probePt_Ratio_Map'):
#        print 'The file contains the ratio map!'
#        rmaps.append(R.Get('Ratio_Map'))
#    else:
#        print 'The file does not contain the ratio map, must create it'
#        heatmap_post_process_log.write('Creating the Ratio Map... \n')
        #num = R.Get('Closure_probeEta_probePt_efficiency')
        #denom = R.Get('Efficiency_TransferMap_vs_probeEta_and_probePt')
#        num = R.Get('Closure_probePt_efficiency')
#        denom = R.Get('probePt_efficiency')
#        cl = num.Clone()
#        cl.Divide(denom)
#        cl.SetTitle('Ratio vs Pt (GeV)')
#        cl.GetYaxis().SetTitle('Ratio')
#        cl.SetName('probePt_Ratio_Map')
#        R.Write()
#        rmaps.append(cl)
#        heatmap_post_process_log.write('Ratio Map has been created. \n')
#
#    print 'Phase 1 complete, now moving to Phase 2'
#    print 'We have this Ratio Map'
#    print(rmaps)
#    return


    heatmap_post_process_log.write('Now retriving rootfile... \n')
    f = root.TFile(ops.rootfile, 'update')
    
    heatmaps = list()
    for h in open(ops.heatmaps, 'r'):
        h = h.rstrip()
        heatmaps.append(h)

    heatmap_post_process_log.write('We have %i heatmaps to process' % len(heatmaps) + ' \n')

    for h in heatmaps:
        hm = f.Get(h)
        heatmap_post_process_log.write('Attempting to clone heatmap... \n')
        clone = hm.Clone()
        numXbins = clone.GetNbinsX()
        numYbins = clone.GetNbinsY()
        heatmap_post_process_log.write('Heatmap has been cloned, it has %i x-axis and %i y-axis bins' % (numXbins, numYbins)+' \n')
        for xbin in range(numXbins):
            for ybin in range(numYbins):
                val = clone.GetBinContent(xbin,ybin)
                if val==0.0:
                    heatmap_post_process_log.write('The original value in bin (%i, %i) is %f' % (xbin, ybin, val)+' \n')
                    total = 0
                    sides = 0
                    errsq = 0
                    if xbin-1 >=0:
                        heatmap_post_process_log.write('The value in bin (%i, %i) is %f' % (xbin-1, ybin, clone.GetBinContent(xbin-1,ybin))+' \n')
                        total += clone.GetBinContent(xbin-1,ybin)
                        errsq += clone.GetBinError(xbin-1,ybin) ** 2
                        sides += 1
                    if xbin+1 < numXbins:
                        heatmap_post_process_log.write('The value in bin (%i, %i) is %f' % (xbin+1, ybin, clone.GetBinContent(xbin+1,ybin))+' \n')
                        total += clone.GetBinContent(xbin+1,ybin)
                        errsq += clone.GetBinError(xbin+1,ybin) ** 2
                        sides += 1
                    if ybin-1 >=0:
                        heatmap_post_process_log.write('The value in bin (%i, %i) is %f' % (xbin, ybin-1, clone.GetBinContent(xbin,ybin-1))+' \n')
                        total += clone.GetBinContent(xbin,ybin-1)
                        errsq += clone.GetBinError(xbin,ybin-1) ** 2
                        sides += 1
                    if ybin+1 < numYbins:
                        heatmap_post_process_log.write('The value in bin (%i, %i) is %f' % (xbin, ybin+1, clone.GetBinContent(xbin,ybin+1))+' \n')
                        total += clone.GetBinContent(xbin,ybin+1)
                        errsq += clone.GetBinError(xbin,ybin+1) ** 2
                        sides += 1
                    if sides > 0:
                        newval = total/sides
                        error = math.sqrt(errsq)
                    else:
                        newval = 0
                        error = 0
                    clone.SetBinContent(xbin,ybin,newval)
                    clone.SetBinError(xbin,ybin,error)
                    heatmap_post_process_log.write('The new value in bin (%i,%i) is %f' % (xbin,ybin,newval)+' \n')
                heatmap_post_process_log.write('The value in bin (%f, %f) is %f' % (clone.GetXaxis().GetBinLowEdge(xbin), clone.GetYaxis().GetBinLowEdge(ybin), clone.GetBinContent(xbin,ybin))+' \n')
                heatmap_post_process_log.write('The uncertainty in bin (%f, %f) is %f' % (clone.GetXaxis().GetBinLowEdge(xbin), clone.GetYaxis().GetBinLowEdge(ybin), clone.GetBinError(xbin,ybin))+' \n')
        clone.SetName('Post_Processed_%s' % h)
        f.Write()

    return
#-------------------------------------------------------------------------------
if __name__ == '__main__': main()

# EOF


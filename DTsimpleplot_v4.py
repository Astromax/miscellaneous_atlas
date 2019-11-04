#!/usr/bin/env python 
"""
NAME
      DTsimpleplot_v4.py: the purpose of this script is to make simple plots from an input root file
OPTIONS
     -h, --help: Print manual & exits
     
2018-01-30

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
#import numpy as np
#from rootpy.plotting import Hist2D, Canvas
#from rootpy.tree import Cut

root.gROOT.SetBatch(root.kFALSE)

## Globals                                                                 
time = time.strftime('%m-%d')
baseplotdir = '/gpfs/slac/atlas/fs1/u/mbaugh/RooFit/muon_study_24Jun'
#--------------------------------------------------------------------------------
### Options
def options():
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--filebase', default='/gpfs/slac/atlas/fs1/d/mbaugh',
            help='Base of the path to the relevant file')
    parser.add_argument('-f', '--filepath', default='analysis/mc16a/ZTagAndProbe_zll_v03/merged',
            help='The path to the file under investigation')
    parser.add_argument('-r', '--rootfile', default='zll_PowhegPythia8_merged_processed.root',
            help='The root file with the desired tree')
    parser.add_argument('-p', '--plotdir', default='muon_study_24Jun',
            help='The directory the plots are dumped into')
    parser.add_argument('-s', '--selection',  default='muon_study_selections.txt',
            help="Selection criteria for this plot")
    parser.add_argument('-sn', '--savename', default='mumu_linear',
            help='Name to append to the plots for easier identification')
    parser.add_argument('-o', '--output',  default='out.txt',
            help="Name of output file.")
    return parser.parse_args()

#--------------------------------------------------------------------------------
### Main
#--------------------------------------------------------------------------------
def main():
    ops = options()

    print 'Beginning program'

    filebase    = ops.filebase
    filepath    = ops.filepath
    rootfile    = ops.rootfile
    plotdir     = ops.plotdir
    selections  = ops.selection
    savename    = ops.savename

    numsels = list()

    kw1 = 'zmass'

    selectionset = open(selections, 'r')
    for line in selectionset:
        selection = ast.literal_eval(line)
        if selection['Category'] == 'Numerator':
            numsel = selection['Criteria']
            numsels.append(numsel)
        elif selection['Category'] == 'Denominator':
            denomsel = selection['Criteria']

    assert filepath

    kw0 = '_NoSys'
    trees = list()
    #Access the root file, get the tree, and get the branches
    fb  = os.path.join(filebase, filepath)
    tgt = os.path.join(fb, rootfile)
    f = root.TFile(tgt)
    keys = f.GetListOfKeys()
    for key in keys:
        if kw0 in key.GetName():
            trees.append(key.GetName())

    for t in trees:
        print 'We have this tree: %s' % t

    for t in trees:
        mytree = f.Get(t)
        mybranches = mytree.GetListOfBranches()
        mt_words = t.split("_")
        sample = mt_words[0]
        svn = '%s_%s' % (sample, savename)
        for branch in mybranches:
            name = branch.GetName()
            if name==kw1:
                print 'This branch is named %s' % name
                #var = name.split("_")[1]
                var = name
                plot_prod(mytree, name, '%s' % var, 'M_{Z} [GeV]', 100, 0, 5000, denomsel, numsels, svn)
            
    print 'The loop has finished!'


#------------------------------------------------------------------------------
### Plot producer
#------------------------------------------------------------------------------
def plot_prod(mytree, vname, vtitle, vxtitle, vbins, vlow, vhigh, denomsel, numsels, savename ):

    print 'Drawing variable name %s' % vname

    #1st section: initialize temp hists according to parameters from the input
    var_denom = root.TH1D('var_denom', 'var_denom', vbins, vlow, vhigh)
    var_bl_eff = root.TH1D('var_bl_eff', 'var_bl_eff', vbins, vlow, vhigh)

    var_nums = list()
    var_effs = list()

    var_num1 = root.TH1D('var_num1', 'var_num1', vbins, vlow, vhigh)
    var_num2 = root.TH1D('var_num2', 'var_num2', vbins, vlow, vhigh)
    var_num3 = root.TH1D('var_num3', 'var_num3', vbins, vlow, vhigh)
    var_num4 = root.TH1D('var_num4', 'var_num4', vbins, vlow, vhigh)
    var_num5 = root.TH1D('var_num5', 'var_num5', vbins, vlow, vhigh)
    var_num6 = root.TH1D('var_num6', 'var_num6', vbins, vlow, vhigh)
    var_num7 = root.TH1D('var_num7', 'var_num7', vbins, vlow, vhigh)

    var_nums.append(var_num1)
    var_nums.append(var_num2)
    var_nums.append(var_num3)
    var_nums.append(var_num4)
    var_nums.append(var_num5)
    var_nums.append(var_num6)
    var_nums.append(var_num7)

    var_eff1 = root.TH1D('var_eff1', 'var_eff1', vbins, vlow, vhigh)
    var_eff2 = root.TH1D('var_eff2', 'var_eff2', vbins, vlow, vhigh)
    var_eff3 = root.TH1D('var_eff3', 'var_eff3', vbins, vlow, vhigh)
    var_eff4 = root.TH1D('var_eff4', 'var_eff4', vbins, vlow, vhigh)
    var_eff5 = root.TH1D('var_eff5', 'var_eff5', vbins, vlow, vhigh)
    var_eff6 = root.TH1D('var_eff6', 'var_eff6', vbins, vlow, vhigh)
    var_eff7 = root.TH1D('var_eff7', 'var_eff7', vbins, vlow, vhigh)

    var_effs.append(var_eff1)
    var_effs.append(var_eff2)
    var_effs.append(var_eff3)
    var_effs.append(var_eff4)
    var_effs.append(var_eff5)
    var_effs.append(var_eff6)
    var_effs.append(var_eff7)

    mytree.Draw('%s>>var_denom' % vname, denomsel)
    var_bl_eff.Divide(var_denom, var_denom, 1, 1, 'B')
   
    dummy = root.TCanvas()
    colors = [root.kGreen, root.kBlue, root.kRed, root.kOrange, root.kYellow, root.kViolet, root.kCyan]
    legends = list()

    for i, sel in enumerate(numsels):
        mytree.Draw('%s>>%s' % (vname, var_nums[i].GetName()), sel)
        try:
            print 'Average efficiency is %f from num %f and den %f ' % ( (var_nums[i].Integral() / var_denom.Integral()), var_nums[i].Integral(), var_denom.Integral() )
        except ZeroDivisionError:
            print 'The denominator plot is empty!'
        var_effs[i].Divide(var_nums[i], var_denom, 1, 1, 'B')
        var_effs[i].SetLineColor(colors[i])
        subsels = sel.split('&&')
        unique_sel = subsels[-1]
        legends.append('%s Efficiency' % unique_sel)
        
    var_canvas = root.TCanvas()
    var_bl_eff.GetYaxis().SetRangeUser(0.0001, 1.4)
    var_bl_eff.SetLineColor(root.kBlack)
    var_bl_eff.SetTitle(vtitle)
    var_bl_eff.SetXTitle(vxtitle)
    var_bl_eff.SetYTitle('Efficiency')
    var_bl_eff.SetStats(0)
    #var_canvas.SetLogy()
    var_bl_eff.Draw()

    tmp_leg = root.TLegend(0.65, 0.70, 0.9, 0.95)
    tmp_leg.AddEntry(var_bl_eff, 'Baseline Efficiency')

    for j, leg in enumerate(legends):
        print 'The legend entry is %s' % leg
        var_effs[j].Draw('same')
        tmp_leg.AddEntry(var_effs[j], leg)

    tmp_leg.Draw()
    var_canvas.SaveAs('%s/%s_%s_%s.png' % (baseplotdir, vname, time, savename))

    del var_denom
    del var_bl_eff
    del var_canvas
    del var_nums
    del var_effs
    del colors
    del tmp_leg

#------------------------------------------------------------------------------
### Significance 
#------------------------------------------------------------------------------
def calculate_significance(s, b):

    significance  = s/math.sqrt(b)
    print 'The significance has been found!'



#-------------------------------------------------------------------------------
### Bayesian Blocks
#-------------------------------------------------------------------------------
'''
def Bayesian_Blocks(histo):
    t = np.sort(histo)
    N = t.size()
    
    edges = np.concatenate([t[:1], 0.5 * (t[1:] + t[:-1]), t[-1:]])
    block_length = t[-1] - edges

    nn_vec = np.ones(N)
    best = np.zeros(N, dtype=float)
    last = np.zeros(N, dtype=int)

    for k in range(N):
        #Compute the width & count of the final bin for all possible
        #locations of the kth changepoint
        width = block_length[:k + 1] - block_length[k + 1]
        count_vec = np.cumsum(nn_vec[:k + 1][::-1])[::-1]

        #Evaluate fitness function for these possibilities
        fit_vec = count_vec * (np.log(count_vec) - np.log(width))
        fit_vec -=4 #4 comes from the prior on the number of changepoints
        fit_vec[1:] += best[:k]

        #Find the max of the fitness: this is the kth changepoint
        i_max = np.argmax(fit_vec)
        last[k] = i_max
        best[k] = fit_vec[i_max]

    #----------------------------------------------------------------
    ## Recover changepoints be iteratively peeling off the last block
    #----------------------------------------------------------------

    change_points = np.zeros(N, dtype=int)
    i_cp = N
    ind  = N
    while True:
        i_cp -= 1
        change_points[i_cp] = ind
        if ind == 0:
            break
        ind = last[ind - 1]
    change_points = change_points[i_cp:]

    return edges[change_points]
'''





#--------------------------------------------------------------------------------
if __name__ == '__main__': main()

# EOF

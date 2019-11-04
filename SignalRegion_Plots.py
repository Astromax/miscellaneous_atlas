#!/usr/bin/env python 
"""
NAME
      SignalRegion_plots.py the purpose of this script is to make simple plots from an input root file
OPTIONS
     -h, --help: Print manual & exits
     
2019-07-15

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
v = 'v00'
time = time.strftime('%m-%d')
path = 'SRS_Plots_%s_%s' % (time, v)
os.mkdir(path)
SRS_log = open('%s/SRS_Plot_log_%s_%s.txt' % (path,time,v), 'w')
#--------------------------------------------------------------------------------
### Options
def options():
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--targets', default='SRS_targets.txt',
            help='The list of rootfiles to run over')
    parser.add_argument('-s', '--selection',  default='SR_signalMC_selections.txt',
            help="Selection criteria for this plot")
    parser.add_argument('-k', '--kinematics', default='SR_kinematics.txt',
            help='The kinematic variables we want to use')
    return parser.parse_args()

#--------------------------------------------------------------------------------
### Main
#--------------------------------------------------------------------------------
def main():
    ops = options()

    print 'Beginning program'
    SRS_log.write('Commencing Signal Region Selection Plots... \n')

    # Particulars: things which will be specific to what you need
    logmode = False

    selection   = ops.selection

    selections = list()
    selectionset = open(selection, 'r')
    for line in selectionset:
        selection = ast.literal_eval(line)
        selections.append(selection)
        SRS_log.write('Adding selection: %s' % selection+' \n')

    kinematics = list()
    kinematicset = open(ops.kinematics, 'r')
    for line in kinematicset:
        kin = ast.literal_eval(line)
        kinematics.append(kin)
        SRS_log.write('Adding kinematic variable: %s' % kin['Branch']+' \n')

    #Access the root file, get the tree, and get the branches
    targets = list()
    targetset = open(ops.targets, 'r')
    for file in targetset:
        target = file.rstrip()
        targets.append(target)
        SRS_log.write('Adding target root file: %s' % target+' \n')
        targetword = target.split('/')[-1]
        target_ID = targetword.split('_')[2]+targetword.split('_')[3]+targetword.split('_')[4]
        SRS_log.write('The target ID is: %s' % target_ID+' \n')

    print 'Attempting to retrieve the tree' 
    newtargets = combine(targets)

    print 'There are %i new targets' % len(newtargets)

    for t in newtargets:
        trees = list()
        specs = specify(t)
        f = root.TFile(t)
        kw0 = specs['kw0']
        keys = f.GetListOfKeys()
        for key in keys:
            if kw0 in key.GetName():
                tree = f.Get(key.GetName())
                print 'The name of the key is %s' % key.GetName()
                trees.append(tree)

        print 'The name for this target is: %s' % specs['name']
        print 'These are the available trees: '
        print(trees)
        for k in kinematics:
            var = k['Branch']
            bins = k['bin_edges']
            plot_prod(trees, var, bins, selections, specs, logmode, v)
        del trees

    print 'The loop has finished!'
    SRS_log.write('The Circle is now Complete \n')
#------------------------------------------------------------------------------                                                                                                                                                                       
### Plot producer                                                                                                                                                                                                                                   
#------------------------------------------------------------------------------                                                                                                                                                                    
def plot_prod(trees, vxtitle, bins, sels, specs, logmode, v ):

    print 'Drawing variable name %s' % vxtitle
    SRS_log.write('Drawing variable name: %s' % vxtitle+' \n')

    vbins = len(bins) - 2
    vlow  = min(bins)
    vhigh = max(bins)

    #vtitle = 'Number of Lepton Scatters vs %s' % vxtitle

    #1st section: initialize temp hists according to parameters from the input                                                                                                                                                                         
    var_nums = create_TH1Ds('var_num', sels, vbins, vlow, vhigh)

    colors = [root.kRed, root.kBlue, root.kGreen, root.kOrange, root.kYellow, root.kViolet, root.kCyan, root.kMagenta, root.kAzure, root.kTeal]
    legends = list()

    ymin = 0
    ymax = 0

    mass_map = {'102400':299.8, '137400':400.3, '172600':500.2, '208000':599.9, '243800':700.1, '279600':799.7}
    legend_map = {'lifeWeight_base_0.20_reweight_0.05':'0.2 --> 0.05 ns', 'lifeWeight_base_0.20_reweight_0.10':'0.2 --> 0.10 ns', 'lifeWeight_base_1.00_reweight_0.05':'1.0 --> 0.05 ns', 'lifeWeight_base_1.00_reweight_0.10':'1.0 --> 0.10 ns', 'lifeWeight_base_1.00_reweight_0.20':'1.0 --> 0.20 ns', 'lifeWeight_base_1.00_reweight_0.40':'1.0 --> 0.40 ns', 'lifeWeight_base_1.00_reweight_0.60':'1.0 --> 0.60 ns', 'lifeWeight_base_4.00_reweight_0.05':'4.0 --> 0.05 ns', 'lifeWeight_base_4.00_reweight_0.10':'4.0 --> 0.10 ns', 'lifeWeight_base_4.00_reweight_0.20':'4.0 --> 0.20 ns', 'lifeWeight_base_4.00_reweight_0.40':'4.0 --> 0.40 ns', 'lifeWeight_base_4.00_reweight_0.60':'4.0 --> 0.60 ns',  'lifeWeight_base_4.00_reweight_1.00':'4.0 --> 1.0 ns', 'lifeWeight_base_4.00_reweight_2.00':'4.0 --> 2.0 ns', 'lifeWeight_base_4.00_reweight_3.00':'4.0 --> 3.0 ns', 'lifeWeight_base_10.00_reweight_0.05':'10.0 --> 0.05 ns', 'lifeWeight_base_10.00_reweight_0.10':'10.0 --> 0.10 ns', 'lifeWeight_base_10.00_reweight_0.20':'10.0 --> 0.20 ns', 'lifeWeight_base_10.00_reweight_0.40':'10.0 --> 0.40 ns', 'lifeWeight_base_10.00_reweight_0.60':'10.0 --> 0.60 ns',  'lifeWeight_base_10.00_reweight_1.00':'10.0 --> 1.0 ns', 'lifeWeight_base_10.00_reweight_2.00':'10.0 --> 2.0 ns', 'lifeWeight_base_10.00_reweight_3.00':'10.0 --> 3.0 ns', 'lifeWeight_base_10.00_reweight_5.00':'10.0 --> 5.0 ns', ' trackletFlavor==3':'default'}


    print 'Analyzing signal sample: %s' % specs['name']
    mass_ID = specs['name'].split('_')[2]
    lifetime_ID = specs['name'].split('_')[3]
    print 'The lifetime ID is: %s' % lifetime_ID
    lifetime = lifetime_ID[2:]
    print 'The lifetime is: %s' % lifetime
    lfwords = lifetime.split('p')
    life = lfwords[0]+'.'+lfwords[1]+'0'
    mass = mass_map[mass_ID]
    vtitle = 'Chargino Mass: %0.1f GeV' % mass
    print 'The mass of this chargino is: %f GeV' % mass
    print 'The lifetime of this chargino is: %s ns' % life

    lifetime_weights = list()
    branches = trees[0].GetListOfBranches()
    for b in branches:
        if 'lifeWeight' in b.GetName():
            base = b.GetName().split('_')[2]
            rew = b.GetName().split('_')[4]
            if base == life and base != rew:
                newweight = '%s*%s' % (b.GetName(), specs['weight'])
                lifetime_weights.append(newweight)

    for i, sel in enumerate(sels):
        var_trees = create_TH1Ds('var_tree', trees, vbins, vlow, vhigh)
        for k,tree in enumerate(trees):
            tree.Draw('%s>>%s' % (vxtitle, var_trees[k].GetName()), '('+sel+')'+'*%f*%s' %  (specs['lumi'], specs['weight']))
            var_nums[i].Add(var_trees[k])
        var_nums[i].SetLineColor(root.kBlack)
        if var_nums[i].GetMaximum() > ymax:
            ymax = var_nums[i].GetMaximum()
        subsels = sel.split('&&')
        unique_sel = subsels[-1]
        unique_sel = unique_sel.split(')')[0]
        legends.append('%s' % legend_map[unique_sel])
        SRS_log.write('Applying unique selection: %s' % unique_sel+' \n')
        del var_trees

        SRS_log.write('Chargino Mass: %f, Base Lifetime: %s' % (mass, life)+' \n')
        for a in range(vbins+2):
            lowbinedge = var_nums[i].GetBinLowEdge(a)
            baseSubInt = var_nums[i].Integral(a, -1)
            SRS_log.write('Bin: %i' % a+' \n')
            SRS_log.write('Bin Low Edge: %f' % lowbinedge+' \n')
            SRS_log.write('Integral Above Low Edge: %f' % baseSubInt+' \n')
            SRS_log.write('---------------------------------------------------------- \n')
        SRS_log.write('Base Chargino above, Lifetime reweights below \n')

        alters = create_TH1Ds('alter', lifetime_weights, vbins, vlow, vhigh)
        for j, lf in enumerate(lifetime_weights):
            var_trees = create_TH1Ds('var_tree', trees, vbins, vlow, vhigh)
            for k,tree in enumerate(trees):
                tree.Draw('%s>>%s' % (vxtitle, var_trees[k].GetName()), '('+sel+')'+'*%f*%s' %  (specs['lumi'], lf))
                alters[j].Add(var_trees[k])
            alters[j].SetLineColor(colors[j])
            if alters[j].GetMaximum() > ymax:
                ymax = alters[j].GetMaximum()
            if 'lifeWeight' in lf:
                lfID = lf.split('*')[0]
            else:
                lfID = lf
            legends.append('%s' % legend_map[lfID])
            SRS_log.write('Applying lifetime: %s' % legend_map[lfID]+' \n')
            for b in range(vbins+2):
                lowbinedge = alters[j].GetBinLowEdge(b)
                alterSubInt = alters[j].Integral(b, -1)
                SRS_log.write('Mass: %f, Lifetime: %s' % (mass, legend_map[lfID])+' \n')
                SRS_log.write('Bin: %i' % b+' \n')
                SRS_log.write('Bin Low Edge: %f' % lowbinedge+' \n')
                SRS_log.write('Integral Above Low Edge: %f' % alterSubInt+' \n')
                SRS_log.write('--------------------------------------------------------- \n')
            del var_trees
        print 'Completed lifetime reweight loop'

        var_canvas = root.TCanvas()
        var_nums[i].GetYaxis().SetRangeUser(ymin, ymax * 1.25)
        var_nums[i].SetTitle(vtitle)
        var_nums[i].SetXTitle(vxtitle)
        var_nums[i].SetYTitle('Number of Tracklets')
        var_nums[i].SetStats(0)
        if logmode:
            var_canvas.SetLogy()
        var_nums[i].Draw()
        alters[0].Draw('same')

        tmp_leg = root.TLegend(0.65, 0.70, 0.90, 0.95)
        tmp_leg.SetTextSize(0.025) 
        for l, leg in enumerate(legends):
            if l==0:
                tmp_leg.AddEntry(var_nums[l], leg)
                continue
            print 'The legend entry is %s' % leg
            tmp_leg.AddEntry(alters[l-1], leg)
            alters[l-1].SetStats(0)
            alters[l-1].Draw('same')

        tmp_leg.Draw()
        var_canvas.SaveAs('%s/%s_Tracklets_vs_%s_%s_selection_%i_%s.pdf' % (path, specs['name'], vxtitle, time, i, v))

        del alters
        del var_canvas
        del tmp_leg
    del var_nums
    del colors


#------------------------------------------------------------------------------
### Significance 
#------------------------------------------------------------------------------
def calculate_significance(s, b):
    significance  = s/math.sqrt(b)
    print 'The significance has been found!'
    
    return significance
#-------------------------------------------------------------------------------
### Combine C1C1 and C1N1 files
#-------------------------------------------------------------------------------
def combine(roots):
    matches = list()
    combined_roots = list()
    pairs = list()
    for r in roots:
        if r in matches:
            continue
        rword = r.split('/')[-1]
        ctype = rword.split('_')[1]
        r_ID = rword.split('_')[2]+rword.split('_')[3]+rword.split('_')[4]
        for r1 in roots:
            if r1 == r or r1 in matches:
                continue
            r1word = r1.split('/')[-1]
            r1_ID = r1word.split('_')[2]+r1word.split('_')[3]+r1word.split('_')[4]
            if r1_ID == r_ID:
                matches.append(r)
                matches.append(r1)
                pair = (r, r1)
                pairs.append(pair)

    for p in pairs:
        pword = p[0].split('/')[-1]
        p_ID = pword.split('_')[0]+'_'+pword.split('_')[2]+'_'+pword.split('_')[3]+'_'+pword.split('_')[4]+'_'+pword.split('_')[5]
        os.system('hadd -f %s/%s.root %s %s' % (path, p_ID, p[0], p[1]))
        combined_roots.append('%s/%s.root' % (path,p_ID))

    singles = [root for root in roots if root not in matches]

    combined_roots.extend(singles)

    return combined_roots
#-------------------------------------------------------------------------------
### Create list of TH1Ds
#-------------------------------------------------------------------------------
def create_TH1Ds(nametype, sels, vbins, vlow, vhigh):
    plot_list = list()
    for i in range(len(sels)):
        temp = root.TH1D("%s_%i" % (nametype, i), "%s_%i" % (nametype, i), vbins, vlow, vhigh)
        plot_list.append(temp)

    return plot_list
#-------------------------------------------------------------------------------
### Specifications
#-------------------------------------------------------------------------------
def specify(target):
    targetwords = target.split('/')

    if 'data' in targetwords[-1]:
        kw0 = 'data'
        weight = 'leptonVetoWeight'
        lumi = 1.0
        name = 'data'
        vaxis = 'Lepton_Scatters'
    else:
        kw0 = '_NoSys'
        weight = 'eventWeight*genWeight'
        lumi = 35000
        name_initial = targetwords[-1].split('.')[0]
        if 'merged_processed' in name_initial:
            name = name_initial[:-17]
        else:
            name = name_initial
        vaxis = 'Tracklets'

    specifications = {'kw0': kw0, 'weight': weight, 'lumi': lumi, 'name': name, 'vaxis': vaxis}

    return specifications
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

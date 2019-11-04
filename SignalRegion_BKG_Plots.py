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
v = 'v01'
time = time.strftime('%m-%d')
path = 'SRS_BKG_Plots_%s_%s' % (time, v)
os.mkdir(path)
SRS_log = open('%s/SRS_BKG_Plot_log_%s_%s.txt' % (path,time,v), 'w')
TFName = 'Post_Processed_Efficiency_TransferMap_vs_probeEta_and_probePt'
TFMaps = {'electron': 'RateMaps/ZTAP_electron_efficiencies_BW+L_10-11_v13.root', 'muon': 'RateMaps/ZTAP_muon_efficiencies_BW+L_10-11_v12.root'}
CRName = 'probePt_Ratio_Map'
CRMaps = {'electron': 'RateMaps/Closure_vs_ZTAP_electron_comparison_NoTM_11Oct.root', 'muon': 'RateMaps/Closure_vs_ZTAP_muon_comparison_NoTM_11Oct.root'}
#--------------------------------------------------------------------------------
### Options
def options():
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--targets', default='SRS_bkg_targets.txt',
            help='The list of rootfiles to run over')
    parser.add_argument('-s', '--selection',  default='SR_bkg_selections.txt',
            help="Selection criteria for this plot")
    parser.add_argument('-k', '--kinematics', default='SR_bkg_kinematics.txt',
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
        print 'Adding target root file: %s' % target
        SRS_log.write('Adding target root file: %s' % target+' \n')
        targetword = target.split('/')[-1]

    print 'Attempting to retrieve the tree' 
    newtargets = targets

    print 'There are %i new targets' % len(newtargets)

    for t in newtargets:
        print 'This is the new target'
        print(t)
        specs = specify(t)
        SRS_log.write('Analyzing target %s, specs as follows' % t+' \n')
        SRS_log.write('Name: %s, Weight: %s, Lumi: %f' % (specs['name'], specs['weight'], specs['lumi'])+' \n')
        f = root.TFile(target)
        trees = list()
        kw0 = specs['kw0']
        keys = f.GetListOfKeys()
        for key in keys:
            if kw0 in key.GetName():
                tree = f.Get(key.GetName())
                print 'The name of the key is %s' % key.GetName()
                SRS_log.write('Tree Name: %s' % key.GetName()+' \n')
                trees.append(tree)

        print 'The name for this target is: %s' % specs['name']
        for k in kinematics:
            var = k['Branch']
            bins = k['bin_edges']
            plot_prod(f, trees, bins, selections, var, specs, logmode, v)

    print 'The loop has finished!'
    SRS_log.write('The Circle is now Complete \n')
#------------------------------------------------------------------------------                                                                                                                                                                   
### Plot producer                                                                                                                                                                                                                                   
#------------------------------------------------------------------------------                                                                                                                                                                    
def plot_prod(file, trees, bins, sels, var, specs, logmode, v ):

    print 'Drawing variable name %s' % var
    SRS_log.write('Drawing variable name: %s' % var+' \n')

    vbins = len(bins) - 2 
    vlow  = min(bins)
    vhigh = max(bins)

    vtitle = 'Number of %s vs %s' % (specs['vaxis'], var)

    #1st section: initialize temp hists according to parameters from the input                                                                                                                                                                 
    var_nums = create_TH1Ds('var_num', sels, vbins, vlow, vhigh)
    minima   = create_TH1Ds('minima', sels, vbins, vlow, vhigh)
    maxima   = create_TH1Ds('maxima', sels, vbins, vlow, vhigh)

    colors = [root.kRed, root.kBlue, root.kGreen, root.kOrange, root.kYellow, root.kViolet, root.kCyan, root.kMagenta, root.kAzure]
    flavor_map = {' trackletFlavor==1': 'Muon', ' trackletFlavor==2': 'Electron'}
    legends = list()

    ymin = 0
    ymax = 0

    alterweights = list()
    altersmears = list()

    maintree = trees[0]

    branches = maintree.GetListOfBranches()
    for B in branches:
        if 'leptonVetoWeight_' in B.GetName():
            alterweights.append(B.GetName())
            if var == 'trackletPt_smeared':
                alterterm = B.GetName().split('_')[1]
                newsmear = 'trackletPt_smeared'+alterterm
                altersmears.append(newsmear)
                SRS_log.write('New altersmear: %s' % newsmear+' \n')
            else:
                altersmears.append(var)
                SRS_log.write('New altersmear: %s' % var+' \n')

    alters = create_TH1Ds('alter', alterweights, vbins, vlow, vhigh)

    for i, sel in enumerate(sels):
        upper_limits = list()
        lower_limits = list()
        upper_weights = list()
        lower_weights = list()
        nominals = list()
        syst_uncs = list()
        upper_integral = 0
        nom_integral   = 0
        lower_integral = 0        
        alter_integrals = list()
        alter_subintegrals = {}

        SRS_log.write('Beginning of loop iteration %i' % i+' \n')

        maintree.Draw('%s>>%s' % (var, var_nums[i].GetName()), '('+sel+')'+'*%s' % specs['weight'])
        integral = var_nums[i].Integral()
        print 'True Integral: %f' % integral

        stat_uncs = statistical_uncertainty(maintree, sel, specs['weight'])
        #print 'The statistical uncertainty is %f' % stat_unc
        file.cd()

        ptstatmax = max(stat_uncs.keys())

        #Get non-closure systematic
        selwords = sel.split('&&')
        sel_ID = selwords[-1]
        sel_ID = sel_ID.split(')')[0]
        if sel_ID == ' trackletFlavor==1':
            mappath = CRMaps['muon']
        elif sel_ID == ' trackletFlavor==2':
            mappath = CRMaps['electron']

        deviations = {}
        ptmax = 0
        fM = root.TFile(mappath)
        TM = fM.Get(CRName)
        temp0 = TM.Clone()
        ptbins = temp0.GetXaxis().GetNbins()
        for p in range(ptbins):
            pt = temp0.GetBinLowEdge(p)
            print 'In bin %i we have low bin edge %f' % (p, pt)
            SRS_log.write('In bin %i we have low bin edge %f' % (p, pt)+' \n')
            ratio = temp0.GetBinContent(p)
            #deviation = abs(ratio - 1.0)
            deviation = 0.2
            #print 'In bin %i we have deviation %f' % (p, deviation)
            SRS_log.write('In bin %i we have ratio %f and deviation %f' % (p,ratio,deviation)+' \n')
            deviations.update({pt: deviation})
            if pt > ptmax:
                ptmax = pt
        file.cd() #Always remember to return the primary!! 

        for b in range(vbins+2):
            nom = var_nums[i].GetBinContent(b)
            nompt = var_nums[i].GetBinLowEdge(b)
            syst_unc = nom
            if nompt in deviations.keys():
                syst_unc *= deviations[nompt]
            elif nompt > ptmax:
                syst_unc *= deviations[ptmax]
            else:
                syst_unc *= 0.2
            syst_uncs.append(syst_unc)
            nominals.append(nom)
            upper_limits.append(nom)
            lower_limits.append(nom)
            upper_weights.append('nom')
            lower_weights.append('nom')
            SRS_log.write('Bin %i, Bin Low Edge %f, Systematic Uncertainty: %f' % (b, nompt, syst_unc)+' \n')
        
        SRS_log.write('True Integral: %f' % integral+' \n')
        max_alter_int = 0
        min_alter_int = integral
        max_alter_ID = 'nom'
        min_alter_ID = 'nom'
        for j,a in enumerate(alterweights):
            print 'Now using alterweight: %s' % a
            print 'Now using altersmear: %s' % altersmears[j]
            SRS_log.write('Now using alterweight: %s' % a+' \n')
            SRS_log.write('Now using altersmear: %s' % altersmears[j]+' \n')
            maintree.Draw('%s>>%s' % (altersmears[j], alters[j].GetName()), '('+sel+')'+'*%s' % a)
            alter_int = alters[j].Integral()
            alter_integrals.append(alter_int)
            SRS_log.write('Alterweight %s Integral: %f' %(a, alters[j].Integral())+' \n')
            if alter_int > max_alter_int:
                max_alter_int = alter_int
                max_alter_ID = a
            elif alter_int < min_alter_int:
                min_alter_int = alter_int
                min_alter_ID = a
            subintegrals = {}
            for b in range(vbins+2):
                binedge = alters[j].GetBinLowEdge(b)
                binSubInt = alters[j].Integral(b, -1)
                subintegrals[binedge] = binSubInt
                val = alters[j].GetBinContent(b)
                if val > upper_limits[b]:
                    upper_limits[b] = val
                    upper_weights[b] = '%s' % a
                elif val < lower_limits[b]:
                    lower_limits[b] = val
                    lower_weights[b] = '%s' % a
                SRS_log.write('Alterweight %s, Integral from bin edge %f upwards: %f' % (a, binedge, binSubInt)+' \n')
            alter_subintegrals[a] = subintegrals
            #print 'The subintegrals for %s are:' % a
            #print(alter_subintegrals[a])

        for b in range(vbins+2):
            alter_skew_up = upper_limits[b] - nominals[b]
            alter_skew_down = nominals[b] - lower_limits[b]
            full_syst_up = math.sqrt(alter_skew_up ** 2 + syst_uncs[b] ** 2)
            full_syst_down = math.sqrt(alter_skew_down ** 2 + syst_uncs[b] ** 2)
            minlevel = nominals[b] - full_syst_down
            maxlevel = nominals[b] + full_syst_up 
            minima[i].SetBinContent(b, minlevel)
            maxima[i].SetBinContent(b, maxlevel)
            upper_integral += upper_limits[b]
            nom_integral   += nominals[b]
            lower_integral += lower_limits[b]
            lowbinedge = alters[0].GetBinLowEdge(b)
            binner = lowbinedge
            if lowbinedge > ptstatmax:
                binner = ptstatmax
            SRS_log.write('Bin: %i' % b+' \n')
            SRS_log.write('Bin Low Edge: %f' % lowbinedge+' \n')
            SRS_log.write('Nominal Value %f' % nominals[b]+' \n')
            SRS_log.write('Statistical uncertainty: %f' % stat_uncs[binner]+' \n')
            SRS_log.write('Non-Closure systematic: %f' % syst_uncs[b]+' \n')
            SRS_log.write('Alterweight upward systematic: %f' % alter_skew_up+' \n')
            SRS_log.write('Alterweight downward systematic: %f' % alter_skew_down+' \n')
            SRS_log.write('Full systematic upward uncertainty: %f' % full_syst_up+' \n')
            SRS_log.write('Full systematic downward uncertainty: %f' % full_syst_down+' \n')
            SRS_log.write('Full systematic upward bin content: %f' % maxlevel+' \n')
            SRS_log.write('Full systematic downward bin content: %f' % minlevel+' \n')
            SRS_log.write('Upper Limit: %f, Alterweight: %s' % (upper_limits[b], upper_weights[b])+' \n')
            SRS_log.write('Lower Limit: %f, Alterweight: %s' % (lower_limits[b], lower_weights[b])+' \n')
            SRS_log.write('-------------------------------------------------------------------- \n')

        for n in range(vbins+2):
            lowbinedge    = alters[0].GetBinLowEdge(n)
            nomSubInt     = var_nums[i].Integral(n,-1)
            SysUpSubInt   = maxima[i].Integral(n, -1)
            SysDownSubInt = minima[i].Integral(n, -1)
            SRS_log.write('Bin: %i' % n+' \n')
            SRS_log.write('Bin Low Edge: %f' % lowbinedge+' \n')
            SRS_log.write('Nominal Integral above Low Edge: %f' % nomSubInt+' \n')
            SRS_log.write('Systematic Up Integral above Low Edge: %f' % SysUpSubInt+' \n')
            SRS_log.write('Systematic Down Integral above Low Edge: %f' % SysDownSubInt+' \n')
            SRS_log.write('==================================================================== \n')

        var_nums[i].SetLineColor(colors[i])
        minima[i].SetLineColor(colors[i])
        minima[i].SetLineStyle(3)
        maxima[i].SetLineColor(colors[i])
        maxima[i].SetLineStyle(7)
        scale = (integral / nom_integral)
        likely_maximum = scale * upper_integral
        likely_minimum = scale * lower_integral
        if max(upper_limits) > ymax:
            ymax = max(upper_limits)
        subsels = sel.split('&&')
        unique_sel = subsels[-1]
        unique_sel = unique_sel.split(')')[0]
        legends.append('%s' % flavor_map[unique_sel])
        SRS_log.write('Applying unique selection: %s' % unique_sel+' \n')
        SRS_log.write('True Integral: %f' % integral+' \n')
        SRS_log.write('Full Maximum Integral: %f' % maxima[i].Integral()+' \n')
        SRS_log.write('Full Minimum Integral: %f' % minima[i].Integral()+' \n')
        SRS_log.write('Maximum Alterweight: %f' % max(alter_integrals)+' \n')
        SRS_log.write('Minimum Alterweight: %f' % min(alter_integrals)+' \n')
        SRS_log.write('Likely Maximum: %f' % likely_maximum+' \n')
        SRS_log.write('Likely Minimum: %f' % likely_minimum+' \n')
        SRS_log.write('Nominal Integral: %f' % nom_integral+' \n')
        SRS_log.write('Upper Limit Integral: %f' % upper_integral+' \n')
        SRS_log.write('Lower Limit Integral: %f' % lower_integral+' \n')
        SRS_log.write('End of loop iteration %i' % i+' \n')
        SRS_log.write('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ \n')

    var_canvas = root.TCanvas()
    var_nums[0].GetYaxis().SetRangeUser(ymin, ymax * 1.25)
    var_nums[0].SetTitle(vtitle)
    var_nums[0].GetXaxis().SetTitle(var)
    var_nums[0].GetYaxis().SetTitle('Number of %s' % specs['vaxis'])
    var_nums[0].SetStats(0)
    if logmode:
        var_canvas.SetLogy()
    var_nums[0].Draw()
    minima[0].Draw('same')
    maxima[0].Draw('same')

    tmp_leg = root.TLegend(0.65, 0.70, 0.90, 0.95)

    for j, leg in enumerate(legends):
        tmp_leg.AddEntry(var_nums[j], leg)
        if j==0:
            continue
        print 'The legend entry is %s' % leg
        var_nums[j].SetStats(0)
        var_nums[j].Draw('same')
        minima[j].Draw('same')
        maxima[j].Draw('same')

    tmp_leg.Draw()
    var_canvas.SaveAs('%s/%s_%s_vs_%s_%s_%s.pdf' % (path, specs['name'], specs['vaxis'], var, time, v))

    del var_canvas
    del var_nums
    del colors
    del tmp_leg
#------------------------------------------------------------------------------
### Statistical Uncertainty
#------------------------------------------------------------------------------
def statistical_uncertainty(tree, sel, weight):
    selwords = sel.split('&&')
    sel_ID = selwords[-1]
    sel_ID = sel_ID.split(')')[0]
    if sel_ID == ' trackletFlavor==1':
            mappath = TFMaps['muon']
    elif sel_ID == ' trackletFlavor==2':
            mappath = TFMaps['electron']
    else:
        print 'Please include trackletFlavor as the final part of the selection'
        return 0

    total = 0
    fM = root.TFile(mappath)
    TM = fM.Get(TFName)
    temp0 = TM.Clone()
    xbins = temp0.GetXaxis().GetNbins() 
    ybins = temp0.GetYaxis().GetNbins() 
    xmin = temp0.GetXaxis().GetBinLowEdge(1)
    xmax = temp0.GetXaxis().GetBinUpEdge(xbins)
    ymin = temp0.GetYaxis().GetBinLowEdge(1)
    ymax = temp0.GetYaxis().GetBinUpEdge(ybins)

    truetemp = root.TH2D("truetemp", "truetemp", xbins, xmin, xmax, ybins, ymin, ymax)
    tree.Draw('trackletPt_smeared:trackletEta>>truetemp', '('+sel+')'+'*%s' % weight)

    uncertainties = {}

    for j in range(ybins):
        ytotal = 0
        for i in range(xbins):
            TF = temp0.GetBinContent(i,j)
            TF_sig = temp0.GetBinContent(i,j)
            NCR_orig = truetemp.GetBinContent(i,j)
            NCR_sig = math.sqrt(NCR_orig)

            bkg_sig_sqr = ((NCR_orig * TF_sig)**2 + (TF*NCR_sig)**2)
            ytotal += bkg_sig_sqr
            total += bkg_sig_sqr
        print 'Completed x loop'
        y_uncertainty = math.sqrt(ytotal)
        ybinedge = truetemp.GetYaxis().GetBinLowEdge(j) 
        print 'The statistical uncertainty in y-bin %i is: %f' % (j, y_uncertainty)
        uncertainties.update({ybinedge:y_uncertainty})
        SRS_log.write('The uncertainty in bin %i with low edge %f is: %f' % (j, ybinedge, y_uncertainty)+' \n')
    del temp0
    del truetemp

    uncertainty = math.sqrt(total)
    SRS_log.write('The total statistical uncertainty is: %f' % uncertainty+' \n')

    return uncertainties
#------------------------------------------------------------------------------
### Significance 
#------------------------------------------------------------------------------
def calculate_significance(s, b):

    significance  = s/math.sqrt(b)
    print 'The significance has been found!'
    
    return significance
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
### Create list of TH1Ds
#-------------------------------------------------------------------------------
def create_TH1Ds(nametype, sels, vbins, vlow, vhigh):
    plot_list = list()
    for i in range(len(sels)):
        temp = root.TH1D("%s_%i" % (nametype, i), "%s_%i" % (nametype, i), vbins, vlow, vhigh)
        plot_list.append(temp)

    return plot_list

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

#!/usr/bin/env python
"""
NAME
      smear_test.py: the purpose of this script is to test variations of the probePt smearing function
OPTIONS
     -h, --help: Print manual & exits

2019-07-20

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
import numpy as np
import matplotlib.pyplot as plt

root.gROOT.SetBatch(root.kFALSE)

## Globals
# Smearing function parameters                                                                                                                                                                                                                                              
g_par_Mean = -0.450114
g_par_Sigma = 13.5062
g_par_Alpha = 1.68961

#Easy constant
m_pi_2 = math.pi/2

v = 'v0'
time = time.strftime('%m-%d')
smear_test_log = open('ST_log_%s_%s.txt' % (time,v), 'w')
#-------------------------------------------------------------------------------
### Options
#-------------------------------------------------------------------------------
def options():
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--rootfile', default='smear_test_dummy.root',
                        help='The rootfile with the relevant heatmap')
    return parser.parse_args()
#-------------------------------------------------------------------------------
### Main
#-------------------------------------------------------------------------------
def main():
    ops = options()

    print 'Commence primary ignition'
    smear_test_log.write('Smearing function parameter g_par_Mean: %f' % g_par_Mean+' \n')
    smear_test_log.write('Smearing function parameter g_par_Sigma: %f' % g_par_Sigma+' \n')
    smear_test_log.write('Smearing function parameter g_par_Alpha: %f' % g_par_Alpha+' \n')

    print 'Initialize uniform, random set of dummy Pts'
    l = np.random.uniform(0,150,1000)

    g_fDqoverpt = root.TF1("FittedSmearingFunction", crystalBall, -1000, 1000, 4)
    g_fDqoverpt.SetParameter(0, 1.)
    g_fDqoverpt.SetParameter(1, g_par_Mean)
    g_fDqoverpt.SetParameter(2, g_par_Sigma)
    g_fDqoverpt.SetParameter(3, g_par_Alpha)

    canvas = root.TCanvas()
    g_fDqoverpt.Draw()
    canvas.SaveAs('Fitted_Smearing_test.png')

    smear_test_log.write('Now retriving rootfile... \n')
    f = root.TFile(ops.rootfile, 'update')
    tree = f.Get(data)
    new_smears = root.TH1D("new_smears", "new_smears",)
    

    print 'Smear the set of dummy Pts'
    smears = [PtSmear(pt, g_fDqoverpt) for pt in l]

    print 'First value in dummy list: %f' % l[0]
    print 'First value in smeared list: %f' % smears[0]

    plt.scatter( l, smears)
    plt.yscale('log')
    plt.ylim(min(smears), max(smears))
    plt.savefig('plottest5.png')

    return
#-------------------------------------------------------------------------------
### Crystal Ball
#-------------------------------------------------------------------------------
def crystalBall( x, par):
    constant = par[0]
    mean = par[1]
    sigma = par[2]
    alpha = par[3]
       
    #evaluate the crystal ball function
    if (sigma < 0.):
        return 0.
    if (alpha < 0.):
        return 0.
    z = (x[0] - mean)/sigma
    alpha = abs(alpha)
    norm1 = sigma*math.sqrt(2*math.pi)*root.TMath.Erf(alpha/math.sqrt(2))
    norm2 = sigma*np.exp(-alpha*alpha/2)/alpha
    norm3 = norm2
    constant /= (norm1 + norm2 + norm3)
    if (z  < - alpha):
        return constant * np.exp( alpha * (z + 0.5 * alpha))
    elif (z  > + alpha):
        return constant * np.exp(-alpha * (z - 0.5 * alpha))
    else:
        return constant * np.exp(- 0.5 * z * z)
                                                                        
#-------------------------------------------------------------------------------
### Crystal Ball Integral
#-------------------------------------------------------------------------------
def crystalBallIntegral( x, par):
    constant = 1
    mean = par[0]
    sigma = par[1]
    alpha = par[2]

    #evaluate the crystal ball function
    if (sigma < 0.):
        return 0.
    if (alpha < 0.):
        return 0.
    z = (x[0] - mean)/sigma
    alpha = abs(alpha)
    norm1 = sigma*math.sqrt(2*math.pi)*root.TMath.Erf(alpha/math.sqrt(2))
    norm2 = sigma*np.exp(-alpha*alpha/2)/alpha
    norm3 = norm2
    constant /= (norm1 + norm2 + norm3)
    if (z  < - alpha):
        return constant * (+1) * sigma / alpha * np.exp( alpha * (z + 0.5 * alpha))
    elif (z  > + alpha):
        add0 = constant * (+1) * sigma / alpha * np.exp( alpha * (- alpha + 0.5 * alpha))
        sub0 = constant * (-1) * math.sqrt(m_pi_2) * sigma * root.TMath.Erf(alpha / math.sqrt(2))
        add1 = constant * (-1) * math.sqrt(m_pi_2) * sigma * root.TMath.Erf(- alpha / math.sqrt(2))
        sub1 = constant * (-1) * sigma / alpha * np.exp(-alpha * (alpha - 0.5 * alpha))
        return constant * (-1) * sigma / alpha * np.exp(-alpha * (z - 0.5 * alpha)) + add0 + add1 - sub0 - sub1
    else:
        add0 = constant * (+1) * sigma / alpha * np.exp( alpha * (- alpha + 0.5 * alpha))
        sub0 = constant * (-1) * math.sqrt(m_pi_2) * sigma * root.TMath.Erf(alpha / math.sqrt(2))
        return constant * (-1) * math.sqrt(m_pi_2) * sigma * root.TMath.Erf(- z / math.sqrt(2)) + add0 - sub0
        
#-------------------------------------------------------------------------------
### Pt Smearing (in GeV)
#-------------------------------------------------------------------------------
def PtSmear(pt, g_fDqoverpt):
    qoverpt = 1./pt

#    if !g_fDqoverpt:
#        print 'g_fDqoverpt does not exist!'
#        return 0

    smearedQoverpt = qoverpt + g_fDqoverpt.GetRandom()*(0.001) #TeV^-1 -> GeV^-1

    return abs(1./smearedQoverpt)
#-------------------------------------------------------------------------------
if __name__ == '__main__': main()

# EOF


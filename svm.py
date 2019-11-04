from __future__ import division
import matplotlib.pyplot as plt
import math
import ROOT, os, sys
from ROOT import gROOT
from ROOT import TPad
#from sklearn.svm import LinearSVC
from sklearn import svm
import numpy as np
import random
from random import sample 


#User Settings
dorandom = 1
nevents = 10000
quicktest = 1
debug = 1
dologpt = 1
ROOT.gStyle.SetOptStat(0)
#npzfilelocation = "/gpfs/slac/atlas/fs1/d/woodsn/myscripts/training_minions/npzfiles/bkg_100k.npz"
#print npzfilelocation
npzfilelocation = "/gpfs/slac/atlas/fs1/d/woodsn/myscripts/training_minions/npzfiles/bkg_all.npz"
print npzfilelocation

def main():
    #create outputfiles
    os.system('rm -r output')
    os.system('mkdir output')
    myfile = open('output/linearpt_slope_intercept.txt','w')
    logmyfile = open('output/linearlogpt_slope_intercept.txt','w')
    rocmyfile = open('output/linearroc.txt','w')
    roclogmyfile = open('output/linearlogroc.txt','w')

    #load npz file
    npzfile=np.load(npzfilelocation)
    oldclassification=npzfile['classification']
    classification = list()
    for i in oldclassification:
        if i==0:
            i=-1
        classification.append(i)
            
    samples=npzfile['samples']

    #input file processing
    if dorandom:
        indices = np.random.choice(len(samples), nevents, replace=False)
        testsamples = list()
        testclassifications = list()
        for i in indices:
            testsamples.append(samples[i])
            testclassifications.append(classification[i])
    if dologpt:
        log_pt = [math.log(s[0]) for s in testsamples]
        ntrk = [s[1] for s in testsamples]
        testsamples = zip(log_pt, ntrk)


    weight=npzfile['weight']

    pt_array = [s[0] for s in testsamples]
    ntrk_array = [s[1] for s in testsamples]

    gluon_pt_array = [a for a, b in zip(pt_array, testclassifications) if b == -1]
    gluon_ntrk_array = [a for a, b in zip(ntrk_array, testclassifications) if b == -1]

    quark_pt_array = [a for a, b in zip(pt_array, testclassifications) if b == 1]
    quark_ntrk_array = [a for a, b in zip(ntrk_array, testclassifications) if b == 1]

    prefit_plots(pt_array, ntrk_array, gluon_pt_array, gluon_ntrk_array, quark_pt_array, quark_ntrk_array)

    #svm
    #lambdaarray=[1, 10, 40, 60, 80,  100, 200, 400, 600, 800, 1000]
    lambdaarray=[10]
    for lambdax in lambdaarray:
        print 'x is: ' + str(lambdax)
        #clf = svm.SVC(kernel='linear', C=lambdax)
        clf = svm.LinearSVC(C=lambdax)
        clf.fit(testsamples, testclassifications)

        w=clf.coef_[0]
        a = -w[0] / w[1]
        b = -clf.intercept_[0] / w[1]
        print 'slope: ' + str(a)
        print 'intercept: ' + str(b)

        #postfit_plots(pt_array, ntrk_array, classification, gluon_pt_array, gluon_ntrk_array, quark_pt_array, quark_ntrk_array, a, b, lambdax)

        #quark tag efficiency = nQuarks tagged quark/ nQuarks
        pass_quarks_pt = [x for x,y in zip(quark_pt_array, quark_ntrk_array) if y <= a*x + b]
        pass_quarks_ntrk = [y for x,y in zip(quark_pt_array, quark_ntrk_array) if y <= a*x + b]

        fail_gluons_pt = [x for x,y in zip(gluon_pt_array, gluon_ntrk_array) if y > a*x + b]
        fail_gluons_ntrk = [y for x,y in zip(gluon_pt_array, gluon_ntrk_array) if y > a*x + b]

        postfit_plots(pt_array, ntrk_array, testclassifications, fail_gluons_pt, fail_gluons_ntrk, pass_quarks_pt, pass_quarks_ntrk, a, b, lambdax)


        print 'nQuarks: ' + str(len(quark_pt_array))
        print 'nQuarks correctly classified: ' + str(len(pass_quarks_pt))
        q_tag_eff = float(len(pass_quarks_pt)/len(quark_pt_array))
        print 'q_tag_eff: %f' % q_tag_eff


        print 'nGluons: ' + str(len(gluon_pt_array))
        print 'nGluons correctly classified: ' + str(len(fail_gluons_pt))
        g_rejection = float(len(fail_gluons_pt)/len(gluon_pt_array))
        print 'g_rejection: %f' % g_rejection

        fitness = q_tag_eff * g_rejection
        print 'Fitness: %f' % fitness

        myfile.write(str(lambdax) + '\t' + str(a) + '\t' + str(b) + '\n')


def postfit_plots(pt_array, ntrk_array, classification, gluon_pt_array, gluon_ntrk_array, quark_pt_array, quark_ntrk_array, a, inter, lambdax):
    if len(gluon_ntrk_array)==0 or len(quark_ntrk_array)==0 or len(gluon_pt_array)==0 or len(quark_pt_array)==0: return
 
    mxx = np.linspace(0,max(pt_array))
    #myy = a * mxx - (clf.intercept_[0]) / w[1]
    myy = a * mxx + inter

    fig = plt.figure()
    plt.scatter(pt_array, ntrk_array, c=classification)
    plt.title('PostFit All background pt vs ntrk') 
    plt.plot(mxx,myy,'y-', linewidth=3)
    fig.savefig('output/result' + str(lambdax) + 'linear.pdf')

    gluon_ntrk_pt_hist = plt.figure()
    plt.hist2d(gluon_pt_array, gluon_ntrk_array, bins=(50,max(gluon_ntrk_array)))
    plt.title('PostFit Gluon pt vs ntrk') 
    plt.plot(mxx,myy,'y-', linewidth=3)
    gluon_ntrk_pt_hist.savefig('output/svm' + str(lambdax) + 'gluon_ntrk_pt_hist.pdf')

    quark_ntrk_pt_hist = plt.figure()
    plt.hist2d(quark_pt_array, quark_ntrk_array, bins=(50,max(quark_ntrk_array)))
    plt.title('PostFit Quark pt vs ntrk') 
    plt.plot(mxx,myy,'y-', linewidth=3)
    quark_ntrk_pt_hist.savefig('output/svm' + str(lambdax) + 'quark_ntrk_pt_hist.pdf')

def prefit_plots(pt_array, ntrk_array, gluon_pt_array, gluon_ntrk_array, quark_pt_array, quark_ntrk_array):
    #diagnostic plots of quarks and gluons from bkgd+signal
    if len(gluon_ntrk_array)==0 or len(quark_ntrk_array)==0 or len(gluon_pt_array)==0 or len(quark_pt_array)==0: return

    bkg_ntrk_pt = plt.figure()
    plt.scatter(pt_array, ntrk_array)
    plt.title('prefit bkg ntrk vs pt')
    print 'statement 0'
    bkg_ntrk_pt.savefig('output/bkg_ntrk_pt.pdf')
    print 'statement 1'

    ntrk_pt_hist = plt.figure()
    plt.hist2d(pt_array, ntrk_array, bins=(50,max(ntrk_array)))
    plt.title('prefit bkg ntrk vs pt')
    ntrk_pt_hist.savefig('output/ntrk_pt_hist.pdf')

    gluon_ntrk_pt_hist = plt.figure()
    plt.hist2d(gluon_pt_array, gluon_ntrk_array, bins=(50,max(gluon_ntrk_array)))
    plt.title('prefit gluon ntrk vs pt')
    gluon_ntrk_pt_hist.savefig('output/gluon_ntrk_pt_hist.pdf')

    quark_ntrk_pt_hist = plt.figure()
    plt.title('prefit quark ntrk vs pt')
    plt.hist2d(quark_pt_array, quark_ntrk_array, bins=(50,max(quark_ntrk_array)))
    quark_ntrk_pt_hist.savefig('output/quark_ntrk_pt_hist.pdf')

if __name__ == '__main__':
    main()


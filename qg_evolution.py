from __future__ import division
import matplotlib.pyplot as plt
import math
import ROOT, os, sys, time
from ROOT import gROOT
from ROOT import TPad
#from sklearn.svm import LinearSVC
from sklearn import svm
import numpy as np
import random
from random import sample 


#Time
time = time.strftime('%m-%d')

#User Settings
quicktest = 1
debug = 1
dologpt = 1
ROOT.gStyle.SetOptStat(0)
npzfilelocation = "/gpfs/slac/atlas/fs1/d/woodsn/myscripts/training_minions/npzfiles/bkg_all.npz"
print npzfilelocation

#Hyperparameters
doPureEvolution = True
selection_fraction = 0.2
n_generations = 100
breeding_mode = 'Random'
nevents = 10000
lambdax = 10
population_limit = 100
initial_population = 100
slope_mutation_scale = 0.2
intercept_mutation_scale = 1.0
doUltimate = True
v = 'v0'

def main():
    #create outputfiles
    #os.system('rm -r evolution_output')
    #os.system('mkdir evolution_output')
    evolutionlog = open('evolution_output/evolution_output_log_%s_%s.txt' % (time,v), 'w')
    #myfile = open('output/linearpt_slope_intercept.txt','w')
    #logmyfile = open('output/linearlogpt_slope_intercept.txt','w')
    #rocmyfile = open('output/linearroc.txt','w')
    #roclogmyfile = open('output/linearlogroc.txt','w')

    #load npz file
    npzfile=np.load(npzfilelocation)
    oldclassification=npzfile['classification']
    classification = list()
    for i in oldclassification:
        if i==0:
            i=-1
        classification.append(i)
            
    samples=npzfile['samples']
    weight=npzfile['weight']

    evolutionlog.write('Total number of samples: %i' % len(samples) + ' \n')

    # Create Generation Zero 
    generation_zero = list()
    for i in range(initial_population):
        #input file processing
        indices = np.random.choice(len(samples), nevents, replace=False)
        testsamples = list()
        testclassifications = list()
        for j in indices:
            testsamples.append(samples[j])
            testclassifications.append(classification[j])
        if dologpt:
            log_pt = [math.log(s[0]) for s in testsamples]
            ntrk = [s[1] for s in testsamples]
            testsamples = zip(log_pt, ntrk)

        pt_array = [s[0] for s in testsamples]
        ntrk_array = [s[1] for s in testsamples]

        gluon_pt_array = [a for a, b in zip(pt_array, testclassifications) if b == -1]
        gluon_ntrk_array = [a for a, b in zip(ntrk_array, testclassifications) if b == -1]
        quark_pt_array = [a for a, b in zip(pt_array, testclassifications) if b == 1]
        quark_ntrk_array = [a for a, b in zip(ntrk_array, testclassifications) if b == 1]

        if doPureEvolution:
            if i==0: 
                print 'Using pure evolution'
                evolutionlog.write('Creating initial population using pure evolution \n')
            a = (8 - 4) * np.random.random_sample() + 4
            b = -((24 - 12) * np.random.random_sample() + 12)
        if not doPureEvolution:
            #svm
            if i==0: 
                print 'Using SVM to create initial population'
                evolutionlog.write('Creating initial population using SVM, lambda value %f' % lambdax + ' \n')
            clf = svm.LinearSVC(C=lambdax)
            clf.fit(testsamples, testclassifications)

            w=clf.coef_[0]
            a = -w[0] / w[1]
            b = -clf.intercept_[0] / w[1]
    
        organism = (a, b)

        #quark tag efficiency = nQuarks tagged quark/ nQuarks
        pass_quarks_pt = [x for x,y in zip(quark_pt_array, quark_ntrk_array) if y <= a*x + b]
        pass_quarks_ntrk = [y for x,y in zip(quark_pt_array, quark_ntrk_array) if y <= a*x + b]

        fail_gluons_pt = [x for x,y in zip(gluon_pt_array, gluon_ntrk_array) if y > a*x + b]
        fail_gluons_ntrk = [y for x,y in zip(gluon_pt_array, gluon_ntrk_array) if y > a*x + b]

        f = fitness(quark_pt_array, pass_quarks_pt, gluon_pt_array, fail_gluons_pt)
        
        scored_organism = (f, organism)
        generation_zero.append(scored_organism)

        del indices
        del testsamples
        del testclassifications

    print 'Generation Zero has been created, population %i' % len(generation_zero)
    evolutionlog.write('Generation Zero has been created with population: %f' % len(generation_zero) + ' \n')
    nchildren = int(2 * (1.0/selection_fraction))
    mutation_scales = [slope_mutation_scale, intercept_mutation_scale]
    generation = generation_zero
    elites = list()
    for g in range(n_generations):
        sorted_generation = sorted(generation, reverse=True)
        elites.append(sorted_generation[0][1])
        print 'In generation %i, the highest fitness is %f with slope %f and intercept %f' % (g, sorted_generation[0][0], sorted_generation[0][1][0], sorted_generation[0][1][1])
        evolutionlog.write('In generation %i, the highest fitness is %f with slope %f and intercept %f' % (g, sorted_generation[0][0], sorted_generation[0][1][0], sorted_generation[0][1][1]) + ' \n')
        sorted_generation = sorted_generation[:population_limit]
        selection_number = int(selection_fraction * len(generation))
        scored_parents = sorted_generation[:selection_number]
        parents = [p for f,p in scored_parents]

        nextgeneration = next_generation(parents, breeding_mode, nchildren, mutation_scales)
        generation = list()
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

        pt_array = [s[0] for s in testsamples]
        ntrk_array = [s[1] for s in testsamples]
        
        gluon_pt_array = [a for a, b in zip(pt_array, testclassifications) if b == -1]
        gluon_ntrk_array = [a for a, b in zip(ntrk_array, testclassifications) if b == -1]
        quark_pt_array = [a for a, b in zip(pt_array, testclassifications) if b == 1]
        quark_ntrk_array = [a for a, b in zip(ntrk_array, testclassifications) if b == 1]

        for n in nextgeneration:
            pass_quarks_pt = [x for x,y in zip(quark_pt_array, quark_ntrk_array) if y <= n[0]*x + n[1]]
            pass_quarks_ntrk = [y for x,y in zip(quark_pt_array, quark_ntrk_array) if y <= n[0]*x + n[1]]
            fail_gluons_pt = [x for x,y in zip(gluon_pt_array, gluon_ntrk_array) if y > n[0]*x + n[1]]
            fail_gluons_ntrk = [y for x,y in zip(gluon_pt_array, gluon_ntrk_array) if y > n[0]*x + n[1]]

            f = fitness(quark_pt_array, pass_quarks_pt, gluon_pt_array, fail_gluons_pt)
            scored_organism = (f, n)
            generation.append(scored_organism)
        if g==(n_generations-1):
            final_generation = sorted(generation, reverse=True)
            elites.append(final_generation[0][1])
            print 'The final fitness level is %f, slope of %f and intercept of %f' % (final_generation[0][0], final_generation[0][1][0], final_generation[0][1][1])
            evolutionlog.write('The final fitness level is %f, slope of %f and intercept of %f' % (final_generation[0][0], final_generation[0][1][0], final_generation[0][1][1]) + ' \n')

        del indices
        del testsamples
        del testclassifications

    print 'All generations have been produced'

    if doUltimate:
        print 'Ultimate Mode Activated'
        print 'The number of elites is %i' % len(elites)
        slopes = [e[0] for e in elites]
        intercepts = [e[1] for e in elites]
        #ultimate_slope = sum(slopes)/len(elites)
        #ultimate_intercept = sum(intercepts)/len(elites)

        ultimate_slope = 5.912559
        ultimate_intercept = -17.286117

        print 'The ultimate slope is %f' % ultimate_slope
        print 'The ultimate intercept is %f' % ultimate_intercept

        pt_array = [math.log(s[0]) for s in samples]
        ntrk_array = [s[1] for s in samples]
        
        gluon_pt_array = [a for a,b in zip(pt_array, classification) if b == -1]
        gluon_ntrk_array = [a for a,b in zip(ntrk_array, classification) if b == -1]
        quark_pt_array = [a for a,b in zip(pt_array, classification) if b == 1]
        quark_ntrk_array = [a for a,b in zip(ntrk_array, classification) if b == 1]

        pass_quarks_pt = [x for x,y in zip(quark_pt_array, quark_ntrk_array) if y <= ultimate_slope*x + ultimate_intercept]
        pass_quarks_ntrk = [y for x,y in zip(quark_pt_array, quark_ntrk_array) if y <= ultimate_slope*x + ultimate_intercept]
        fail_gluons_pt = [x for x,y in zip(gluon_pt_array, gluon_ntrk_array) if y > ultimate_slope*x + ultimate_intercept]
        fail_gluons_ntrk = [y for x,y in zip(gluon_pt_array, gluon_ntrk_array) if y > ultimate_slope*x + ultimate_intercept]

        f = fitness(quark_pt_array, pass_quarks_pt, gluon_pt_array, fail_gluons_pt)

        print 'The ultimate q_tag_eff is %f' % float(len(pass_quarks_pt)/len(quark_pt_array))
        print 'The ultimate fitness is %f' % f

        evolutionlog.write('The number of elites is %i' % len(elites) + ' \n')
        evolutionlog.write('The ultimate slope is %f' % ultimate_slope + ' \n')
        evolutionlog.write('The ultimate intercept is %f' % ultimate_intercept + ' \n')
        evolutionlog.write('The q_tag_eff is %f' % float(len(pass_quarks_pt)/len(quark_pt_array)) + ' \n')
        evolutionlog.write('The g_rejection is %f' % float(len(fail_gluons_pt)/len(gluon_pt_array)) + ' \n')
        evolutionlog.write('The ultimate fitness is %f' % f + ' \n')


        #myfile.write(str(lambdax) + '\t' + str(a) + '\t' + str(inter) + '\n')
# Next Generation function-------------------------------------------------
def next_generation(parents, breeding_mode, nchildren, mutation_scales):
    next_generation = list()
    next_generation.append(parents[0])

    pairs = pairings(parents, breeding_mode)
    for pair in pairs:
        children = breeding(pair[0], pair[1], nchildren, mutation_scales)
        for child in children:
            next_generation.append(child)

    return next_generation
# --------------------------------------------------------------------------
# Pairing function----------------------------------------------------------
def pairings(parents, breeding_mode):
    N = len(parents)
    if breeding_mode == 'Egalitarian':
        pairs = list()
        for i in range(N/2):
            pair = (parents[i], parents[N-i-1])
            pairs.append(pair)
    elif breeding_mode == 'Ruthless':
        pairs = list()
        for j in range(0, N-1, 2):
            pair = (parents[j], parents[j+1])
            pairs.append(pair)
    elif breeding_mode == 'Random':
        random.shuffle(parents)
        pairs = zip(*[iter(parents)]*2)

    return pairs
# --------------------------------------------------------------------------
# Fitness function----------------------------------------------------------
def fitness(quark_pt_array, pass_quarks_pt, gluon_pt_array, fail_gluons_pt):
    q_tag_eff = float(len(pass_quarks_pt)/len(quark_pt_array))
    g_rejection = float(len(fail_gluons_pt)/len(gluon_pt_array))
    fitness = q_tag_eff * g_rejection

    return fitness
# --------------------------------------------------------------------------
# Mutation function---------------------------------------------------------
def mutation(organism, mutation_scales):
    new_slope = organism[0] + np.random.normal(0, mutation_scales[0]) 
    new_intercept = organism[1] + np.random.normal(0, mutation_scales[1])
    new_organism = (new_slope, new_intercept)

    return new_organism
# --------------------------------------------------------------------------
# Breeding function---------------------------------------------------------
def breeding(organism1, organism2, nchildren, mutation_scales):
    children = list()
    average_slope = float((organism1[0] + organism2[0])/2)
    average_intercept = float((organism1[1] + organism2[1])/2)
    base = (average_slope, average_intercept)
    children.append(base)

    for i in range(nchildren-1):
        child = mutation(base, mutation_scales)
        children.append(child)

    return children
# --------------------------------------------------------------------------
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


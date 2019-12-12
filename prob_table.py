#! /usr/bin/python
import os
import random
from math import factorial
from math import exp
def prob_table():
    N_M = [] #Expected number based on the particular model
    N_D = [] #Actual number of detections

    for i in range(23, 51):
        N_M.append(i)
        N_D.append(i)

    model_count = len(N_M)
    detect_count = len(N_D)
    Probs = []

    for i in xrange(model_count):
        if i == 0:
            Probs.append([0])
        else:
            Probs.append([N_M[i-1]])
        for j in xrange(detect_count): 
            if i == 0:
                Probs[i].append(N_D[j])
            else:
                P = N_M[i-1] ** N_D[j] * exp(-N_M[i-1])/float(factorial(N_D[j]))
                Probs[i].append(P)

    for i in Probs:
        for j in range(0,len(i)):
            print str(i[j]) + " "
        print "\n"

prob_table()
    ###Poisson distribution: P(N_D, N_M) = N_M^N_D * exp{-N_M}/N_D!

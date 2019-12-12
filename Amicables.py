#! /usr/bin/python
import math
import os
def Amicables():
    N = 10000 #The upper limit under consideration
    primes = [] #List of prime numbers up to N
    amicables = [] #will fill with the amicable numbers
    perfects = [] #fill with the perfect numbers

   ###The idea here is to create a list of the sums of the proper divisors for the given coordinate within the array
    divisorsum = [0]

    for j in range(1,N):
        divisorsum.append(0)
        divisors_j = [1]
        cut  = math.sqrt(j) + 1
        for k in range(1,int(cut)):
            if (j % k == 0) and (k != 1):
                y = j/k
                divisors_j.append(y)
                divisors_j.append(k)
        divisorsum[j] = sum(divisors_j)
        


    ###Now we want to look at the divisorsum list & find all pairs s.t. the indice of one is the value of the other & vice-versa
    for d in range(1,N):
        pair = divisorsum[d]
        if pair <= len(divisorsum):
            if (divisorsum[pair] == d) and (pair not in amicables) and (pair != d):
                amicables.append(d)
                amicables.append(pair)
        if (d == pair) and (d != 1):
            perfects.append(d)
        if (pair == 1) and (d != 1):
            primes.append(d)

    for f in range(1,N):
        #if f in primes:
          #  print "The number %i is prime" % f
        if f in amicables:
            print "The number %i is half of an amicable pair" % f
            f +=1
        if f in perfects:
            print "The number %i is a perfect number" %f

    amicable_total = sum(amicables)
    print "The sum of amicable pairs under %i is %i" % (N, amicable_total)

Amicables()


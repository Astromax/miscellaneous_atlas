#! /usr/bin/python
import math
import os
def sum_of_primes():
    N = 2000000 #The upper limit under consideration
    primes = []

    divisorsum = [0]

    for j in range(1,N):
        divisorsum.append(0)
        divisors_j = [1]
        cut = math.sqrt(j) + 1
        for k in range(1,int(cut)):
            if (j % k == 0) and (k != 1):
                y = j/k
                divisors_j.append(y)
                divisors_j.append(k)
        divisorsum[j] = sum(divisors_j)
        if (divisorsum[j] == 1) and (j != 1):
            primes.append(j)

    prime_sum = sum(primes)
    print "The sum of primes below %i is %i" % (N, prime_sum)

sum_of_primes()

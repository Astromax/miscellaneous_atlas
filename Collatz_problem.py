#! /usr/bin/python
import os
import math
def Collatz_problem():
    N = 1000000

    k = math.floor(math.log(N,2)) #find the largest integer power of 2 in the set
    Collatz_max = k+1  #The longest chain will be at least this long

    troublemaker = 2**k #The number with the longest chain so far

    sixth = math.ceil(N/6) #smallest integer above 1/6 of the start value
    if sixth % 2 == 0:
        sixth_odd = sixth + 1
    else:
        sixth_odd = sixth

    half = math.ceil(N/2) #smallest integer above the halfway point

    special_even = 3*sixth_odd + 1 #This gives the smallest large even that is linked to a small #

    Evens = [] #will fill with the large evens that can't have the longest chains

    for x in range(int(sixth)):
        num = special_even + 6*x
        Evens.append(num)

    special_evens = set(Evens)
    
    b = range(int(half), N) #initially the top half of the full set (bottom half is provably weaker)
    Bases = set(b)
    
    Bases.difference_update(special_evens) #remove the special evens, they are provably weaker

    for i in Bases:
        root = i
        count = 1
        rem = [] #list of numbers to be removed from the base sequence
        while i != 1:
            if i % 2 == 0:
                i *= 0.5
                if i in Bases: #if this number is ahead in base list, add to remove list
                    rem.append(i) 
            else:
                i = 3*i + 1
                if i in Bases: #same as above
                    rem.append(i) 
            count += 1
        removals = set(rem)
        #print "There are %i elements in removals" % len(removals)
        #Bases = Bases.difference_update(removals)
        if count > Collatz_max:
            Collatz_max = count
            troublemaker = root
        

    print "The number below %i with the longest Collatz sequence is %i with length %i" % (N,troublemaker, Collatz_max)


Collatz_problem()

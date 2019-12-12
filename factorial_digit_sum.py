#! /usr/bin/python
import os
import math
def factorial_digit_sum():
    N = 99 #The factorial whose digits we want
    m = 1  #The multiplier 
    base = 10 #Which base representation we are using
    digits = [1] #The digits of the target number in reverse order


    for x in range(N):
        if m < base:
            digits = [m * z for z in digits]
            size = len(digits)

        if (m >= base) and (m < base**2):
            a = m % base                     #This gets the "ones place" multiplier
            b = math.floor(m/base)           #This gets the "tens place" multiplier
            dig_a = [a * z for z in digits]  
            dig_a.append(0)                  #Add a 0 to the end of this list
            digits.insert(0,0)               #Adds a 0 at the start of the list
            dig_b = [b * z for z in digits]
            digits = [sum(q) for q in zip(dig_a, dig_b)] #Add element by element
            size = len(digits)

        for y in range(size):
            if (digits[y] > (base - 1)) and (y != size - 1):
                carry = math.floor(digits[y]/base)
                digits[y] = digits[y] % base
                digits[y+1] += carry
            if (digits[y] > (base - 1)) and (y == size - 1):
                carry = math.floor(digits[y]/base)
                digits[y] = digits[y] % base
                digits.append(carry)

        m += 1
        

    digitsum = sum(digits)
    print "The sum of the base %i digits of %i factorial is %i" % (base,N,digitsum)

factorial_digit_sum()
###Somewhere in 20<N<30 range this thing borks up.  No idea why, topically
###it seems identical to long_digit_sum.py & 2*1000 >> 100!
###The problem is that the jumps in size start to cover several digits
###For 10 =< m < 100, find a way to "left append" a zero & then do the ordinary stuff
###Ex. m = 11: make a new list w/ an extra 0 on the left, then add the original to it
###and apply the carry function.
###For m = 123: 2 new lists, one w/ 2 extra left 0s, the other is 2x & has 1 extra 0,
###then combined with 3x the first list
###I think I've got the idea for generalization, but I'm not sure how to efficiently 
###implement it

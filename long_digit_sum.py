#! /usr/bin/python
import os
import math
def long_digit_sum():
    N = 1000 #The exponent we are using
    b = 2  #The base 
    digits = [1] #The digits of the target number in reverse order


    for x in range(N):
        digits = [b * z for z in digits]
        size = len(digits)
        for y in range(size):
            if (digits[y] > 9) and (y != size - 1):
                carry = math.floor(digits[y]/10)
                digits[y] = digits[y] % 10
                digits[y+1] += carry
            if (digits[y] > 9) and (y == size - 1):
                carry = math.floor(digits[y]/10)
                digits[y] = digits[y] % 10
                digits.append(carry)

    digitsum = sum(digits)
    print "The sum of the base-10 digits of %i raised to the %i is %i" % (b,N,digitsum)

long_digit_sum()

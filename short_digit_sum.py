#! /usr/bin/python
import os
import math
def short_digit_sum():
    N = 10 #The number in the exponent
    B = math.factorial(100) #The number whose digits we are summing
    
    total = 0
    while B > 0:
        total += B % 10
        B = math.floor(B/10)

    print "The sum of the digits is %i" %  total

#it seems to have trouble with the big target number, I'm guessing because the 
#computer can't hold 1000 digits in memory so it's rolling over
#some sort of modular math trick should solve this fairly easily
#For sure, modular math is the way to go
#The last digit of 2^N: 2,4,8,6,2 for N = 1,2,3,4,5 respectively.  
#There's a group-like structure, the last digit is easily determined by
#finding Y = N mod 4, N=1,2,3,4 --> Y = 1,2,3,0 --> LD = 2,4,8,6
#So, 2^1000 has N = 1000 --> Y = 0 --> LD = 6
#What about the shift part? Shift amounts to divide-by-ten & floor
#This is equivalent to N -= 3 & floor((5/6)*B').  There should be a cool way
#to do this, but what is it?  
#Last digit cycles with periodicity 4, last pair of digits cycles with
#periodicity of 20 starting with 04 (2^22 ends with 04 but 2^21 doesn't end with 02)
#1000 % 20 = 0 --> last pair should be 76
#There must be a general prescription, otherwise this method is too tedious to work
#2^1000 has ~300 digits in it
#At best the computer is 64-bit --> 16 chunks, but how to do the breaking???
#Montegomery reduction seems like the way to go

short_digit_sum()

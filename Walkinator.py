#!/usr/bin/env python
import os
import math
def Walkinator():

    July_steps = [10156, 10348, 10467, 8489, 9889, 13185, 13135, 12643, 10340, 17668, 9638, 14728, 14015, 10848, 11753, 16085, 11503, 10078, 12510, 22763, 13379, 12406, 11012, 14061, 11554]

    August_steps = [6626, 19069, 17266, 13865, 13605, 12403, 6559, 13201, 12696, 25168, 13923, 10882, 10945, 10077, 12352, 11433, 8322, 14419, 10570, 13306, 13549, 9469, 16588, 15196, 12712, 16934, 12565, 18195, 22298, 13010, 14447]

    August_miles = [3.3, 9.5, 8.7, 6.9, 6.9, 6.2, 3.3, 6.7, 6.5, 12.8, 7.0, 5.3, 5.8, 5.2, 6.5, 5.8, 4.2, 7.1, 5.4, 6.6, 6.7, 4.8, 8.6, 7.7, 6.3, 8.6, 6.3, 9.0, 11.1, 6.7, 7.4]

    August_Cal = [298, 829, 749, 605, 596, 535, 289, 580, 551, 1099, 603, 470, 515, 439, 551, 497, 363, 623, 464, 568, 593, 409, 741, 670, 546, 766, 554, 810, 988, 562, 619]

    stuff = August_steps

    L = len(stuff)
    total = sum(stuff)
    mean = total/L
    print "The total number of steps is " +repr(total)
    print "The mean number of steps is " +repr(mean)
    
    total_miles = sum(August_miles)
    mean_miles = total_miles/L
    print "The total number of miles is " +repr(total_miles)
    print "The mean number of miles is " +repr(mean_miles)
    

    total_cal = sum(August_Cal)
    mean_cal = total_cal/L
    print "The total number of Calories is " +repr(total_cal)
    print "The mean number of Calories is " +repr(mean_cal)
Walkinator()

#How to modify this: the stats for each month are lists as seen above, set some variable "stuff  = month", and have it look at each list that has "stuff" in it's name and calculate the relevant quantities, then print out "The total number of steps in <month> is " +repr(total_steps)




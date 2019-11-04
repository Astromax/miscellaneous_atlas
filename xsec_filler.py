#!/usr/bin/env python 
"""
NAME
    xsec_filler.py: reads in csv file & updates the relevant list
OPTIONS
    
DESCRIPTION

"""
#------------------------------------------------------------------------------                                                                                                                                  
### Imports                                                                                                                                                                                          
#------------------------------------------------------------------------------                                                                                                                        
## std                                                                                                                                                                                                  
import os, argparse, sys, glob, time

#special                                                                                                                                                                                                        
import ROOT as root
import math
import csv
import ast

root.gROOT.SetBatch(root.kFALSE)

## Globals
time = time.strftime('%m-%d')
#--------------------------------------------------------------------------------                                                                                                            
### Options                                                                                                                                                                                  
#--------------------------------------------------------------------------------             
def options():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--csv', default='MC16a_Chargino_mAMSB_Samples.csv',
            help='csv file with the list of xsec')
    parser.add_argument('-s', '--signals', default='dummy_signal.txt',
            help='List of signals which need xsec info')
    return parser.parse_args()

#------------------------------------------------------------------------------
### Main
#------------------------------------------------------------------------------
def main():
    ops = options()
    
    threshold = 0.001
    ##First section: scroll through the existing file to the end 
  
    existing_sigs = list()
    new_sigs = list()

    existing_signals = open(ops.signals, 'r+')
    for line in existing_signals:
        line = line.rstrip()
        exDSID = line.split(' ')[0]
        existing_sigs.append(exDSID)

    xsec_index = -1
    fileff_index = -1
    mcterm_index = -1
    with open(ops.csv, 'r') as csvfile:
        csv_reader = csv.reader(csvfile, delimiter=',')
        line_count = 0    
        for row in csv_reader:
            if line_count == 0:
                print 'The header row is: '
                print(row)
                line_count += 1
                for i,r  in enumerate(row):
                    if r == 'Cross section [pb]':
                        xsec_index = i
                    elif r == 'Filter efficiency':
                        fileff_index = i
                    elif r == 'JobOptions':
                        mcterm_index = i
            else:
                if float(row[xsec_index]) > threshold:
                    DSID = row[mcterm_index].split('.')[1]
                    if DSID not in existing_sigs:
                        existing_signals.write('%s      %i      %s       %0.1f       %s    %0.1f' % (DSID, 0, row[xsec_index], 1.0, row[fileff_index], 1.0) + '\n')
                        new_sigs.append(DSID)

    for ds in new_sigs:
        print 'Use rucio to find sample with DSID: %s' % ds

    return
#--------------------------------------------------------------------------------                                                                                                                               
if __name__ == '__main__': main()

# EOF

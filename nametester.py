#!/usr/bin/env python
#--------------------------------------------------------------------------------------                                                                                                                                                           
### Imports                                                                                                                                                                                                                                 
#--------------------------------------------------------------------------------------                                                                                                                                                     
## std                                                                                                                                                                                                                                
import os, argparse, sys

#--------------------------------------------------------------------------------                                                                                                                                                
### Options                                                                                                                                                                                                                        
#--------------------------------------------------------------------------------                                                                                                                                          
def options():
    parser = argparse.ArgumentParser()
    parser.add_argument('-sg', '--signals', default='decoy_signals.txt',
            help='List of signal files')
    return parser.parse_args()

#--------------------------------------------------------------------------------
### Main
#--------------------------------------------------------------------------------
def main():

    ops = options()
    uniquefiles = list()
    filepairs = list()
    truepairs = list()
    truepaths = list()
    rawfiles = open(ops.signals, 'r')

    for i, file in enumerate(rawfiles):
        file = file.rstrip()
        if file not in uniquefiles and os.path.isfile(file):
            uniquefiles.append(file)
            bname = file.split('/')[-1]
            bname = bname[:-5] #drop the .root portion
            filepair = (file, bname)
            filepairs.append(filepair)

    print 'The number of filepairs is: %i' % len(filepairs)
    for fp in filepairs:
        if fp[0] in truepaths:
            continue
        path = fp[0]
        bn = fp[1]
        print 'This is the path: %s' % path
        print 'This is the basename: %s' % bn
        dups = [f for f in filepairs if f[1] == bn and f[0] != path]
        print 'There are %i duplicates' % len(dups)
        if len(dups) > 0:
            initpathwords = path.split('/')
            for i,d in enumerate(dups):
                print 'This is the first duplicate path: %s' % d[0]
                pathwords = d[0].split('/')
                for j,pw in enumerate(pathwords):
                    if pw == initpathwords[j]:
                        continue
                    else:
                        distinct = pw
                        distinct_index = j
                        continue
                new_bname = distinct+'_'+d[1]
                newpair = (d[0], new_bname)
                truepairs.append(newpair)
            distinct = initpathwords[distinct_index]
            new_fname = distinct+'_'+bn
            new_pair = (path, new_fname)
            truepairs.append(new_pair)
        else:
            truepairs.append(fp)
        truepaths = [t[0] for t in truepairs]
                
    print 'The number of truepairs is: %i' % len(truepairs)
    for i, tp in enumerate(truepairs):
        print 'This is the %ith truepair' % i
        print(tp)

    return
#--------------------------------------------------------------------------------                                                                                                                                     
if __name__ == '__main__': main()

# EOF
    

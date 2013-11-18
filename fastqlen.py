#! /usr/bin/env python

import sys
from collections import Counter
from Bio.SeqIO.QualityIO import FastqGeneralIterator

def calc_lengths(fh):
    lengths = [ len(seq) for title, seq, qual in FastqGeneralIterator(fh) ]
    counts = Counter()
    for l in lengths:
        counts[l] += 1
    return counts

def format_counts(counts):
    outlines = [ '{}\t{}'.format( key, counts[key]) 
                     for key in sorted(counts.keys()) ]
        
    sys.stdout.write('length\tcount\n')
    sys.stdout.write('\n'.join(outlines))
    sys.stdout.write('\n')
    return

def main(args):
    if args:
        for fname in args:
            print fname
            with open(fname, 'r') as fh:
                format_counts(calc_lengths(fh))
    else:
        with sys.stdin as fh:
            format_counts(calc_lengths(fh))
    
if __name__ == '__main__':
    main(sys.argv[1:])

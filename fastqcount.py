#! /usr/bin/env python

import sys
from collections import Counter
from Bio.SeqIO.QualityIO import FastqGeneralIterator

def count_fastq(fh):
    counts = Counter()
    for n, data in enumerate(FastqGeneralIterator(fh)):
        seq = data[2]
        counts[seq] += 1
    unique = len(counts)
    total  = n + 1     # enumerate is base-0!
    return total, unique

def main(args):
    if args:
        sys.stdout.write('file\ttotal\tunique\t%unique\n')
        for fname in args:
            with open(fname, 'r') as fh:
                total, unique = count_fastq(fh)
                line = '{}\t{}\t{}\t{}\n'.format(
                         fname, total,
                         unique, round(unique / total * 100, 1) )
                sys.stdout.write(line)

    else:
        sys.stdout.write('file\ttotal\tunique\t%unique\n')
        with sys.stdin as fh:
            total, unique = count_fastq(fh)
            line = '{}\t{}\t{}\n'.format(
                total,
                unique, round(unique / total * 100, 1) )
            sys.stdout.write(line)
    
if __name__ == '__main__':
    main(sys.argv[1:])

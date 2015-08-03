#!/usr/bin/python
# by Aaron (ams@bio.aau.dk)

from Bio import SeqIO
from toolbox import biofileformat
import matplotlib.pyplot as plt

def get_sizes(fname, sub=False):
    """Parses a sequence file and returns a list of sequence lengths"""

    with biofileformat.FileType('rb')(fname) as fh:
        seqformat = biofileformat.from_handle(fh)

        okformats = [ "fasta", "fastq" ]
        if seqformat not in okformats:
            print "takes only fasta/fastq w/wo compression"
            return

        if sub:
            sizes = []
            for n, rec in enumerate(SeqIO.parse(fh, seqformat)):
                sizes.append(len(rec))
                if n == (sub - 1):
                    break
        else:
            sizes = [len(rec) for rec in SeqIO.parse(fh, seqformat)]

    return sizes

def plot_length(fname, sub=False):
    """
    Parses a sequence file and returns a plot of sequence lengths.
    Optional argument to subset the file.
    """
    sizes = get_sizes(fname, sub)

    plt.hist(sizes)
    plt.title("%s\n%i sequences, range %i to %i bp" \
                % (fname, len(sizes), min(sizes),max(sizes)))
    plt.xlabel("Sequence length (bp)")
    plt.ylabel("Count")
    plt.show()

    return

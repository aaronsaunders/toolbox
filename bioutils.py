#!/usr/bin/python
# by Aaron (ams@bio.aau.dk)

DNA_complements = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
                   'a': 't', 't': 'a', 'c': 'g', 'g': 'c',}

RNA_complements = {'A': 'U', 'U': 'A', 'C': 'G', 'G': 'C',
                   'a': 'u', 'u': 'a', 'c': 'g', 'g': 'c',}

DNA_codon_table = {
    'TTT': 'F', 'TCT': 'S', 'TAT': 'Y', 'TGT': 'C',
    'TTC': 'F', 'TCC': 'S', 'TAC': 'Y', 'TGC': 'C',
    'TTA': 'L', 'TCA': 'S', 'TAA': '-', 'TGA': '-',
    'TTG': 'L', 'TCG': 'S', 'TAG': '-', 'TGG': 'W',
    'CTT': 'L', 'CCT': 'P', 'CAT': 'H', 'CGT': 'R',
    'CTC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',
    'CTA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',
    'CTG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',
    'ATT': 'I', 'ACT': 'T', 'AAT': 'N', 'AGT': 'S',
    'ATC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',
    'ATA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',
    'ATG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',
    'GTT': 'V', 'GCT': 'A', 'GAT': 'D', 'GGT': 'G',
    'GTC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',
    'GTA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',
    'GTG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G',
    }

RNA_codon_table = {
    'UUU': 'F', 'UCU': 'S', 'UAU': 'Y', 'UGU': 'C',
    'UUC': 'F', 'UCC': 'S', 'UAC': 'Y', 'UGC': 'C',
    'UUA': 'L', 'UCA': 'S', 'UAA': '-', 'UGA': '-',
    'UUG': 'L', 'UCG': 'S', 'UAG': '-', 'UGG': 'W',
    'CUU': 'L', 'CCU': 'P', 'CAU': 'H', 'CGU': 'R',
    'CUC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',
    'CUA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',
    'CUG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',
    'AUU': 'I', 'ACU': 'T', 'AAU': 'N', 'AGU': 'S',
    'AUC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',
    'AUA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',
    'AUG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',
    'GUU': 'V', 'GCU': 'A', 'GAU': 'D', 'GGU': 'G',
    'GUC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',
    'GUA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',
    'GUG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G',
    }

def regex_from_string(seq):
    '''Makes an IUPAC DNA sequence string into a regular expression'''
    bases = { 'A' : 'A', 'T': 'T', 'G': 'G', 'C': 'C',
	      'R' : '[AG]', 'Y': '[CT]', 'S': '[GC]',
	      'W' : '[AT]', 'K': '[GT]', 'M': '[AC]',
	      'B': '[CGT]', 'D': '[AGT]', 'H': '[ACT]',
	      'V': '[ACG]', 'N': '[AGCT]'}
    seq = seq.upper()

    return ''.join([ bases[character] for character in seq ])

def regex_to_string(seq):
    '''Makes an IUPAC DNA sequence string into a regular expression'''
    bases = {'A': 'A',
             'AC': 'M',
             'ACG': 'V',
             'ACGT': 'N',
             'ACT': 'H',
             'AG': 'R',
             'AGT': 'D',
             'AT': 'W',
             'C': 'C',
             'CG': 'S',
             'CGT': 'B',
             'CT': 'Y',
             'G': 'G',
             'GT': 'K',
             'T': 'T'}
    seq = seq.upper()

    newseq = []
    regex  = []
    collect = False

    for base in list(seq):
        if base == '[':
            collect = True
            continue
        if base == ']':
            collect = False
            newseq.append(bases[''.join(sorted(regex))])
            regex = []
            continue
        if collect == True:
            regex.append(base)
            continue
        newseq .append(base)

    return ''.join(newseq)

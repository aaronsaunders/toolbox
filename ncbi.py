#!/usr/bin/python
# by Aaron (ams@bio.aau.dk)

########## Import functions ###########

import sys 
from Bio import Entrez
from  Bio.Blast  import  NCBIWWW
from  Bio.Blast  import  NCBIXML


"""
A library of functions for efetch scripting
"""


############## Functions   ###############

def  efetch_nt_fasta(acc, start, stop, strand):
""" wrapper for Bio.efetch"""
  
    if strand == 'PLUS':
        strand = 1
    if strand == 'MINUS':
        strand = -1
    handle = Entrez.efetch(db='nuccore', id=acc, seq_start=start,
            seq_stop=stop, strand=strand, retmode='fasta')
    record = SeqIO.read(handle, 'fasta')
    
    return  record

def  efetch_nt_gb(acc, start, stop, strand):
""" wrapper for Bio.efetch"""
    
    if strand == 'PLUS':
        strand = 1
    if strand == 'MINUS':
        strand = -1
    handle = Entrez.efetch(db='nuccore', id=acc, seq_start=start,
            seq_stop=stop, strand=strand, retmode='gbwithparts')
    record = SeqIO.read(handle, 'genbank')  
  
  return  record
  
def ncbi_blast_fastafile(filehandle):
    """ """
    records  =  SeqIO.parse(filehandle,  format="fasta")
    blast_records = []

    for record in records:
        result_handle  =  NCBIWWW.qblast('blastn', 'nr', record.format('fasta'))
        blast_records.append(NCBIXML.read(result_handle))
   

#!/usr/bin/python
# by Aaron (ams@bio.aau.dk)

########## Import functions ###########

import sys 
from Bio import Entrez
from  Bio.Blast  import  NCBIWWW
from  Bio.Blast  import  NCBIXML
import re
from Bio import SeqIO


Entrez.email = 'ams@bio.aau.dk'


"""
A library of functions for efetch scripting
"""


############## Functions   ###############

def parse_cds(record):
      """
      Takes Biopython SeqRecord as input
      returns CDS acc, start, stop, strand
      """

      # parse out coded_by qualifier from biopython SeqRecord
      coded_by = ''
      CDS_count = 0
      for gp_feature in record.features:
          if gp_feature.type == 'CDS':
              coded_by = gp_feature.qualifiers['coded_by']
              CDS_count =+ 1
      if not coded_by:
          logfile.write('Error: CDS not found for %s\n'%record.id)
          cds_not_found.append(record.id)
          return      
      if CDS_count > 1:
          logfile.write('Error: Multiple CDS found for %s\n'%record.id)
          cds_not_found.append(record.id)
          return
      
      # accepts
      # NM_010510.1:21..569 
      # NM_010510.1<21..569
      # NM_010510.1:21..>569  
      # NM_010510.1:<1..>569)
      # complement(NM_010510.1:21..569)
      
      mobj = re.search('^complement\(([^)]+)\)', str(coded_by))
      if mobj:
          cds_strand = 'MINUS'
          coded_by = mobj.group(1)
      else:
          cds_strand = 'PLUS'
      
      pattern = '(\w+\.\d):<?(\d+)\.\.>?(\d+)'
      mobj = re.search(pattern, str(coded_by))
      if mobj:
          cds_acc = mobj.group(1) 
          cds_start = mobj.group(2) 
          cds_stop = mobj.group(3) 

      return cds_acc, cds_start, cds_stop, cds_strand


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
  


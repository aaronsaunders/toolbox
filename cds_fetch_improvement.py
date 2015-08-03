#Script to take a file of proteins in GenBank/GenPept format, examine
#their annotation, and use this to download their CDS from the NCBI.
#Written and tested on Biopython 1.49, on 2008/01/19
# by Peter biopython at maubp.freeserve.co.uk

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio import Entrez
#Edit this next line (and read the NCBI Entrez usage guidelines)
#Entrez.email = "Your.Name.Here at example.com"

def get_nuc_by_name(name, start=None, end=None) :
   """Fetches the sequence from the NCBI given an identifier name.

   Note start and end should be given using one based counting!

   Returns a Seq object."""
   record = SeqIO.read(Entrez.efetch("nucleotide",
                                     id=name.strip(),
                                     seq_start=start,
                                     seq_stop=end,
                                     retmode="text",
                                     rettype="fasta"), "fasta")
   return record.seq

def get_nuc_from_coded_by_string(source) :
   """Fetches the sequence from the NCBI for a "coded_by" string.

   e.g. "NM_010510.1:21..569" or "AF376133.1:<1..>553"
   or "join(AB061020.1:1..184,AB061020.1:300..1300)"
   or "complement(NC_001713.1:67323..68795)"

   Note - joins and complements are handled by recursion.

   Returns a Seq object."""
   
   if source.startswith("complement(") :
       assert source.endswith(")")
       #For simplicity this works by recursion
       return get_nuc_from_coded_by_string(source[11:-1]).reverse_complement()
   
   if source.startswith("join(") :
       assert source.endswith(")")
       #For simplicity this works by recursion.
       #Note that the Seq object (currently) does not have a join
       #method, so convert to strings and join them, then go back
       #to a Seq object:
       return Seq("".join(str(get_nuc_from_coded_by_string(s)) \
                          for s in source[5:-1].split(",")))
   
   if "(" in source or ")" in source \
   or source.count(":") != 1 or source.count("..") != 1 :
       raise ValueError("Don't understand %s" % repr(source))
   name, loc = source.split(":")
   #Remove and ignore any leading < or > for fuzzy locations which
   #indicate the full CDS extends beyond the region sequenced.
   start, end = [int(x.lstrip("<>")) for x in loc.split("..")]
   #We could now download the full sequence, and crop it locally:
   #return get_nuc_by_name(name)[start-1:end]
   #However, we can ask the NCBI to crop it and then download
   #just the bit we need!
   return get_nuc_by_name(name,start,end)

def find_protein_within_nuc(protein_seq, nuc_seq, table) :
   """Search all six frames to find a protein's CDS."""
   for frame in [0,1,2] :
       start = nuc_seq[frame:].translate(table).find(protein_seq)
       if start != -1 :
          return nuc_seq[frame+3*start:frame+3*(start+len(protein_seq))]
   rev_seq = nuc_seq.reverse_complement()
   for frame in [0,1,2] :
       start = rev_seq[frame:].translate(table).find(protein_seq)
       if start != -1 :
          return rev_seq[frame+3*start:frame+3*(start+len(protein_seq))]
   raise ValueError("Could not find the protein sequence "
                    "in any of the six translation frames.")


def get_nuc_record(protein_record, table="Standard") :
   """Given a protein record, returns a record with the CDS nucleotides.

   The protein's annotation is used to determine the CDS sequence(s)
   which are downloaded from the NCBI using Entrez.

   The translation table specified is used to check the nucleotides
   actually do give the expected protein sequence.

   Tries to get the CDS information from a "coded_by" qualifier,
   failing that it falls back on a DB_SOURCE xref entry (which
   does not specify which bit of the nucleotide sequence referenced
   is required - this is deduced from the expected translation).  This
   could get the wrong region if the happens to be two genes with
   different nucleotides encoding the same protein sequence!
   """
   if not isinstance(protein_record, SeqRecord) :
       raise TypeError("Expect a SeqRecord as the protein_record.")
   feature = None
   for f in protein_record.features :
       if f.type == "CDS" and "coded_by" in f.qualifiers :
           feature = f
           break
   if feature :
       #This is the good situation, there is a precise "coded_by" string
       #Check this CDS feature is for the whole protein:
       assert feature.location.start.position == 0
       assert feature.location.end.position == len(protein_record)
       source = feature.qualifiers["coded_by"][0]
       print "Using %s" % source
       return SeqRecord(Seq(""))
       nuc = get_nuc_from_coded_by_string(source)
       #See if this included the stop codon - they don't always!
       if str(nuc[-3:].translate(table)) == "*" :
           nuc = nuc[:-3]
   elif "db_source" in protein_record.annotations :
       #Note the current parsing of the DBSOURCE lines in GenPept
       #files is non-optimal (as of Biopython 1.49).  If the
       #parsing is changed then the following code will need
       #updating to pull out the first xrefs entry.
       parts = protein_record.annotations["db_source"].split()
       source = parts[parts.index("xrefs:")+1].strip(",;")
       print "Using %s" % source
       nuc_all = get_nuc_by_name(source)
       nuc = find_protein_within_nuc(protein_record.seq, nuc_all, table)
   else :
       raise ValueError("Could not determine CDS source from record.")
   assert str(nuc.translate(table)) == str(protein_record.seq), \
          "Translation:\n%s\nExpected:\n%s" \
          % (translate(nuc,table), protein_record.seq)
   return SeqRecord(nuc, id=protein_record.id,
                    description="(the CDS for this protein)")

#Now use the above functions to fetch the CDS sequence for some proteins...
gbk_input = "Diatoms_in.gp" #any proteins in GenBank/GenPept format.
nucs = (get_nuc_record(p, table="Standard") for p \
         in SeqIO.parse(open(gbk_input),"genbank"))

handle = open("nucleotide.fasta","w")
SeqIO.write(nucs, handle, "fasta")
handle.close()
print "Done"
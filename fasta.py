
#!/usr/bin/python
# fasta.py
# ams@bio.aau.dk

## NOTE: Use the biopython simple parsers
## to parse fasta to strings (fast!) use
# from Bio.SeqIO.FastaIO import SimpleFastaParser
# for name, seq in SimpleFastaParser(fh):

## to parse fastq to strings (fast!) use
# from Bio.SeqIO.QualityIO import FastqGeneralIterator
# for title, seq, qual in FastqGeneralIterator(fh)

def parse_rnammer_fasta_header(header):
    """take a header and returns 
    gi, acc, start, stop, strand(PLUS/MINUS), molecule
    """
    main, molecule, score = header.strip().split()
    fields = main.split('|')
    gi = fields[1]
    acc, version = fields[3].split('.')
    position = fields[4][1:]
    coords, strand = position.split('_')
    start, stop = coords.split('-')
    if strand == 'DIR-':
            strand = 'MINUS'
    if strand == 'DIR+':
            strand = 'PLUS'
    title, molecule = molecule[1:].split('=')
    molecule = molecule.replace('_', ' ')

    return gi, acc, start, stop, strand, molecule



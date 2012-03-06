
#!/usr/bin/python
# fasta.py
# ams@bio.aau.dk


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


def read_FASTA_strings(filename):
    with open(filename) as file:
        return file.read().split('>')[1:]
    
def read_FASTA_entries(filename):
    return [seq.partition('\n') for seq in read_FASTA_strings(filename)]

def get_items_from_file(filename, testfn=None):
    """Return all the items in the file named filename; if testfn
    then include only those items for which testfn is true"""
    with open(filename) as file:
        return get_items(file, testfn)

def find_item_in_file(filename, testfn=None):
    """Return the first item in the file named filename; if testfn
    then return the first item for which testfn is true"""
    with open(filename) as file:
        return find_item(file, testfn)

def find_item(src, testfn):
    """Return the first item in src; if testfn then return the first
    item for which testfn is true"""
    gen = item_generator(src)
    item = next(gen)
    if not testfn:
        return item
    else:
        try:
            while not testfn(item):
                item = next(gen)
            return item
        except StopIteration:
            return None

def get_items(src, testfn=None):
    """Return all the items in src; if testfn then include only those
    items for which testfn is true"""
    return [item for item in item_generator(src)
            if not testfn or testfn(item)]

def item_generator(src):
    """Return a generator that produces a FASTA sequence from src
    each time it is called"""
    skip_intro(src)
    seq = ''
    description = src.readline().split('|')
    line = src.readline()
    while line:
        while line and line[0] != '>':
            seq += line
            line = src.readline()
        yield (description, seq)
        seq = ''
        description = line[1:-1].split('|')
        line = src.readline()

def skip_intro(src):
    """Skip introductory text that appears in src before the first
    item"""
    pass          # no introduction in a FASTA file

def test():
    filename = '../data/aa010.fasta'
    for item in get_items_from_file(filename):
        print(item[0], item[1], sep='\n')
    print('\nKlebsiella sequences:\n')
    for item in get_items_from_file(filename,
                                    lambda itm:
                                        -1 < itm[0][-1].find('Klebsiella')):
        print(item[0], item[1], sep='\n')
    print('\nGI 6693805:\n')
    print(find_item_in_file(filename,
                            lambda itm: '6693805' == itm[0][1]))
    

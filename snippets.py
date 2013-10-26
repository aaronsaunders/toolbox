#Python code snippets
# ams@bio.aau.dk

# boiler plate
if __name__ == '__main__':
    main()

# to iterate and index (NOTE: base-zero index)
for i, item in enumerate(iterator)):
    pass

# index starts at 1!!!
for i, item in enumerate(iterator, 1)):
    pass

# take the first string in a delimited sting
first, rest = item.split(';', 1)

seqs = [ seq.strip() for seq in fasta.split('\n')[1:] ]


# To write out a list of lists
with open('out.txt', 'w') as fh:
    fh.write('\n'.join( '\t'.join(rec) for rec in list_of_lists ))


####################################
## Path manipulations
# These functions use the "os" module.
    import os
# To get or change the current working directory
    os.getcwd()
    # returns the full path of the pwd
    os.chdir(path)
# To check if a file or directory exists:
    os.listdir(path) 
    # returns a list of the names in the pwd (not .files)
    os.path.exists(path) # returns True or False
    os.path.isdir(path) # returns True or False
    os.path.isfile(path) # returns True or False
# To work with path and filename:
    os.path.split(path) # returns (head, tail)
    os.path.dirname(path) # returns head
    os.path.basename(path) # returns tail
    os.path.join(path1[, path2[, ...]]) # joins intelligently, eg. // = /
# Get the filenames in a directory:
    filenames = os.listdir(inputdir)
# to create (make) a directory:
    os.mkdir()
    # if exists OSError is raised
   os.mkdir('outdirname')
   cwd = os.getcwd()
   savepath = os.path.join(cwd, 'outdirname')
# To make a unique output file/folder name with time/date stamp
    import datetime
    now = datetime.datetime.now()
    now.strftime("%Y%m%d_%H%M")

# http://docs.python.org/library/datetime.html#strftime-and-strptime-behavior
 



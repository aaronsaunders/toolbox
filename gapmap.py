#!/usr/bin/env python
usage="""
-----------------------------------------------------------------------
Insert gaps into DNA sequences given a corresponding protein alignment.

USAGE:

	gapmap.py AAfile.fta DNAfile.fta


veraion 1.1  5 Mar 09 - added argument parsing to read in file names
version 1.0 17 Jun 08
-----------------------------------------------------------------------

"""

import sys
import re
from sets import Set as set # useful for version 2.3?

# INITIALIZE STUFF
# CHANGE ME!!!
# fileroot='seq'
# 
# inputfilename=fileroot+'AA.fta'
# inputDNA=fileroot+'DNA.fta'
# 

if len(sys.argv)<3:
	print usage	

else:
	# load the amino-acid file
	inputfilename = sys.argv[1]

	# load the DNA file
	inputDNA=sys.argv[2]		

	outname=inputDNA
	outfilename="gap_" + outname

	### Change settings here
	keepall=False # don't chop out conserved sites? just redundant seqs
	printtoscreen=True # false means it will save to the outfile
	savetofile=True



	#other definitinos
	firstseq = False
	seqname= ''
	dataset={}
	nucs={}
	skipnext=False
	nbrf=False

	# First, load in the seqs into a dictionary

	try: 
		inFile = file(inputfilename, 'r') 
		alllines = inFile.readlines() # trying a different approach
	# this is easier for the \r detection, but for very large files,
	# it is probably better to use the original formulation
	# plus, does readlines() work in 2.4??

	except IOError:
		print"Can't find the file %s. Are you in the right directory?" % inputfilename
		print 
		print
		printtoscreen=False # false means it will save to the outfile
		printhtml=False
		savetofile=False
		alllines=[]

	if alllines[0].find('\r')>0:
		print
		print "MacFormat (use unix file when possible)."
		lines=alllines[0].split('\r')
		# print "Found %d lines" % len(lines)	
	
	else:
		lines=alllines
	
	# print len(lines)

	for inline in lines:
		line=inline.strip('\r\n') 
		if skipnext: # needed for NBRF format
			skipnext=False
		else:
			if line and line[0]=='>':
				if line.startswith('>DL;') or line.startswith('>P1;'):
					nbrf=True
					# print 'NBRF format'
					seqname=line.split()[1]  # defaults to space
					skipnext=True
				else:
					seqname=line[1:]
				firstseq=True
			elif firstseq: # we know it's not the first line
				# print  'we have sequence'
				try:
					dataset[seqname] +=  line
				except:
					dataset[seqname]  =  line
				# this should be modified to append

	inFile.close()

	# so dataset[] is the AA sequence with gaps
	# do the same for the Nuc list
	# First, load in the seqs into a dictionary

	try: 
		inFile = file(inputDNA, 'r') 
		alllines = inFile.readlines() # trying a different approach
	# this is easier for the \r detection, but for very large files,
	# it is probably better to use the original formulation
	# plus, does readlines() work in 2.4??

	except IOError:
		print"Can't find the file %s. Are you in the right directory?" % inputfilename
		print 
		print
		printtoscreen=False # false means it will save to the outfile
		printhtml=False
		savetofile=False
		alllines=[]

	if alllines[0].find('\r')>0:
		print
		print "MacFormat (use unix file when possible)."
		lines=alllines[0].split('\r')
		# print "Found %d lines" % len(lines)	
	
	else:
		lines=alllines
	

	for inline in lines:
	# print len(lines)
		line=inline.strip('\r\n') 
		if skipnext: # needed for NBRF format
			skipnext=False
		else:
			if line and line[0]=='>':
				if line.startswith('>DL;') or line.startswith('>P1;'):
					nbrf=True
					# print 'NBRF format'
					seqname=line.split()[1]  # defaults to space
					skipnext=True
				else:
					seqname=line[1:]
				firstseq=True
			elif firstseq: # we know it's not the first line
				# print  'we have sequence'
				try:
					nucs[seqname] +=  line
				except:
					nucs[seqname]  =  line
				# this should be modified to append

	inFile.close()




	if nbrf: print "NBRF format" 

	# print dataset.keys()
	# print dataset.values()
	# print '\n'

	# first get rid of all the invariable rows
	# make a reverse dictionary

	backseqs={}
	backseqtemp={}
	dataname=[]
	# removes identical sequences -- make another way
	# to pupulate catacol_list if you want to keep these
	# ->it would just be dataset.items()

	for key, value in dataset.items():
		newseq=''
		# print "HERE IS THE " + key
		findex=0
		af=value.find('-')
		while af>-1 and (af < len(value)+1):
			# print "Found dash at %d" % af
			gf=3*(af)
			newseq=nucs[key][:gf]+'---'+nucs[key][gf:]
			nucs[key]=newseq	
			af=value.find('-',af+1)  
		# print nucs[key]
	
	######### need to check that we aren't messing things up with all this sorting!

	sortkeys=nucs.keys()
	firstout= "Found %d sequences...\r\r" % (len(sortkeys))
	print firstout
	sortkeys.sort()
	if printtoscreen:
		for pri in sortkeys:  # this happens to be the keys for outdict
			print ">" + pri
			print nucs[pri]
			# if  indexdict[pri] in labelindex.keys():
			# 	print 'yes', str(labelindex[indexdict[pri]])

	if savetofile: # not printtoscreen, i.e., print to file
		outFile = file(outfilename, 'w')
		for pri in sortkeys:
			labelstring='>'+ pri +'\n'
			outFile.write(labelstring)
			outFile.write(nucs[pri])
			outFile.write('\n')
		outFile.close()
		print "Saved to file \'" + outfilename + "\' ...I think"
	
	

	# [len(row) for row in myCols]

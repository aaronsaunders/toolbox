#!/usr/bin/python
# ams@bio.aau.dk

import subprocess

# make an alignment 
# muscle -in strep.mod.seqs.txt -out strep.align.txt

# degap the alignment

# calculate a tree
# FastTree -nt strep.align.txt > strep.tree.txt

# print the tree
# in R
# library(ape)
# tree = read.tree('tree.nwk')
# plot(tree)



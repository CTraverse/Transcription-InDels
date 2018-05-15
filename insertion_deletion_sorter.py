# -*- coding: utf-8 -*-
"""
Created on Wed Aug 24 13:02:47 2016

@author: chuck
"""


from sys import argv
import re

script, infile = argv

#Just a short script to separate insertions and deletions into two separate files
indels = list(open(infile, "r"))

insertions = open(re.sub(".txt", "_insertions.txt", infile), "w")
deletions = open(re.sub(".txt", "_deletions.txt", infile), "w")

for indel in indels:
    
    indel_split = indel.split("\t")
    if indel_split[5].count("I") == 1:
        insertions.write(indel)
    else:
        deletions.write(indel)






















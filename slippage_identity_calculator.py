# -*- coding: utf-8 -*-
"""
Created on Fri May 20 11:25:05 2016

@author: chuck
"""

"""
This script will calculate the %identity after a slippage event
"""

from sys import argv
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

script, infile = argv

deletions = list(open(infile, "r"))

slipped_seq = str()
del_seq = str()

#List of the different types of mispairings that can occur
mispair_list = ["rA_dA", "rA_dC", "rA_dG", 
                "rT_dT", "rT_dC", "rT_dG",
                "rC_dA", "rC_dC", "rC_dT",
                "rG_dA", "rG_dT", "rG_dG"]

#Build a dictionary of dictionaries for each position in the RNA:DNA hybrid
#Each position in the hybrid contains a dictionary that contains the counts of each type of mispair
hybrid_position = {}

#Populate dictionary with mismatches and mismatch types, all set to 0
for i in range(9):
    hybrid_position [i] = {}
    hybrid_position [i] ["matches"] = 0
    hybrid_position [i] ["mismatches"] = 0
    
    for mispair in mispair_list:
        hybrid_position[i][mispair] = 0

#Loop through each deletion and find how many mismatches there are with the original sequence
for deletion in deletions:
    #Process the deletion into the sequence it slipped over and the sequence right before the deletion
    del_split = deletion.split("\t")
    slipped_seq = del_split[3]
    del_seq = del_split[3] + del_split[4]
    del_len = len(del_split[4])    
    
    #Position where the deletion eneed
    slipped_to_seq = del_seq[ len(del_seq) - 9 : len(del_seq) ]
        
      
    #Loop through each base in the deletion and determine:
    #If there as a mismatch, and what type of mismatch it was
    same_count = 0    
    j = 0
    for base_1, base_2 in zip(slipped_seq, slipped_to_seq):
        #If they're identical, count it as a match
        if base_1 == base_2:
            hybrid_position[j]["matches"] += 1
            j += 1
        
        #If they aren't identical, count it as a mismatch and count the type of mismatch
        else:
            hybrid_position[j]["r%s_d%s" % (base_1, str(Seq(base_2, IUPAC.unambiguous_dna).reverse_complement()))] += 1
            hybrid_position[j]["mismatches"] += 1
            j += 1

#Print out the mismatches
for i in range(9):
    print("Position %s of the RNA:DNA hybrid" % (i))
    print("Matches: %s" % (hybrid_position [i]["matches"]))
    print("Mismatches: %s" % (hybrid_position [i]["mismatches"]))
    
    for mispair in mispair_list:
        print("%s\t%s" % (mispair, hybrid_position [i] [mispair]))



















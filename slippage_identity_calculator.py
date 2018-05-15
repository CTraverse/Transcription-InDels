# -*- coding: utf-8 -*-
"""
Created on Fri May 20 11:25:05 2016

@author: chuck
"""


from sys import argv
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

script, infile = argv

# This script was a test script for calculating the %identity after a slippage event

deletions = list(open(infile, "r"))

slipped_seq = str()
del_seq = str()

mispair_list = ["rA_dA", "rA_dC", "rA_dG", 
                "rT_dT", "rT_dC", "rT_dG",
                "rC_dA", "rC_dC", "rC_dT",
                "rG_dA", "rG_dT", "rG_dG"]

mispairs = {}

mispairs["rA_dA"] = 0
mispairs["rA_dC"] = 0
mispairs["rA_dG"] = 0

mispairs["rT_dT"] = 0
mispairs["rT_dC"] = 0
mispairs["rT_dG"] = 0

mispairs["rC_dA"] = 0
mispairs["rC_dC"] = 0
mispairs["rC_dT"] = 0

mispairs["rG_dA"] = 0
mispairs["rG_dT"] = 0
mispairs["rG_dG"] = 0

hybrid_position = {}
for i in range(9):
    
    hybrid_position [i] = {}
    hybrid_position [i] ["matches"] = 0
    hybrid_position [i] ["mismatches"] = 0
    
    for mispair in mispair_list:
        
        hybrid_position[i][mispair] = 0


for deletion in deletions:

    del_split = deletion.split("\t")
    slipped_seq = del_split[3]
    del_seq = del_split[3] + del_split[4]
    del_len = len(del_split[4])    
        
    slipped_to_seq = del_seq[ len(del_seq) - 9 : len(del_seq) ]
        
    same_count = 0    
    j = 0
    for base_1, base_2 in zip(slipped_seq, slipped_to_seq):
        if base_1 == base_2:
            hybrid_position[j]["matches"] += 1
            j += 1
        
        else:
            
            hybrid_position[j]["r%s_d%s" % (base_1, str(Seq(base_2, IUPAC.unambiguous_dna).reverse_complement()))] += 1
            hybrid_position[j]["mismatches"] += 1
            j += 1

        
for i in range(9):
    print("Position %s of the RNA:DNA hybrid" % (i))
    print("Matches: %s" % (hybrid_position [i]["matches"]))
    print("Mismatches: %s" % (hybrid_position [i]["mismatches"]))
    
    for mispair in mispair_list:
        print("%s\t%s" % (mispair, hybrid_position [i] [mispair]))



















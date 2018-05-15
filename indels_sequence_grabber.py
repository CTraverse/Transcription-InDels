# -*- coding: utf-8 -*-
"""
Created on Mon Aug 29 16:29:56 2016

@author: chuck
"""
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from sys import argv
#### This script will sort through each indel and determine the deleted/inserted sequence
##### along with the 9 bases on either side of the indel
#### Determine strand and grab from raw genome sequence
#### Reverse complement if on reverse strand
#### Also, I built in automatic adjustment for ambiguous sequence



script, infile, reference_genome, outfile = argv

indels = list(open(infile, 'r'))
ref = str(list(SeqIO.parse(reference_genome, 'fasta'))[0].seq)


### Determine strand

deletion_list = []
insertion_list = []

for indel in indels:
    indel_split = indel.split("\t")
    
    position = int(indel_split[1])
    CIGAR = indel_split[2]
    complement = indel_split[3]
    read_seq = indel_split[4]
    
    if indel_split[2].count("D") == 1: #If it's a deletion
        
        #Start and end positions of deletion
        del_start = position +  int(CIGAR.split("M")[0]) - 1
        del_end = del_start + int(CIGAR.split("M")[1].split("D")[0])
        
        #DNA sequence of deletion
        del_seq = ref[del_start:del_end]
        upstream_seq = ref[del_start - 9:del_start]
        downstream_seq = ref[del_end:del_end + 9]

        if indel_split[3] == "False": #Not complemented
            if len(del_seq) >= 0:
                
                while del_seq[0] == downstream_seq[0]:

                    del_start += 1
                    del_end += 1
                    
                    del_seq = ref[del_start:del_end]
                    upstream_seq = ref[del_start - 9:del_start]
                    downstream_seq = ref[del_end:del_end + 9]
                    
            deletion_list.append("%s\t%s\t%s\t%s\t%s\t%s\n" % (indel_split[0], del_start, indel_split[2], upstream_seq, del_seq, downstream_seq))
            
        else: #Otherwise, the read aligned to a complemented gene
           
            temp_upstream_seq = str(Seq(downstream_seq, IUPAC.unambiguous_dna).reverse_complement())
            temp_downstream_seq = str(Seq(upstream_seq, IUPAC.unambiguous_dna).reverse_complement())
            
            del_seq = str(Seq(del_seq, IUPAC.unambiguous_dna).reverse_complement())
            downstream_seq = temp_downstream_seq
            upstream_seq = temp_upstream_seq
            if len(del_seq) >= 0:
                
                while del_seq[0] == downstream_seq[0]:

                    del_start -= 1
                    del_end -= 1
                    
                    del_seq = str(Seq(ref[del_start:del_end], IUPAC.unambiguous_dna).reverse_complement())
                    downstream_seq = str(Seq(ref[del_start - 9:del_start], IUPAC.unambiguous_dna).reverse_complement())
                    upstream_seq = str(Seq(ref[del_end:del_end + 9], IUPAC.unambiguous_dna).reverse_complement())

            deletion_list.append("%s\t%s\t%s\t%s\t%s\t%s\n" % (indel_split[0], del_start, indel_split[2], upstream_seq, del_seq, downstream_seq))
                    

    
    else: #Otherwise, it is an insertion
        
        ins_start = position +  int(CIGAR.split("M")[0]) - 1
        ins_end = ins_start + int(CIGAR.split("M")[1].split("I")[0])
        
        ins_seq = read_seq[int(CIGAR.split("M")[0]): int(CIGAR.split("M")[0]) + int(CIGAR.split("M")[1].split("I")[0])]
        upstream_seq = ref[ins_start - 9:ins_start]
        downstream_seq = ref[ins_start:ins_start + 9]

        if indel_split[3] == "False": #Not complemented

            read_seq_mod = 0
            while ins_seq[0] == downstream_seq[0]:

                ins_start += 1
                read_seq_mod += 1
                
                ins_seq = read_seq[int(CIGAR.split("M")[0]) + read_seq_mod: int(CIGAR.split("M")[0]) + int(CIGAR.split("M")[1].split("I")[0]) + read_seq_mod]
                upstream_seq = ref[ins_start - 9:ins_start]
                downstream_seq = ref[ins_start:ins_start + 9]
                
            insertion_list.append("%s\t%s\t%s\t%s\t%s\t%s\n" % (indel_split[0], ins_start, indel_split[2], upstream_seq, ins_seq, downstream_seq))                
                
        else:

            ins_seq = str(Seq(ins_seq, IUPAC.unambiguous_dna).reverse_complement())
            temp_downstream_seq = str(Seq(upstream_seq, IUPAC.unambiguous_dna).reverse_complement())
            temp_upstream_seq = str(Seq(downstream_seq, IUPAC.unambiguous_dna).reverse_complement())
            
            upstream_seq = temp_upstream_seq
            downstream_seq = temp_downstream_seq
            
            while ins_seq[0] == downstream_seq[0]:

                ins_start -= 1
                read_seq_mod -= 1
                
                ins_seq = str(Seq(read_seq[int(CIGAR.split("M")[0]) + read_seq_mod: int(CIGAR.split("M")[0]) + int(CIGAR.split("M")[1].split("I")[0]) + read_seq_mod], IUPAC.unambiguous_dna).reverse_complement())
                downstream_seq = str(Seq(ref[ins_start - 9:ins_start], IUPAC.unambiguous_dna).reverse_complement())
                upstream_seq = str(Seq(ref[ins_start:ins_start + 9], IUPAC.unambiguous_dna).reverse_complement())
                
            insertion_list.append("%s\t%s\t%s\t%s\t%s\t%s\n" % (indel_split[0], ins_start, indel_split[2], upstream_seq, ins_seq, downstream_seq))
            
del_out = open(re.sub(".txt", "", outfile) + "_deletions.txt", "w")
ins_out = open(re.sub(".txt", "", outfile) + "_insertions.txt", "w")

for insertion in insertion_list:
    ins_out.write(insertion)

slipped_seq = str()
del_seq = str()

for deletion in deletion_list:

    del_split = deletion.split("\t")
    
    slipped_seq = del_split[3]
    del_seq = del_split[3] + del_split[4]
    slipped_to_seq = del_seq[ len(del_seq) - 9 : len(del_seq) ]
    
    
    same_count = 0    
        
    for base_1, base_2 in zip(slipped_seq, slipped_to_seq):
        if base_1 == base_2:
            same_count += 1
    
    del_out.write(deletion.strip("\n") + "\t" + str(same_count) + "\n")


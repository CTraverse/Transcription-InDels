# -*- coding: utf-8 -*-
"""
Created on Wed Aug 24 13:32:50 2016

@author: chuck
"""

"""
This script determines location of indel (insertion/deletion) in the genome.
infile is 9_alignment and 10_alignment concatenated into one sam file. Data from CirSeqDistribution_v3 pipeline
"""

from Bio import SeqIO
from sys import argv

script, infile, replicate_name, reference_genome, outfile = argv

genome = list(SeqIO.parse(reference_genome, "fasta"))

#This function will find where in the genome the deletion is located
def genome_loc(position_of_indel, reference_genome):
    flag = False # Flag to tell if indel position is within a protein coding gene. Default False
    complement = str("False")
    
    for gene in reference_genome:
        if "complement" in str(gene.description): #Checks for complement because the position annotation is different in NCBI reference genomes for complement vs non-complement genes
            start = int(str(gene.description).split("complement(")[1].split("..")[0])
            end = int(str(gene.description).split("complement(")[1].split("..")[1].split(")")[0])
            complement = str("True")
        else:
            start = int(str(gene.description).split("location=")[1].split("..")[0])
            end = int(str(gene.description).split("location=")[1].split("..")[1].split("]")[0])
            complement = str("False")
        if start <= position_of_indel <= end: #Asks if the given position is between the start and end of the current gene
            flag = True
            break #Breaks the loop if position is between gene because we don't want intergenic indels
            
    return flag, complement

sam_file = list(open(infile, "r"))


#Creates dictionary of all rotated reads
## keys are read identifiers and each key calls a list of all rotated reads with that read identifier 
read_identifier = ""
sam_dict = {}
for read in sam_file:
    read_split = read.split("\t")   
    if read_split[0] != read_identifier:
        sam_dict.update( { read_split[0] : [] } )
        sam_dict[read_split[0]].append(read)
        read_identifier = read_split[0]
    elif read_split[0] == read_identifier:
        sam_dict[read_split[0]].append(read)

# Loop through the rotated reads and find all insertions/deletions
master_list = []
for read in sam_dict:
    rotated_list = sam_dict[read]
    # Finds maximum alignment score of all rotated reads for that read identifier
    max_AS = max(rotated_list, key = lambda rotate : int(rotate.split("\t")[11].split(":")[2])).split("\t")[11] 
    
    #Grab the full read information associated with the maximum alignment score
    for rotate in rotated_list:
        rotate_split = rotate.split("\t")
        if rotate_split[11] == max_AS:
            rotate_split
            #the CIGAR string tells us if there is a misalignment. I or D indicates a deletion
            CIGAR = rotate_split[5]
            #Get rid of the numbers, we're only interested in the string designations
            CIGAR_replace = filter(lambda x: x.isalpha(), CIGAR)
            read_seq = rotate_split[9]
            
            ### Find Insertions and Deletions within the lists
            if CIGAR_replace == "MDM": ### Determines if max_AS for list is a deletion with no soft clipping
                #Grab the sequences directly before and after the deletion using the numeric parts of the CIGAR indicator
                left_del_loc = int(CIGAR.split("M")[0]) - 1
                right_del_loc = left_del_loc + int(CIGAR.split("M")[1].split("D")[0])
                #Determine the quality of the bases before and after the deletion
                position = int(rotate_split[3]) + left_del_loc                
                left_Q = ord(rotate_split[10][left_del_loc]) - 33 # Quality score left of deletion
                right_Q = ord(rotate_split[10][right_del_loc]) - 33 # Quality score right of deletion  
                
                #Only keep deletions if the surrounding bases had high quality scores
                if left_Q >= 20 and right_Q >= 20:
                    #Grab the location in the genome for the current deletion
                    in_gene, gen_complement = genome_loc(position, genome)
                    #Only take deletions if they fall within a gene, and not between two genes
                    if in_gene == True:                       
                        if len(master_list) > 0:
                            #Check if there's a duplicate already written. Sometimes two rotations can have identical 
                            #alignment scores. If this is the case, just keep the first one and discard any others
                            if master_list[len(master_list) - 1].split("\t")[0] == rotate_split[0]:
                                pass #don't write
                            else:
                                #If the last entry is not of the same read, write the deletion 
                                master_list.append(str(rotate) + "\t%s" % gen_complement)                       
                        else:
                            #Write the first entry
                            master_list.append(str(rotate) + "\t%s" % gen_complement)                         

            ### Determines if max_AS for list is an insertion with no soft clipping                 
            elif CIGAR_replace == "MIM":
                #Grab the inserted sequences using the CIGAR indicator
                left_ins_loc = int(CIGAR.split("M")[0]) - 1
                right_ins_loc = left_ins_loc + int(CIGAR.split("M")[1].split("I")[0]) + 1
                #Determine the quality of the bases before and after the insertion
                position = int(rotate_split[3]) + left_ins_loc
                left_Q = ord(rotate_split[10][left_ins_loc]) - 33
                right_Q = ord(rotate_split[10][right_ins_loc]) - 33
                ins_Q = rotate_split[10][left_ins_loc + 1:right_ins_loc]
                Q_list = []
                
                #Populate list of quality scores of each inserted base
                for Q in ins_Q: 
                    Q_list.append(ord(Q) - 33)
                # Check if all of the inserted bases are at least 20
                if min(Q_list) >= 20: 
                    #Get location of insertion in the genome
                    in_gene, gen_complement = genome_loc(position, genome)
                    #Only take insertions if they fall within a gene, and not between two genes
                    if in_gene == True:                   
                        if len(master_list) > 0:
                            #Check if there's a duplicate already written. Sometimes two rotations can have identical 
                            #alignment scores. If this is the case, just keep the first one and discard any others
                            if master_list[len(master_list) - 1].split("\t")[0] == rotate_split[0]:
                                pass #don't write
                            else:
                                #If the last entry is not of the same read, write the deletion 
                                master_list.append(str(rotate) + "\t%s" % gen_complement)                
                        else:
                            master_list.append(str(rotate) + "\t%s" % gen_complement)       
        
#Write the master list into a file of all of the insertions and deletions that we found
out = open(outfile, "w")
for item in master_list:
    item_split = item.split("\t")
    out.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (replicate_name, item_split[3], item_split[5], item_split[len(item_split) - 1], item_split[9], item_split[10])) 



# -*- coding: utf-8 -*-
"""
Created on Fri Jun 24 10:14:51 2016

@author: chuck
"""



from sys import argv
from Bio import SeqIO


script, coding_sequence, coverage_file, indels, outfile = argv

#This script find all of the mononucleotide repeats in the coding regions of a genome

genes = list(SeqIO.parse(coding_sequence, "fasta"))
Q_20 = list(open(coverage_file, "r"))
error_list = list(open(indels, "r"))


repeat_dict = {}
repeat_dict[4], repeat_dict[5], repeat_dict[6], repeat_dict[7], repeat_dict[8]  = [], [], [], [], []
repeat_dict[9], repeat_dict[10], repeat_dict[11], repeat_dict[12], repeat_dict[13] = [], [], [], [], []


gene_num = 0
for gene in genes:
    gene_num += 1
    gene_len = len(gene)
    position = 0
    gene_seq = str(gene.seq)
    
    #Determine if gene is complemented or not to find the gene start and end locations
    if "complement" in str(gene.description):
        correct_orientation = False     
        gene_start = int(gene.description.split("complement(")[1].split("..")[0].strip("]").strip(">").strip("<"))
        gene_end =  int(gene.description.split("complement(")[1].split("..")[1].strip(")]").strip(">").strip("<"))
    
    else:
        gene_start = int(gene.description.split("location=")[1].split("..")[0].strip("]").strip(">").strip("<"))
        gene_end =  int(gene.description.split("location=")[1].split("..")[1].strip(")]").strip(">").strip("<"))
        correct_orientation = True
    
    #Walk across the gene up until your two bases away from the end
    while position < gene_len - 2:
        
        base_identity = str(gene.seq[position])
        repeat_len = 1
        position_list = [position]
        
        #Find consecutive identical bases
        while gene_seq[position] == gene_seq[position + repeat_len]:
            position_list.append(position + repeat_len)
            repeat_len += 1
            print(position + repeat_len, gene_len, gene_num)
            if position + repeat_len >= len(gene_seq):
                print("Greater!")
                break
            print(repeat_len)
        
        #Only continue if 4 or more consecutive bases were identical
        if repeat_len >= 4:

            errors = 0
            repeat_cov = 0
            #Calculate the coverage for each mononucleotide reeat
            if correct_orientation == True:
                
                for i in position_list:
                    repeat_cov += int(Q_20[gene_start + i - 1].split("\t")[2]) + int(Q_20[gene_start + i - 1].split("\t")[3]) + int(Q_20[gene_start + i - 1].split("\t")[4]) + int(Q_20[gene_start + i - 1].split("\t")[5].strip("\n"))
                    
                repeat_end = gene_start + position + repeat_len - 1
                
                for error in error_list:
                    error_position = int(error.split("\t")[1]) + 1
                    
                    if error_position - 1 <= repeat_end <= error_position + repeat_len:
                        errors += 1
                        
                repeat_dict[repeat_len].append("%s\t%s\t%s\t%s" % (str(repeat_end), base_identity, repeat_cov, errors))
                
            elif correct_orientation == False:
                
                for i in position_list:
                    repeat_cov += int(Q_20[gene_end - i - 1].split("\t")[2]) + int(Q_20[gene_end - i -1].split("\t")[3]) + int(Q_20[gene_end - i -1].split("\t")[4]) + int(Q_20[gene_end - i -1].split("\t")[5].strip("\n"))

                repeat_end = gene_end - position - repeat_len 
                
                for error in error_list:
                    error_position = int(error.split("\t")[1])

                    if error_position - 2 <= repeat_end <= error_position +2:
                        errors += 1

                repeat_dict[repeat_len].append( "%s\t%s\t%s\t%s" % (str(repeat_end), base_identity, repeat_cov, errors))

        else:
            pass

        position += repeat_len


out = open(outfile, "w")
i = 4
while i <= 13:
    positions = repeat_dict[i]
    for position in positions:
        out.write(str(position) + "\t" + str(i) + "\n")
    i += 1
out.close()




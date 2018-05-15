# -*- coding: utf-8 -*-
"""
Created on Thu May 26 14:03:13 2016

@author: chuck
"""
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import SeqIO
from random import randint
from sys import argv

script, coding_seq, ref, deletions, threshold, outfile = argv

### This script generates a simulated set of deletions sampled across the transcriptome.
### It first tabulates the overall transcript expression level of each gene across the genome
### It then randomly samples the genome, weighted by its expression level and generates a deletion randomly within the gene

deletion_list = list(open(deletions, "r"))
genome = list(SeqIO.parse(coding_seq, "fasta"))
threshold = list(open(threshold, "r"))
out = open(outfile, "w")
ref =str(list(SeqIO.parse(ref, "fasta"))[0].seq)

per_gene_cov = 0
genome_cov = 0

master_list = []

## Tabulate the overall per gene coverage 
for gene in genome:
    
    #Determine the start and end of the current gene and determine if current gene is reverse complemented relative to reference  
    correct_orientation = True
    if "complement" in str(gene.description):
        correct_orientation = False     
        gene_start = int(gene.description.split("complement(")[1].split("..")[0].strip("]").strip(">").strip("<"))
        gene_end = int(gene.description.split("complement(")[1].split("..")[1].strip(")]").strip(">").strip("<"))
    else:

        gene_start = int(gene.description.split("location=")[1].split("..")[0].strip("]").strip(">").strip("<"))
        gene_end = int(gene.description.split("location=")[1].split("..")[1].strip(")]").strip(">").strip("<"))
        
    location = gene_start - 1
    
    gene_coverage = 0
    
    #Loop through each position in the current gene
    while location <= gene_end - 1:
        
        current_line = threshold[location].strip("\n").split("\t") #Split current line into each column
        base_position = int(current_line[0]) #Current position in reference
        coverage = int(current_line[2]) + int(current_line[3]) + int(current_line[4]) + int(current_line[5]) #calculate total coverage for position
        current_base = current_line[1] #Identity of the current base
        current_gene =  str(gene.seq) #Sequence of the current gene
        
        gene_coverage += coverage
        genome_cov += coverage
        location += 1
    
    per_base_cov = int(float(float(gene_coverage) / float(len(gene.seq))) * 100)
    per_gene_cov += per_base_cov
    
    
    if per_base_cov == 0:
        pass
    else:
        master_list.append([per_gene_cov, str(gene.seq), gene_start, correct_orientation])
    
#Randomly sample the genome, weighted by gene coverage and generate deletions
random_deletions = []
j = len(deletion_list)
for i in range(2):
    j = len(deletion_list)
    print(i)
    for deletion in deletion_list:
        del_split = deletion.split("\t")
        j -= 1
        random_gene = randint(1, per_gene_cov)
        for gene in master_list:
            
            ## This is a special case for the first gene in the genome list
            if int(gene[0]) == int(master_list[0][0]):
                if 0 < random_gene <= gene[0]:
    
                    gene_seq = gene[1]
    
                    random_deletion = randint(14, len(gene_seq) - len(del_split[4]) - 15)
                    
                    del_start = int(gene[2]) + random_deletion
                    del_end = int(gene[2]) + random_deletion + len(del_split[4])                    
                    
                    left_seq = ref[del_start - 15: del_start]
                    del_seq = ref[del_start: del_end]
                    
                    right_seq = ref[del_end : del_end + 15]
    

                    while True:
                        if del_seq[0] == right_seq[0]:
                            del_start += 1
                            del_end += 1
                            
                            left_seq = ref[del_start - 15: del_start]
                            del_seq = ref[del_start: del_end]
                            right_seq = ref[del_end : del_end + 15]
                        else:
                            break
                    random_deletions.append("%s\t%s\t%s\t%s\t%s\t%s\n" % (del_split[0] + "_sim", del_split[1], del_split[2], left_seq, del_seq, right_seq))

                else:
                    pass
               
            else:
    
                if int(master_list[master_list.index(gene) - 1][0]) < random_gene > int(gene[0]):
                    pass
                
                elif int(master_list[master_list.index(gene) - 1][0]) < random_gene <= int(gene[0]):
                    gene_seq = gene[1]
    
                    random_deletion = randint(14, len(gene_seq) - len(del_split[4]) - 15)
                    
                    del_start = int(gene[2]) + random_deletion
                    del_end = int(gene[2]) + random_deletion + len(del_split[4])                    
                    
                    left_seq = ref[del_start - 15: del_start]
                    del_seq = ref[del_start: del_end]
                    right_seq = ref[del_end : del_end + 15]

                    if gene[3] == False:
                        while True:
                            if del_seq[0] == right_seq[0]:
                                del_start += 1
                                del_end += 1
                                
                                left_seq = ref[del_start - 15: del_start]
                                del_seq = ref[del_start: del_end]
                                right_seq = ref[del_end : del_end + 15]
                            else:
                                break
                    
                    else:
                        temp_left = right_seq ## Swap left and right because the sequences are being reverse complemented
                        temp_right = left_seq
                        
                        left_seq = str(Seq(str(temp_left), IUPAC.unambiguous_dna).reverse_complement())
                        right_seq = str(Seq(str(temp_right), IUPAC.unambiguous_dna).reverse_complement())
                        del_seq = str(Seq(str(del_seq), IUPAC.unambiguous_dna).reverse_complement())
                        

                        while True:
                            if del_seq[0] == right_seq[0]:
                                del_start -= 1
                                del_end -= 1
                                
                                del_seq = str(Seq(str(ref[del_start:del_end]), IUPAC.unambiguous_dna).reverse_complement())
                                right_seq = str(Seq(str(ref[del_start - 15:del_start]), IUPAC.unambiguous_dna).reverse_complement())
                                left_seq = str(Seq(str(ref[del_end:del_end + 15]), IUPAC.unambiguous_dna).reverse_complement())

                            else:
                                break

                    random_deletions.append("%s\t%s\t%s\t%s\t%s\t%s\n" % (del_split[0] + "_sim", del_split[1], del_split[2], left_seq, del_seq, right_seq))
                
                else:
                    pass



for rand in random_deletions:
    del_split = rand.split("\t")
    
    slipped_seq = del_split[3]
    del_seq = del_split[3] + del_split[4]
    slipped_to_seq = del_seq[ len(del_seq) - 15 : len(del_seq) ]
    
    
    same_count = 0    
        
    for base_1, base_2 in zip(slipped_seq, slipped_to_seq):
        if base_1 == base_2:
            same_count += 1
    
    
    out.write(rand.strip("\n") + "\t" + str(same_count) + "\n")
    
    
out.close()






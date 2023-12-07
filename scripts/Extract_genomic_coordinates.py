#!/usr/bin/python3

######################################################################
#ALEJANDRO VALENZUELA - RETRIEVAL GENOMIC COORDINATES TEMPLATE SCRIPT#
######################################################################


#import pandas as pd
#import numpy as np
import sys
#import requests
from itertools import combinations
#import seaborn as sns
#from scipy import stats
#import pickle
from collections import Counter
#import copy
#from scipy.stats import sem, t, chi2_contingency, norm
#from scipy import mean
import re
#import os
import gzip
import tarfile
#import time
#from Bio import AlignIO


"""
Here we create the function for checking the input parameters and saving
them in different variables, error if the usage is not good
"""

if len(sys.argv) == 9:
    alignment_path = sys.argv[1]
    fasta_path = sys.argv[2]
    tracking_file = sys.argv[3]
    gff_file = sys.argv[4]
    caas_file = sys.argv[5]
    output_coordinates = sys.argv[6]
    gene_equivalences = sys.argv[7]
    gene_name_file = sys.argv[8]
else:
    sys.exit("The usage shoud be: length 8 files")


#align = AlignIO.read(alignment_path, "fasta")
selected_positions = []
removed_positions = []


#FUNCTION A) --> GET INDEX OF FILTERED POSITIONS IN ALIGNMENTS (html files)
# Read gene names from the TSV file
#gene_names = []
with open(gene_name_file, 'r') as gene_file:
    for line in gene_file:
        # Assuming gene names are separated by tabs
        # Extract the gene name as a string
        gene_name = line.rstrip().split('\t')[0]
        print(gene_name)
        #gene_names.extend(gene_in_line)
        # Iterate over each gene name
        #for gene_name in gene_names:
        first_char = gene_name[0]  # First character
        second_char = gene_name[:2]  # First two characters
        
        with tarfile.open(tracking_file, "r:gz") as tar:
            # List all files and directories in the zip
            for member in tar.getmembers():
                # Check if the member name matches the pattern
                if member.name.startswith(f"{first_char}/{second_char}/{gene_name}."):
                    with tar.extractfile(member) as file:
                        for line in file:
                            line = line.decode('utf-8').rstrip()  # Decode line
                            if "selected:" in line:
                                print(line)
                                intervals = line.split(" ")
                                for i in range(4,len(intervals)):
                                    if "-" in intervals[i]:
                                        init = int(intervals[i].split("-")[0])
                                        end = int(intervals[i].split("-")[1])
                                        nums = list(range(init, end+1))
                                        for num in nums:
                                            selected_positions.append(num)
                                    else:
                                        selected_positions.append(int(intervals[i]))
                            elif "removed:" in line:
                                print(line)
                                intervals = line.split(" ")
                                for i in range(5,len(intervals)):
                                    if "-" in intervals[i]:
                                        init = int(intervals[i].split("-")[0])
                                        end = int(intervals[i].split("-")[1])
                                        nums = list(range(init, end+1))
                                        for num in nums:
                                            removed_positions.append(num)
                                    else:
                                        removed_positions.append(int(intervals[i]))


    #FUNCTION B.1) --> GET INFORMATION FROM EQUIVALENCE BETWEEN GENE NAMES AND PROTEINS

        # Reset variables for each gene
        protein_id = None
        target_sequence = False
        aligned_CDS = ""

        # Existing code modified to be inside the loop
        with open(gene_equivalences, "r") as in_fh2a:
            for line in in_fh2a:
                line = line.rstrip()
                if line.split()[0] == gene_name:
                    protein_id = line.split()[1]
                    break  # Exit the loop once the matching gene is found

        # Check if protein_id was found
        if protein_id is None:
            print(f"Protein ID for gene {gene_name} not found.")
            continue  # Skip to the next gene_name

        # Rest of your code that uses protein_id
        with open(gene_equivalences, "r") as in_fh2a:
            target_sequence = False
            aligned_CDS = ""
            for line in in_fh2a:
                line = line.rstrip()
                if line.split()[0] == gene_name:
                    #print(line.split()[1])
                    protein_id = line.split()[1]
        if protein_id is None:
            print(f"Protein ID for gene {gene_name} not found.")
            continue  # Skip to the next gene_name
        #print(transcript_name)

        #FUNCTION B) --> GET CORRESPONDING HUMAN CDS POSITIONS ALIGNED
        with tarfile.open(alignment_path, "r:gz") as in_fh2:
            # List all files and directories in the zip
            for member in in_fh2.getmembers():
                # Check if the member name matches the pattern
                if member.name.startswith(f"{first_char}/{second_char}/{gene_name}."):
                    with in_fh2.extractfile(member) as file:
                        target_sequence = False
                        aligned_CDS = ""
                        for line in file:
                            line = line.decode('utf-8').rstrip()  # Decode line
                        #for line in in_fh2:
                            if line.startswith(">Homo"):
                                target_sequence = True
                            elif line.startswith(">") and not line.startswith(">Homo"):
                                target_sequence = False
                            elif not line.startswith(">") and target_sequence:
                                aligned_CDS += line

        #FUNCTION C) read fasta sequence

        #print(fasta_path)
        with gzip.open(fasta_path, "rt") as in_fh3:
            fasta_header = protein_id
            #print(fasta_header)
            real_CDS = ""
            for line in in_fh3:
                line = line.rstrip()
                if line.startswith(">"+fasta_header):
                    target_sequence = True
                elif line.startswith(">") and not line.startswith(">"+fasta_header):
                    target_sequence = False
                elif not line.startswith(">") and target_sequence:
                    real_CDS += line

        #print(aligned_CDS)
        #print(real_CDS)



        #FUNCTION D) Get indexes of human sequence based on alignment

        counter_nt = 3
        counter_codons = 0

        prev_ambiguity=False
        missmatch = False
        n_ambiguities=0
        ref_aligned_position = 0
        indexes = {}
        for real_position in range(0,len(real_CDS),3):
            codon = real_CDS[real_position:real_position+3]
            #print(codon, "real")
            for aligned_position in range(ref_aligned_position,len(aligned_CDS),3):
                aligned_codon = aligned_CDS[aligned_position:aligned_position +3]
                #print(aligned_codon, "aligned")
                if aligned_CDS[aligned_position] == "-":
                    counter_nt += 0
                    counter_codons += 1
                    indexes[counter_codons] = 0
                elif aligned_codon == codon:
                    counter_nt += 3
                    counter_codons += 1
                    break
                elif aligned_codon != codon:
                    missmatch=True
                    break
            ref_aligned_position = 3*counter_codons
            indexes[counter_codons] = counter_nt-3
            if missmatch:
                break

        #Selected_positions_filtered is a list (with the 1-based set of codons filtered in your alignments)
        selected_positions_filt = selected_positions
        #Selected indexes is a dictionary (it accumulates the count of CDS nt along the MSA) 
        selected_indexes = indexes
        print(selected_positions_filt)
        print(selected_indexes)

        #FUNCTION F) GET coordinate info from gff gff_file

        starts = []
        ends = []
        CDS_positions = []
        with open(gff_file, "r") as in_fh4:
            fasta_header = protein_id
            real_CDS = ""
            for line in in_fh4:
                line = line.rstrip()
                if fasta_header in line and line.split("\t")[2] == "CDS":
                    #print(line)
                    fields = line.split("\t")
                    chrom = fields[0]
                    start = int(fields[3])
                    starts.append(start)
                    end = int(fields[4])
                    ends.append(end)
                    strand = fields[6]
                    print(strand)
            if strand == "+":
                init_coord = starts[0]
                for i in range(0, len(starts)):
                    for j in range(starts[i], ends[i]+1):
                        #print(j)
                        CDS_positions.append(j)
            elif strand == "-":
                init_coord = ends[-1]
                for i in range(0, len(starts)):
                    curr_end = ends[-(i+1)]
                    curr_start = starts[-(i+1)]
                    for j in range(curr_end, curr_start-1, -1):
                        #print(j)
                        CDS_positions.append(j)

        #####################################################
        #Translate selected positions to output_coordinates##
        #####################################################

        #print(CDS_positions)
        #print(strand)

        with open(caas_file, "r") as in_fh5:
            for line in in_fh5:
                fields=line.rstrip().split("\t")
                if fields[0] == gene_name and int(fields[2]) < len(selected_positions_filt):
                    caas_position = int(fields[2])
                    prefiltered_position = selected_positions_filt[caas_position]
                    if prefiltered_position <= len(selected_indexes):
                        prefiltered_nt = selected_indexes[prefiltered_position]
        #LAST ELEMENT ADDING
                        if prefiltered_nt == len(CDS_positions):
                            if strand == "+":
                                current_coord = CDS_positions[-1]+1
                                real_coord = CDS_positions[-1]
                            elif strand == "-":
                                current_coord = CDS_positions[-1]-1
                                real_coord = CDS_positions[-1]
                            with open(output_coordinates, "a") as in_fh6:
                                print("{}\t{}\t{}\t{}\t{}\t{}\t{}".format(fields[0], fields[1], fields[2], fields[3], chrom, current_coord, real_coord), file=in_fh6)
        #NON-GAPPY filtered positions in reference
                        elif prefiltered_nt != 0:
                            #print(prefiltered_nt, type(prefiltered_nt), "NON-O")
                            #Current coord does not correspond to the real coordinate from the codon presenting the CAAS (it is used as a python index)
                            current_coord = CDS_positions[prefiltered_nt]
                            #Real coord is the actual coordinate corresponding to the actual genomic position targeted by the codon presenting CAAS
                            real_coord = CDS_positions[prefiltered_nt-1]
                            with open(output_coordinates, "a") as in_fh6:
                                print("{}\t{}\t{}\t{}\t{}\t{}\t{}".format(fields[0], fields[1], fields[2], fields[3], chrom, current_coord, real_coord), file=in_fh6)

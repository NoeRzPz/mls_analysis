#!/usr/bin/python3

###########################################################
###
## ALEJANDRO VAALENZUELA SEBA - GENOMIC PRIMATE COORD TO MAMMAL AA ZOONOMIA 
###########################################################################

import pandas as pd
import numpy as np
from itertools import combinations
#import seaborn as sns
from scipy import stats
from collections import Counter
from scipy.stats import sem, t, chi2_contingency, norm
from scipy import mean
from Bio import AlignIO


"""
FUNCTION to retrieve MAMMALIAN AA positions correponding from primate coordinates with CAAS
"""

CDS_indexes = {}
with open(primate_coordinates, "r") as in_fh8, open(mammal_aa, "w") as out_fh:
    print(CDS_positions, len(CDS_positions))    # No estan definidas en el script
        #CDS_positions.reverse()
    for line in in_fh8:
        fields=line.rstrip().split("\t")
        print(fields)
        chr = fields[4]
        pos = int(fields[5])
        #print(chr, pos, strand)
        print("{}\t{}\t{}\t{}\t{}\t{}\t".format(fields[0], fields[1], fields[2], fields[3], fields[4], fields[6]), end="", file=out_fh)
        if pos == CDS_positions[-1]+1 and strand == "+":
            CDS_index = len(CDS_positions)
        elif pos == CDS_positions[-1]-1 and strand == "-":
            CDS_index = len(CDS_positions)
        else:
            CDS_index = CDS_positions.index(pos)
        print(str(CDS_index)+"cdsindex")
        #Get CDS distance from start
        if CDS_index in list(selected_indexes.values()):
            print(str(CDS_index)+"cdsindex")
            print(selected_indexes)
            prefiltered_position = (list(selected_indexes.keys())[list(selected_indexes.values()).index(CDS_index)])
            print(list(selected_indexes.keys()))
            print(list(selected_indexes.values()))
            #print(CDS_index)
            print(prefiltered_position)
            if prefiltered_position in selected_positions_filt and prefiltered_position != 0:
                print("hay posicion para ese CAAS")
                print(prefiltered_position) 
                caas_position = selected_positions_filt.index(prefiltered_position)
                print(caas_position)
                with open(output_coordinates, "a") as in_fh6:
                    print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(fields[0], fields[1], fields[2], fields[3], caas_position, chrom, current_coord, real_coord), file=out_fh)

                #alignment = AlignIO.read(protein_alignment_path, "phylip-relaxed")
                #for record in alignment:
                 #       #print(record.id, record.seq[caas_position])
                 #       print("{}|{} ".format(record.id, record.seq[caas_position]), file=out_fh, end="")
                #print(file=out_fh)
            else:
                print("filtered_position")
                with open(output_coordinates, "a") as in_fh6:
                    print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(fields[0], fields[1], fields[2], fields[3], "filtered_position", chrom, current_coord, real_coord), file=out_fh)
                #print("filtered_position", file=out_fh)
        else:
            #print("filtered_position", file=out_fh)
            with open(output_coordinates, "a") as in_fh6:
                    print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(fields[0], fields[1], fields[2], fields[3], "filtered_position", chrom, current_coord, real_coord), file=out_fh)
                

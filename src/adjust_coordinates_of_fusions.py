#!/usr/bin/env python3
from distutils.command.sdist import sdist
import sys
import argparse
#%%
def parse_index(index_file):
    chr2length = {}
    with open(index_file, "r") as index:
        for line in index:
            cols = line.strip().split('\t')
            chr, length = cols[0], cols[1]
            chr2length[chr] = length
    return(chr2length)

def parse_LFFF_output(fusion_info):
    chr2max_busco = {}
    with open(fusion_info, "r") as fusion_file:
        for line in fusion_file:
            cols = line.strip().split('\t')
            chr, busco_pos = cols[0], round(float(cols[3]), 1)
            if chr in chr2max_busco.keys():
                if busco_pos > chr2max_busco[chr]:
                    chr2max_busco[chr]= busco_pos # update value
            else:
                chr2max_busco[chr] = busco_pos
    return(chr2max_busco)


def adjust_last_coordinate_per_chr(chr2length, chr2max_busco, fusion_info, prefix):
    with open(prefix + 'adjusted_fusion_positions.tsv', 'w') as output_file:
        with open(fusion_info, 'r') as fusion_file:
            for line in fusion_file:
                cols = line.strip().split('\t')
                chr, busco_pos = cols[0], round(float(cols[3]), 1)
                max_pos = chr2max_busco[chr]
                if max_pos == busco_pos:
                    length = chr2length[chr]
                    cols[3] = length
                output_file.write(("%s\t%s\t%s\t%s\n") % (cols[0], cols[1], cols[2], cols[3]))
    return()
#%%

# The coordinates of fusion chromosomes output by LFFF are based on ortholog positions, 
# therefore the end of the last segment of chromosome is the position of the last ortholog. 
# This script adjusts this value to the precise chromosome end (i.e. chr length)

if __name__ == "__main__":
    SCRIPT = "adjust_fusion_coordinates.py"
    # argument set up
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--fusions", type=str, help = "Fusion position output from LFSF", required=True)
    parser.add_argument("-i", "--index", type=str, help = "index file for genome", required=True)
    parser.add_argument("-p", "--prefix", type=str, help = "prefix for output file", required=True)
    args = parser.parse_args()
    fusion_info = args.fusions
    index_file = args.index
    prefix = args.prefix
    chr2length = parse_index(index_file)
    chr2max_busco = parse_LFFF_output(fusion_info)
    chr2max_busco
    adjust_last_coordinate_per_chr(chr2length, chr2max_busco, fusion_info, prefix)
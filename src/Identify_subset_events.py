#!/usr/bin/env python3
#%%
import os, sys, glob, argparse
import warnings # suppresses warning from importing ete3 (due to codeml)
warnings.filterwarnings("ignore", category=SyntaxWarning)
from ete3 import Tree
from pathlib import Path

# %%
def parse_tree(tree_file):
	tip_list = []
	t = Tree(tree_file, format=1) # Format_1 means node names are present and read - we have node names thanks to the replacement of branch support values with unique numbers done above :)
	# First get the list of tips in the tree
	for node in t.traverse("preorder"):
		if node.is_leaf():
			tip_list.append(node.name)	
	return(t, tip_list)

def parse_table(table_file):
    with open(table_file, 'r') as table:
        Merians2Node_fusion = []
        for line in table:
            if not line.startswith("Merian"):
                cols = line.rstrip("\n").split("\t")
                merians, tips, node, event = cols[0], cols[1], cols[2], cols[3]
                merian_node = str(merians + '.' + node)
                if '(' in merians:
                    continue # deal with these later - bracketed events are inherently subsets
                else:   
                    if str(event) == 'fusion':
                        Merians2Node_fusion.append(merian_node)
    return(Merians2Node_fusion)

def identify_potential_subsets(Merians2Node_fusion): # Identify subsets of larger fusions - potential fusion subsets
    potential_subsets = {}
    for i in Merians2Node_fusion:
        Merians = i.split('.')[0]
        Merians = Merians.split(',')
        # check if these merians are a subset of another set of merians:
        copy_list = Merians2Node_fusion
        for merian in Merians:
            possible_substrings = [str(merian + '.'), str(merian + ',')] # match exactly e.g. only M1 not M10 as well
            copy_list = [str for str in copy_list if any(sub in str for sub in possible_substrings)]
        copy_list.remove(i)
        if len(copy_list) != 0: # i.e. if more than just the original event (i) has matched
            to_remove = []
            for entry in copy_list:
                entry_merians = entry.split('.')[0].split(',')
                matches = set(entry_merians) & set(Merians)
                #print(i, entry, matches)
                if len(entry_merians) == len(Merians):
                    to_remove.append(entry)
            final_list = [ele for ele in copy_list if ele not in to_remove]
            if len(final_list) != 0:
                potential_subsets[i] = final_list
    return(potential_subsets)

def infer_subsets(t, potential_subsets): # now need to check if phylogenetically each potential subset is actually a subset
    inferred_subsets = []
    node_content = t.get_cached_content(store_attr='name') # first get all node names
    for subset in potential_subsets:
        queries = potential_subsets[subset]
        subset_node = subset.split('.')[1]
        for node in t.traverse():
            if node.name == subset_node:
                subset_tips = node_content[node]
        for query in queries:
            query_node = query.split('.')[1] # now check if node is upstream of query node 
            for node in t.traverse():
                if node.name == query_node:
                    query_tips = node_content[node]
            matched = set(subset_tips) & set(query_tips)
            if len(matched) == len(query_tips): # i.e. if all spp with the "larger" fusion also have a subset of the fusion
#                print(subset, 'full_fusion node:', query_node, 'subset node is:', subset_node)
                inferred_subsets.append(str(subset + '.' + query))
    return(inferred_subsets)

def print_subsets_table(inferred_subsets, subsets_table_file): # make a final dataframe with same format as original mapped dataframe, just with subsets present
    with open(subsets_table_file, 'w') as subset_table:
        subset_table.write("%s\t%s\t%s\t%s" % ("Merians", "Node", "Event", "Parent_fusion") + "\n")
        for entry in inferred_subsets:
            inferred_subsets = entry.split('.')
            subset_merians = inferred_subsets[0]
            subset_node = inferred_subsets[1]
            full_merians = inferred_subsets[2]
            full_node = inferred_subsets[3]
            subset_table.write("%s\t%s\t%s\t%s" % (subset_merians, subset_node, full_merians, full_node) + "\n")

#%%
if __name__ == "__main__":
    SCRIPT = "Identify_subset_events.py"
    # argument set up
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_data", type=str, help = "path to lep_fusion_fission_finder output", required=True)
    parser.add_argument("-tree", "--tree", type=str, help = "Phylogenetic tree", default="fsf")
    parser.add_argument("-f", "--prefix", type=str, help = "Prefix for all output files", default="fsf")
    args = parser.parse_args()
    table_file = args.input_data
    tree_file = args.tree
    prefix = args.prefix    

    # Run functions
    print("\t[+] Identifying subsets of fusion events.")
    t, tip_list = parse_tree(tree_file)
    Merians2Node_fusion = parse_table(table_file) # makes list of Merians 2 node
    potential_subsets = identify_potential_subsets(Merians2Node_fusion)
    inferred_subsets = infer_subsets(t, potential_subsets)
    subsets_table_file = prefix + '_subset_fusions.tsv'
    print("\t[+] Writing subset events to " + prefix + "_subset_fusion.tsv")
    print_subsets_table(inferred_subsets, subsets_table_file)

# %%
# prefix = 'test'
# tree_file = '/lustre/scratch123/tol/teams/blaxter/projects/lepidoptera_genomics/cw22/Leps_200/Analysis/LFSF/final_analysis/r2_m17.newick.txt'
# table_file = '/lustre/scratch123/tol/teams/blaxter/projects/lepidoptera_genomics/cw22/Leps_200/Analysis/LFSF/final_analysis/noComplex_281022/mapped_fusions_fissions_noComplex_181022.tsv'
# t, tip_list = parse_tree(tree_file)
# Merians2Node_fusion = parse_table(table_file) # makes list of Merians 2 node
# potential_subsets = identify_potential_subsets(Merians2Node_fusion)
# inferred_subsets = infer_subsets(t, potential_subsets)
# subsets_table_file = prefix + '_summary.tsv'
# print_subsets_table(inferred_subsets, subsets_table_file)
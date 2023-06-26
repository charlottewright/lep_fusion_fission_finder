#!/usr/bin/env python3
#%%
import os
import warnings # suppresses warning from importing ete3 (due to codeml)
warnings.filterwarnings("ignore", category=SyntaxWarning)
from ete3 import Tree
import merian_tools
import importlib
importlib.reload(merian_tools)
from merian_tools import *
#%%
if __name__ == "__main__":
	SCRIPT = "map_fusion_fissions_client.py"
	# argument set up
	parser = argparse.ArgumentParser()
	parser.add_argument("-i", "--input_data", type=str, help = "path to lep_fusion_fission_finder output", required=True)
	parser.add_argument("-tree", "--tree", type=str, help = "Phylogenetic tree", default="fsf")
	parser.add_argument("-o", "--output", type=str, help = "output location relative to working directory", required=True)
	parser.add_argument("-f", "--prefix", type=str, help = "Prefix for all output files", default="fsf")
	parser.add_argument("-t", "--threshold", type=int, help = "Threshold for rearrangement to be shared between tips", default=0.75)
	parser.add_argument("-l", "--label_status", type=str, help = "Specify if tree already contains internal node labels", default='False')
	args = parser.parse_args()
	input_data = args.input_data
	tree_file = args.tree
	output_location = args.output
	prefix = args.prefix
	threshold = args.threshold
	label_status = args.label_status

# Read in arguments manually
#input_data = 'species_files/'
#tree_file = 'r2_m17.newick.txt'
#output_location = 'test_080123'
#prefix = 'test'
#threshold = 1
#label_status = 'True'

# Run functions
print("[+] Running map_fusion_fissions.py with a threshold of " + str(threshold))
print("\t[+] Parsing chromosome assignment files")
file_list = parse_info(input_data)
print("\t[+] Parsing the tree")
t, tip_list = parse_tree(tree_file, label_status)
print("\t[+] Summarising status of each chromosome")
df_combined, spp_list = gather_stats_and_make_table(file_list, input_data)
if len(spp_list) != len(tip_list):
	sys.exit("\t ERROR: Not all species in tree in are in the tsv set or vice versa.")
List_unique_fusions, list_of_fusion_objects = assign_spp_to_fusions(df_combined)
print("\t[+] Identifying unique fusion and fission events and assiging to species")
List_unique_fusions.sort(key=sort_list_high2low)   
List_unique_splits, list_of_split_objects = assign_spp_to_splits(df_combined)
print("\t[+] Mapping fusions and splits onto the tree")	
list_of_split_objects = map_splits(list_of_split_objects, t, threshold)
list_of_fusion_objects = map_fusions(list_of_fusion_objects, t, threshold)
fusions_yet_dealt_with = find_remaining_fusions_to_map(list_of_fusion_objects) # has Thera & Chloro as expected
list_of_fusion_objects = map_subsets_of_fusions_where_possible(fusions_yet_dealt_with, list_of_fusion_objects, t, get_possible_subsets, 1)
ambiguous_fusions = find_remaining_fusions_to_map(list_of_fusion_objects) # re-run function
assert len(fusions_yet_dealt_with) >= len(ambiguous_fusions) # sanity check
ambiguous_fusions, list_of_ambiguous_fusion_objects = map_ambiguous_fusions(ambiguous_fusions, t)
# resolved_fusions = set(fusions_yet_dealt_with).difference(set(ambiguous_fusions))
print("\t[+] Writing results to " + str(output_location) + "/mapped_fusions_fissions_" + prefix + ".tsv")
Directory = str(output_location) 
if not os.path.exists(Directory):
  os.mkdir(Directory)
write_results(output_location, prefix, list_of_fusion_objects, list_of_ambiguous_fusion_objects, list_of_split_objects, write_mapped_events)

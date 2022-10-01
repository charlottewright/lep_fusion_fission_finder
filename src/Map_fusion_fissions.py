#!/usr/bin/python
#%%
# Code taken from Map_fusion_splits_v2.6.py (dated 27/8/21) on 1/10/22 
#from lustre.scratch116.tol.teams.team301.projects.lepidoptera_genomics.cw22.Datafreeze_080621.Analysis.Features.Ancestral_assignments.lep_fusion_split_finder_v1.Split_chr.Map_all_splits import List_to_plot
import os
import glob
import pandas as pd
import sys
from ete3 import Tree, TreeStyle, AttrFace, faces
import re
from pathlib import Path
#conda install -c conda-forge biopython
import argparse


parser = argparse.ArgumentParser(description='This maps each fusion and split event onto the phylogeny')
parser.add_argument("--output", help="This is the location for the output file")
parser.add_argument("--prefix", help="This is the prefix to all files")

args = parser.parse_args()

output_location = args.output
prefix = args.prefix

#%%
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# GET LIST OF FUSIONS READY
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#df_combined = pd.read_csv('results/Overall_fusions.tsv', sep=',')
#print(df_combined)

def parse_info(path):
        file_list = []
        for filename in Path(path).glob('*.tsv'):
                file_list.append(str(filename))
        return file_list

def is_non_zero_file(fpath):
        return os.path.isfile(fpath) and os.path.getsize(fpath) > 0

PATH= "/lustre/scratch116/tol/teams/team301/projects/lepidoptera_genomics/cw22/Datafreeze_080621/Analysis/Features/Ancestral_assignments/lep_fusion_split_finder_v2/"
file_list = parse_info(PATH)

#%%
# Gather stats
max_number_fused_chrs = 0
total_number_fusions = 0
total_number_splits = 0

chr_list = []
d = []

for FILE in file_list:
        result = is_non_zero_file(FILE)
#       print(result)
        if str(result) == "True":
                file = pd.read_csv(FILE, sep='\t', index_col = False)
#               print(file)
                # Header = "query_chr", "status", "assigned_ref_chr"
                Spp_name = FILE.split('/')
                Spp_name = Spp_name[14]
                Spp_name = Spp_name.split('_')
                Genus = Spp_name[0]
                Species = Spp_name[1]
                Full_name = Genus + '_' + Species
                chr_list = []
                for value in file.query_chr:
                        chr_list.append(value)
                for QUERY in chr_list:
                        filtered_data = file[file['query_chr'] == QUERY]
                        STATUS = str(filtered_data.iloc[0].status[::])
                        ASSIGNED = str(filtered_data.iloc[0].assigned_ref_chr[::])
                        List_of_fusions = ASSIGNED.split(',')
                        List_of_fusions = sorted(List_of_fusions)
                        if len(List_of_fusions) > max_number_fused_chrs:
                                max_number_fused_chrs = len(List_of_fusions)
                        if STATUS == "fusion":
                                total_number_fusions = total_number_fusions + 1
                        if STATUS == "split":
                                total_number_splits = total_number_splits + 1
                        entry = {'Chrom_ID':QUERY, 'Status':STATUS, 'Assigned_ref_chr':ASSIGNED, 'Spp':Full_name }
                        d.append(entry)

df_combined = pd.DataFrame(d)

List_unique_fusions = []
List_unique_splits = []

filtered_fusions = df_combined[df_combined['Status'] == "fusion"]
filtered_splits = df_combined[df_combined['Status'] == "split"]
unique_combos_fusions = filtered_fusions.Assigned_ref_chr.astype('str').unique()
unique_merians_splits = filtered_splits.Assigned_ref_chr.astype('str').unique()

# FOR FUSIONS

threshold = 1
for i in unique_combos_fusions:
 #       print("The first unique Merian combo is", i)
  #      print(type(i))
#		value = str(i)
	spp_with_fusion = []
	for index, row in filtered_fusions.iterrows():
			Merians = str(row['Assigned_ref_chr'])
			Merians =  Merians.split(',')
			value = i.split(",")
			matched_merians = set(Merians)&set(value)
		#	print("The merians in the row are:", Merians, "The matched merians are:", matched_merians, "The total merians in the unique combo are:", value)
			number_matched = len(matched_merians)
			number_in_unique_combo = len(value)
			total_number_merians = len(Merians)
			percent_match = number_matched / total_number_merians
			if number_matched == number_in_unique_combo:
			#if percent_match >= threshold:
				Species = str(row['Spp'])
				spp_with_fusion.append(Species)
	entry = {'Unique_Merian_combo':i, 'spp_with_fusion':spp_with_fusion}
	List_unique_fusions.append(entry)

# Sort list_unique_fusions based on number of merian units per fusion (high to low)
def sort_list_high2low(d):
    return len(d['Unique_Merian_combo'])

List_unique_fusions.sort(key=sort_list_high2low)


# FOR SPLITS
for i in unique_merians_splits:
#        print("The first unique Merian is", i)
#        print(type(i))
        i = str(i)
        subset = filtered_splits[filtered_splits['Assigned_ref_chr'] == i]
     #   print(subset)
        spp_with_split = subset["Spp"].to_list()
      #  print("Spp with split are:", spp_with_split)
        entry = {'Unique_Merian_split':i, 'spp_with_split':spp_with_split}
        List_unique_splits.append(entry)

#print("The list of unique fusions is:", List_unique_fusions)
#print("The list of unique splits is:", List_unique_splits)

#%%
# TEMPORARY
print(List_unique_fusions)
#%%
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# GET TREE READY
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Format_1 means node names are present and read - we have node names thanks to the replacement of branch support values with unique numbers done above :)
t = Tree("supermatrix_120721.treefile", format=1)

# Root the tree then label the internal nodes
t.set_outgroup(t&"Hydropsyche_tenuis")

count = 0
for node in t.traverse("preorder"):
	if node.is_leaf():
		continue
	else:
		count = count + 1
		node.name = str(count)

#%%
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# LETS GO: MAP FUSIONS
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def search(Node, mapped_fusions_dict):
    return [element for element in mapped_fusions_dict if element['Node'] == Node]

mapped_fusions_dict = []
tracking_fusion_number_dict = []
# Added on 18/8:
lost_fusions_dict = []
nodes_with_fusions_list = {}
total_tally_fusions = 0
# We can also use the add_feature() method to dynamically add new features.
#label = input("Fusion")
#	Need to make this as a modifiable number
threshold = 0.75

for key in List_unique_fusions:
	x = key["Unique_Merian_combo"]
#	print(x)
	Fusion_event = str(x)
	spp_list = key['spp_with_fusion']
#	print(isinstance(spp_list, list))
	for node in t.traverse("preorder"):
		tip_list = []
		for leaf in node:
			tip = leaf.name
			tip_list.append(tip)
	# Get descendents
#		if sorted(tip_list) == sorted(spp_list):
#                        print("Success!")
# New code here - parametre:
#		if all(x in spp_list for x in tip_list):
		total_spp_list = len(spp_list)
	#		print("The spp list is:", spp_list)
	#		print("The total number spp is:", total_spp_list)
	#		print("The tips are:", tip_list)
		matches = set(spp_list)&set(tip_list)
		spp_with_loss = [d for d in tip_list if d not in spp_list]
	#		print("The matches are:", matches)
		number_matches = len(matches)
		number_tips = len(tip_list)
		proportion_match = number_matches / number_tips
	#		print("The number of matches is:", number_matches)
		if proportion_match >= threshold:
			#print("The node that matches is:", node.name)
		#	print("The matched list is:", matches, "while the total tips at the node is:", tip_list)
			entry = {'Fusion':Fusion_event, 'Tips':str(matches), 'Node':node.name}
			mapped_fusions_dict.append(entry)
			if len(spp_with_loss) != 0:
				loss_entry = {'Fusion':Fusion_event, 'Tips':spp_with_loss, 'Node':node.name}
				lost_fusions_dict.append(loss_entry)
				spp1 = spp_with_loss[0]
				spp1 = str(spp1)
				# the shorcut to finding nodes by name
				spp1 = t&spp1
				spp1.add_feature("Loss1", Fusion_event)
			if len(spp_with_loss) == 2:
				#loss_entry = {'Fusion':Fusion_event, 'Tips':str(spp_with_loss), 'Node':node.name}
				#lost_fusions_dict.append(loss_entry)
				spp1 = spp_with_loss[0]
				spp2 = spp_with_loss[1]
				spp1 = str(spp1)
				spp2 = str(spp2)
				spp1 = t&spp1
				spp2 = t&spp2
				spp1.add_feature("Loss1", Fusion_event)
				spp2.add_feature("Loss1", Fusion_event)
			if len(spp_with_loss) == 3:
				#loss_entry = {'Fusion':Fusion_event, 'Tips':str(spp_with_loss), 'Node':node.name}
				#lost_fusions_dict.append(loss_entry)
				spp1 = spp_with_loss[0]
				spp2 = spp_with_loss[1]
				spp3 = spp_with_loss[2]
				spp1 = str(spp1)
				spp2 = str(spp2)
				spp3 = str(spp3)
				spp1 = t&spp1
				spp2 = t&spp2
				spp3 = t&spp3
				spp1.add_feature("Loss1", Fusion_event)
				spp2.add_feature("Loss1", Fusion_event)
				spp3.add_feature("Loss1", Fusion_event)
			value = str(Fusion_event)
			NODE = node.name
			NODE = str(NODE)
			total_tally_fusions = total_tally_fusions + 1
			if NODE in nodes_with_fusions_list:
				nodes_with_fusions_list[NODE] += 1
				Fusion_x = nodes_with_fusions_list[NODE]
				Fusion_x = "Fusion" + str(Fusion_x)
				node.add_feature(Fusion_x, value)
				# These lines are all to do with identifying if a node already has a subset of a Merian combo (e.g. M11,M2 when M11,M2,MZ is now being processed)
			#	Fusion_Merians = x.split(",")
			#	Number_Merian_in_Fusion = len(Fusion_Merians)
			#	results = search(node.name, mapped_fusions_dict)
			#	entry2 = {'Fusion':Fusion_event, 'Tips':str(matches), 'Node':node.name, 'Fusion_x':Fusion_x}
			#	tracking_fusion_number_dict.append(entry2)
				# Get values of particular key in list of dictionaries
			#	res = [ sub['Fusion'] for sub in results ]
		#		print("The results are:", res)
			#	for entry in res:
			#		entry_str = entry
			#		dups = set(entry_list)&set(Fusion_Merians)
			#		entry_list = entry.split(",")
			#		number_dups = len(dups)
			#		if number_dups >= 2:
			#			if Number_Merian_in_Fusion > number_dups:
			#				if all(x in Fusion_Merians for x in entry_list):
						# Find out what "Fusion_x" the smaller entry is associated with
						# By searching  tracking_fusion_number_dict for entry
			#					results2 = search(node.name, tracking_fusion_number_dict)
			#					corresponding_fusion_number = [d['Fusion_x'] for d in results2 if d['Fusion'] == entry]
			#					corresponding_fusion_number = str(corresponding_fusion_number)
			#					corresponding_fusion_number = re.sub('[\W_]+', '', corresponding_fusion_number)
			#					node.del_feature(corresponding_fusion_number)
			#					for i in range(len(mapped_fusions_dict)):
			#						if mapped_fusions_dict[i]['Node'] == node.name and mapped_fusions_dict[i]['Fusion'] == entry_list:
			#					if node.name == "Pieris_napi":
			#							del(mapped_fusions_dict[i])
			#						print("The node is:", node.name, "The Fusion_Merians are:", Fusion_Merians, "The entry list is", entry_list, "The corresponding fusion number is:", corresponding_fusion_number)
								# This line added 23/8 updates the mapped_fusions_dict to remove sub-fusions
							#	mapped_fusions_dict =  [item for item in mapped_fusions_dict if item['Node'] != node.name and item['Fusion'] != str(entry)]
			else:
				node.add_feature("Fusion1", value)
				nodes_with_fusions_list[NODE] = 1
#				node.add_feature("Fusion_2", value)
#				print("THE NEW LABEL IS:", node.Fusion_1)
			for element in matches:
				if element in spp_list:
					spp_list.remove(element)
#			print("After removing these tips, the total list now is:", spp_list)
#print("The final mapped fusions are:", mapped_fusions_dict)
# So far the max number of (sensible) fusions per node is 8
print(t.get_ascii(show_internal=True, attributes = ["name", 'Fusion1', "Fusion2", "Fusion3", "Fusion4", "Fusion5", "Fusion6", "Fusion7", "Fusion8", "Fusion9", "Fusion10", "Fusion11", "Fusion12", "Fusion13", "Fusion14", "Fusion15", "Fusion16", "Fusion17", "Fusion18", "Fusion19", "Fusion20", "Fusion21", "Fusion22", "Fusion23", "Fusion24", "Fusion25", "Fusion26", "Fusion27", "Fusion28", "Fusion29", "Fusion30", "Fusion31", "Fusion32", "Fusion33", "Fusion34", "Fusion35", "Fusion36", "Fusion37", "Fusion38", "Fusion39", "Fusion40", "Fusion41"]))

#%%
#print(lost_fusions_dict)
print(mapped_fusions_dict)

#%%
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# LETS GO - MAP SPLITS
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
mapped_splits_dict = []
lost_splits_dict = []
nodes_with_splits_list = {}
total_tally_splits = 0
#	Need to make this as a modifiable number
threshold = 0.75

for key in List_unique_splits:
	x = key["Unique_Merian_split"]
#	print(x)
	Split_event = str(x)
	spp_list = key['spp_with_split']
#	print(isinstance(spp_list, list))
	for node in t.traverse("preorder"):
#	print("The node name is:")
#	print(node.name)
		tip_list = []
#	DESCENDENTS = t.get_children()
#	print("The descendents of:", node.name, "are", DESCENDENTS)
		for leaf in node:
			tip = leaf.name
			tip_list.append(tip)
#	print("The leaf nodes are:")
#	print(tip_list)
#		if sorted(tip_list) == sorted(spp_list):
#			print("Success!")
#		if all(x in spp_list for x in tip_list):
		total_spp_list = len(spp_list)
	#		print("The spp list is:", spp_list)
	#		print("The total number spp is:", total_spp_list)
	#		print("The tips are:", tip_list)
		matches = set(spp_list)&set(tip_list)
		spp_with_loss = [d for d in tip_list if d not in spp_list]
	#		print("The matches are:", matches)
		number_matches = len(matches)
		number_tips = len(tip_list)
		proportion_match = number_matches / number_tips
	#		print("The number of matches is:", number_matches)
		if proportion_match >= threshold:
#			print("The node that matches is:", node.name)
#			print("The matched list is:", matches, "while the total tips at the node is:", tip_list)
			entry = {'Split':Split_event, 'Tips':str(matches), 'Node':node.name}
			mapped_splits_dict.append(entry)
			if len(spp_with_loss) != 0:
				loss_entry = {'Split':Split_event, 'Tips':str(spp_with_loss), 'Node':node.name}
				lost_splits_dict.append(loss_entry)
			value = str(Split_event)
			NODE = node.name
			NODE = str(NODE)
			total_tally_splits = total_tally_splits + 1
			if NODE in nodes_with_splits_list:
				nodes_with_splits_list[NODE] += 1
				Split_x = nodes_with_splits_list[NODE]
				Split_x = "Split" + str(Split_x)
				node.add_feature(Split_x, value)
			else:
				node.add_feature("Split1", value)
				nodes_with_splits_list[NODE] = 1
			print("THE NEW LABEL IS:", node.Split1)
			for element in tip_list:
				if element in spp_list:
					spp_list.remove(element)
			print("After removing these tips, the total list now is:", spp_list)
print("The final mapped splits are:", mapped_splits_dict)
# So far the max number of (sensible) fusions per node is 8
print(t.get_ascii(show_internal=True, attributes = ["name", "Split1", "Split2", "Split3", "Split4", "Split5", "Split6", "Split7", "Split8", "Split9", "Split10", "Split11", "Split12", "Split13", "Split14", "Split15"]))
#%%
# Convert lost splits & lost fusios to dataframes
lost_fusions_df = pd.DataFrame(lost_fusions_dict)
lost_splits_df = pd.DataFrame(lost_splits_dict)

#%%
lost_fusions_df
#%%
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# PREPARE LIST OF UNIQUE "LOSSES"
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

unique_combos_losses= lost_fusions_df.Fusion.astype('str').unique()
List_unique_lost_fusions = []
mapped_lost_fusions_dict = []
nodes_with_lost_fusions_list = {}

threshold = 1
for i in unique_combos_losses:
 #       print("The first unique Merian combo is", i)
  #      print(type(i))
#		value = str(i)
	spp_with_lost_fusion = []
	for index, row in lost_fusions_df.iterrows():
			Merians = str(row['Fusion'])
			Merians =  Merians.split(',')
			value = i.split(",")
			matched_merians = set(Merians)&set(value)
		#	print("The merians in the row are:", Merians, "The matched merians are:", matched_merians, "The total merians in the unique combo are:", value)
			number_matched = len(matched_merians)
			number_in_unique_combo = len(value)
			total_number_merians = len(Merians)
			percent_match = number_matched / total_number_merians
			if number_matched == number_in_unique_combo:
				print("SUCESS", Merians)
			#if percent_match >= threshold:
				Species = str(row['Tips'])
				print(Species)
				spp_with_lost_fusion.append(Species)
	entry = {'Unique_Merian_combo':i, 'Spp_with_lost_fusion':spp_with_lost_fusion}
	List_unique_lost_fusions.append(entry)
#%%
# Now map the lost_fusions
for key in List_unique_lost_fusions:
	x = key["Unique_Merian_combo"]
#	print(x)
	Lost_fusion = str(x)
	spp_list = str(key['Spp_with_lost_fusion'])
	print(spp_list)
	spp_list = re.sub('[!*)@#%(&$?\'\".\]\[^]', '', spp_list)
	print(spp_list)
	spp_list = spp_list.split(",")
	print(spp_list[0])
	#spp_list = str(spp_list)
	#spp_list = spp_list.split(",")
#	check_list = isinstance(spp_list, list)
	print(spp_list)
#	print(isinstance(spp_list, list))
	for node in t.traverse("preorder"):
		tip_list = []
		for leaf in node:
			tip = leaf.name
			tip_list.append(tip)
		if all(x in spp_list for x in tip_list):
			#print("The node that matches is:", node.name)
		#	print("The matched list is:", matches, "while the total tips at the node is:", tip_list)
			entry = {'Lost_fusion':Lost_fusion, 'Tips':str(tip_list), 'Node':node.name}
			mapped_lost_fusions_dict.append(entry)
			value = str(Lost_fusion)
			NODE = node.name
			NODE = str(NODE)
			if NODE in nodes_with_lost_fusions_list:
				nodes_with_lost_fusions_list[NODE] += 1
				Loss_x = nodes_with_lost_fusions_list[NODE]
				Loss_x = "Loss" + str(Loss_x)
				node.add_feature(Loss_x, value)
			else:
				node.add_feature("Loss1", value)
				nodes_with_lost_fusions_list[NODE] = 1
			for element in matches:
				if element in spp_list:
					spp_list.remove(element)
#			print("After removing these tips, the total list now is:", spp_list)
#print("The final mapped fusions are:", mapped_fusions_dict)
#print(t.get_ascii(show_internal=True, attributes = ["name", "Loss1", 'Loss2', "Loss3", "Loss4", "Loss5", "Loss6", "Loss7", "Loss8", "Loss9", "Loss10", "Loss11", "Loss12", "Loss13", "Loss14", "Loss14", "Loss15"]))
#%%

print(t.get_ascii(show_internal=True, attributes = ["name", "Loss1", 'Loss2', "Loss3", "Loss4", "Loss5", "Loss6", "Loss7", "Loss8", "Loss9", "Loss10", "Loss11", "Loss12", "Loss13", "Loss14", "Loss14", "Loss15"]))

#%%	
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Render tree with mapped fusions
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# layouts are functions that allow you to make modifications on each node before they've been drawn. 
# They accept a given node as input and allow you to decide what to do with that node 

def split_layout(node):
	for i in range(50):
		x = str(i)
		Fusion_x = ("Fusion" + x)
		Fx = ("F" + x)
		if Fusion_x in node.features:
			Fx = AttrFace(Fusion_x, fgcolor="blue")
			faces.add_face_to_node(Fx, node, column=0, position="branch-top")
	for i in range(20):
		x = str(i)
		Split_x = ("Split" + x)
		Sx = ("s" + x)
	#	print(Split_x)
		if Split_x in node.features:
			Sx = AttrFace(Split_x, fgcolor="red")
			faces.add_face_to_node(Sx, node, column=0, position="branch-top")
	if node.is_leaf():
		name_face = AttrFace("name")
	else:
		# If internal noded, draws label with smaller font size
		name_face = AttrFace("name", fsize=10)
# Adds the name face to the image at the preferred position
		faces.add_face_to_node(name_face, node, column=0, position="branch-right")

ts = TreeStyle()
# Use my custom layout
ts.layout_fn = split_layout
# Do not add leaf names automatically
ts.show_leaf_name = False
#ts.show_branch_length = True
name = str(prefix) + "_Supermatrix_tree_annotated_splits_v2.6.png"
t.render(name, w=183, units="mm", tree_style=ts)
#t.show(tree_style=ts)

#%%
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Save your results
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Directory = str(output_location) 
if not os.path.exists(Directory):
	os.mkdir(Directory)
print("Directory '% s' created" % Directory)

# Save "overall assignments" as a csv
output = str(output_location) + '/' + str(prefix) + '_Overall_assignments_v2.6.tsv'
df_combined.to_csv(output, header=True, index=False, sep = '\t')

# Save "mapped fusions" as a csv
mapped_fusions_df = pd.DataFrame(mapped_fusions_dict)
output = str(output_location)  + '/' + str(prefix) +'_Mapped_fusions_v2.6.tsv'
mapped_fusions_df.to_csv(output, header=True, index = False, sep = '\t')

# Save "lost fusions" as a csv - list of spp that lost a fusion
lost_fusions_df = pd.DataFrame(lost_fusions_dict)
output = str(output_location)  + '/' + str(prefix) + '_Lost_fusions_v2.6.tsv'
lost_fusions_df.to_csv(output, header=True, index = False, sep = '\t')

# Save "lost splits" as a csv - list of spp that lost a split
lost_splits_df = pd.DataFrame(lost_splits_dict)
output = str(output_location)  + '/' + str(prefix) + '_Lost_splits_v2.6.tsv'
lost_splits_df.to_csv(output, header=True, index = False, sep = '\t')

# Save "mapped splits" as a csv
mapped_splits_df = pd.DataFrame(mapped_splits_dict)
output = str(output_location)  + '/' + str(prefix) + '_Mapped_splits_v2.6.tsv'
mapped_splits_df.to_csv(output, header=True, index = False, sep = '\t')

# Save annotated trees as nwk files
fusion_tree = str(prefix) + "_annotated_fusions_tree_v2.6.nw"
split_tree = str(prefix) + "_annotated_splits_tree_v2.6.nw"
t.write(features=["name", "Fusion1", "Fusion2", "Fusion3", "Fusion4", "Fusion5", "Fusion6", "Fusion7", "Fusion8", "Fusion9", "Fusion10", "Fusion11"], outfile=fusion_tree)
t.write(features=["name", "Split1", "Split2", "Split3", "Split4", "Split5", "Split6", "Split7", "Split8", "Split9", "Split10", "Split11", "Split12", "Split13", "Split14", "Split15"], outfile=split_tree)
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  1. For each node in tree from root node:
#                a.  Get descendants
#                b. If all descendants in trait list - then print node ID and delete descendents from list of spp that have the trait
#                c. Else pass
#Once list is empty - stop

exit ()

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Check number of fusions & splits that map to each node
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
FusionsCount = {}
for key in mapped_fusions_dict:
	Node = key["Node"]
	Node = str(Node)
	if Node in FusionsCount:
		FusionsCount[Node] += 1
	else:
		FusionsCount[Node] = 1
#%%
SplitsCount = {}
for key in mapped_splits_dict:
	Node = key["Node"]
	Node = str(Node)
	if Node in SplitsCount:
		SplitsCount[Node] += 1
	else:
		SplitsCount[Node] = 1

print("Total number of mapped fusions per spp:", FusionsCount)
print("total number of mapped fusion events is:", total_tally_fusions)
print("Total number of mapped splits per spp:", SplitsCount)
print("total number of mapped split events is:", total_tally_splits)

print("Total number of fusions:", total_number_fusions)
print("Total number of splits:", total_number_splits)
print("The maximum number of fused chrs is:", max_number_fused_chrs)

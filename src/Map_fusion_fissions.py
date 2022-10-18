#!/usr/bin/env python3
#%%
# Code taken from Map_fusion_splits_v2.6.py (dated 27/8/21) on 1/10/22 
import os, sys, glob, argparse, re
import pandas as pd
from ete3 import Tree, TreeStyle, AttrFace, faces
from pathlib import Path
import itertools
#%%
### Define functions ###
def parse_info(input_data):  # get list of fusions ready
	file_list = []
	for filename in Path(input_data).glob('*.tsv'):
		file_list.append(str(filename))
	return file_list

def is_non_zero_file(fpath):
        return os.path.isfile(fpath) and os.path.getsize(fpath) > 0

def gather_stats_and_make_table(file_list):
	max_number_fused_chrs, total_number_fusions, total_number_splits = 0, 0, 0
	chr_list, d = [], []
	for file_name in file_list:
		result = is_non_zero_file(file_name)
		if str(result) == "True":
			file = pd.read_csv(file_name, sep='\t', index_col = False) # Header = "query_chr", "status", "assigned_ref_chr"
			Spp_name = file_name.split('/')
			Spp_name = Spp_name[int(input_data.count('/'))] # altered depending on file path to get just file name
			Spp_name = Spp_name.replace('_chromosome_assignments.tsv','')
			chr_list = []
			for value in file.query_chr:
					chr_list.append(value)
			for query_chr in chr_list:
					filtered_data = file[file['query_chr'] == query_chr]
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
					entry = {'Chrom_ID':query_chr, 'Status':STATUS, 'Assigned_ref_chr':ASSIGNED, 'Spp':Spp_name }
					d.append(entry)
	df_combined = pd.DataFrame(d)
	return(df_combined)

def get_fusions_splits_and_unique_events(df_combined):
	filtered_fusions = df_combined[df_combined['Status'] == "fusion"]
	filtered_splits = df_combined[df_combined['Status'] == "split"]
	unique_combos_fusions = filtered_fusions.Assigned_ref_chr.astype('str').unique()
	unique_merians_splits = filtered_splits.Assigned_ref_chr.astype('str').unique()
	return(filtered_fusions, filtered_splits, unique_combos_fusions, unique_merians_splits)

def assign_spp_to_fusions(unique_combos_fusions):
	List_unique_fusions = []
	for i in unique_combos_fusions: # for each unique combination of Merians (query)
		value, spp_with_fusion = str(i), []
		for index, row in filtered_fusions.iterrows():
			Merians = str(row['Assigned_ref_chr']).split(',')
			value = i.split(",")
			matched_merians = set(Merians)&set(value)
			if (len(matched_merians) == len(Merians)) & (len(matched_merians) == len(value)): # if all Merians in query fusion are in a fusion chr of a given spp (i.e. not a subset, all Merians in fusion are there)
				#print("The merians in the row are:", Merians, "The matched merians are:", matched_merians, "The total merians in the unique combo are:", value)
				Species = str(row['Spp']) # i.e. if query fusion is M2,M4, don't matcch M2,4,M6 in spp1
				spp_with_fusion.append(Species)
		entry = {'Unique_Merian_combo':i, 'spp_with_fusion':spp_with_fusion}
		List_unique_fusions.append(entry)
	return(List_unique_fusions)

def sort_list_high2low(d): # # Sort list_unique_fusions based on number of merian units per fusion (high to low)
    return len(d['Unique_Merian_combo'])

def assign_spp_to_splits(unique_merians_splits):
	List_unique_splits = []
	for i in unique_merians_splits: # for each unique Merian 
			subset = filtered_splits[filtered_splits['Assigned_ref_chr'] == str(i)]
			spp_with_split = subset["Spp"].to_list() # get the spp with the split
			entry = {'Unique_Merian_split':i, 'spp_with_split':spp_with_split}
			List_unique_splits.append(entry)
	return(List_unique_splits)

def parse_tree(tree_file, label_status):
	t = Tree(tree_file, format=1) # Format_1 means node names are present and read - we have node names thanks to the replacement of branch support values with unique numbers done above :)
	if label_status == 'False':
		count = -1 # Label each internal node with a number starting from 0
		for node in t.traverse("preorder"):
			if node.is_leaf():
				continue
			else:
				count = count + 1
				node.name = str('n') + str(count)
	return(t)

def search(Node, mapped_fusions_dict): # I don't think this function is ever used???
    return [element for element in mapped_fusions_dict if element['Node'] == Node]

def map_fusions(List_unique_fusions, t, threshold):
	copy_list_unique_fusions = List_unique_fusions.copy()
	mapped_subsets = []
	mapped_fusions_dict = []
	lost_fusions_dict = []
	for key in List_unique_fusions:
		Fusion_event = key["Unique_Merian_combo"].split(',')
		spp_list = list(key['spp_with_fusion'])
		for node in t.traverse("preorder"):
			temp_spp_list = list(spp_list) # need to make a copy of list like so to prevent lists being linked
			tip_list = []
			new_spp = 0 # initalise at zero for each node
			new_spp_details = {}
			for leaf in node:
				tip = leaf.name
				tip_list.append(tip)
			matches = sorted(set(spp_list)&set(tip_list))
			if len(matches) >= 1: # if at least one matched merian..
				# check if any other spp has this fusion as part of a bigger fusion (involving more merians)
				for entry in List_unique_fusions:
					unique_ms = entry['Unique_Merian_combo'].split(',')
					matched_ms = set(unique_ms) & set(Fusion_event)
					if (len(unique_ms) > len(Fusion_event)) & (len(matched_ms) == len(Fusion_event)): # if all Merians in query fusion event are in this fusion, and this fusion has more merians than the query (i.e. the query is a subset of it). Need the greater than or else spp would already be presen in spp list.
					#	print(Fusion_event, len(unique_ms), len(Fusion_event), len(matched_ms))
						spp_matched = entry['spp_with_fusion']
						new_spp = set(spp_matched).difference(spp_list)
						new_spp = list(new_spp) # convert to list
						if len(new_spp) >= 1: # if at least 1 new spp
							for spp in new_spp: # add new spp to spp_list if found at this node
								if spp in tip_list:
									temp_spp_list.append(spp)
									new_spp_details[spp] = unique_ms
							matches = sorted(set(temp_spp_list)&set(tip_list)) # Now update matches again..
				spp_with_loss = [d for d in tip_list if d not in temp_spp_list]
				proportion_match = len(matches) / len(tip_list)
				if proportion_match >= threshold:
					matched = re.sub(r"[\{\}\[\]']", '', str(matches))
					fusion_merians = re.sub(r"[\{\}\[\]']", '', str(Fusion_event))
					fusion_merians = fusion_merians.replace(' ','')
					entry = {'Merians':fusion_merians, 'Tips':str(matched), 'Node':node.name}
					mapped_fusions_dict.append(entry)
					if len(spp_with_loss) != 0:
						loss_entry = {'Merians':Fusion_event, 'Tips':spp_with_loss, 'Node':node.name}
						lost_fusions_dict.append(loss_entry)
					for element in matches:
						if element in spp_list:
							spp_list.remove(element)
					if len(new_spp_details) >= 1: # if at least 1 new spp
						for spp in new_spp_details:
							spp_as_list = []
							unique_ms = new_spp_details[spp]
							reform_unique_ms = re.sub(r"[\{\}\[\]']", '', str(unique_ms))
							reform_unique_ms = reform_unique_ms.replace(' ','')
							spp_as_list.append(spp)
							entry = {'Unique_Merian_combo':reform_unique_ms, 'spp_with_fusion':spp_as_list}
							mapped_subsets.append(entry)
							if len(unique_ms) >3: # i.e. if there are still parts of the larger fusion to map
								remaining_merians = sorted(set(unique_ms).difference(matched_ms))
								remaining_merians = re.sub(r"[\{\}\[\]']", '', str(remaining_merians))
								remaining_merians = remaining_merians.replace(' ','')
								reform_spp = []
								reform_spp.append(spp)
								if not any(d['Unique_Merian_combo'] == remaining_merians for d in copy_list_unique_fusions):
									entry = {'Unique_Merian_combo':remaining_merians, 'spp_with_fusion':reform_spp}
									copy_list_unique_fusions.append(entry)
									print(entry)
								else:
									for d in copy_list_unique_fusions: # find right dict within list
										current_spp = d['spp_with_fusion']
										merians = d['Unique_Merian_combo']
										if merians == remaining_merians:
											current_plus_new_spp = current_spp.append(spp)
											d.update((k, current_plus_new_spp) for k, v in d.items() if k == remaining_merians)
	lost_fusions_df = pd.DataFrame(lost_fusions_dict) # Convert lost fusions to dataframe
	return(mapped_fusions_dict, lost_fusions_df, t, mapped_subsets, copy_list_unique_fusions)

def combinantorial(lst):
    index, pairs = 1, []
    for element1 in lst:
        for element2 in lst[index:]:
            pairs.append((element1, element2))
        index += 1
    return pairs


def account_for_subsets_of_fusions(copy_list_unique_fusions, t, mapped_subsets): # this function adds in fusions which are subsets of larger fusions (e.g. sp1 may have M1,M2,M3 and sp2 may have M1,M2,M4 --> M1,M2 shared event)
	copy_list_unique_fusions = [x for x in copy_list_unique_fusions if x not in mapped_subsets] # remove fusions involving >2 Merians that have subsets that have already been mapped so that don't remap subsets here
	for key in copy_list_unique_fusions:
		Fusion_event = key["Unique_Merian_combo"].split(',')
		spp_list = key['spp_with_fusion']
		if len(Fusion_event) > 2:
			spp_remaining = list(spp_list)
			mapped_subsets_counter = 0
			fusion_pairs = combinantorial(Fusion_event) # make tuple of all possible pairs of merians to fuse
			for Subset_fusion_event in fusion_pairs:
		#		if Subset_fusion_event in previously_mapped_fusions:
		#			print(Subset_fusion_event, 'yas')
				Subset_fusion_event =  list(Subset_fusion_event) # convert to list so same format as 'Fusion_event'
				for node in t.traverse("preorder"):
					temp_spp_list = list(spp_remaining) # need to make a copy of list like so to prevent lists being linked
					total_new_spp, tip_list = 0, []
					for leaf in node:
						tip = leaf.name
						tip_list.append(tip)
					matches = sorted(set(spp_list)&set(tip_list))
					if len(matches) >= 1: # if at least one matched merian..
						for entry in copy_list_unique_fusions: # check if any other spp has this fusion as part of a bigger fusion (involving more merians)
							unique_ms = entry['Unique_Merian_combo'].split(',')
							matched_ms = set(unique_ms) & set(Subset_fusion_event)
							if (len(unique_ms) > len(Subset_fusion_event)) & (len(matched_ms) == len(Subset_fusion_event)): # if all Merians in query fusion event are in this fusion, and this fusion has more merians than the query (i.e. the query is a subset of it). Need the greater than or else spp would already be presen in spp list.
								spp_matched = entry['spp_with_fusion']
								new_spp = set(spp_matched).difference(list(spp_list))
								new_spp = list(new_spp) # convert to list
								if len(new_spp) >= 1: # if at least 1 new spp
									total_new_spp += 1
									for spp in new_spp:  # add new spp to spp_list if found at this node
										if spp in tip_list:
											temp_spp_list.append(spp)
									matches = sorted(set(temp_spp_list)&set(tip_list)) # Now update matches again..
									proportion_match = len(matches) / len(tip_list) # only want to do following code if have a new spp
									if (proportion_match >= threshold) & (len(matches) >= 2):
										matched = re.sub(r"[\{\}\[\]']", '', str(matches))
										fusion_merians = re.sub(r"[\{\}\[\]']", '', str(Subset_fusion_event))
										fusion_merians = fusion_merians.replace(' ','')
										entry = {'Merians':fusion_merians, 'Tips':str(matched), 'Node':node.name}
										for element in matches:
											if element in spp_remaining:
												spp_remaining.remove(element)
										if entry not in mapped_fusions_dict:
											mapped_fusions_dict.append(entry)
											mapped_subsets_counter =+1
										else:
											mapped_subsets_counter =+ 1
			if mapped_subsets_counter == 0:  # i.e. no subsets exist at any node, means just have [A,B,C] fusion. So need to add [A,B | A,C | A,B] to dict. But first need to know node where this happened e..g. two species may have [A,B,C]
				spp_remaining = list(spp_list)  # need to make a copy of list like so to prevent lists being linked
				print(Fusion_event, spp_list)
				for node in t.traverse("preorder"): # so re-traverse the tree..
					temp_spp_list = list(spp_list) 
					tip_list = []
					for leaf in node:
						tip = leaf.name
						tip_list.append(tip)
					matches = sorted(set(spp_remaining)&set(tip_list))
					if len(matches) >= 1: # if at least one matched merian..
						proportion_match = len(matches) / len(tip_list)
						if proportion_match >= threshold:
							matched = re.sub(r"[\{\}\[\]']", '', str(matches))
							fusion_merians = re.sub(r"[\{\}\[\]']", '', str(Fusion_event))
							fusion_merians = fusion_merians.replace(' ','')
							if len(Fusion_event) == 3:
								entry = {'Merians':fusion_pairs, 'Tips':str(matched), 'Node':node.name}
								mapped_fusions_dict.append(entry)
								for element in matches:
									if element in spp_remaining:
										spp_remaining.remove(element)
							else: # if >3 Merians in fusion event then we need to add more than one entry to the list of events
								# e.g. if 4 Merians, then 3 events must have occured. Already have "final" event in list so need to add two "subset" fusions
								events_to_add = len(Fusion_event) -2
								for i in range(1,events_to_add+1, 1):
									entry = {'Merians':fusion_pairs, 'Tips':str(matched), 'Node':node.name}
									mapped_fusions_dict.append(entry)
									for element in matches:
										if element in spp_remaining:
											spp_remaining.remove(element)
	return(mapped_fusions_dict, t)

def map_splits(List_unique_splits, t, threshold):
	mapped_splits_dict = []
	lost_splits_dict = []
	for key in List_unique_splits:
		x = key["Unique_Merian_split"]
		Split_event = str(x)
		spp_list = key['spp_with_split']
		for node in t.traverse("preorder"):
			tip_list = []
			for leaf in node:
				tip = leaf.name
				tip_list.append(tip)
			matches = set(spp_list)&set(tip_list)
			spp_with_loss = [d for d in tip_list if d not in spp_list]
			number_matches = len(matches)
			number_tips = len(tip_list)
			proportion_match = number_matches / number_tips
			if proportion_match >= threshold:
				matched = re.sub(r"[\{\}']", '', str(matches))
				entry = {'Merians':Split_event, 'Tips':str(matched), 'Node':node.name}
				mapped_splits_dict.append(entry)
				if len(spp_with_loss) != 0:
					loss_entry = {'Merians':Split_event, 'Tips':str(spp_with_loss), 'Node':node.name}
					lost_splits_dict.append(loss_entry)
				for element in matches:
					if element in spp_list:
						spp_list.remove(element)
	lost_splits_df = pd.DataFrame(lost_splits_dict) # Convert lost splits to dataframes
	return(mapped_splits_dict, lost_splits_df, t)

def make_list_unique_losses(lost_fusions_df):
	unique_combos_losses= lost_fusions_df.Merians.astype('str').unique()
	List_unique_lost_fusions = []
#	threshold = 1
	for i in unique_combos_losses:
		spp_with_lost_fusion = []
		for index, row in lost_fusions_df.iterrows():
				Merians = str(row['Merians']).split(',')
				value = i.split(",")
				matched_merians = set(Merians)&set(value)
				number_matched = len(matched_merians)
				number_in_unique_combo = len(value)
			#	total_number_merians = len(Merians)
			#	percent_match = number_matched / total_number_merians
				if number_matched == number_in_unique_combo:
				#if percent_match >= threshold:
					Species = str(row['Tips'])
					spp_with_lost_fusion.append(Species)
		entry = {'Unique_Merian_combo':i, 'Spp_with_lost_fusion':spp_with_lost_fusion}
		List_unique_lost_fusions.append(entry)
	return(List_unique_lost_fusions)

def map_lost_fusions(List_unique_lost_fusions, t):
	mapped_lost_fusions_dict = []
	for key in List_unique_lost_fusions:
		x = key["Unique_Merian_combo"]
		Lost_fusion = str(x)
		spp_list = str(key['Spp_with_lost_fusion'])
		spp_list = re.sub('[!*)@#%(&$?\'\".\]\[^]', '', spp_list)
		spp_list = spp_list.split(",")
		for node in t.traverse("preorder"):
			tip_list = []
			for leaf in node:
				tip = leaf.name
				tip_list.append(tip)
			if all(x in spp_list for x in tip_list):
				tips = re.sub(r"[\[\]']", '', str(tip_list))
				entry = {'Merians':Lost_fusion, 'Tips':str(tips), 'Node':node.name}
				mapped_lost_fusions_dict.append(entry)
				for element in tip_list: 
					if element in spp_list:
						spp_list.remove(element)
	return(mapped_lost_fusions_dict, t)

def split_layout(node): # layouts are functions that allow you to make modifications on each node before they've been drawn. 
	for i in range(50): # They accept a given node as input and allow you to decide what to do with that node 
		x = str(i)
		Fusion_x = ("Fusion" + x)
		Fx = ("F" + x)
		if Fusion_x in node.features:
			Fx = AttrFace(Fusion_x, fgcolor="blue")
			faces.add_face_to_node(Fx, node, column=0, position="branch-top")
	for i in range(20):
		Split_x = ("Split" + str(i))
		Sx = ("s" + str(i))
		if Split_x in node.features:
			Sx = AttrFace(Split_x, fgcolor="red")
			faces.add_face_to_node(Sx, node, column=0, position="branch-top")
	if node.is_leaf():
		name_face = AttrFace("name")
	else: 		# If internal noded, draws label with smaller font size
		name_face = AttrFace("name", fsize=10)
		faces.add_face_to_node(name_face, node, column=0, position="branch-right") # Adds the name face to the image at the preferred position

def annotate_tree_with_fusions(mapped_fusions_dict, t):
	nodes_with_fusions_dict = {}
	for entry in mapped_fusions_dict:
		fusion_event, tips, node_name = entry['Merians'], entry['Tips'], entry['Node'] # all fusions, regardless of if its a subset of another fusion or not, are plotted in the same way
		node= t&node_name # the shorcut to finding nodes by name
		if node_name in nodes_with_fusions_dict:
			nodes_with_fusions_dict[node_name] += 1
			Fusion_x = nodes_with_fusions_dict[node_name]
			Fusion_x = "Fusion" + str(Fusion_x)
			node.add_feature(Fusion_x, fusion_event)
		else:
			node.add_feature("Fusion1", fusion_event)
			nodes_with_fusions_dict[node_name] = 1
	return(t)

def annotate_tree_with_lost_fusions(mapped_lost_fusions_dict, t): # only need to annotated tree with lost fusions if there are any
	for entry in mapped_lost_fusions_dict:
		fusion_loss_event, spp_with_loss = entry['Merians'], entry['Tips']
		spp_with_loss = spp_with_loss.split(' ')
		spp1 = str(spp_with_loss[0])
		spp1 = t&spp1 # the shorcut to finding nodes by name
		spp1.add_feature("Loss1", fusion_loss_event)
		if len(spp_with_loss) == 2:
			spp1, spp2 = str(spp_with_loss[0]), str(spp_with_loss[1])
			spp1 = t&spp1
			spp2 = t&spp2
			spp1.add_feature("Loss1", fusion_loss_event)
			spp2.add_feature("Loss1", fusion_loss_event)
		if len(spp_with_loss) == 3:
			spp1, spp2, spp3 = str(spp_with_loss[0]), str(spp_with_loss[1]), str(spp_with_loss[2])
			spp1 = t&spp1
			spp2 = t&spp2
			spp3 = t&spp3
			spp1.add_feature("Loss1", fusion_loss_event)
			spp2.add_feature("Loss1", fusion_loss_event)
			spp3.add_feature("Loss1", fusion_loss_event)
	return(t)

def annotate_tree_with_splits(mapped_splits_dict, t):
	nodes_with_splits_dict = {}
	for entry in mapped_splits_dict:
		split_event, tips, node_name = entry['Merians'], entry['Tips'], entry['Node'] # all fusions, regardless of if its a subset of another fusion or not, are plotted in the same way
		node= t&node_name # the shorcut to finding nodes by name
		if node_name in nodes_with_splits_dict:
			nodes_with_splits_dict[node_name] += 1
			Event_x = nodes_with_splits_dict[node_name]
			Event_x = "Fission" + str(Event_x)
			node.add_feature(Event_x, split_event)
		else:
			node.add_feature("Fission1", split_event)
			nodes_with_splits_dict[node_name] = 1
	return(t)

def write_results(output_location, t, df_combined, mapped_fusions_dict, mapped_splits_dict, lost_fusions_df, lost_splits_df, prefix): # save results
	Directory = str(output_location) 
	if not os.path.exists(Directory):
		os.mkdir(Directory)
#	print("Directory '% s' created" % Directory)
	output_location = str(output_location) + '/'
	mapped_fusions_df = pd.DataFrame(mapped_fusions_dict)
	mapped_splits_df = pd.DataFrame(mapped_splits_dict)
	mapped_fusions_df['Event']= 'fusion'
	mapped_splits_df['Event']= 'fission'
	mapped_events_df = pd.concat([mapped_fusions_df, mapped_splits_df], axis=0)
	df_combined.to_csv(str(output_location) + 'overall_assignments_' + prefix + '.tsv', header=True, index=False, sep = '\t') # Save "overall assignments" as a csv 
	mapped_events_df.to_csv(str(output_location) + 'mapped_fusions_fissions_' + prefix + '.tsv', header=True, index = False, sep = '\t') # Save mapped fusions and fissions as a csv
	if len(lost_fusions_df) != 0:
		lost_fusions_df.to_csv(str(output_location) + 'lost_fusions_' + prefix + '.tsv' , header=True, index = False, sep = '\t') # Save "lost fusions" as a csv - list of spp that lost a fusion
	if len(lost_splits_df) != 0:
		lost_splits_df.to_csv(str(output_location) + 'lost_splits_' + prefix + '.tsv', header=True, index = False, sep = '\t')   # Save "lost splits" as a csv - list of spp that lost a split
	fusion_tree = str(output_location) + 'annotated_fusions_tree_' + str(prefix) + '.nw' # Save annotated trees as nwk files
	fission_tree = str(output_location) + 'annotated_fissions_tree_' + str(prefix) + '.nw'
	t.write(features=["name", "Fusion1", "Fusion2", "Fusion3", "Fusion4", "Fusion5", "Fusion6", "Fusion7", "Fusion8", "Fusion9", "Fusion10", "Fusion11"], outfile=fusion_tree)
	t.write(features=["name", "Fission1", "Fission2", "Fission3", "Fission4", "Fission5", "Fission6", "Fission7", "Fission8"], outfile=fission_tree)
	return

def get_event_stats(mapped_fusions_dict, mapped_splits_dict):
	FusionsCount = {}
	for key in mapped_fusions_dict:
		Node = str(key["Node"])
		if Node in FusionsCount:
			FusionsCount[Node] += 1
		else:
			FusionsCount[Node] = 1
	SplitsCount = {}
	for key in mapped_splits_dict:
		Node = str(key["Node"])
		if Node in SplitsCount:
			SplitsCount[Node] += 1
		else:
			SplitsCount[Node] = 1
	#print("Total number of mapped fusions per node:", FusionsCount)
	print("Total number of mapped fusion events:", len(mapped_fusions_dict))
	print("Total number of mapped splits per node:", SplitsCount)
	print("Total number of mapped split events:", len(mapped_splits_dict))
	print("Total number of fusion chromosomes:", sum(FusionsCount.values()))
	print("Total number of split chromosomes:", sum(SplitsCount.values()))
	print("Maximum number of fused chromosomes:", max(FusionsCount.values()))
	if len(mapped_splits_dict) != 0:
		print("Maximum number of split chromosomes:", max(SplitsCount.values()))
	else:
		print("Maximum number of split chromosomes: 0")
#%%
if __name__ == "__main__":
	SCRIPT = "Map_fusion_fissions.py"
	# argument set up
	parser = argparse.ArgumentParser()
	parser.add_argument("-i", "--input_data", type=str, help = "path to lep_fusion_fission_finder output", required=True)
	parser.add_argument("-tree", "--tree", type=str, help = "Phylogenetic tree", default="fsf")
	parser.add_argument("-o", "--output", type=str, help = "output location relative to working directory", required=True)
	parser.add_argument("-f", "--prefix", type=str, help = "Prefix for all output files", default="fsf")
	parser.add_argument("-t", "--threshold", type=int, help = "Threshold for rearrangement to be shared between tips", default=0.75)
	parser.add_argument("-l", "--label", type=str, help = "Specify if tree already contains internal node labels", default='False')
	args = parser.parse_args()
	input_data = args.input_data
	tree_file = args.tree
	output_location = args.output
	prefix = args.prefix
	threshold = args.threshold
	label_status = args.label

	# Run functions
	print("[+] Running map_fusion_fissions.py with a threshold of " + str(threshold))
	print("\t[+] Parsing chromosome assignment files")
	file_list = parse_info(input_data)
	print("\t[+] Summarising status of each chromosome")
	df_combined = gather_stats_and_make_table(file_list)
	print("\t[+] Identifying unique fusion and fission events")
	filtered_fusions, filtered_splits, unique_combos_fusions, unique_merians_splits = get_fusions_splits_and_unique_events(df_combined)
	print("\t[+] Assigning species to each fusion and fission event")
	List_unique_fusions = assign_spp_to_fusions(unique_combos_fusions)
	List_unique_fusions.sort(key=sort_list_high2low)		
	List_unique_splits = assign_spp_to_splits(unique_merians_splits)
	print("\t[+] Parsing the tree")
	t = parse_tree(tree_file, label_status)
	print("\t[+] Mapping fusions and splits onto the tree")
	mapped_fusions_dict, lost_fusions_df, t, mapped_subsets, copy_list_unique_fusions = map_fusions(List_unique_fusions, t, threshold)
	mapped_fusions_dict, t = account_for_subsets_of_fusions(copy_list_unique_fusions, t, mapped_subsets)
	mapped_splits_dict, lost_splits_df, t = map_splits(List_unique_splits, t, threshold)
	if not lost_fusions_df.empty: # if there are any lost fusions in the dataframe
		print("\t[+] At least one loss of a fusion found. Mapping lost fusions")
		List_unique_lost_fusions = make_list_unique_losses(lost_fusions_df) 
		mapped_lost_fusions_dict, t = map_lost_fusions(List_unique_lost_fusions, t)
	print("\t[+] Annotating tree with fusion and fission events")
	t = annotate_tree_with_fusions(mapped_fusions_dict, t) # Optional - if want to draw annotated trees
	if not lost_fusions_df.empty: # only annotate with lost fusions if there are any
		t = annotate_tree_with_lost_fusions(mapped_lost_fusions_dict, t)
	t = annotate_tree_with_splits(mapped_splits_dict, t)
	print("\t[+] Writing results to " + str(output_location) + "overall_assignments_" + prefix + ".tsv")
	write_results(output_location, t, df_combined, mapped_fusions_dict, mapped_splits_dict, lost_fusions_df, lost_splits_df, prefix)
	get_event_stats(mapped_fusions_dict, mapped_splits_dict) # Check number of fusions & splits that map to each node
#%%
quit()

# print(t.get_ascii(show_internal=True, attributes = ["name", "Fission1", "Fission2", "Fission3", "Fission4", "Fission5", "Fission6", "Fission7", "Fission8"]))
#	So far the max number of (sensible) fusions per node is 8
# print(t.get_ascii(show_internal=True, attributes = ["name", 'Fusion1', "Fusion2", "Fusion3", "Fusion4", "Fusion5", "Fusion6", "Fusion7", "Fusion8", "Fusion9", "Fusion10", "Fusion11", "Fusion12", "Fusion13", "Fusion14", "Fusion15", "Fusion16", "Fusion17", "Fusion18", "Fusion19", "Fusion20", "Fusion21", "Fusion22", "Fusion23", "Fusion24", "Fusion25", "Fusion26", "Fusion27", "Fusion28", "Fusion29", "Fusion30", "Fusion31", "Fusion32", "Fusion33", "Fusion34", "Fusion35", "Fusion36", "Fusion37", "Fusion38", "Fusion39", "Fusion40", "Fusion41"]))
#print("The final mapped splits are:", mapped_splits_dict) 	# So far the max number of (sensible) fusions per node is 8
#print(t.get_ascii(show_internal=True, attributes = ["name", "Split1", "Split2", "Split3", "Split4", "Split5", "Split6", "Split7", "Split8", "Split9", "Split10", "Split11", "Split12", "Split13", "Split14", "Split15"]))

# This code currently causes jupyter to crash? Not sure if needed anyway?
ts = TreeStyle() # Render tree with mapped fusions
ts.layout_fn = split_layout # Use my custom layout
ts.show_leaf_name = False # Do not add leaf names automatically
#ts.show_branch_length = True
name = str(prefix) + "_Supermatrix_tree_annotated_splits_v2.6.png"
t.render(name, w=183, units="mm", tree_style=ts)
#t.show(tree_style=ts)

#%%
# input_data = '../data/subset_data/'
# tree_file = 'r2_m17_noCrazies.newick.txt'
# output_location = 'test'
# prefix = 'test'
# threshold = 1
# label_status = 'True'
# file_list = parse_info(input_data)
# df_combined = gather_stats_and_make_table(file_list)
# filtered_fusions, filtered_splits, unique_combos_fusions, unique_merians_splits = get_fusions_splits_and_unique_events(df_combined)
# List_unique_fusions = assign_spp_to_fusions(unique_combos_fusions)
# List_unique_fusions.sort(key=sort_list_high2low)
# List_unique_splits = assign_spp_to_splits(unique_merians_splits)
# t = parse_tree(tree_file, label_status)
# mapped_fusions_dict, lost_fusions_df, t, mapped_subsets, List_unique_fusions = map_fusions(List_unique_fusions, t, threshold)
# mapped_fusions_dict, t = account_for_subsets_of_fusions(List_unique_fusions, t, mapped_subsets)

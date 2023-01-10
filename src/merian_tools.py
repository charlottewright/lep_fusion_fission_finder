from cgi import test
from socket import TIPC_DEST_DROPPABLE
from termios import TIOCPKT_FLUSHWRITE
import os, sys, glob, argparse, re
import pandas as pd
from ete3 import Tree
from pathlib import Path
import itertools

class FusionEvent:
  
  def __init__(self, merian_elements, species):
    self.merian_elements = merian_elements
    self.species = species
    self.node_mapping = {} # this will contain the mapped nodes
    self.node_lost = {} # this will contain nodes which have "lost" a fusion
    self.subsets_mapped = {} # this tracks subsets of the fusion that have beeen mapped
    self.remaining_to_map = {}
    self.full_fusion = {} # ADDED on 9/1/23 
    for sp in species:
      self.remaining_to_map[sp] = merian_elements # initalise each species as having all merian elements still to map

  def merian_contains(self, query):
    test = query in self.merian_elements
    return(test)

  def __repr__(self): # just how you would like it to be represented
    if isinstance(self.merian_elements, (list)): # if list (which is the case for ambiguous fusions)
      return str(['|'.join(i) for i in self.merian_elements])
    else:
      return '|'.join(sorted(list(self.merian_elements)))  # if set (which is the case for mappable fusions)

  def update_mapped_nodes(self, node, tip_list): # dictionary of node: tips with node. 
    node_mapping_dict = self.node_mapping # update each time get a mapped node with species at node
    node_mapping_dict[node] = tip_list
    self.node_mapping = node_mapping_dict # reassign value of updated_dict

  def update_lost_nodes(self, node, tip_list):
    node_lost_dict = self.node_lost
    node_lost_dict[node] = tip_list
    self.node_lost = node_lost_dict # reassign value of updated_dict

  def update_mapped_subsets(self, partial_fusion, spp):
    node_mapping_subset_dict = self.subsets_mapped
    if frozenset(partial_fusion) in node_mapping_subset_dict.keys():
      current_spp_list = node_mapping_subset_dict[frozenset(partial_fusion)]
      current_spp_list.append(spp)
      node_mapping_subset_dict[frozenset(partial_fusion)] = current_spp_list
    else:
      node_mapping_subset_dict[frozenset(partial_fusion)] = [spp]
    self.subsets_mapped = node_mapping_subset_dict # reassign value of subset_dict

  def update_remaining_to_map(self, species):
    mapped_merians = set() # use this to track which merians have been mapped
    for partial_fusion, species_list in self.subsets_mapped.items():
      if species in species_list:
          mapped_merians = mapped_merians.union(partial_fusion) # track total merians mapped for a spp
    remaining_merians = set(self.merian_elements.difference(mapped_merians))
    current_remaining_to_map = self.remaining_to_map
    current_remaining_to_map[species] = remaining_merians
    self.remaining_to_map = current_remaining_to_map # reassign 'remaining_to_map' as updated set of merians yet to map

  def update_full_fusion(self, node, full_fusion): # use this to keep track of the full fusion to which ambiguous events belong to
    full_merian_mapping_dict = self.full_fusion
    full_merian_mapping_dict[node] = full_fusion
    self.full_fusion = full_merian_mapping_dict # reassign value of updated_dict

class SplitEvent:

  def __init__(self, merian_elements, species):
    self.merian_elements = merian_elements
    self.species = species
    self.node_mapping = {} # this will contain the mapped nodes
    self.node_lost = {} # this will contain nodes which have "lost" a fusion

  def update_mapped_nodes(self, node, tip_list): # dictionary of node: tips with node. 
    node_mapping_dict = self.node_mapping # update each time get a mapped node with species at node
    node_mapping_dict[node] = tip_list
    self.node_mapping = node_mapping_dict # reassign value of updated_dict

  def update_lost_nodes(self, node, tip_list):
    node_lost_dict = self.node_lost
    node_lost_dict[node] = tip_list
    self.node_lost = node_lost_dict # reassign value of updated_dict

### Functions ###

# Parse input chr assignment files 
def parse_info(input_data):  # get list of fusions ready
  file_list = []
  for filename in Path(input_data).glob('*.tsv'):
    file_list.append(str(filename))
  return file_list

# Parse tree
def parse_tree(tree_file, label_status):
  t = Tree(tree_file, format=1)
  tip_list = t.get_tree_root().get_leaf_names()
  if label_status == 'False':
    count = -1 # Label each internal node with a number starting from 0
    for node in t.traverse("preorder"):
      if node.is_leaf():
        continue
      else:
        count = count + 1
        node.name = str('n') + str(count)
  return(t, tip_list)

# Gather stats 
def gather_stats_and_make_table(file_list, input_data):
  spp_list, d = [], []
  for file_name in file_list:
    if os.path.isfile(file_name) and os.path.getsize(file_name) > 0:
      file = pd.read_csv(file_name, sep='\t', index_col = False) # Header = "query_chr", "status", "assigned_ref_chr"
      Spp_name = file_name.split('/')[int(input_data.count('/'))].replace('_chromosome_assignments.tsv','')
      spp_list.append(Spp_name)
      for query_chr in file.query_chr:
        filtered_data = file[file['query_chr'] == query_chr]
        STATUS = str(filtered_data.iloc[0].status[::])
        ASSIGNED = str(filtered_data.iloc[0].assigned_ref_chr[::])
        entry = {'Chrom_ID':query_chr, 'Status':STATUS, 'Assigned_ref_chr':ASSIGNED, 'Spp':Spp_name } 
        d.append(entry)
  df_combined = pd.DataFrame(d)
  return(df_combined, spp_list)

# Helper function - add a new fusion object to list of fusion objects
def add_fusion_object_to_list(merian_elements, species, list_to_add_object_to): # add a new fusion object to an existsing list of fusion objects
  list_to_add_object_to.append(
  FusionEvent(
    merian_elements, 
    species
  )
  )
  return(list_to_add_object_to)

# Assign species to fusions
def assign_spp_to_fusions(df_combined):
  filtered_fusions = df_combined[df_combined['Status'] == "fusion"] # filter dataframe to only include fused chromosomes
  unique_combos_fusions = filtered_fusions.Assigned_ref_chr.astype('str').unique() # get all unique combinations of Merians
  List_unique_fusions = [] 
  list_of_fusion_objects = []
  for i in unique_combos_fusions: # for each unique combination of Merians
    spp_with_fusion = [] 
    for index, row in filtered_fusions.iterrows(): # for every row in the dataframe 
      Merians = row['Assigned_ref_chr'].split(',') # get list of Merians
      unique_combo = i.split(",")
      if sorted(Merians) == sorted(unique_combo):
        Species = row['Spp']
        spp_with_fusion.append(Species)
    List_unique_fusions.append({'Unique_Merian_combo':i, 'spp_with_fusion':spp_with_fusion})
    list_of_fusion_objects = add_fusion_object_to_list(set(unique_combo), set(spp_with_fusion), list_of_fusion_objects)
  return(List_unique_fusions, list_of_fusion_objects)

#  Sort list high to low
def sort_list_high2low(d): # # Sort list_unique_fusions based on number of merian units per fusion (high to low)
    return len(d['Unique_Merian_combo'])

# Assign species to splits
def assign_spp_to_splits(df_combined):
  List_unique_splits = []
  list_of_split_objects = []
  filtered_splits = df_combined[df_combined['Status'] == "split"]
  unique_merians_splits = filtered_splits.Assigned_ref_chr.astype('str').unique()
  for i in unique_merians_splits: # for each unique Merian 
    subset = filtered_splits[filtered_splits['Assigned_ref_chr'] == str(i)]
    spp_with_split = subset["Spp"].to_list() # get the spp with the split
    entry = {'Unique_Merian_split':i, 'spp_with_split':spp_with_split}
    List_unique_splits.append(entry)
    list_of_split_objects.append(
        SplitEvent(
          i, 
          set(spp_with_split)
        )
    )
  return(List_unique_splits, list_of_split_objects)

# Helper function - find matches at a node relative to a list of target species
def get_matched_species_at_node(node, target_species):
  tip_set = set(node.get_leaf_names()) # get set of tips at node
  matches = sorted(set(target_species) & tip_set) # find species at node which have the fusion exactly
  return(matches, tip_set)

# Map splits
def map_splits(list_of_split_objects, t, threshold):
  for split in list_of_split_objects:
    spp_set = set(split.species) # set of spp with the fusion
    for node in t.traverse("preorder"): 
      matches, tip_set = get_matched_species_at_node(node, spp_set)
      if len(matches) >= 1: # if at least one matched spp check if any other spp has this fusion as part of a bigger fusion (involving more merians)
        spp_with_loss = [d for d in tip_set if d not in spp_set]
        proportion_match = len(matches) / len(tip_set)
        if proportion_match >= threshold:
          split.update_mapped_nodes(node.name, list(matches)) # add entry to dict of matched node and species with fusion
          if len(spp_with_loss) != 0:
            split.update_lost_nodes(node.name, spp_with_loss) # add entry to dict of matched node and species with fusion
          [spp_set.remove(spp) for spp in matches if spp in spp_set] # remove spp that have now matched
  return(list_of_split_objects)

# Make_possible_subsets (Helper function)
def get_possible_subsets(set_of_merian_elements):
  possible_subsets = []
  for i in range(2,(len(set_of_merian_elements))):
    all_tuples = list(itertools.combinations(set_of_merian_elements, i))
    for j in all_tuples:
      possible_subsets.append(j)
  return(possible_subsets)

# Map fusions
def map_fusions(list_of_fusion_objects, t, threshold):
  for fusion in list_of_fusion_objects:
    spp_set = set(fusion.species) # set of spp with the fusion
    for node in t.traverse("preorder"): 
      larger_fusions_at_node, new_spp = {}, [] # initalise as empty list and dict at each node
      temp_spp_list = list(fusion.species) # need to make a copy of list like so to prevent lists being linked
      matches, tip_set = get_matched_species_at_node(node, spp_set) # find species at node which have the fusion exactly and get set of tips at node
      if len(matches) >= 1: # if at least one matched spp check if any other spp has this fusion as part of a bigger fusion (involving more merians)
        for other_fusion in list_of_fusion_objects: 
          matched_merians = (fusion.merian_elements & other_fusion.merian_elements)
          if (len(other_fusion.merian_elements) > len(fusion.merian_elements)) & (len(matched_merians) == len(fusion.merian_elements)): # if all Merians in query fusion event are in this fusion, and this fusion has more merians than the query (i.e. the query is a subset of it). Need the greater than or else spp would already be presen in spp list.
            new_spp = list(other_fusion.species.difference(spp_set)) # convert to list
           # print('The fusion is:', fusion, 'The other fusion is:', other_fusion, 'the new spp are:', new_spp)
            for spp in new_spp: # add new spp to spp_list if found at this node
              if spp in tip_set:
                temp_spp_list.append(spp)
                larger_fusions_at_node[spp] = other_fusion # add to dict the species and larger fusion
            matches = sorted(set(temp_spp_list)&tip_set) # Now update matches again..
        spp_with_loss = [d for d in tip_set if d not in temp_spp_list]
        proportion_match = len(matches) / len(tip_set)
        if proportion_match >= threshold:
          fusion.update_mapped_nodes(node.name, list(matches)) # add entry to dict of matched node and species with fusion
          if len(spp_with_loss) != 0:
            fusion.update_lost_nodes(node.name, spp_with_loss) # add entry to dict of matched node and species with fusion
          [spp_set.remove(spp) for spp in matches if spp in spp_set] # remove spp that have now matched
          if larger_fusions_at_node: # if at least one new spp
            for spp, other_fusion in larger_fusions_at_node.items(): 
             # print('The fusion is:', fusion, 'The other fusion is:', other_fusion, 'the new spp is:', spp)
              other_fusion.update_mapped_subsets(fusion.merian_elements, spp) # add to the fusion object the fact that part of it has been mapped (i.e. a subset has been mapped). From this we can work out (diff) how much "merian" in the fusion is still to be mapped
          #    if (len(other_fusion.merian_elements) - len(fusion.merian_elements)) >= 2: # if two or more merians still to map
              other_fusion.update_remaining_to_map(spp) # update 'remaining_to_map' merians for this species, regardless of number of merians still to map
  return(list_of_fusion_objects)

# Fing remaining fusions
def find_remaining_fusions_to_map(list_of_fusion_objects): # First lets just make a list of fusions where at least 1 spp is yet to have a subset mapped
  fusions_yet_dealt_with = [] # is a mix of subsets yet to be mapped as didn't previously exist in dataset and those that truly cannot be resolved
  for fusion in list_of_fusion_objects:
    if len(fusion.merian_elements) > 2: # this code is only relevent for fusions greater than two :)
      for spp in fusion.species: # for each spp, check whether all subsets of the fusion have mapped already
        mapped_subsets = fusion.subsets_mapped
        partial_fusions_for_spp = [partial_fusion for partial_fusion, list_spp in mapped_subsets.items() if spp in list_spp] # len(partial_fusions_for_spp) # holds the number of partial fusions
        if (len(fusion.merian_elements) - 1 - len(partial_fusions_for_spp)) > 1: # if there are still subsets of a fusion to be mapped for a given spp..
          fusions_yet_dealt_with.append(fusion)
  return(fusions_yet_dealt_with)

# Helper function - make list of all possible tuples from list
def combinantorial(lst):
    index, pairs = 1, []
    for element1 in lst:
        for element2 in lst[index:]:
            pairs.append((element1, element2))
        index += 1
    return pairs

# Map subsets of remaining fusions where possible
def map_subsets_of_fusions_where_possible(fusions_yet_dealt_with, list_of_fusion_objects, t, get_possible_subsets, threshold):
  for fusion in fusions_yet_dealt_with:
    possible_subsets = get_possible_subsets(fusion.merian_elements)
    mapped_subsets = fusion.subsets_mapped
    for possible_subset in possible_subsets:
      updatable_spp_list = list(fusion.species)
      for mapped_merians, mapped_subset_spp in mapped_subsets.items():
          mapped_merians = tuple(list(mapped_merians))
          if sorted(mapped_merians) == sorted(possible_subset):
            [updatable_spp_list.remove(spp) for spp in fusion.species if spp in mapped_subset_spp]
      if len(updatable_spp_list) != 0: # i.e. if there is at least one species which does not have this subset then do the algorithm, else no need!
        for node in t.traverse("preorder"): 
          larger_fusions_at_node, new_spp = {}, [] # initalise as empty list at each node
          temp_spp_list = list(updatable_spp_list) # need to make a copy of list like so to prevent lists being linked
          matches, tip_set = get_matched_species_at_node(node, updatable_spp_list) # find species at node which have the fusion exactly and get set of tips at node
          if len(matches) >= 1: # if at least one matched spp check if any other spp has this fusion as part of a bigger fusion (involving more merians)
            for other_fusion in list_of_fusion_objects: 
              matched_merians = (set(possible_subset) & other_fusion.merian_elements)
              if (len(other_fusion.merian_elements) > len(possible_subset)) & (len(matched_merians) == len(possible_subset)): # if all Merians in query fusion event are in this fusion, and this fusion has more merians than the query (i.e. the query is a subset of it). Need the greater than or else spp would already be presen in spp list.
                new_spp = list(other_fusion.species.difference(updatable_spp_list)) # convert to list
            #   print('the fusion is:', fusion.merian_elements, 'the fusion spp are:', fusion.species, 'matched_merians are:', matched_merians, 'other fusion is:', other_fusion.merian_elements, 'tmp spp list:', temp_spp_list, 'new spp:', new_spp, 'other fusion spp:', other_fusion.species)
                for spp in new_spp: # add new spp to spp_list if found at this node
                  if spp in tip_set:
                    temp_spp_list.append(spp)
                    larger_fusions_at_node[spp] = other_fusion # add to dict the species and larger fusion
                matches = sorted(set(temp_spp_list)&tip_set) # Now update matches again..
            spp_with_loss = [d for d in tip_set if d not in temp_spp_list]
            proportion_match = len(matches) / len(tip_set)
            if (proportion_match >= threshold) & (len(matches) >1): # need to match more than just one species for the subset mapping to mean anything
              existing_fusion = [obj for obj in list_of_fusion_objects if obj.merian_elements == set(possible_subset)]
          #    print('matched_merians are:', matched_merians, 'other fusion is:', other_fusion.merian_elements, 'node is:', node.name)
              if len(existing_fusion)!= 0: # if fusion doesn't already exist in list of unique fusions
                existing_fusion[0].update_mapped_nodes(node.name, list(matches)) # add entry to dict of matched node and species with fusion
              else:
                list_of_fusion_objects = add_fusion_object_to_list(set(possible_subset), set(matches), list_of_fusion_objects) # make a new FusionEvent where subset is the merian_elements and spp_list is the species matched
                list_of_fusion_objects[-1].update_mapped_nodes(node.name, list(matches)) # add entry to dict of matched node and species with fusion
              for spp in updatable_spp_list:
                if spp in matches:
                  fusion.update_mapped_subsets(set(possible_subset), spp) # add to the fusion object the fact that part of it has been mapped (i.e. a subset has been mapped). From this we can work out (diff) how much "merian" in the fusion is still to be mapped
              if len(spp_with_loss) != 0:
                fusion.update_lost_nodes(node.name, spp_with_loss) # add entry to dict of matched node and species with fusion
              [updatable_spp_list.remove(spp) for spp in matches if spp in updatable_spp_list] # remove spp that have now matched
              for spp, other_fusion in larger_fusions_at_node.items(): 
                other_fusion.update_mapped_subsets(set(possible_subset), spp) # add to the fusion object the fact that part of it has been mapped (i.e. a subset has been mapped). From this we can work out (diff) how much "merian" in the fusion is still to be mapped
                #  other_fusion.update_remaining_to_map(spp) # update 'remaining_to_map' merians for this species, regardless of number of merians still to map
  return(list_of_fusion_objects)

def map_ambiguous_fusions(ambiguous_fusions, t): # currently doesn't record losses of ambiguous subsets
  threshold = 1 # hard coded for now
  list_of_ambiguous_fusion_objects = []
  for fusion in ambiguous_fusions: 
    remaining_to_map = fusion.remaining_to_map
    for spp_set, merians in remaining_to_map.items(): 
      remaining_to_map_tuples = combinantorial(list(merians)) # get all possible tuples
      for node in t.traverse("preorder"): 
        matches, tip_set = get_matched_species_at_node(node, [spp_set]) # find species at node which have the fusion exactly and get set of tips at node
        proportion_match = len(matches) / len(tip_set)
        if proportion_match >= threshold:
          for spp in matches: 
            fusion.update_mapped_subsets(set(remaining_to_map_tuples), spp) # add to the fusion ob
            fusion.update_remaining_to_map(spp) # update 'remaining_to_map' merians for this species, regardless of number of merians still to map
            existing_amgiguous_fusion = [obj for obj in list_of_ambiguous_fusion_objects if obj.merian_elements == set(remaining_to_map_tuples)]
            if len(existing_amgiguous_fusion)!= 0: 
              existing_amgiguous_fusion[0].update_mapped_nodes(node.name, list(matches)) # add entry to dict of matched node and species with fusion
            else:
              list_of_ambiguous_fusion_objects = add_fusion_object_to_list(set(remaining_to_map_tuples), spp_set, list_of_ambiguous_fusion_objects) # make a new object for ambiguous fusion
              list_of_ambiguous_fusion_objects[-1].update_mapped_nodes(spp, spp)
          list_of_ambiguous_fusion_objects[-1].update_full_fusion(spp, fusion.merian_elements) # add the full_fusion info to the new subset fusion object
  return(ambiguous_fusions, list_of_ambiguous_fusion_objects)

# Helper function - i.merian_elemebnts, i.full_fusion
def calculate_number_events_to_report(ambiguous_merians, full_fusion, node):
  set_ambiguous_merians =  set([item for sublist in ambiguous_merians for item in sublist]) # works
  total_merians = len(full_fusion[node]) # get length of fusion i.e. number merians involved
  num_merians_fully_mapped = (total_merians - len(set_ambiguous_merians)) # total num merians - number of ambiguous merians (if two numbers aren't same, fusions must have been mapped)
  if num_merians_fully_mapped != 0:
    num_events_to_report = total_merians - num_merians_fully_mapped -1
  else:
    num_events_to_report = total_merians - 2
  return(num_events_to_report) # should be num_events_to_report

def write_mapped_events(event, output_file, event_type, event_loss):
  if event_loss == "true":
    mapped_nodes = event.node_lost
  else:
    mapped_nodes = event.node_mapping
  for node in mapped_nodes:
    mapped_node = node # mapped node
    tip_spp = str(mapped_nodes[node]) # get tip spp which were at that mapped node with event
    if (event_type == "fusion") | (event_type == "ambiguous"):
      merian_elements = str(sorted(list(event.merian_elements))) # first convert set to list so can sort, then convert to string
      merian_elements =re.sub("[^A-Z,0-9\(\)]", "", merian_elements,0,re.IGNORECASE)
    else:
      merian_elements = re.sub("[^A-Z,0-9\(\)]", "", event.merian_elements,0,re.IGNORECASE)
    tip_spp =re.sub("[^A-Z,_ ]", "", tip_spp,0,re.IGNORECASE)
    if event_type == "ambiguous":
      num_events_to_report = calculate_number_events_to_report(event.merian_elements, event.full_fusion, node)
      for __ in range(0,num_events_to_report): 
        output_file.write(("%s\t%s\t%s\t%s\n") % (merian_elements, tip_spp, node, "fusion")) # write one line per event (i.e. so if 2 events occured, write 2 lines)
    else:
      output_file.write(("%s\t%s\t%s\t%s\n") % (merian_elements, tip_spp, mapped_node, event_type)) # write to outputfile

def write_results(output_location, prefix, list_of_fusion_objects, list_of_ambiguous_fusion_objects, list_of_split_objects, write_mapped_events): # save results
  with open(output_location + '/' + 'mapped_fusions_fissions_' + prefix + '.tsv', "w") as output_file:
    output_file.write(("%s\t%s\t%s\t%s\n") % ("Merian_elements", "Tips", "Node", "Event")) # write header
    for i in list_of_fusion_objects:
      write_mapped_events(i, output_file, "fusion", "false")  # write fusions to outputfile
      write_mapped_events(i, output_file, "fission", "true")  # write lost fusions i.e. fissions to outputfile
    for i in list_of_split_objects:
      write_mapped_events(i, output_file, "fission", "false") # write to outputfile
      write_mapped_events(i, output_file, "fusion", "true") # write to outputfile - fission lost therefore fusion has occured
    for i in list_of_ambiguous_fusion_objects:
      write_mapped_events(i, output_file, "ambiguous", "false")
  return
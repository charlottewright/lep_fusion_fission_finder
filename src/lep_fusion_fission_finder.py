#!/usr/bin/env python3
import sys
import argparse

def parse_reference_table(reference_table_file):
	with open(reference_table_file, 'r') as reference_table: 
		buscoID2merian = {} 
		for line in reference_table:
			if not line.startswith("#"):
				cols = line.rstrip("\n").split()
				if cols[1] == "Complete":
					buscoID, merian = cols[0], cols[2].split(":")[0]
					buscoID2merian[buscoID] = merian # key = buscoID, value = Merian
	return buscoID2merian

def parse_query_table(query_table_file, buscoID2merian):
	with open(query_table_file, 'r') as query_table:
		chr2pos, pos2buscoID = {}, {} 
		for line in query_table:
			if not line.startswith("#"): # ignoring comments at top
				cols = line.rstrip("\n").split() # get columns
				if cols[1] == "Complete":
					buscoID, chr, start, end = cols[0], cols[2].split(":")[0], int(cols[3]), int(cols[4])
					pos = (start + end)/2 # get midpoint coordinate of busco (position)
					if buscoID in buscoID2merian: # only get BUSCOs that are in the assignment table (i.e. filter out 'unassigned')
						try: 
							chr2pos[chr].append(pos) # try and add the position to the list of BUSCO positions on that chromosome
						except KeyError: 
							chr2pos[chr] = [pos] # or create a new list if chromosome hasn't be seen yet
						pos2buscoID[str(pos) + "_" + chr] = buscoID # key = position_chr, value = buscoID
	return chr2pos, pos2buscoID


def get_max_merians(chr2pos, pos2buscoID, window_size, warnings_list):
	max_merian_dict = {}
	for chr, pos_list in chr2pos.items(): # for every chromosome
		max_merian_list = []
		if len(pos_list) >= window_size: # providing there is at least one window's worth of BUSCOs in the chromosome
			pos_list = sorted(pos_list) # get a sorted list of all the BUSCO positions 
			# loop through windows
			for window_stop in range(window_size, len(pos_list)+1, window_size): # for every window of specified size
				window_list = pos_list[window_stop-window_size:window_stop] # get a list of buscos in that window
				merian_list = []
				for pos in window_list:
					try: # Added a 'try' in to account for buscos in query table that may be absent in reference table
						merian_list.append(buscoID2merian[pos2buscoID[str(pos) + "_" + chr]]) # add which Merian they come from to a list
					except KeyError:
						continue
				max_merian_list.append(max(set(merian_list), key=merian_list.count)) # find most common Merian in merian_list
			# deal with last window
			if len(pos_list[window_stop:len(pos_list)]) >= 3 and len(pos_list[window_stop:len(pos_list)])/window_size > 0.5: # if the last window has at least three BUSCOs and is at least 50% of the window_size
				window_list = pos_list[window_stop:len(pos_list)]
				merian_list = []
				for pos in window_list:
					try: # Added a 'try' in to account for buscos in query table that may be absent in reference table
						merian_list.append(buscoID2merian[pos2buscoID[str(pos) + "_" + chr]])
					except KeyError:
						continue
				max_merian_list.append(max(set(merian_list), key=merian_list.count)) # find most common Merian 
			# store max merian list for each chromosome in dict
			max_merian_dict[chr] = max_merian_list
		else: # if there is NOT one window's worth of BUSCOs in the chromosome, print warning that chromosome will be ignored
			warnings_list.append("Ignoring " + chr + " as it has fewer BUSCOs (" + str(len(pos_list)) + ") than the window size (" + str(window_size) + ")\n")
	return max_merian_dict, warnings_list

def get_assignments(max_merian_dict, prefix, warnings_list):
	with open(prefix + "_chromosome_assignments.tsv", "w") as chromosome_assignment_file:
		chromosome_assignment_file.write(("%s\t%s\t%s\n") % ("query_chr", "status", "assigned_ref_chr"))
		chromosome_assignment_dict = {} # used later to seperate ancestral chromosomes from splits
		potential_complex_fusions_dict = {} # used later to check for complex chromosomes
		observed_merian_count = 0 # will be used to raise a warning if less than expected Merian counts are found
		observed_merian_list = [] # will be used to store each Merian found
		# loop through Merian dict
		for chr, max_merian_list in max_merian_dict.items():
			if len(sorted(set(max_merian_list))) > 1: # if the chromosome has windows with > 1 Merian, it's a fusion 
				chromosome_assignment_file.write(("%s\t%s\t%s\n") % (chr, "fusion", ",".join(sorted(set(max_merian_list))))) # write to output file
				for merian in set(max_merian_list):
					observed_merian_list.append(merian)
					try:
						potential_complex_fusions_dict[merian].append(chr) # use this to check for complex chromosomes later
					except KeyError:
						potential_complex_fusions_dict[merian] = [chr]
				observed_merian_count += len(sorted(set(max_merian_list))) # increment observed Merian count by number of Merians in fusion
			else: # otherwise it's either a split or an ancestral chromosome, so store in dict where key=Merian, value=list of chromosomes
				try:
					chromosome_assignment_dict[max_merian_list[0]].append(chr)
				except KeyError:
					chromosome_assignment_dict[max_merian_list[0]] = [chr]
				try:
					observed_merian_list.append(max_merian_list[0])
				except KeyError:
					observed_merian_list.append(max_merian_list[0])
		# loop through list of Merians that are either split or ancestral 
		for merian, chr_list in chromosome_assignment_dict.items():
			if len(chr_list) == 1: # if the Merian is only associated with one chromosome, it's ancestral
				chromosome_assignment_file.write(("%s\t%s\t%s\n") % (chr_list[0], "ancestral", merian)) # write to output file
				observed_merian_count += 1 
				try:
					potential_complex_fusions_dict[merian].append(chr) # use this to check for complex chromosomes later
				except KeyError:
					potential_complex_fusions_dict[merian] = [chr]
			else: # otherwise its a split
				for chr in chr_list: # for every query chromosome the split Merian is associated with 
					chromosome_assignment_file.write(("%s\t%s\t%s\n") % (chr, "split", merian)) # write to outputfile
				observed_merian_count += 1
		# check if Merian count is as expected
		if len(sorted(set(observed_merian_list))) != 31:
			warnings_list.append("Number unique Merians found " + str(len(sorted(set(observed_merian_list)))) + " Merians; expected 31 Merians\n")
		else:
			print("All 31 Merians found!")
		if observed_merian_count != 31:
			warnings_list.append("Genome is composed of  " + str(observed_merian_count) + " of Merians; expected 31 blocks\n")
		complex_chr = []
		for merian, chr_list in potential_complex_fusions_dict.items():
			if len(chr_list) != 1: # if a Merian is found on >1 chromosome and those chr aren't clear split fragments:
				complex_chr.extend(chr_list)
	return warnings_list, potential_complex_fusions_dict

def write_complex_events(prefix, potential_complex_fusions_dict, max_merian_dict):
		complex_chr_list = [] # to keep track of which chr have already been declared
		for merian, chr_list in potential_complex_fusions_dict.items():
			if len(chr_list) != 1: # if a Merian is found on >1 chromosome and those chr aren't clear split fragments:
				for chr in chr_list:
						complex_chr_list.append(chr)
		if len(complex_chr_list) > 1:
			with open(prefix + "_complex_chromosomes.tsv", "w") as complex_chromosome_assignments_file:
				complex_chromosome_assignments_file.write(("%s\n") % ("complex_chr")) # write to output file
				for chr in set(complex_chr_list):
					complex_chromosome_assignments_file.write(("%s\n") % (chr)) # write to output file
	#all_chr_list = list(max_merian_dict.keys()) # uncommment if want to write but I don't think necessaary.
	#	simple_chr = list(set(all_chr_list).difference(set(complex_chr_list))) 
	#	for chr in simple_chr:
	#		complex_chromosome_assignments_file.write(("%s\t%s\n") % (chr, "simple")) # write to output file	

def write_warnings(warnings_list):
	with open(prefix + "_warnings.txt", "w") as warnings_file:
		for warning in warnings_list:
			warnings_file.write(warning)

if __name__ == "__main__":
	SCRIPT = "fusion_split_finder2.py"
	# argument set up
	parser = argparse.ArgumentParser()
	parser.add_argument("-r", "--reference_table", type=str, help = "full_table.tsv file for reference species", required=True)
	parser.add_argument("-q", "--query_table", type=str, help = "full_table.tsv for query species", required=True)
	parser.add_argument("-f", "--prefix", type=str, help = "Prefix for all output files", default="fsf")
	parser.add_argument("-w", "--window_size", type=int, help = "Number of BUSCOs to be used per window (must be odd)", default=17)
	args = parser.parse_args()
	reference_table_file = args.reference_table
	query_table_file = args.query_table
	prefix = args.prefix
	window_size = args.window_size
	# make warnings list and check window size
	warnings_list = []
	if window_size <= 0 or (window_size % 2) == 0:
		sys.exit("\t ERROR: Window size must be an postive and odd.")
	print("[+] Running fusion_split_finder2.py with a window size of " + str(window_size))
	# run fuctions
	print("\t[+] Parsing full table files")
	buscoID2merian = parse_reference_table(reference_table_file)
	chr2pos, pos2buscoID = parse_query_table(query_table_file, buscoID2merian)
	print("\t[+] Finding most common Merian in each window")
	max_merian_dict, warnings_list = get_max_merians(chr2pos, pos2buscoID, window_size, warnings_list)
	print("\t[+] Writing assignments to " + prefix + "_chromosome_assignments.txt")
	warnings_list, potential_complex_fusions_dict = get_assignments(max_merian_dict, prefix, warnings_list)
	print("\t[+] Finding complex events and writing to " + prefix + "_complex_chromosome_assignments.tsv")
	write_complex_events(prefix, potential_complex_fusions_dict, max_merian_dict)
	if len(warnings_list) > 0:
		print("\t[+] Writing " + str(len(warnings_list)) + " warnings to " + prefix + "_warnings.txt")
		write_warnings(warnings_list)

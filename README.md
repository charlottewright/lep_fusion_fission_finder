# lep_fusion_fission_finder
A tool to assign ancestral linkage units and/or identify fusion/fission events in Lepidopteran chromosomes based on a set of reference BUSCO genes as markers.

### Running the scripts

## 1.) Find fusions/fissions
`fusion_split_finder.py` takes the full_table.tsv output file for two species, along with an optional prefix (specified with -f, default "fsf"). The default window size for lepidoptera is 17 BUSCOs but this can be changed with the `-w` flag e.g.:

```
python3 fusion_split_finder.py -q test_data/Aglais_io_full_table.tsv -r test_data/Melitaea_cinxia_full_table.tsv -f Aglais-w 17`
```

This will write two files:

- `Aglais_chromosome_assignments.tsv`: a summary of the assignments for each scaffold in the query genome. For fused/fission chromosomes, their putative origins are listed. 

- `Aglais_warnings.tsv`: list of warnings - lists any contigs with under the threshold of BUSCOs specified (default: 17). Also records the number linkage units found if not 31 as expected. Also records total number of units if not the expected 31.

### Full usage: 

```
usage: fusion_split_finder2.py [-h] -r REFERENCE_TABLE -q QUERY_TABLE [-f PREFIX] [-w WINDOW_SIZE]

optional arguments:
  -h, --help            show this help message and exit
  -r REFERENCE_TABLE, --reference_table REFERENCE_TABLE
                        full_table.tsv file for reference species
  -q QUERY_TABLE, --query_table QUERY_TABLE
                        full_table.tsv for query species
  -f PREFIX, --prefix PREFIX
                        Prefix for all output files
  -w WINDOW_SIZE, --window_size WINDOW_SIZE
                        Number of BUSCOs to be used per window (must be odd)
 ```
 
 ## 2.) Place fusions/fissions in a phylogenetic context
 
 `Map_fusion_fissions.py` takes the output of `fusion_split_finder.py` and infers where fusion/fission occured in a given tree.
 
 `./Map_fusion_fissions.py -i chr_assignments/ -tree spp.treefile -t 1 -o output/ -f test_run`
 
 This will write five files:
 
- annotated_fissions_tree.nw  
- annotated_fusions_tree.nw  
- lost_fusions.tsv  
- lost_splits.tsv  
- mapped_fusions_fissions.tsv  
- overall_assignments.tsv
 
 
 ### Full uages:
 
 ```
usage: Map_fusion_fissions.py [-h] -i INPUT_DATA [-tree TREE] -o OUTPUT [-f PREFIX] [-t THRESHOLD] [-l LABEL]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_DATA, --input_data INPUT_DATA
                        path to lep_fusion_fission_finder output
  -tree TREE, --tree TREE
                        Phylogenetic tree
  -o OUTPUT, --output OUTPUT
                        output location relative to working directory
  -f PREFIX, --prefix PREFIX
                        Prefix for all output files
  -t THRESHOLD, --threshold THRESHOLD
                        Threshold for rearrangement to be shared between tips
  -l LABEL, --label LABEL
                        Specify if tree already contains internal node labels
 ```


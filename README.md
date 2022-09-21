# lep_fusion_fission_finder
A tool to assign ancestral linkage units and/or identify fusion/fission events in Lepidopteran chromosomes based on a set of reference BUSCO genes as markers.

### Running the script
`fusion_split_finder.py` takes the full_table.tsv output file for two species, along with an optional prefix (specified with -f, default "fsf"). The default window size for lepidoptera is 17 BUSCOs but this can be changed with the `-w` flag e.g.:

`python3 fusion_split_finder.py -q test_data/Aglais_io_full_table.tsv -r test_data/Melitaea_cinxia_full_table.tsv -f Aglais-w 17`

### The output
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

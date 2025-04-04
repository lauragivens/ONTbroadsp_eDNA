# ONTbroadsp_eDNA

Data and scripts needed to process environmental DNA metabarcoding output from Oxford Nanopore Technology sequencing, used in the manuscript by Givens et al, 2025. These scripts were run on a high-performance computing cluster.

This workflow starts with the output from [n] GridION R9.4 flow cell runs. Reads were demultiplexed locally with MinKNOW v.[??] and the indexed barcodes were removed. 

# Dependencies  

You will need the following: 
NanoFilt v.2.8.0 (https://github.com/wdecoster/nanofilt)
CDHIT v.4.8.1 (https://github.com/weizhongli/cdhit)
Minimap2 v.2.15 (https://github.com/lh3/minimap2)
samtools v.1.10 (https://github.com/samtools/samtools)
Racon v.1.4.20 (https://github.com/isovic/racon)
Qiime2 v.2021.8 (https://github.com/qiime2/qiime2)

You will also need a BLAST database formatted for use with qiime2

# Step 1: Filter reads by quality and size
While indexing barcodes were removed during demultiplexing, the primer sequences were retained. This step uses Nanofilt to trim primer sequences from both ends and selects for sequences within +-100 bp of the target size. 

# Step 2: Cluster reads into consensus sequences  

## Cluster by similarity  
Reads are first grouped using CD-HIT based on a 95% similarity threshold. Groups of 10 or more sequences are selected for error correction. 

## Map reads onto each other  
Minimap2 is used to align sequences within a cluster and generate a file of alignments that will be used for error correction. 
We use map-ont as the mapping algorithm, which was developed for use with noisy Oxford Nanopore reads. Then, each of the sequences within a cluster are aligned against each other and the resulting mapping alignment file is output. 

## Error correct consensus sequence  
Racon is used to polish the consensus sequences and return a single 'consensus sequence'.

All intermediate files are then zipped.

## Second mapping
Consensus sequences are then run through cd-hit again with a higher percent similarity threshold. 

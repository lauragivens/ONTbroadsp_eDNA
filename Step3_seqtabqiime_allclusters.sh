#!/bin/bash
  
#SBATCH --job-name=Template_seqtabqiime_min
#SBATCH --account=schultzlab -p schultzlab,nsoe-it,common
#SBATCH --output=/path_to_sequences/outfile/Template_seqtabqiime.out
#SBATCH --error=/path_to_sequences/outfile/Template_seqtabqiime.err
#SBATCH --ntasks=10
#SBATCH --mem=15G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lag66@duke.edu

# the fastest way to rename paths in this script is to Find "/work/lag66/AugAG" and Replace it with what corresponds to your "/work/NETID/RUN1" # 

# ---------------------------------------------------------------------------------------------------------------------
################ step 6 ## otu table generation and qiime tax assignment ################
## FASTQ files demultiplexed and barcodes removed by GridION MinKNOW ##
## transferred to DCC /work/lag66/fastq_pass
## Barcode files were concatenated and filtered by size in preprocessing step
## then used cd-hit to cluster sequences in each barcode by 90% similarity
## then only clusters with >mincluster sequences in them are selected for mapping and polishing
## then consensus sequences were polished with racon and concatenated into one fasta file per barcode

## in this step, you are building a table that reports what sequences occur in each barcode
source /hpc/group/schultzlab/lag66/miniconda3/etc/profile.d/conda.sh
# --------------------------------------- To Do ---------------------------------------
# rename the following 

d='/path_to_sequences'
mincluster=10
secondmincluster=3

primer1=COI
primer2=18S
primer3=12S
# --------------------------------------- file processing ---------------------------------------
cd $d/cdhit

echo "Merge all samples"
# to make the files that are going to be BLAST'ed, we are concatenating:
# 1. files that were clustered a second time 
# 2. files that did not pass the second clustering threshold
# we are NOT concatenating the file that passed the first clustering threshold because all of the sequences that were in that file are encompassed by the two above

# However, we are using the sequences that did not pass the first clustering threshold to make the OTU table. There are no additional sequences from the second round of clustering to add to those files

cat barcode0*/*min"$secondmincluster"_oldhead.fasta barcode0*/*belowmin"$secondmincluster"seqs_oldhead.fasta barcode10/*min"$secondmincluster"_oldhead.fasta barcode10/*belowmin"$secondmincluster"seqs_oldhead.fasta barcode11/*min"$secondmincluster"_oldhead.fasta barcode11/*belowmin"$secondmincluster"seqs_oldhead.fasta barcode12/*min"$secondmincluster"_oldhead.fasta barcode12/*belowmin"$secondmincluster"seqs_oldhead.fasta barcode13/*min"$secondmincluster"_oldhead.fasta barcode13/*belowmin"$secondmincluster"seqs_oldhead.fasta barcode14/*min"$secondmincluster"_oldhead.fasta barcode14/*belowmin"$secondmincluster"seqs_oldhead.fasta barcode15/*min"$secondmincluster"_oldhead.fasta barcode15/*belowmin"$secondmincluster"seqs_oldhead.fasta barcode16/*min"$secondmincluster"_oldhead.fasta barcode16/*belowmin"$secondmincluster"seqs_oldhead.fasta > "$primer2"_mapped_polished.fasta #fasta file with barcodeID_seqID; size=N

cat barcode17/*min"$secondmincluster"_oldhead.fasta barcode17/*belowmin"$secondmincluster"seqs_oldhead.fasta barcode18/*min"$secondmincluster"_oldhead.fasta barcode18/*belowmin"$secondmincluster"seqs_oldhead.fasta barcode19/*min"$secondmincluster"_oldhead.fasta barcode19/*belowmin"$secondmincluster"seqs_oldhead.fasta barcode2*/*min"$secondmincluster"_oldhead.fasta barcode2*/*belowmin"$secondmincluster"seqs_oldhead.fasta barcode30/*min"$secondmincluster"_oldhead.fasta barcode30/*belowmin"$secondmincluster"seqs_oldhead.fasta barcode31/*min"$secondmincluster"_oldhead.fasta barcode31/*belowmin"$secondmincluster"seqs_oldhead.fasta barcode32/*min"$secondmincluster"_oldhead.fasta barcode32/*belowmin"$secondmincluster"seqs_oldhead.fasta > "$primer1"_mapped_polished.fasta #fasta file with barcodeID_seqID; size=N

cat barcode33/*min"$secondmincluster"_oldhead.fasta barcode33/*belowmin"$secondmincluster"seqs_oldhead.fasta barcode34/*min"$secondmincluster"_oldhead.fasta barcode34/*belowmin"$secondmincluster"seqs_oldhead.fasta barcode35/*min"$secondmincluster"_oldhead.fasta barcode35/*belowmin"$secondmincluster"seqs_oldhead.fasta barcode36/*min"$secondmincluster"_oldhead.fasta barcode36/*belowmin"$secondmincluster"seqs_oldhead.fasta barcode37/*min"$secondmincluster"_oldhead.fasta barcode37/*belowmin"$secondmincluster"seqs_oldhead.fasta barcode38/*min"$secondmincluster"_oldhead.fasta barcode38/*belowmin"$secondmincluster"seqs_oldhead.fasta barcode39/*min"$secondmincluster"_oldhead.fasta barcode39/*belowmin"$secondmincluster"seqs_oldhead.fasta barcode4*/*min"$secondmincluster"_oldhead.fasta barcode4*/*belowmin"$secondmincluster"seqs_oldhead.fasta > "$primer3"_mapped_polished.fasta #fasta file with barcodeID_seqID; size=N 

# sequences below threshold  
cat barcode*/*belowmin"$mincluster"seqs_oldhead.fasta > clustersbelowmin"$mincluster"threshold.fasta

echo Number of "$primer1" OTUs after merging: $(grep -c ">" "$primer1"_mapped_polished.fasta)
echo Number of "$primer2" OTUs after merging: $(grep -c ">" "$primer2"_mapped_polished.fasta)
echo Number of "$primer3" OTUs after merging: $(grep -c ">" "$primer3"_mapped_polished.fasta)
echo Number of OTUs below threshold after merging: $(grep -c ">" clustersbelowmin"$mincluster"threshold.fasta) 

inputprimer1=""$primer1"_mapped_polished.fasta"
inputprimer2=""$primer2"_mapped_polished.fasta"
inputprimer3=""$primer3"_mapped_polished.fasta"

outputprimer1=""$primer1"_mapped_polished_min"$mincluster""
outputprimer2=""$primer2"_mapped_polished_min"$mincluster""
outputprimer3=""$primer3"_mapped_polished_min"$mincluster""

cd $d 
# --------------------------------------- make otu table ---------------------------------------
module load R/4.0.0
R --no-save --file=makeSeqTab.R

module unload R/4.0.0


# --------------------------------------- prep fasta file for qiime ---------------------------------------

cd $d/cdhit 

awk '/^>/ {$0 = $1 "_" (++n) " " $2 " " $3} 1' "$inputprimer1" | sed -e " /^>/s/barcodelabel=//; /^>/s/otu=//; /^>/s/size.*//; /^>/s/$/ orig_bc=AGCT new_bc=AGCT bc_diffs=0/" > "$outputprimer1".qiimeready.fasta
# final product: >barcodeXX_i barcodeXX_clusterXX orig_bc=AGCT new_bc=AGCT bc_diffs=0
awk '/^>/ {print($0)}; /^[^>]/ {print(toupper($0))}' "$outputprimer1".qiimeready.fasta > "$outputprimer1".qiimeready.subs.fasta
#converts any lowercase nucleotides to uppercase without changing the header


awk '/^>/ {$0 = $1 "_" (++n) " " $2 " " $3} 1' "$inputprimer2" | sed -e " /^>/s/barcodelabel=//; /^>/s/otu=//; /^>/s/size.*//; /^>/s/$/ orig_bc=AGCT new_bc=AGCT bc_diffs=0/" > "$outputprimer2".qiimeready.fasta
# final product: >barcodeXX_i barcodeXX_clusterXX orig_bc=AGCT new_bc=AGCT bc_diffs=0
awk '/^>/ {print($0)}; /^[^>]/ {print(toupper($0))}' "$outputprimer2".qiimeready.fasta > "$outputprimer2".qiimeready.subs.fasta
#converts any lowercase nucleotides to uppercase without changing the header

awk '/^>/ {$0 = $1 "_" (++n) " " $2 " " $3} 1' "$inputprimer3" | sed -e " /^>/s/barcodelabel=//; /^>/s/otu=//; /^>/s/size.*//; /^>/s/$/ orig_bc=AGCT new_bc=AGCT bc_diffs=0/" > "$outputprimer3".qiimeready.fasta
# final product: >barcodeXX_i barcodeXX_clusterXX orig_bc=AGCT new_bc=AGCT bc_diffs=0
awk '/^>/ {print($0)}; /^[^>]/ {print(toupper($0))}' "$outputprimer3".qiimeready.fasta > "$outputprimer3".qiimeready.subs.fasta
#converts any lowercase nucleotides to uppercase without changing the header


conda activate qiime2-2021.8
## import quality controlled and demultiplexed sequences as SampleData[Sequences] 
## SampleData[Sequences] indicates collection of sequences associated with one or more samples
qiime tools import \
       --input-path "$outputprimer1".qiimeready.subs.fasta \
       --output-path "$outputprimer1".qiimeready.qza \
       --type 'FeatureData[Sequence]'


qiime tools import \
       --input-path "$outputprimer2".qiimeready.subs.fasta \
       --output-path "$outputprimer2".qiimeready.qza \
       --type 'FeatureData[Sequence]'

qiime tools import \
       --input-path "$outputprimer3".qiimeready.subs.fasta \
       --output-path "$outputprimer3".qiimeready.qza \
       --type 'FeatureData[Sequence]'

# --------------------------------------------------------------------------------------------- #
########### step 5 classify reads ###########

## qiime feature-classifier classify-sklearn
### this is a naive bayes

## qiime feature-classifier classify-consensus-blast
### assigns taxonomy based on BLAST+ local alignment between query and reference_reads

 #qiime feature-classifier classify-sklearn --i-reads "$output".qiimeready.qza \
  #      --i-classifier /hpc/group/schultzlab/lag66/TrainedClassifiers/qiimeNCBI-12S-refseqs-classifier.qza \
   #             --o-classification "$output".cluster10.taxonomy.qza #\
                        #--p-n-jobs -1 --verbose

  qiime feature-classifier classify-consensus-blast --i-query "$outputprimer1".qiimeready.qza \
         --i-reference-reads /hpc/group/schultzlab/lag66/qiimeNCBI_"$primer1"_derep_seqs.qza \
         --i-reference-taxonomy /hpc/group/schultzlab/lag66/qiimeNCBI_"$primer1"_derep_tax.qza \
         --o-classification "$outputprimer1".blasttaxonomy.qza --p-maxaccepts 5 --p-perc-identity 0.9

 qiime feature-classifier classify-consensus-blast --i-query "$outputprimer2".qiimeready.qza \
         --i-reference-reads /hpc/group/schultzlab/lag66/qiimeNCBI_"$primer2"_derep_seqs.qza \
         --i-reference-taxonomy /hpc/group/schultzlab/lag66/qiimeNCBI_"$primer2"_derep_tax.qza \
         --o-classification "$outputprimer2".blasttaxonomy.qza --p-maxaccepts 5 --p-perc-identity 0.9 
 
 qiime feature-classifier classify-consensus-blast --i-query "$outputprimer3".qiimeready.qza \
         --i-reference-reads /hpc/group/schultzlab/lag66/qiimeNCBI_"$primer3"_derep_seqs.qza \
         --i-reference-taxonomy /hpc/group/schultzlab/lag66/qiimeNCBI_"$primer3"_derep_tax.qza \
         --o-classification "$outputprimer3".blasttaxonomy.qza --p-maxaccepts 5 --p-perc-identity 0.9 
 # --------------------------------------------------------------------------------------------- #
 ########### step 6 export classification ###########
 qiime tools export \
   --input-path "$outputprimer1".blasttaxonomy.qza \
   --output-path "$primer1"blasttaxonomy
 # files will be converted from .qza object into whatever normal file they are

 mv "$primer1"blasttaxonomy/taxonomy.tsv "$primer1".qiimeblast.taxonomy.tsv
rm -r "$primer1"blasttaxonomy


qiime tools export \
   --input-path "$outputprimer2".blasttaxonomy.qza \
   --output-path "$primer2"blasttaxonomy
 # files will be converted from .qza object into whatever normal file they are

 mv "$primer2"blasttaxonomy/taxonomy.tsv "$primer2".qiimeblast.taxonomy.tsv
rm -r "$primer2"blasttaxonomy


qiime tools export \
   --input-path "$outputprimer3".blasttaxonomy.qza \
   --output-path "$primer3"blasttaxonomy
 # files will be converted from .qza object into whatever normal file they are

 mv "$primer3"blasttaxonomy/taxonomy.tsv "$primer3".qiimeblast.taxonomy.tsv
rm -r "$primer3"blasttaxonomy



echo "SLURM_JOBID: " $SLURM_JOBID
echo $(date)


## edited 12/18/2023 to adjust job name and add partitions



#!/bin/bash

# --------------------------------------------------------------------------------------------------------------------- #
################ Preprocessing ## Filter for quality and size using Nanofilt (https://github.com/wdecoster/nanofilt) ################

### This script is built with the expectation that you are working in a unique folder for each run 
### There should be a file_list.txt file of barcode names 

################ the following should have happened before this script runs ################
## FASTQ files demultiplexed and barcodes removed by GridION MinKNOW ##

# --------------------------------------- To Do --------------------------------------- #
# change the following
d="/path_to_sequences" #direct computer to where the fastq_pass folder is stored
minlength=600 #minimum length of desired sequences
maxlength=900 #maximum length of desired sequences
headcrop=15 #trim n nucleotides from start of reads
tailcrop=15 #trim n nucleotides from end of reads
# --------------------------------------- prepare conda environment --------------------------------------- #

source /hpc/group/schultzlab/lag66/miniconda3/etc/profile.d/conda.sh ##tells the computer where conda is saved

# --------------------------------------- unzip files and QC with Nanofilt --------------------------------------- #
# note: If you're trying to merge two fastq_pass folders (for example, if you used two runs to sequence your dataset), the following code will recursively add the folders in fastq_pass_X3 to the folders in fastq_pass
# so fastq_pass_X3/barcode01/fastqfile1.gz will be ADDED to fastq_pass/barcode01 without overwriting the folders already in fastq_pass/barcode01
# cp -r fastq_pass_X3/* fastq_pass/
cd $d 

[ -d "outfile" ] && echo "Directory outfile exists, proceed." || (echo "Directory outfile does not exist, making directory." && mkdir outfile && mkdir outfile/outfilePre && mkdir outfile/outfileCDHIT)
[ -d "Barcodes" ] && echo "Directory Barcodes exists, proceed." || (echo "Directory Barcodes does not exist, making directory." && mkdir Barcodes)
[ -d "Nanofilt" ] && echo "Directory Nanofilt exists, proceed." || (echo "Directory Nanofilt does not exist, making directory." && mkdir Nanofilt)
[ -d "cdhit" ] && echo "Directory cdhit exists, proceed." || (echo "Directory cdhit does not exist, making directory." && mkdir cdhit)

conda activate Nanofilt
barcodeID=$(awk "NR==${SLURM_ARRAY_TASK_ID}" file_list.txt)

echo "MiFish primers estimated 135 bp without indices etc; filtered between 100-400"
echo "Forward primer length: 21 bp"
echo "Reverse primer length: 27 bp"
echo 
echo "Folmer primers estimated 700 bp; should be filtered between 600-900"
echo "Forward primer length: 25 bp"
echo "Reverse primer length: 26 bp"
echo
echo "Leray, 18S 1F/400R, and UPA primers estimated 300-400 bp; should be filtered between 300-600"
echo "Leray Forward primer length: 26 bp"
echo "Leray Reverse primer length: 26 bp"
echo "18S Forward primer length: 22 bp"
echo "18S Reverse primer length: 16 bp"

cd $d/fastq_pass

if [ -d "$barcodeID" ]; then
        gunzip -r "$barcodeID"/*
        cat "$barcodeID"/* > "$d"/Barcodes/"$barcodeID".fastq
        NanoFilt -l $minlength --maxlength $maxlength --headcrop $headcrop --tailcrop $tailcrop "$d"/Barcodes/"$barcodeID".fastq > "$d"/Nanofilt/"$barcodeID".fastq
		echo "Nanofilt trimmed $headcrop bp from start and $tailcrop bp from end of reads; retained sequences between $minlength and $maxlength bp in length and at Q9 and above"
else
        cat "$barcodeID"/*.fastq > $d/Barcodes/"$barcodeID".fastq
        echo "fastq_pass files were already unzipped"
        NanoFilt -l $minlength --maxlength $maxlength --headcrop $headcrop --tailcrop $tailcrop "$d"/Barcodes/"$barcodeID".fastq > "$d"/Nanofilt/"$barcodeID".fastq
		echo "Nanofilt trimmed $headcrop bp from start and $tailcrop bp from end of reads; retained sequences between $minlength and $maxlength bp in length and at Q9 and above"
fi

conda deactivate

echo $(date)

# ----------------------------------------------------------------------------------
## edited 07/01/22 to add SLURM output echo and IF/ELSE test
## edited 08/25/22 to add directory commands, length variables and report desired lengths of each primer set
## edited 09/05/22 to add head/tail primer trimming commands
## edited 04/24/23 to expand MiFish parameters  
## edited 12/18/23 to adjust job name and add more partitions

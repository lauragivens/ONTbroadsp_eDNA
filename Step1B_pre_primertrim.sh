#!/bin/bash

#SBATCH --job-name=Template_Primer_Pre_job
#SBATCH --account=schultzlab -p schultzlab,nsoe-it,common
#SBATCH --array=01-48
#SBATCH --output=/path_to_sequences/Template_Primer_Pre_array_job.%A_%a.out
#SBATCH --error=/path_to_sequences/Template_Primer_Pre_array_job.%A_%a.err
#SBATCH --ntasks=10
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lag66@duke.edu

# the fastest way to rename paths in this script is to Find "/path_to_sequences" and Replace it with what corresponds to your "/work/NETID/RUN1" # 

# --------------------------------------------------------------------------------------------------------------------- #
################ Preprocessing ## Filter for quality and size using Nanofilt (https://github.com/wdecoster/nanofilt) ################

### This script is built with the expectation that you are working in a unique DCC folder for each run (/work/lag66/RUN1)
### make a folder to hold the zipped fastq files within the run folder (/work/lag66/RUN1/fastq_pass)
### then upload the fastq_pass barcode files that correspond to your samples into that folder (will upload as /work/lag66/RUN1/fastq_pass/barcodeXX/X.fastqz)
### then upload the file_list.txt file to run folder

################ the following should have happened before this script runs ################
## FASTQ files demultiplexed and barcodes removed by GridION MinKNOW ##
## transferred to DCC /path_to_sequences/fastq_pass
## file_list.txt transferred to DCC /path_to_sequences

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

echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID
echo $(date)

# ----------------------------------------------------------------------------------
## edited 07/01/22 to add SLURM output echo and IF/ELSE test
## edited 08/25/22 to add directory commands, length variables and report desired lengths of each primer set
## edited 09/05/22 to add head/tail primer trimming commands
## edited 04/24/23 to expand MiFish parameters  
## edited 12/18/23 to adjust job name and add more partitions
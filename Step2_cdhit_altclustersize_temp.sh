#!/bin/bash

#SBATCH --job-name=CDHITTemplate
#SBATCH --account=schultzlab -p schultzlab
#SBATCH --array=01-48%5
#SBATCH --output=/path_to_sequences/outfile/outfileCDHIT/CDHITTemplate_array.%A_%a.out
#SBATCH --error=/path_to_sequences/outfile/outfileCDHIT/CDHITTemplate_array.%A_%a.err
#SBATCH --ntasks=10
#SBATCH --mem=15G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lag66@duke.edu

# the fastest way to rename paths in this script is to Find "path_to_sequences" and Replace it with what corresponds to your "/work/NETID/RUN1" # 
# in array, the %5 restricts the number of jobs running at a time to 5
# --------------------------------------------------------------------------------------------------------------------- #
## FASTQ files demultiplexed and barcodes removed by GridION MinKNOW ##
## transferred to DCC /work/lag66/fastq_pass
## Barcode files were concatenated and filtered by size in preprocessing step

## this script uses cd-hit to cluster sequences in each barcode by 90% similarity
## then only clusters with >mincluster sequences in them are selected for mapping and polishing
## we end with racon

# --------------------------------------- load modules --------------------------------------- #
source /hpc/group/schultzlab/lag66/miniconda3/etc/profile.d/conda.sh
echo "activate cd-hit"
echo
conda activate cd-hit #this conda environment contains the packages seqtk and cd-hit
echo "load modules"
module load Minimap2/2.15 
module load samtools/1.10
module load Racon/1.4.20-rhel8
module load compatbin
echo
# --------------------------------------- To Do --------------------------------------- #
# change the following 

d="/path_to_sequences"
clustering=0.9 #max=1
mincluster=10
maxseq=75000
outfolder="clustermin"$mincluster""


# --------------------------------------- subsample --------------------------------------- #
NanofiltFILE=""$d"/Nanofilt/"$barcodeID".fastq"
echo "test if files exist to continue"
[ -f "$NanofiltFILE" ] && echo " "$NanofiltFILE" exists" || (echo " "$NanofiltFILE" does not exist, exiting" && exit)
echo
echo
NUMfastq="$(echo $(cat "$NanofiltFILE" |wc -l)/4|bc)" 
echo
echo "test to subsample reads"
echo
if (($NUMfastq > $maxseq)) ; then
	echo " "$NanofiltFILE" contains "$NUMfastq" reads"
	echo "randomly subsample Nanofilt reads to max 50k of "$barcodeID".fastq"
	echo "set random seed = 123"
	seqtk sample -s123 $d/Nanofilt/"$barcodeID".fastq 50000 > $d/Nanofilt/"$barcodeID"_sub50k.fastq

	echo "converting subsampled Nanofilt fastq files to fasta"
	seqtk seq -a $d/Nanofilt/"$barcodeID"_sub50k.fastq > $d/Nanofilt/"$barcodeID"_sub50k.fasta
	
	fasta=""$d"/Nanofilt/"$barcodeID"_sub50k.fasta"
else
	echo " "$NanofiltFILE" contains "$NUMfastq" reads"
	echo "will not subsample"
	echo "converting Nanofilt fastq files to fasta"
	seqtk seq -a $d/Nanofilt/"$barcodeID".fastq > $d/Nanofilt/"$barcodeID".fasta

	fasta=""$d"/Nanofilt/"$barcodeID".fasta"
fi
echo
echo "fasta file to cluster is $fasta "
echo
# --------------------------------------- cluster --------------------------------------- #
echo "test if files exist to continue"
[ -f "$fasta" ] && echo " "$fasta" exists" || (echo " "$fasta" does not exist, exiting" && exit)
echo
echo 
echo "generate separate fasta file for each cluster over specifed size"
echo "minimum number of sequences in clusters is "$mincluster" "
cd $d/cdhit/"$barcodeID"
make_multi_seq.pl "$fasta" "$barcodeID"_cluster90.clstr "$outfolder" "$mincluster"
# basic command is make_multi_seq.pl input.fasta input.clstr outputfolder clustersize
echo
conda deactivate

echo "test if files exist to continue"
echo
#z is true if the length of the string is 0; n is true if the length is nonzero
echo 
[ -z "$(ls -A "$d/cdhit/"$barcodeID"/"$outfolder"")" ] && (echo "no clusters with more than "$mincluster" sequences exist for "$barcodeID" " && exit) || echo " there are $(ls $d/cdhit/"$barcodeID"/"$outfolder" | wc -l) clusters for "$barcodeID" "
echo
echo "rename cluster fasta files"
echo
cd $d/cdhit/"$barcodeID"/"$outfolder"

for file in * ; do
        if [[ -f $file ]]; then
                ext=`ls $file | rev | cut -d '.' -f 1 | rev`
                file_name=`ls $item | rev | cut -d '.' -f 2 | rev`
                if [ "$ext" == "fasta" ]; then
                        echo "files already renamed"
                else
                        if [[ ("$ext" == "clstr" || "$ext" == "sh" || "$ext" == "txt" || "$ext" == "err" || "$ext" == "out") ]]; then
                                exit
                        else
                        ClusterName=$(echo "$(basename ${file%.fasta})")
                        mv "$file" "$barcodeID"_cluster"$ClusterName".fasta
                        echo " cluster "$ClusterName" file has "$(grep -c ">" "$barcodeID"_cluster"$ClusterName".fasta)" sequences"
                        fi
                fi
        fi
done

echo
echo
clusterfolder=""$d"/cdhit/"$barcodeID"/"$outfolder""
[ -d ""$clusterfolder"/mapped" ] && echo "Directory $clusterfolder/mapped exists, proceed." || echo "Directory $clusterfolder/mapped does not exist, making directory." && mkdir $clusterfolder/mapped
[ -d ""$clusterfolder"/polished" ] && echo "Directory $clusterfolder/polished exists, proceed." || echo "Directory $clusterfolder/polished does not exist, making directory." && mkdir $clusterfolder/polished && mkdir $clusterfolder/polished/renamed
[ -d ""$clusterfolder"/referenceFASTA" ] && echo "Directory $clusterfolder/referenceFASTA exists, proceed." || echo "$clusterfolder/referenceFASTA does not exist, making directory." && mkdir $clusterfolder/referenceFASTA
echo "at this point, you should have a SequencingRunDirectory. Inside the directory is a Barcodes folder and a Nanofilt folder with one fastq file per barcode"
echo "Inside the directory is a cdhit directory with folders for each barcode"
echo "Inside the cdhit/barcode folder, there is a text file of clusters and the reads assigned to them, a fasta file of representative clusters, and the cluster directory"
echo "Inside the cluster directory is one fasta file per cluster"
echo 
echo 
echo "this path is built like so"
echo "/SequencingRunDirectory"
echo "/SequencingRunDirectory/Nanofilt"
echo "/SequencingRunDirectory/cdhit"
echo "/SequencingRunDirectory/cdhit/barcode01"
echo "/SequencingRunDirectory/cdhit/barcode01/metadata.files"
echo "/SequencingRunDirectory/cdhit/barcode01/"$outfolder""
echo "/SequencingRunDirectory/cdhit/barcode01/"$outfolder"/barcode01_cluster1.fasta"
echo "/SequencingRunDirectory/cdhit/barcode01/"$outfolder"/barcode01_etc.fasta"
echo
echo "if you were to look inside the fasta file, you would have only the raw sequences from the Nanofilt file that are assigned to that individual cluster"
echo
echo "the remaining actions are performed inside the cdhit/barcode/cluster file"
echo "thus each array (which defines the barcode the script is run on) will also have loops (which will ensure the script is run on each cluster within the barcode)"
echo
# --------------------------------------- pipeline --------------------------------------- #
echo "change into cluster directory"
cd $clusterfolder
echo "reminder that clusterfolder directory is $clusterfolder"
echo
echo 
echo "choose one of the fasta sequences in each cluster as a reference"
# we were using AVA to generate consensus sequences, so since with cd-hit we've already got the consensus sequences grouped we can go to RVC
echo "map sequences in each cluster to the reference with minimap"
echo "polish mapped alignments with racon"


for file in barcode*; do
clusterID=$(echo "$(basename ${file%.fasta})")
echo "$(cat "$file" | head -2)" > ref_"$file"

echo "making minimap alignments"
        echo ref_"$file"
        minimap2 -x map-ont -a ref_"$file" "$file" > $clusterfolder/mapped/"$clusterID"_mapped.sam
        [ -s ""$clusterfolder"/mapped/"$clusterID"_mapped.sam" ] && echo " "$file" mapping complete; continuing to racon polish" || (echo " "$file" mapping did not complete, exiting" && exit)
echo "performing racon polish"
        racon $clusterfolder/"$outfolder"/"$clusterID".fasta $clusterfolder/mapped/"$clusterID"_mapped.sam ref_"$file" > $clusterfolder/polished/"$clusterID"_raconpolish.fasta
        [ -s ""$clusterfolder"/polished/"$clusterID"_raconpolish.fasta" ] && echo " "$file" polishing complete" || (echo " "$file" polishing did not complete, exiting" && exit)
        # s argument says if the file exists and is not empty then continue
        mv ref_"$file" referenceFASTA
echo "renaming headers"
echo "adding quantities into file header"
        quant="$(echo "$(grep -c ">" "$file")")"
        clusterName="$(echo "$(echo "$clusterID" | cut -d '_' -f 2)" )"
        newclusterID="$(echo "$clusterID" | sed 's/_.*/;otu=/')"
        echo "cluster is made of "$quant" reads"
        echo "clusterName is "$clusterName""
                echo "clusterID is "$clusterID""
        echo "newclusterID is "$newclusterID""
        awk -v clusterName="$clusterID" -v quant="$quant" '{print ">"clusterName";size="quant; getline; print}' $clusterfolder/polished/"$clusterID"_raconpolish.fasta > $clusterfolder/polished/renamed/"$clusterID"_quant.fasta
echo "addition of quantities finished" # results in >barcodexx_clusterxx;size=xx
echo "adding barcodelabel to identify which barcode 01-48 the sequence belongs to and a label to identify the cluster"
        sed '/^>/s/>*/>barcodelabel='"$barcodeID"';otu=/' $clusterfolder/polished/renamed/"$clusterID"_quant.fasta > $clusterfolder/polished/renamed/"$clusterID"_renamed.fasta
                echo "addition of barcode label and unique cluster identifiers finished"
                echo


done
echo "Concatenate all clusters into one fasta file"
cat $clusterfolder/polished/renamed/*_renamed.fasta > $clusterfolder/"$barcodeID"_min"$mincluster".fasta
echo "Number of fasta seqs in $clusterfolder/"$barcodeID"_min"$mincluster".fasta: $(grep -c ">" $clusterfolder/"$barcodeID"_min"$mincluster".fasta)"
# note if your minimum cluster size returns too many sequences
# awk -F "=" '/^>/ {if ($4>100) {print; getline; print}}' barcode*/*all.fasta >> newclustersize.fasta
# change whatever $4 is greater than to your chosen cutoff

echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID
echo $(date)

######## change names in 'assembly' ava minimap.fasta file so that racon will iterate later #########
#FILERENAME=""$d"/1_AVA/"$barcodeID"_sub50pct_minimap.fasta" # path to input file for renaming
#awk '/^>/{print ">_seq" ++i; next}{print}' $FILERENAME | sed '/^>/s/>*/>'"$barcodeID"'/' > "$barcodeID"_sub50pct_minimap_renamed.fasta # awk will print '>_seq' and a sequential number after every time it finds '>'; then print the next line; this is passed to sed which will substitute > with >barcodeXX

## understanding if operators, conditional statements, and various bracketing ##
# https://acloudguru.com/blog/engineering/conditions-in-bash-scripting-if-statements#h-the-basic-rules-of-bash-conditions
# https://developer.ibm.com/tutorials/l-bash-test/
# ---------------------------------------------------------------------------------------------------------------------
## written 08/10/2022 
## edited 08/26/2022 to update some file tests 
## edited 1/23/2024 to add bc module (for counting fastq sequences)
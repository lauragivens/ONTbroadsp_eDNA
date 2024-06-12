#!/bin/bash

#SBATCH --job-name=Template_2clus
#SBATCH --account=schultzlab -p schultzlab,nsoe-it,common
#SBATCH --array=01-48%5
#SBATCH --output=/path_to_sequences/outfile/Template_2clus.%A_%a.out
#SBATCH --error=/path_to_sequences/outfile/Template_2clus.%A_%a.err
#SBATCH --mem=15G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lag66@duke.edu

#set -e 
#set -e will Exit immediately if a pipeline (which may consist of a single simple command), a list, or a compound command (see SHELL GRAMMAR above), exits with a non-zero status. The shell does not exit if the command that fails is part of the command list immediately following a while or until keyword, part of the test following the if or elif reserved words, part of any command executed in a && or || list except the command following the final && or ||, any command in a pipeline but the last, or if the command's return value is being inverted with !. If a compound command other than a subshell returns a non-zero status because a command failed while -e was being ignored, the shell does not exit. A trap on ERR, if set, is executed before the shell exits. This option applies to the shell environment and each subshell environment separately (see COMMAND EXECUTION ENVIRONMENT above), and may cause subshells to exit before executing all the commands in the subshell. 

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

d="/path_to_sequences"
clustering=0.98
clusteringlong=98
mincluster=10 #original cluster size
secondmincluster=3
outfolder="clustersize"$clusteringlong"_secondary"

cd $d
barcodeID=$(awk "NR==${SLURM_ARRAY_TASK_ID}" file_list.txt)
clusterfolder=""$d"/cdhit/"$barcodeID""
cd $clusterfolder
#fasta="$clusterfolder/"$barcodeID"_belowmin"$secondmincluster"seqs_oldhead.fasta"
fasta="$clusterfolder/"$barcodeID"_min"$mincluster"_oldhead.fasta"
cd-hit-est -i "$fasta" -o $d/cdhit/"$barcodeID"/"$barcodeID"_cluster"$clusteringlong"_secondary -c "$clustering" -n 10 -d 0 -M 0 -T 0 -g 1
# c denotes clustering percentage
# n indicates word length (10 is default)
# d indicates length of description in .clstr file, set to 0 takes the fasta defline and stops at first space
	## d needs to be 0 for make_multi_seq script to work
# M indicates memory limit in MB, set to 0 indicates no memory limit
# T indicates number of threads, set to 0 means all CPUs will be used
echo
echo
echo
echo "output of cd-hit-est is file of representative sequences and a text file of lists of clusters"
echo "sequence file needs fasta extension added"
mv $d/cdhit/"$barcodeID"/"$barcodeID"_cluster"$clusteringlong"_secondary $d/cdhit/"$barcodeID"/"$barcodeID"_cluster"$clusteringlong"_secondary.fasta
echo
echo "fasta file extension added"
echo
echo "plot the cluster distribution using clstr file"
plot_len1.pl $d/cdhit/"$barcodeID"/"$barcodeID"_cluster"$clusteringlong"_secondary.clstr \
  1,2-4,5-9,10-19,20-49,50-99,100-299,500-99999 \
  100-349,350-699,700-999,999-999999 > $d/cdhit/"$barcodeID"/"$barcodeID"_secondarycluster_distribution.txt
#second line is size of clusters (standard is 1,2-4,5-9,10-19,20-49,50-99,100-299,500-99999)
#third line is lengths of sequences (standard is 10-59,60-149,150-499,500-1999,2000-999999)
echo
echo
echo "generate separate fasta file for each cluster over specifed size"
echo "minimum number of sequences in clusters is "$secondmincluster" "
cd $d/cdhit/"$barcodeID"

make_multi_seq.pl "$fasta" "$barcodeID"_cluster"$clusteringlong"_secondary.clstr secondary_clustermin0 0 # make a folder with fasta files for ALL clusters

make_multi_seq.pl "$fasta" "$barcodeID"_cluster"$clusteringlong"_secondary.clstr "$outfolder" "$secondmincluster"
# basic command is make_multi_seq.pl input.fasta input.clstr outputfolder clustersize

echo
conda deactivate

echo "test if files exist to continue"
echo
#z is true if the length of the string is 0; n is true if the length is nonzero
echo
[ -z "$(ls -A "$d/cdhit/"$barcodeID"/"$outfolder"")" ] && echo "no clusters with more than "$secondmincluster" sequences exist for "$barcodeID" " && exit || echo " there are $(ls $d/cdhit/"$barcodeID"/"$outfolder" | wc -l) clusters for "$barcodeID" "
#When you run a command in () you're spawning a subshell. So when you call exit within that subshell, you're just exiting it, and not your top level script 
#https://unix.stackexchange.com/questions/296526/set-e-in-a-subshell
#[ -z "$(ls -A "$d/cdhit/"$barcodeID"/"$outfolder"")" ] && {echo "no clusters with more than "$secondmincluster" sequences exist for "$barcodeID" " && exit} 
#using curly braces instead of parentheses did not work
# note: need to have a space after the opening curly bracket and before echo because apparently bash gets confused if there isn't a space 
echo
echo "rename cluster fasta files"
echo
cd $d/cdhit/"$barcodeID"/"$outfolder"

for f in *; do case "$f" in *.*) echo skipped $f;; *) mv "$f" "$barcodeID"_cluster"$f".fasta; esac; done # looks at all of the files within the outfolder directory, any that have a '.' are skipped (and you are told they're skipped), the rest have the fasta extension added

cd $clusterfolder/"$outfolder"
[ -d ""$clusterfolder"/"$outfolder"/mapped" ] && echo "Directory $clusterfolder/"$outfolder"/mapped exists, proceed." || echo "Directory $clusterfolder/"$outfolder"/mapped does not exist, making directory." && mkdir $clusterfolder/"$outfolder"/mapped
[ -d ""$clusterfolder"/"$outfolder"/polished" ] && echo "Directory $clusterfolder/"$outfolder"/polished exists, proceed." || echo "Directory $clusterfolder/"$outfolder"/polished does not exist, making directory." && mkdir $clusterfolder/"$outfolder"/polished && mkdir $clusterfolder/"$outfolder"/polished/renamed
[ -d ""$clusterfolder"/"$outfolder"/referenceFASTA" ] && echo "Directory $clusterfolder/"$outfolder"/referenceFASTA exists, proceed." || echo "$clusterfolder/"$outfolder"/referenceFASTA does not exist, making directory." && mkdir $clusterfolder/"$outfolder"/referenceFASTA

echo "choose one of the fasta sequences in each cluster as a reference"
# we were using AVA to generate consensus sequences, so since with cd-hit we've already got the consensus sequences grouped we can go to RVC
echo "map sequences in each cluster to the reference with minimap"
echo "polish mapped alignments with racon"

#  ----------------------------------------------------------------------------
# Polish and concatenate the sequences above the minimum threshold
#  ----------------------------------------------------------------------------
echo $PWD
for file in barcode*; do #for all files that start with barcode*
clusterID=$(echo "$(basename ${file%.fasta})")
echo "$(cat "$file" | head -2)" > ref_"$file"

echo "making minimap alignments"
        echo ref_"$file"
        minimap2 -x map-ont -a ref_"$file" "$file" > $clusterfolder/"$outfolder"/mapped/"$clusterID"_mapped.sam
        [ -s ""$clusterfolder"/"$outfolder"/mapped/"$clusterID"_mapped.sam" ] && echo " "$file" mapping complete; continuing to racon polish" || { echo >&2 " "$file" mapping did not complete, exiting"; exit 1; }
echo "performing racon polish"
        racon $clusterfolder/"$outfolder"/"$clusterID".fasta $clusterfolder/"$outfolder"/mapped/"$clusterID"_mapped.sam ref_"$file" > $clusterfolder/"$outfolder"/polished/"$clusterID"_raconpolish.fasta
        [ -s ""$clusterfolder"/"$outfolder"/polished/"$clusterID"_raconpolish.fasta" ] && echo " "$file" polishing complete" || { echo >&2 " "$file" polishing did not complete, exiting"; exit 1; }
        # s argument says if the file exists and is not empty then continue
        mv ref_"$file" referenceFASTA
echo "renaming headers"
echo "adding quantities into file header"
	quant="$(echo "$(grep ';size=' "$file" | cut -d "=" -f4 | awk '{sum += $1} END {print sum}')")"
	#quant="$(echo "$(grep -c ">" "$file")")"
        clusterName="$(echo "$(echo "$clusterID" | cut -d '_' -f 2)" )"
        newclusterID="$(echo "$clusterID" | sed 's/_.*/;otu=/')"
        echo "cluster is made of "$quant" reads"
        echo "clusterName is "$clusterName""
                echo "clusterID is "$clusterID""
        echo "newclusterID is "$newclusterID""

		#awk -v barcodeID="$barcodeID" -v quant="$quant" '{print ">"barcodeID";size="quant; getline; print}' "$clusterfolder"/"$outfolder"/polished/"$clusterID"_raconpolish.fasta > "$clusterfolder"/"$outfolder"/polished/renamed/"$clusterID"_quant_newhead.fasta
		awk -v clusterName="$clusterID" -v quant="$quant" '{print ">"clusterName";size="quant; getline; print}' $clusterfolder/"$outfolder"/polished/"$clusterID"_raconpolish.fasta > $clusterfolder/"$outfolder"/polished/renamed/"$clusterID"_quant_oldhead.fasta

echo "addition of quantities finished" # results in >barcodexx_clusterxx;size=xx
echo "adding barcodelabel to identify which barcode 01-48 the sequence belongs to and a label to identify the cluster"

		#sed '/^>/s/>*/>otu='"$clusterName"';barcodelabel=/' "$clusterfolder"/"$outfolder"/polished/renamed/"$clusterID"_quant_newhead.fasta > "$clusterfolder"/"$outfolder"/polished/renamed/"$clusterID"_renamed_newhead.fasta
        sed '/^>/s/>*/>barcodelabel='"$barcodeID"';otu=clustered/' $clusterfolder/"$outfolder"/polished/renamed/"$clusterID"_quant_oldhead.fasta > $clusterfolder/"$outfolder"/polished/renamed/"$clusterID"_renamed_oldhead.fasta

                echo "addition of barcode label and unique cluster identifiers finished"
                echo

done
echo "Concatenate all clusters into one fasta file"
cat $clusterfolder/"$outfolder"/polished/renamed/*_renamed_oldhead.fasta > $clusterfolder/"$barcodeID"_min"$secondmincluster"_oldhead.fasta
#cat $clusterfolder/"$outfolder"/polished/renamed/*_renamed_newhead.fasta > $clusterfolder/"$barcodeID"_min"$secondmincluster"_newhead.fasta
echo "Number of fasta seqs in $clusterfolder/"$barcodeID"_min"$secondmincluster"_oldhead.fasta: $(grep -c ">" $clusterfolder/"$barcodeID"_min"$secondmincluster"_oldhead.fasta)"
#echo "Number of fasta seqs in $clusterfolder/"$barcodeID"_min"$secondmincluster"_newhead.fasta: $(grep -c ">" $clusterfolder/"$barcodeID"_min"$secondmincluster"_newhead.fasta)"


#  ----------------------------------------------------------------------------
# now transition to making a file for the sequences below the minimum threshold
#  ----------------------------------------------------------------------------

#  ----------------------------------------------------------------------------
# for loop to make sure we're in the right directory to write renamed files 
#  ----------------------------------------------------------------------------
cd $clusterfolder 

# is there a clustersbelow directory or zipped clustersbelow directory? 
if [ -d secondary_clustersbelow"$secondmincluster" ] ; then 
	echo "secondary_clustersbelow"$secondmincluster" directory already unzipped" 
	# first check whether there is a renamed folder within directory 
	if [ -f $clusterfolder/secondary_clustersbelow"$secondmincluster"/renamed.zip ] ; then
		# if there is a zipped renamed directory within clustersbelow folder, remove it 
                 rm $clusterfolder/secondary_clustersbelow"$secondmincluster"/renamed.zip
		echo "removed $clusterfolder/secondary_clustersbelow"$secondmincluster"/renamed.zip" 
         elif [ -d $clusterfolder/secondary_clustersbelow"$secondmincluster"/renamed ] ; then
		 # if there is a renamed directory within the clustersbelowfolder, remove it 
                 rm -r $clusterfolder/secondary_clustersbelow"$secondmincluster"/renamed
		 echo "removed $clusterfolder/secondary_clustersbelow"$secondmincluster"/renamed"
         else    
		 # if there is no renamed directory within the clsutersbelow folder, say so and move on 
                 echo " no renamed folder or zipped file within secondary_clustersbelow"$secondmincluster" "
         fi
		 echo "working directory is "$PWD" " 
         mkdir $clusterfolder/secondary_clustersbelow"$secondmincluster"/renamed
         
	 # next check whether the fasta files are within another clustersbelow folder 
	 if [[ $(find $clusterfolder/secondary_clustersbelow"$secondmincluster" -maxdepth 1 -type f -name "*.fasta" | wc -l) -gt 1 ]] ; then
                 # -gt 1 sees if the output of find has more than 1 fasta file 
				 # if there are loose fasta files within the first clustersbelow folder
		 echo "fasta files detected within $clusterfolder/secondary_clustersbelow"$secondmincluster" directory"
         
		 # if there are loose fasta  files we then check to see if there's a clustersbelow subdirectory within clustersbelow and remove if there is        
		 if [ -d $clusterfolder/secondary_clustersbelow"$secondmincluster"/secondary_clustersbelow"$secondmincluster" ] ; then
			 rm -r $clusterfolder/secondary_clustersbelow"$secondmincluster"/secondary_clustersbelow"$secondmincluster"
			 echo "removed $clusterfolder/secondary_clustersbelow"$secondmincluster"/secondary_clustersbelow"$secondmincluster" directory because fasta files are outside the directory" 
			 
		 elif [ -f $clusterfolder/secondary_clustersbelow"$secondmincluster"/secondary_clustersbelow"$secondmincluster".zip ] ; then 
			 # if there is a zipped clustersbelow directory within the clustersbelow directory, remove it 
                        rm $clusterfolder/secondary_clustersbelow"$secondmincluster"/secondary_clustersbelow"$secondmincluster".zip
                        echo "removed $clusterfolder/secondary_clustersbelow"$secondmincluster"/secondary_clustersbelow"$secondmincluster" zipped directory because fasta files are outside the directory"
		else 
			# if there is no clustersbelow file within the clustersbelow directory, say so 
		       echo "no secondary_clustersbelow"$secondmincluster" files within $clusterfolder/secondary_clustersbelow"$secondmincluster" directory"	
		fi 
		# now change directory and move on 
                 cd $clusterfolder/secondary_clustersbelow"$secondmincluster" || exit

	# if there are not loose fasta files in the main clustersbelow directory, check if there is a folder they are in
	# there are three options
 	# fasta files can be inside a zipped clustersbelow folder
	# fasta files can be inside an unzipped clustersbelow folder
	# the code didn't work and the files are actually inside the main clustersbelow directory 	
	 elif [ -f $clusterfolder/secondary_clustersbelow"$secondmincluster"/secondary_clustersbelow"$secondmincluster".zip ] ; then
                 # if the fasta files are in a zipped clustersbelow folder 
                 unzip $clusterfolder/secondary_clustersbelow"$secondmincluster"/secondary_clustersbelow"$secondmincluster".zip
		 echo "unzipped $clusterfolder/secondary_clustersbelow"$secondmincluster"/secondary_clustersbelow"$secondmincluster".zip"
                cd $clusterfolder/secondary_clustersbelow"$secondmincluster"/secondary_clustersbelow"$secondmincluster"  || exit
	elif [ -d $clusterfolder/secondary_clustersbelow"$secondmincluster"/secondary_clustersbelow"$secondmincluster" ] ; then
		# if the fasta files are in an unzipped clustersbelow folder
		cd $clusterfolder/secondary_clustersbelow"$secondmincluster"/secondary_clustersbelow"$secondmincluster"
	else    
                 cd $clusterfolder/secondary_clustersbelow"$secondmincluster" || exit
                 echo "no zipped secondary_clustersbelow"$secondmincluster" file within secondary_clustersbelow"$secondmincluster" " 
	 fi
	 echo $PWD

	 # if the clustersbelow folder is instead zipped
elif [ -f secondary_clustersbelow"$secondmincluster".zip ] ; then 
	unzip *.zip	
	echo "unzipped secondary_clustersbelow"$secondmincluster".zip " 
	# first check if there is a renamed folder within the clustersbelow directory
	# there are three options
	# there is a zipped renamed folder
	# there is an unzipped renamed folder
	# there is no renamed folder 
	if [ -f $clusterfolder/secondary_clustersbelow"$secondmincluster"/renamed.zip ] ; then 
		rm $clusterfolder/secondary_clustersbelow"$secondmincluster"/renamed.zip
	        echo "removed $clusterfolder/secondary_clustersbelow"$secondmincluster"/renamed.zip"
	elif [ -d $clusterfolder/secondary_clustersbelow"$secondmincluster"/renamed ] ; then 
		rm -r $clusterfolder/secondary_clustersbelow"$secondmincluster"/renamed
		echo "removed $clusterfolder/secondary_clustersbelow"$secondmincluster"/renamed"
		#mkdir $clusterfolder/clustersbelow"$secondmincluster"/renamed
	else 
		echo " no renamed folder or zipped file within secondary_clustersbelow"$secondmincluster" "
	fi
	echo $PWD
	mkdir $clusterfolder/secondary_clustersbelow"$secondmincluster"/renamed 
	
	# next check whether the fasta files are within another clustersbelow folder
	# if there are loose fasta files remove any zipped or unzipped nested clustersbelow directories
	if [[ $(find $clusterfolder/secondary_clustersbelow"$secondmincluster" -maxdepth 1 -type f -name "*.fasta" | wc -l) -gt 1 ]] ; then

		echo "fasta files detected within $clusterfolder/secondary_clustersbelow"$secondmincluster" directory"
		if [ -f $clusterfolder/secondary_clustersbelow"$secondmincluster"/secondary_clustersbelow"$secondmincluster".zip ] ; then
		       # if there is a zipped clustersbelow directory within the clustersbelow directory, remove it	
			rm $clusterfolder/secondary_clustersbelow"$secondmincluster"/secondary_clustersbelow"$secondmincluster".zip
			echo "removed $clusterfolder/secondary_clustersbelow"$secondmincluster"/secondary_clustersbelow"$secondmincluster" zipped directory because fasta files are outside the directory"
		elif [ -d $clusterfolder/secondary_clustersbelow"$secondmincluster"/secondary_clustersbelow"$secondmincluster" ] ; then 
			# if there is an unzipped clustersbelow directory, remove it  
			rm -r $clusterfolder/secondary_clustersbelow"$secondmincluster"/secondary_clustersbelow"$secondmincluster"
			echo "removed $clusterfolder/secondary_clustersbelow"$secondmincluster"/secondary_clustersbelow"$secondmincluster" directory because fasta files are outside the directory"
		else     
                       echo "no secondary_clustersbelow"$secondmincluster" files within $clusterfolder/secondary_clustersbelow"$secondmincluster" directory"
		fi 
		cd $clusterfolder/secondary_clustersbelow"$secondmincluster" || exit

        # if there are not loose fasta files in the main clustersbelow directory, check if there is a folder they are in
        # there are three options
        # fasta files can be inside a zipped clustersbelow folder
        # fasta files can be inside an unzipped clustersbelow folder
        # the code didn't work and the files are actually inside the main clustersbelow directory  
	elif [ -f $clusterfolder/secondary_clustersbelow"$secondmincluster"/secondary_clustersbelow"$secondmincluster".zip ] ; then
                 # if the fasta files are in a zipped clustersbelow folder 
                 unzip $clusterfolder/secondary_clustersbelow"$secondmincluster"/secondary_clustersbelow"$secondmincluster".zip
                 echo "unzipped $clusterfolder/secondary_clustersbelow"$secondmincluster"/secondary_clustersbelow"$secondmincluster".zip"
                cd $clusterfolder/secondary_clustersbelow"$secondmincluster"/secondary_clustersbelow"$secondmincluster"  || exit
        elif [ -d $clusterfolder/secondary_clustersbelow"$secondmincluster"/secondary_clustersbelow"$secondmincluster" ] ; then
                # if the fasta files are in an unzipped clustersbelow folder
                cd $clusterfolder/secondary_clustersbelow"$secondmincluster"/secondary_clustersbelow"$secondmincluster"
        else
                 cd $clusterfolder/secondary_clustersbelow"$secondmincluster" || exit
                 echo "no zipped secondary_clustersbelow"$secondmincluster" file within secondary_clustersbelow"$secondmincluster" " 
	fi
	echo $PWD
	
else
	mkdir $clusterfolder/secondary_clustersbelow"$secondmincluster"
	cd $clusterfolder/secondary_clustersbelow"$secondmincluster" || exit
		if [ -d renamed ] ; then 
			rm -r renamed 
			mkdir renamed
		elif [ -f renamed.zip ] ; then
		       rm renamed.zip
		       mkdir renamed 	       
		fi
		
fi
echo $PWD
#  ----------------------------------------------------------------------------
# end of for loop to make sure we're in the right directory to write renamed files 
#  ----------------------------------------------------------------------------


cd $clusterfolder/secondary_clustermin0
for f in *; do case "$f" in *.*) echo skipped $f;; *) mv "$f" "$barcodeID"_cluster"$f".fasta; esac; done

cd $clusterfolder # go back to the 'barcode' folder inside cdhit  
diff --brief --from-file=secondary_clustermin0 $outfolder | grep secondary_clustermin0: | awk '{print $4}' > secondary_clustersbelow"$secondmincluster".txt #find all of the secondary_clusters that were below the minimum threshold and write their names to a txt folder
rsync -av --files-from secondary_clustersbelow"$secondmincluster".txt --remove-source-files secondary_clustermin0/ secondary_clustersbelow"$secondmincluster"/ #move all the secondary_clusters below the minimum threshold to a new folder

cd $clusterfolder/secondary_clustersbelow"$secondmincluster"
 
echo "Concatenate all secondary clusters below threshold into one fasta file"

cd $clusterfolder/secondary_clustersbelow"$secondmincluster"
find ./ -type f -name "*.fasta" -exec cat {} + > $clusterfolder/"$barcodeID"_belowmin"$secondmincluster"seqs_oldhead.fasta
echo "Number of fasta seqs in $clusterfolder/"$barcodeID"_belowmin"$secondmincluster"seqs_oldhead.fasta: $(grep -c ">" $clusterfolder/"$barcodeID"_belowmin"$secondmincluster"seqs_oldhead.fasta)"
#find ./ -type f -name "*_newhead.fasta" -exec cat {} + > $clusterfolder/"$barcodeID"_belowmin"$secondmincluster"seqs_newhead.fasta
#echo "Number of fasta seqs in $clusterfolder/"$barcodeID"_belowmin"$secondmincluster"seqs_newhead.fasta: $(grep -c ">" $clusterfolder/"$barcodeID"_belowmin"$secondmincluster"seqs_newhead.fasta)"


#  ----------------------------------------------------------------------------
# Clean up the folders 
#  ----------------------------------------------------------------------------

cd $clusterfolder 
#mkdir $clusterfolder/clustersbelow"$secondmincluster"/clustersbelow"$secondmincluster" # too many files for this to work on  
#mv clustersbelow"$secondmincluster"/*fasta clustersbelow"$secondmincluster"/clustersbelow"$secondmincluster" 

rm -r $clusterfolder/secondary_clustermin0

if [ -f $clusterfolder/secondary_clustersbelow"$secondmincluster".zip ] ; then 
rm secondary_clustersbelow"$secondmincluster".zip 
fi 
zip -r secondary_clustersbelow"$secondmincluster".zip secondary_clustersbelow"$secondmincluster"
rm -r secondary_clustersbelow"$secondmincluster"

cd $clusterfolder/"$outfolder"
zip -r mapped.zip mapped
rm -r mapped

zip -r polished.zip polished
rm -r polished

zip -r referenceFASTA.zip referenceFASTA
rm -r referenceFASTA

echo "Number of fasta seqs in $clusterfolder/"$barcodeID"_belowmin"$secondmincluster"seqs_oldhead.fasta: $(grep -c ">" $clusterfolder/"$barcodeID"_belowmin"$secondmincluster"seqs_oldhead.fasta)"
#echo "Number of fasta seqs in $clusterfolder/"$barcodeID"_belowmin"$secondmincluster"seqs_newhead.fasta: $(grep -c ">" $clusterfolder/"$barcodeID"_belowmin"$secondmincluster"seqs_newhead.fasta)"

echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID
echo $(date)

## edited 12/18/2023 to adjust job name and add more partitions 
## edited 1/23/2024 to add module load compatbin (for bc for counting fastq seqs)

#!/bin/bash

#SBATCH --job-name=Template_2clus
#SBATCH --array=01-48%5
echo
echo
echo "activate cd-hit"
echo
#source /conda/location/path #add path to conda source if necessary
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
fasta="$clusterfolder/"$barcodeID"_min"$mincluster".fasta"
cd-hit-est -i "$fasta" -o $d/cdhit/"$barcodeID"/"$barcodeID"_cluster"$clusteringlong"_secondary -c "$clustering" -n 10 -d 0 -M 0 -T 0 -g 1
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
echo
[ -z "$(ls -A "$d/cdhit/"$barcodeID"/"$outfolder"")" ] && echo "no clusters with more than "$secondmincluster" sequences exist for "$barcodeID" " && exit || echo " there are $(ls $d/cdhit/"$barcodeID"/"$outfolder" | wc -l) clusters for "$barcodeID" "
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
        clusterName="$(echo "$(echo "$clusterID" | cut -d '_' -f 2)" )"
        newclusterID="$(echo "$clusterID" | sed 's/_.*/;otu=/')"
        echo "cluster is made of "$quant" reads"
        echo "clusterName is "$clusterName""
                echo "clusterID is "$clusterID""
        echo "newclusterID is "$newclusterID""
		echo
		awk -v clusterName="$clusterID" -v quant="$quant" '{print ">"clusterName";size="quant; getline; print}' $clusterfolder/"$outfolder"/polished/"$clusterID"_raconpolish.fasta > $clusterfolder/"$outfolder"/polished/renamed/"$clusterID"_quant.fasta

echo "addition of quantities finished" # results in >barcodexx_clusterxx;size=xx
echo "adding barcodelabel to identify which barcode 01-48 the sequence belongs to and a label to identify the cluster"
echo
sed '/^>/s/>*/>barcodelabel='"$barcodeID"';otu=clustered/' $clusterfolder/"$outfolder"/polished/renamed/"$clusterID"_quant.fasta > $clusterfolder/"$outfolder"/polished/renamed/"$clusterID"_renamed.fasta
echo "addition of barcode label and unique cluster identifiers finished"
echo

done
echo "Concatenate all clusters into one fasta file"
cat $clusterfolder/"$outfolder"/polished/renamed/*_renamed.fasta > $clusterfolder/"$barcodeID"_min"$secondmincluster".fasta
echo "Number of fasta seqs in $clusterfolder/"$barcodeID"_min"$secondmincluster".fasta: $(grep -c ">" $clusterfolder/"$barcodeID"_min"$secondmincluster".fasta)"
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
echo
cd $clusterfolder/secondary_clustermin0
for f in *; do case "$f" in *.*) echo skipped $f;; *) mv "$f" "$barcodeID"_cluster"$f".fasta; esac; done

cd $clusterfolder # go back to the 'barcode' folder inside cdhit  
diff --brief --from-file=secondary_clustermin0 $outfolder | grep secondary_clustermin0: | awk '{print $4}' > secondary_clustersbelow"$secondmincluster".txt #find all of the secondary_clusters that were below the minimum threshold and write their names to a txt folder
rsync -av --files-from secondary_clustersbelow"$secondmincluster".txt --remove-source-files secondary_clustermin0/ secondary_clustersbelow"$secondmincluster"/ #move all the secondary_clusters below the minimum threshold to a new folder

cd $clusterfolder/secondary_clustersbelow"$secondmincluster"
 
echo "Concatenate all secondary clusters below threshold into one fasta file"

cd $clusterfolder/secondary_clustersbelow"$secondmincluster"
find ./ -type f -name "*.fasta" -exec cat {} + > $clusterfolder/"$barcodeID"_belowmin"$secondmincluster"seqs.fasta
echo "Number of fasta seqs in $clusterfolder/"$barcodeID"_belowmin"$secondmincluster"seqs.fasta: $(grep -c ">" $clusterfolder/"$barcodeID"_belowmin"$secondmincluster"seqs.fasta)"

#  ----------------------------------------------------------------------------
# Clean up the folders 
#  ----------------------------------------------------------------------------

cd $clusterfolder 
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

echo "Number of fasta seqs in $clusterfolder/"$barcodeID"_belowmin"$secondmincluster"seqs.fasta: $(grep -c ">" $clusterfolder/"$barcodeID"_belowmin"$secondmincluster"seqs.fasta)"

echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID
echo $(date)

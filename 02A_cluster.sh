#!/bin/bash

#SBATCH --job-name=Template_CDHIT
#SBATCH --array=01-48%5

# --------------------------------------- load modules --------------------------------------- #
#source /conda/location/path #add path to conda source if necessary
echo "activate cd-hit"
conda activate cd-hit #this conda environment contains the packages seqtk and cd-hit
echo "load modules"
module load Minimap2/2.15 
module load samtools/1.10
module load Racon/1.4.20-rhel8
module load compatbin
echo

d="/path_to_sequences"
clustering=0.95 #max=1
clusteringlong=95
mincluster=10
maxseq=300000
outfolder="clustermin"$mincluster""

cd $d
barcodeID=$(awk "NR==${SLURM_ARRAY_TASK_ID}" file_list.txt)
# --------------------------------------- move out and err files ----------------------------------- #
if [ "$(find -maxdepth 2 -name "Pre*_[1-48].err")" != "" ] ; then
        mv Pre*.* outfile/outfilePre
        echo "moved preprocessing output files" 
else
        echo "preprocessing output files already moved"
fi

[ -d ""$d"/cdhit" ] && echo "Directory $d/cdhit exists, proceed." || (echo "Directory $d/cdhit does not exist, making directory." && mkdir $d/cdhit)
mkdir $d/cdhit/"$barcodeID"

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
	echo "randomly subsample Nanofilt reads of "$barcodeID".fastq"
	echo "set random seed = 123"
	seqtk sample -s123 $d/Nanofilt/"$barcodeID".fastq $maxseq > $d/Nanofilt/"$barcodeID"_sub.fastq

	echo "converting subsampled Nanofilt fastq files to fasta"
	seqtk seq -a $d/Nanofilt/"$barcodeID"_sub.fastq > $d/Nanofilt/"$barcodeID"_sub.fasta
	
	fasta=""$d"/Nanofilt/"$barcodeID"_sub.fasta"
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
echo "fastq file was converted to fasta"
echo
echo "begin clustering with cd-hit-est"
echo "cd-hit-est is specifically used for nucleotide data"
echo "generation of "$clustering" similarity clusters (max = 1)"
echo
cd-hit-est -i "$fasta" -o $d/cdhit/"$barcodeID"/"$barcodeID"_cluster"$clusteringlong" -c "$clustering" -n 10 -d 0 -M 0 -T 0 -g 1
echo
echo "output of cd-hit-est is file of representative sequences and a text file of lists of clusters"
echo "sequence file needs fasta extension added"
mv $d/cdhit/"$barcodeID"/"$barcodeID"_cluster"$clusteringlong" $d/cdhit/"$barcodeID"/"$barcodeID"_cluster"$clusteringlong".fasta
echo
echo "fasta file extension added"
echo
echo "plot the cluster distribution using clstr file"
plot_len1.pl $d/cdhit/"$barcodeID"/"$barcodeID"_cluster"$clusteringlong".clstr \
  1,2-4,5-9,10-19,20-49,50-99,100-299,500-99999 \
  100-349,350-699,700-999,999-999999 > $d/cdhit/"$barcodeID"/"$barcodeID"_cluster_distribution.txt
#second line is size of clusters (standard is 1,2-4,5-9,10-19,20-49,50-99,100-299,500-99999)
#third line is lengths of sequences (standard is 10-59,60-149,150-499,500-1999,2000-999999)
echo
echo
echo "generate separate fasta file for each cluster over specifed size"
echo "minimum number of sequences in clusters is "$mincluster" "
cd $d/cdhit/"$barcodeID"
make_multi_seq.pl "$fasta" "$barcodeID"_cluster"$clusteringlong".clstr clustermin0 0 # make a folder with fasta files for ALL clusters
make_multi_seq.pl "$fasta" "$barcodeID"_cluster"$clusteringlong".clstr "$outfolder" "$mincluster"
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

for f in *; do case "$f" in *.*) echo skipped $f;; *) mv "$f" "$barcodeID"_cluster"$f".fasta; esac; done # looks at all of the files within the outfolder directory, any that have a '.' are skipped (and you are told they're skipped), the rest have the fasta extension added
echo

clusterfolder=""$d"/cdhit/"$barcodeID""
[ -d ""$clusterfolder"/"$outfolder"/mapped" ] && echo "Directory $clusterfolder/"$outfolder"/mapped exists, proceed." || echo "Directory $clusterfolder/"$outfolder"/mapped does not exist, making directory." && mkdir $clusterfolder/"$outfolder"/mapped
[ -d ""$clusterfolder"/"$outfolder"/polished" ] && echo "Directory $clusterfolder/"$outfolder"/polished exists, proceed." || echo "Directory $clusterfolder/"$outfolder"/polished does not exist, making directory." && mkdir $clusterfolder/"$outfolder"/polished && mkdir $clusterfolder/"$outfolder"/polished/renamed
[ -d ""$clusterfolder"/"$outfolder"/referenceFASTA" ] && echo "Directory $clusterfolder/"$outfolder"/referenceFASTA exists, proceed." || echo "$clusterfolder/"$outfolder"/referenceFASTA does not exist, making directory." && mkdir $clusterfolder/"$outfolder"/referenceFASTA
echo "change into cluster directory"
cd $clusterfolder/"$outfolder"
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
        [ -s ""$clusterfolder"/"$outfolder"/mapped/"$clusterID"_mapped.sam" ] && echo " "$file" mapping complete; continuing to racon polish" || (echo " "$file" mapping did not complete, exiting" && exit)
echo "performing racon polish"
        racon $clusterfolder/"$outfolder"/"$clusterID".fasta $clusterfolder/"$outfolder"/mapped/"$clusterID"_mapped.sam ref_"$file" > $clusterfolder/"$outfolder"/polished/"$clusterID"_raconpolish.fasta
        [ -s ""$clusterfolder"/"$outfolder"/polished/"$clusterID"_raconpolish.fasta" ] && echo " "$file" polishing complete" || (echo " "$file" polishing did not complete, exiting" && exit)
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
        echo
		awk -v clusterName="$clusterID" -v quant="$quant" '{print ">"clusterName";size="quant; getline; print}' $clusterfolder/"$outfolder"/polished/"$clusterID"_raconpolish.fasta > $clusterfolder/"$outfolder"/polished/renamed/"$clusterID"_quant.fasta

echo "addition of quantities finished" # results in >barcodexx_clusterxx;size=xx
echo "adding barcodelabel to identify which barcode 01-48 the sequence belongs to and a label to identify the cluster"
echo
sed '/^>/s/>*/>barcodelabel='"$barcodeID"';otu=/' $clusterfolder/"$outfolder"/polished/renamed/"$clusterID"_quant.fasta > $clusterfolder/"$outfolder"/polished/renamed/"$clusterID"_renamed.fasta
echo "addition of barcode label and unique cluster identifiers finished"
echo

done
echo "Concatenate all clusters into one fasta file"
cat $clusterfolder/"$outfolder"/polished/renamed/*_renamed.fasta > $clusterfolder/"$barcodeID"_min"$mincluster".fasta
echo "Number of fasta seqs in $clusterfolder/"$barcodeID"_min"$mincluster".fasta: $(grep -c ">" $clusterfolder/"$barcodeID"_min"$mincluster".fasta)"
#  ----------------------------------------------------------------------------
# now transition to making a file for the sequences below the minimum threshold
#  ----------------------------------------------------------------------------
#
#  ----------------------------------------------------------------------------
# for loop to make sure we're in the right directory to write renamed files 
#  ----------------------------------------------------------------------------
cd $clusterfolder 
echo
# is there a clustersbelow directory or zipped clustersbelow directory? 
if [ -d clustersbelow"$mincluster" ] ; then 
	echo "clustersbelow"$mincluster" directory already unzipped" 
	# first check whether there is a renamed folder within directory 
	if [ -f $clusterfolder/clustersbelow"$mincluster"/renamed.zip ] ; then
		# if there is a zipped renamed directory within clustersbelow folder, remove it 
                 rm $clusterfolder/clustersbelow"$mincluster"/renamed.zip
		echo "removed $clusterfolder/clustersbelow"$mincluster"/renamed.zip" 
         elif [ -d $clusterfolder/clustersbelow"$mincluster"/renamed ] ; then
		 # if there is a renamed directory within the clustersbelowfolder, remove it 
                 rm -r $clusterfolder/clustersbelow"$mincluster"/renamed
		 echo "removed $clusterfolder/clustersbelow"$mincluster"/renamed"
         else    
		 # if there is no renamed directory within the clsutersbelow folder, say so and move on 
                 echo " no renamed folder or zipped file within clustersbelow"$mincluster" "
         fi
		 echo "working directory is "$PWD" " 
         mkdir $clusterfolder/clustersbelow"$mincluster"/renamed
         
	 # next check whether the fasta files are within another clustersbelow folder 
	 if [[ $(find $clusterfolder/clustersbelow"$mincluster" -maxdepth 1 -type f -name "*.fasta" | wc -l) -gt 1 ]] ; then
                 # -gt 1 sees if the output of find has more than 1 fasta file 
				 # if there are loose fasta files within the first clustersbelow folder
		 echo "fasta files detected within $clusterfolder/clustersbelow"$mincluster" directory"
         
		 # if there are loose fasta  files we then check to see if there's a clustersbelow subdirectory within clustersbelow and remove if there is        
		 if [ -d $clusterfolder/clustersbelow"$mincluster"/clustersbelow"$mincluster" ] ; then
			 rm -r $clusterfolder/clustersbelow"$mincluster"/clustersbelow"$mincluster"
			 echo "removed $clusterfolder/clustersbelow"$mincluster"/clustersbelow"$mincluster" directory because fasta files are outside the directory" 
			 
		 elif [ -f $clusterfolder/clustersbelow"$mincluster"/clustersbelow"$mincluster".zip ] ; then 
			 # if there is a zipped clustersbelow directory within the clustersbelow directory, remove it 
                        rm $clusterfolder/clustersbelow"$mincluster"/clustersbelow"$mincluster".zip
                        echo "removed $clusterfolder/clustersbelow"$mincluster"/clustersbelow"$mincluster" zipped directory because fasta files are outside the directory"
		else 
			# if there is no clustersbelow file within the clustersbelow directory, say so 
		       echo "no clustersbelow"$mincluster" files within $clusterfolder/clustersbelow"$mincluster" directory"	
		fi 
		# now change directory and move on 
                 cd $clusterfolder/clustersbelow"$mincluster" || exit

	# if there are not loose fasta files in the main clustersbelow directory, check if there is a folder they are in
	# there are three options
 	# fasta files can be inside a zipped clustersbelow folder
	# fasta files can be inside an unzipped clustersbelow folder
	# the code didn't work and the files are actually inside the main clustersbelow directory 	
	 elif [ -f $clusterfolder/clustersbelow"$mincluster"/clustersbelow"$mincluster".zip ] ; then
                 # if the fasta files are in a zipped clustersbelow folder 
                 unzip $clusterfolder/clustersbelow"$mincluster"/clustersbelow"$mincluster".zip
		 echo "unzipped $clusterfolder/clustersbelow"$mincluster"/clustersbelow"$mincluster".zip"
                cd $clusterfolder/clustersbelow"$mincluster"/clustersbelow"$mincluster"  || exit
	elif [ -d $clusterfolder/clustersbelow"$mincluster"/clustersbelow"$mincluster" ] ; then
		# if the fasta files are in an unzipped clustersbelow folder
		cd $clusterfolder/clustersbelow"$mincluster"/clustersbelow"$mincluster"
	else    
                 cd $clusterfolder/clustersbelow"$mincluster" || exit
                 echo "no zipped clustersbelow"$mincluster" file within clustersbelow"$mincluster" " 
	 fi
	 echo $PWD

	 # if the clustersbelow folder is instead zipped
elif [ -f clustersbelow"$mincluster".zip ] ; then 
	unzip clustersbelow"$mincluster".zip	
	echo "unzipped clustersbelow"$mincluster".zip " 
	# first check if there is a renamed folder within the clustersbelow directory
	# there are three options
	# there is a zipped renamed folder
	# there is an unzipped renamed folder
	# there is no renamed folder 
	if [ -f $clusterfolder/clustersbelow"$mincluster"/renamed.zip ] ; then 
		rm $clusterfolder/clustersbelow"$mincluster"/renamed.zip
	        echo "removed $clusterfolder/clustersbelow"$mincluster"/renamed.zip"
	elif [ -d $clusterfolder/clustersbelow"$mincluster"/renamed ] ; then 
		rm -r $clusterfolder/clustersbelow"$mincluster"/renamed
		echo "removed $clusterfolder/clustersbelow"$mincluster"/renamed"
		#mkdir $clusterfolder/clustersbelow"$mincluster"/renamed
	else 
		echo " no renamed folder or zipped file within clustersbelow"$mincluster" "
	fi
	echo $PWD
	mkdir $clusterfolder/clustersbelow"$mincluster"/renamed 
	
	# next check whether the fasta files are within another clustersbelow folder
	# if there are loose fasta files remove any zipped or unzipped nested clustersbelow directories
	if [[ $(find $clusterfolder/clustersbelow"$mincluster" -maxdepth 1 -type f -name "*.fasta" | wc -l) -gt 1 ]] ; then

		echo "fasta files detected within $clusterfolder/clustersbelow"$mincluster" directory"
		if [ -f $clusterfolder/clustersbelow"$mincluster"/clustersbelow"$mincluster".zip ] ; then
		       # if there is a zipped clustersbelow directory within the clustersbelow directory, remove it	
			rm $clusterfolder/clustersbelow"$mincluster"/clustersbelow"$mincluster".zip
			echo "removed $clusterfolder/clustersbelow"$mincluster"/clustersbelow"$mincluster" zipped directory because fasta files are outside the directory"
		elif [ -d $clusterfolder/clustersbelow"$mincluster"/clustersbelow"$mincluster" ] ; then 
			# if there is an unzipped clustersbelow directory, remove it  
			rm -r $clusterfolder/clustersbelow"$mincluster"/clustersbelow"$mincluster"
			echo "removed $clusterfolder/clustersbelow"$mincluster"/clustersbelow"$mincluster" directory because fasta files are outside the directory"
		else     
                       echo "no clustersbelow"$mincluster" files within $clusterfolder/clustersbelow"$mincluster" directory"
		fi 
		cd $clusterfolder/clustersbelow"$mincluster" || exit

        # if there are not loose fasta files in the main clustersbelow directory, check if there is a folder they are in
        # there are three options
        # fasta files can be inside a zipped clustersbelow folder
        # fasta files can be inside an unzipped clustersbelow folder
        # the code didn't work and the files are actually inside the main clustersbelow directory  
	elif [ -f $clusterfolder/clustersbelow"$mincluster"/clustersbelow"$mincluster".zip ] ; then
                 # if the fasta files are in a zipped clustersbelow folder 
                 unzip $clusterfolder/clustersbelow"$mincluster"/clustersbelow"$mincluster".zip
                 echo "unzipped $clusterfolder/clustersbelow"$mincluster"/clustersbelow"$mincluster".zip"
                cd $clusterfolder/clustersbelow"$mincluster"/clustersbelow"$mincluster"  || exit
        elif [ -d $clusterfolder/clustersbelow"$mincluster"/clustersbelow"$mincluster" ] ; then
                # if the fasta files are in an unzipped clustersbelow folder
                cd $clusterfolder/clustersbelow"$mincluster"/clustersbelow"$mincluster"
        else
                 cd $clusterfolder/clustersbelow"$mincluster" || exit
                 echo "no zipped clustersbelow"$mincluster" file within clustersbelow"$mincluster" " 
	fi
	echo $PWD
	
else
	mkdir $clusterfolder/clustersbelow"$mincluster"
	cd $clusterfolder/clustersbelow"$mincluster" || exit
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
cd $clusterfolder/clustermin0
for f in *; do case "$f" in *.*) echo skipped $f;; *) mv "$f" "$barcodeID"_cluster"$f".fasta; esac; done

cd $clusterfolder # go back to the 'barcode' folder inside cdhit  
diff --brief --from-file=clustermin0 $outfolder | grep clustermin0: | awk '{print $4}' > clustersbelow"$mincluster".txt #find all of the clusters that were below the minimum threshold and write their names to a txt folder
rsync -av --files-from clustersbelow"$mincluster".txt --remove-source-files clustermin0/ clustersbelow"$mincluster"/ #move all the clusters below the minimum threshold to a new folder

cd $clusterfolder/clustersbelow"$mincluster"
mkdir renamed
for i in barcode*; do
	clusterID=$(echo "$(basename ${i%.fasta})")
	        quant="$(echo "$(grep -c ">" "$i")")"
	        clusterName="$(echo "$(echo "$clusterID" | cut -d '_' -f 2)" )"
	        newclusterID="$(echo "$clusterID" | sed 's/_.*/;otu=/')"
	        echo "renaming headers; "$clusterID" is made of "$quant" reads"
	        
	       awk -v clusterName="$clusterID" -v quant="$quant" '{print ">"clusterName";size="quant; getline; print}NR==2{exit}' "$clusterID".fasta > $clusterfolder/clustersbelow"$mincluster"/renamed/"$clusterID"_quant.fasta
	# results in >barcodexx_clusterxx;size=xx
	# adding barcodelabel to identify which barcode 01-48 the sequence belongs to and a label to identify the cluster
    		sed '/^>/s/>*/>barcodelabel='"$barcodeID"';otu=/' $clusterfolder/clustersbelow"$mincluster"/renamed/"$clusterID"_quant.fasta > $clusterfolder/clustersbelow"$mincluster"/renamed/"$clusterID"_renamed.fasta
	        
done 
echo "Concatenate all clusters below threshold into one fasta file"

cd $clusterfolder/clustersbelow"$mincluster"/renamed
find ./ -type f -name "*renamed.fasta" -exec cat {} + > $clusterfolder/"$barcodeID"_belowmin"$mincluster"seqs.fasta
echo "Number of fasta seqs in $clusterfolder/"$barcodeID"_belowmin"$mincluster"seqs.fasta: $(grep -c ">" $clusterfolder/"$barcodeID"_belowmin"$mincluster"seqs.fasta)"

#  ----------------------------------------------------------------------------
# Clean up the folders 
#  ----------------------------------------------------------------------------

cd $clusterfolder 
rm -r $clusterfolder/clustermin0

if [ -f $clusterfolder/clustersbelow"$mincluster".zip ] ; then 
rm clustersbelow"$mincluster".zip 
fi 
zip -r clustersbelow"$mincluster".zip clustersbelow"$mincluster"
rm -r clustersbelow"$mincluster"

cd $clusterfolder/"$outfolder"
zip -r mapped.zip mapped
rm -r mapped

zip -r polished.zip polished
rm -r polished

zip -r referenceFASTA.zip referenceFASTA
rm -r referenceFASTA

echo "Number of fasta seqs in $clusterfolder/"$barcodeID"_belowmin"$mincluster"seqs.fasta: $(grep -c ">" $clusterfolder/"$barcodeID"_belowmin"$mincluster"seqs.fasta)"

echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID
echo $(date)

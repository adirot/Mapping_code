#!/bin/bash


echo " "
echo " "
echo " -------------------------------------------------"
echo "| Read cleanning, trimming, mapping and counting |"
echo " -------------------------------------------------"
echo " "
echo " Look at the code comments before running, you need to choose your file locations and some more settings. The lines you need to change are marked by ##"

code_dir="$PWD"
## change this to your genome dir: 
genome_dir="${code_dir}/Genome_files"

## set a path to the input FASTA files:
inputFilesFolder="${code_dir}/raw_data"

## set a path to the output merged fastq files:
mergedFastqFolder="${code_dir}/output/merged_fastq"

## set a path to the output fastQC files:
fastqcFolder="${code_dir}/output/fastQC"

## give the path to the referance fasta file for alignment
referenceFile="${genome_dir}/KLY_ERCCtrimmed_plasmidpza.fasta" # KLY + ERCC + plasmid + fixed annotations using blast
referenceFileFolder="${genome_dir}" # needs to be thesame as the folder in "referenceFile"

## choose a base name for the referance index
referenceIndexBaseName=KLYmoreannotations

## choose clean floder name
cleanOutputFolder="${code_dir}/output/clean_fastq"

## choose name for alignment results folder
alignOutputFolder="${code_dir}/output/aligned_KLYmoreannotations"

## choose name for count results folder
countOutputFolder="${code_dir}/output/count_KLYmoreannotations"

## gtf file:
gtfFile="${genome_dir}/kly_more_names_ERCC_pzaplasmid.gtf"

## gtf file name
gtfFileName="kly_more_names_ERCC_pzaplasmid.gtf" # needs to be the same as the file name in "gtfFile"

# Create the output folders
mkdir -p "${fastqcFolder}"
mkdir -p "${alignOutputFolder}/log"
mkdir -p "${countOutputFolder}"
mkdir -p "${cleanOutputFolder}/log"
mkdir -p "${originalFastqFolder}"



echo " "
echo " "
echo "clean ployA and low quality and adptors"
echo "---------------------------------------"
cd "${originalFastqFolder}"

# -m: minimum length filtering
# -O 1 : trim end even if its of length 1. (since we want 1/0 missmach in the mapping this is what Yuval Nevo recomanded). 
# -q: quality edge trimming (3 end)
# -a: 3 end -g: 5 end
# -e: allowed error rate

# Set the library adaptors and indexies for cleaning
nebnext_5adaptor="GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC"
nebnext_3adaptor="GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC" # this is also the illumina read 1 adaptor


indexies=(
	ATCACG
	CGATGT
	TTAGGC
	TGACCA
	ACAGTG
	GCCAAT
	CAGATC
	ACTTGA
	GATCAG
	TAGCTT
	GGCTAC
	CTTGTA
#	ACATCG # index 2
#	TGGTCA # index 4
#	CACTGT # index 5
	)

cd $originalFastqFolder
# unzip the files
for FILE in *.fastq.gz
do
	gzip -d $FILE
done

# clean the reads using cutadapt
for FILE in *.fastq
do
	echo " "
	echo $FILE
	echo "3 end:"
	end3=$(reverse_compliment $nebnext_3adaptor)     
	echo $end3
	cutadapt -e 0.2 -m 15 -O 1 -q 20 -a "${end3}" -o "${FILE%.fastq}_clean_3end.fastq" $FILE > "${cleanOutputFolder}/log/${FILE%.fastq}_clean_log_3end.txt"
	
	#echo " "
	#echo "5 end (with index):"
	#end5=$(reverse_compliment "${indexies[$i]}")$(reverse_compliment $nebnext_5adaptor)
	#echo $end5

	echo "5 end (with wildcard index):"
	end5=$(reverse_compliment "NNNNNN")$(reverse_compliment $nebnext_5adaptor)
	echo $end5


	cutadapt -e 0.2 -m 15 -O 1 -q 20 -g "${end5}" -o "${FILE%.fastq}_clean.fastq" "${FILE%.fastq}_clean_3end.fastq" > "${cleanOutputFolder}/log/${FILE%.fastq}_clean_log_5end.txt"
	
	# move the files to the "cleanOutputFolder"
	mv *_clean*.fastq "${cleanOutputFolder}"
done

echo " "
echo " "
echo "Map to reference"
echo "----------------"

 
cd "${referenceFileFolder}"
echo "build index file for the referance sequence"
bowtie2-build "${referenceFile}" $referenceIndexBaseName 
mv "${referenceIndexBaseName}"*".bt2" "${cleanOutputFolder}"

cd "${cleanOutputFolder}"

for FILE in *_clean.fastq
do
	  FILEBaseName=$(basename ${FILE%_clean.fastq})
  	echo $FILEBaseName
  	bowtie2 --local --very-sensitive-local -a -x $referenceIndexBaseName -U $FILE -S "${FILEBaseName}_local_aligned.sam" --un "${FILEBaseName}_local_unmapped.fastq"
  	
  	mv "${FILEBaseName}_local_aligned.sam" "${alignOutputFolder}"
  	mv "${FILEBaseName}_local_unmapped.fastq" "${alignOutputFolder}"
done
mv "${referenceIndexBaseName}"*".bt2" "${referenceFileFolder}"


echo " "
echo " "
echo "Calculate expresion levels"
echo "--------------------------"
cp "${gtfFile}" "${alignOutputFolder}"
cd "${alignOutputFolder}"


for FILE in *_aligned.sam
do
  	echo $FILE
  	
    htseq-count $FILE $gtfFileName --stranded=reverse -i transcript_id -a 0 --additional-attr gene_name -t exon > "${FILE%.sam}_count_stranded_reverse_a_0.tsv"
    mv "${FILE%.sam}_count_stranded_reverse_a_0.tsv" "${countOutputFolder}"
  	
  	
done

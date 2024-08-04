#!/bin/bash


echo " "
echo " "
echo " -------------------------------------------------"
echo "| Read cleanning, trimming, mapping and counting |"
echo " -------------------------------------------------"
echo " "
echo " Look at the code comments before running, you need to choose your file locations and some more settings. The lines you need to change are marked by ##"

### Inputs. Sets the path to input data and Output folders###
code_dir="$PWD"
## change this to your genome dir: 
genome_dir="${code_dir}/Genome_files"

## set a path to the input FASTQ files:
#inputFilesFolder="${code_dir}/FASTQ_raw_data/Experiment1" # run on the data from Experiment1
#inputFilesFolder="${code_dir}/FASTQ_raw_data/Experiment2" # run on the data from Experiment2
#inputFilesFolder="${code_dir}/FASTQ_raw_data/Exponential"
inputFilesFolder="${code_dir}/FASTQ_raw_data/Time_in_SHX"
#inputFilesFolder="${code_dir}/FASTQ_raw_data/subsample" # for debugging. delete before upload

EXPname=$(basename "$inputFilesFolder") # folder name of the input file (Experiment1 or 2)

## Path to output folder:
outputFolder="${code_dir}/map_count_outputs/${EXPname}"

## set a path to the output fastQC files:
fastqcFolder="${outputFolder}/fastQC"

## give the path to the referance fasta file for alignment
#referenceFile="${genome_dir}/KLY_ERCCtrimmed_plasmidpza.fasta"
referenceFile="${genome_dir}/KLY.fasta"


## gtf file path and name:
#gtfFile="${genome_dir}/KLY_ERCCtrimmed_plasmidpza.gtf"
gffFile="${genome_dir}/NZ_CP008801.1_corrected.gff"

## ERCC Fasta and gff file - these files will be concatenated to the other fasta and gff files
ERCC_fasta="${genome_dir}/ERCC.fasta"
ERCC_gff="${genome_dir}/ERCC.gff"

## plasmid Fasta and gff files - these files will be concatenated to the other fasta and gff files
plasmid_fasta="${genome_dir}/pza21mCherry.fasta"
plasmid_gff="${genome_dir}/pza21mCherry.gff"

## Extra genes Fasta and gff files - these files will be concatenated to the other fasta and gff files
extra_fasta="${genome_dir}/Extra_genes.fasta"
extra_gff="${genome_dir}/Extra_genes.gff"

## choose a base name for the referance index
referenceIndexBaseName=KLY_ERCC

## choose clean fastq floder name
cleanOutputFolder="${outputFolder}/clean_fastq"

## choose name for alignmed results folder
alignOutputFolder="${outputFolder}/aligned"

## choose name for count results folder
countOutputFolder="${outputFolder}/counts"

# Concatenate ERCC and plasmid fasta files to the genome fasta
awk 'FNR==1{print ""}1' $referenceFile $ERCC_fasta $plasmid_fasta $extra_fasta > "${genome_dir}/referenceFile_combined.fasta"
referenceFile="${genome_dir}/referenceFile_combined.fasta"
 

# Concatenate ERCC and plasmid gff files to the genome gff
awk 'FNR==1{print ""}1' $gffFile $ERCC_gff $plasmid_gff $extra_gff > "${genome_dir}/gffFile_combined.gff"
gffFile="${genome_dir}/gffFile_combined.gff"



gffFileName=$(basename "$gffFile") # needs to be the same as the file name in "gffFile"



# Create the output folders
mkdir -p "${alignOutputFolder}/log"
mkdir -p "${countOutputFolder}"
mkdir -p "${cleanOutputFolder}/log"

: <<CLEAN
echo " "
echo " "
echo "clean ployA tail, low quality, and adptors"
echo "---------------------------------------"
cd "${inputFilesFolder}"


# Set the library adaptors and indexies for cleaning
nebnext_5adaptor="GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC"
nebnext_3adaptor="GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC" # this is also the illumina read 1 adaptor

function reverse_compliment {

	echo $1 | tr ACGTacgtN TGCAtgcaN | rev

}

# clean the reads using cutadapt
for FILE in *.fastq
do
	# Cutadapt parameters:	
	# -m: minimum length filtering
	# -O 1 : trim end even if its of length 1. (since we want 1/0 missmach in the mapping this is what Yuval Nevo recomanded). 
	# -q: quality edge trimming (3 end)
	# -a: 3 end -g: 5 end
	# -e: allowed error rate
	echo " "
	echo $FILE
	echo "3 end:"
	end3=$(reverse_compliment $nebnext_3adaptor)     
	echo $end3
	cutadapt -e 0.2 -m 15 -O 1 -q 20 -a "${end3}" -o "${FILE%.fastq}_clean_3end.fastq" $FILE > "${cleanOutputFolder}/log/${FILE%.fastq}_clean_log_3end.txt"
	
	echo "5 end (with wildcard index):"
	end5=$(reverse_compliment "NNNNNN")$(reverse_compliment $nebnext_5adaptor)
	echo $end5

	cutadapt -e 0.2 -m 15 -O 1 -q 20 -g "${end5}" -o "${FILE%.fastq}_clean.fastq" "${FILE%.fastq}_clean_3end.fastq" > "${cleanOutputFolder}/log/${FILE%.fastq}_clean_log_5end.txt"
	
	# clean up temporary files
	rm *_clean_3end.fastq
	# move the files to the "cleanOutputFolder"
	mv *_clean.fastq "${cleanOutputFolder}"
	
	
done
CLEAN

echo " "
echo " "
echo "Map to reference"
echo "----------------"

 
cd "${genome_dir}"
echo "build index file for the referance sequence"
bowtie2-build "${referenceFile}" $referenceIndexBaseName 
mv "${referenceIndexBaseName}"*".bt2" "${cleanOutputFolder}"

cd "${cleanOutputFolder}"

for FILE in *_clean.fastq
do
	FILEBaseName=$(basename ${FILE%_clean.fastq})
  	echo $FILEBaseName
  	bowtie2 --local --very-sensitive-local -a -x $referenceIndexBaseName -U $FILE -S "${FILEBaseName}.sam" --un "${FILEBaseName}_unmapped.fastq"
  	
  	mv "${FILEBaseName}.sam" "${alignOutputFolder}"
  	mv "${FILEBaseName}_unmapped.fastq" "${alignOutputFolder}"
done

echo " "
echo " "
echo "alignment QC"
echo "------------"
cd "${alignOutputFolder}"
for FILE in *.sam
do
	echo $FILE
	#htseq-qa $FILE 
	
	samtools stats $FILE > "log/${FILE%.sam}_samtools_stat.txt"
	samtools flagstat $FILE > "log/${FILE%.sam}_samtools_flagstat.txt"
		
	
done

echo " "
echo " "
echo "Calculate expresion levels"
echo "--------------------------"
cp "${gffFile}" "${alignOutputFolder}"
cd "${alignOutputFolder}"


for FILE in *.sam
do
  	echo $FILE
  	
    htseq-count $FILE $gffFileName --stranded=reverse -i Name -a 0 -t gene > "${FILE%.sam}.tsv"
	htseq-count $FILE $gffFileName --stranded=yes -i Name -a 0 -t gene > "${FILE%.sam}_reverse.tsv"
    mv "${FILE%.sam}.tsv" "${countOutputFolder}"
  	mv "${FILE%.sam}_reverse.tsv" "${countOutputFolder}"
  	
done

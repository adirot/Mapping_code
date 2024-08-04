## Read_cleanning_mapping_and_counting.bash
This is the code used to get from raw FASTQ files to normalized read counts.

## FASQ_raw_data
The raw FASTQ files are not provided in this github since they are large, but you can download the whole github, and manually add them (see read me files in folder). 

## Experiment metadata
Metadata for the samples. The code uses this metadata for calculations of the normalizations. 

## Normalization of raw counts
This folder contains a jupyter notebook that preforms the different normalizations (downsampling, ERCC spike-in nomalization, and RPK after down-sampling) of the raw read counts. It needs as an input the raw counts that are the output of HTseq-count. If you ran the code in "Read_cleanning_mapping_and_counting" you should have a folder with these files under "map_count_outputs/Experiment#/counts".
After you run this code, new folders will be created containg csv files with the normalized read counts. 
For DEseq2 normalization, see the subfolder "DEseq normalization code".

### DEseq normalization code
See the read me file in this subfolder. This folder contains the R code used to normalize the reads using DESeq2. Note that this code only normalizes the reads, but does not compute the log2 fold changes. 

## Genome Files
The Genome files folder contains the genome data necessary for mapping of the raw reads to the genome. 
The two files are 

KLY_ERCCtrimmed_plasmidpza.fasta: The sequence to map to. It includes the KLY sequence (CP008801.1)*, the PZA21mcherry plasmid, the ERCC sequnces.
                                    
                                    * A few genes are added, which are in the KLYR but not the KLY: spc, tetR
  
  
kly_more_names.gtf: The annotation file is based on the annotation file of the KLY sequence (CP008801.1)*, and it also annotates the PZA21mcherry plasmid, the ERCC     
                         sequnces
                        
                        * A few genes are added to this annotation, which are in the KLYR but not the KLY: spc, tetR. Also, some genes names are changed in the       
                          "kly_more_names" to thier common names. 
  


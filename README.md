## Read_cleanning_mapping_and_counting
This is the code used to get from raw FASTQ files to normalized read counts.

## FASQ_raw_data
The raw FASTQ files are not provided in this github since they are large, but you can download the whole github, and manually add them (see read me files in folder). 

## Experiment metadata
Metadata for the samples. The code uses this metadata for calculations of the normalizations. 

## Genome Files
The Genome files folder contains the genome data necessary for mapping of the raw reads to the genome. 
The two files are 

KLY_ERCCtrimmed_plasmidpza.fasta: The sequence to map to. It includes the KLY sequence (CP008801.1)*, the PZA21mcherry plasmid, the ERCC sequnces.
                                    
                                    * A few genes are added, which are in the KLYR but not the KLY: spc, tetR
  
  
kly_more_names.gtf: The annotation file is based on the annotation file of the KLY sequence (CP008801.1)*, and it also annotates the PZA21mcherry plasmid, the ERCC     
                         sequnces
                        
                        * A few genes are added to this annotation, which are in the KLYR but not the KLY: spc, tetR. Also, some genes names are changed in the       
                          "kly_more_names" to thier common names. 
  


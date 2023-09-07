# Mapping_code
This is the code used to get from raw FASTQ files to normalized read counts.


## Genome Files
The Genome files folder contains the genome data necessary for mapping of the raw reads to the genome. 
The two files are 
  KLY_ERCCtrimmed_plasmidpza.fasta: The sequence to map to. It includes the KLY sequence (CP008801.1)*, the PZA21mcherry plasmid, the ERCC sequnces.
                                    * A few genes are added, which are in the KLYR but not the KLY: spc, tetR
  kly_more_names.gff3: The annotation file is based on the annotation file of the KLY sequence (CP008801.1)*, and it also annotates the PZA21mcherry plasmid, the ERCC     
                         sequnces
                        * A few genes are added to this annotation, which are in the KLYR but not the KLY: spc, tetR. Also, some genes names are changed in the       
                          "kly_more_names" to thier common names. 
  

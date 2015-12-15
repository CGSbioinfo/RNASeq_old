# RNASeq
RNASeq analysis pipeline
Requirments: trim galore, FASTQC, STAR, HTSEQ, RNASeqQC, Picard tools
Pipeline for the analysis of RNASeq data. It is divided in several python scripts.
Each script perform one step:
1) Prepare the folder for the rest of the analysis
2) Trim reads and perform fastqc on the data 
3) Map the reads
4) count the read

The "small" version is for small rna seq library.
In order to run the script a parameter file with the following parameters is required:

*Working directory = #path to your analysis folder
GTF File = #path to your gtf file
Reference Genome = #path to the reference genome FOLDER **needs to be formated for STAR**
Number of samples = #number of sample in the analysis
BedFile = #path to Bed file for mapping qc (position of all gene for the genome)  for junction saturation (RNASeqQC)
BedFile10K = # subset of 10k for RNA body coverage (RNASeqQC)
refFlat = #Reformated bed file used for picard tools to check strand and reads mapping location qc 
rRNA_interval_list = # list simillar to bed of rRNA intervals for QC
strand = #Strand setting for HTSEQ*


####################################
##### BS_analysis_workflow #########
####################################

1. check fastqc (fastqc.sh <file.fastq>)
#base content
#base quality
#filter length

2. trimm quality (trimmbase.sh <file.fastq>; trimmquality.sh <file.fastq>; [trimmfull.sh <filename.fastq>])

#firstCut=45
#lastCut=152
#quality=20
#length=30
	
3.fastqc.sh cut

4.bismark_genome_preparation (bismarkprep.sh <genome_folder_name>)

5.bismark align (bismark.sh <genome_folder_name> <query.fastq>)

6.bismark_methylation_extractor -s <results.sam> -o <output_folder>

#bismark (genome preparation; align; report; methylation extractor)
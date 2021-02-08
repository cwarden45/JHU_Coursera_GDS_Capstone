import sys
import re
import os

#NOTE: This code is based upon scripts in the following locations:
#https://storage.googleapis.com/rna_deg/bucket_file_list.txt (from https://sourceforge.net/projects/rnaseq-deg-methodlimit/files/)
#https://github.com/cwarden45/RNAseq_templates/tree/master/TopHat_Workflow

##required samples
unfinishedSamples = ("SRR1554534","SRR1554535","SRR1554536","SRR1554539","SRR1554556","SRR1554561","SRR1554537","SRR1554538","SRR1554541","SRR1554566","SRR1554567","SRR1554568")

##extra samples
#unfinishedSamples = ("SRR1554544","SRR1554546","SRR1554549","SRR1554551","SRR1554553","SRR1554554")

threads = "4"
ref = "/home/cwarden/Ref/STAR/hg19_Bioconductor_UCSC_GTF_50bp"
alignmentFolder = "hg19_STAR_Alignment"
readsFolder = "."

command = "mkdir " + alignmentFolder
os.system(command)

fileResults = os.listdir(readsFolder)

for file in fileResults:
	result = re.search("(.*)_1.fastq.gz$",file)
	
	if result:
		sample = result.group(1)
		
		if (sample in unfinishedSamples):
			print sample
				
			outputSubfolder = alignmentFolder +"/" + sample
			command = "mkdir " + outputSubfolder
			os.system(command)
										
			read1 = re.sub(".gz$","",readsFolder + "/" + file)
			command = "gunzip -c " + read1 + ".gz > " + read1
			os.system(command)

			read2 = re.sub("_1.fastq","_2.fastq",read1)
			command = "gunzip -c " + read2 + ".gz > " + read2
			os.system(command)
			
			starPrefix = outputSubfolder +"/" + sample + "_"
			command = "/opt/STAR-2.7.2d/bin/Linux_x86_64_static/STAR --genomeDir " + ref+ " --readFilesIn " + read1 + " " + read2 + " --runThreadN " +threads+ " --outFileNamePrefix " + starPrefix + " --twopassMode Basic --outSAMstrandField intronMotif"
			os.system(command)

			starSam = starPrefix + "Aligned.out.sam"
			alnBam = outputSubfolder + "/aligned.bam"
			command = "/opt/samtools/samtools view -bS " + starSam + " > " + alnBam
			os.system(command)
			
			userBam = alignmentFolder + "/" + sample + ".bam"
			command = "/opt/samtools/samtools sort -o " + userBam + " " + alnBam
			os.system(command)
			
			command = "rm " + starSam
			os.system(command)

			command = "rm " + alnBam
			os.system(command)

			command = "/opt/samtools/samtools index " + userBam
			os.system(command)
			
			command = "rm " + read1
			os.system(command)
			
			command = "rm " + read2
			os.system(command)
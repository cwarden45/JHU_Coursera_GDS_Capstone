import sys
import re
import os
import subprocess

#NOTE: This code is based upon scripts in the following locations:
#https://storage.googleapis.com/rna_deg/bucket_file_list.txt (from https://sourceforge.net/projects/rnaseq-deg-methodlimit/files/)
#https://github.com/cwarden45/RNAseq_templates/tree/master/TopHat_Workflow

finishedSamples = ()
threads = "4"
ref = "/home/cwarden/Ref/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome"
alignmentFolder = "hg19_TopHat2_Alignment"
readsFolder = "."
strandType = "no"

command = "mkdir " + alignmentFolder
os.system(command)

fileResults = os.listdir(readsFolder)

for file in fileResults:
	result = re.search("(.*)_1.fastq.gz$",file)
	
	if result:
		sample = result.group(1)
		
		if (sample not in finishedSamples):
			print sample
				
			outputSubfolder = alignmentFolder +"/" + sample
			command = "mkdir " + outputSubfolder
			os.system(command)
										
			read1 = readsFolder + "/" + file
			read2 = re.sub("_1.fastq","_2.fastq",read1)

			tophatStrand = ""
			if strandType == "no":
				tophatStrand = "fr-unstranded"
			elif strandType == "reverse":
				tophatStrand = "fr-firststrand"
			elif strandType == "yes":
				tophatStrand = "fr-secondstrand"
			else:
				print "Need to provide TopHat mapping for strand: " + strandType
				sys.exit()
			
			command = "/opt/tophat-2.1.1.Linux_x86_64/tophat2 -o " + outputSubfolder + " -p " + threads + " --no-coverage-search --library-type " + tophatStrand + " " + ref + " " + read1 + " " + read2
			os.system(command)
									
			topHatBam = outputSubfolder + "/accepted_hits.bam"																			
			userBam = alignmentFolder + "/" + sample + ".bam"
			
			command = "mv " + topHatBam + " " + userBam
			os.system(command)

			command = "/opt/samtools/samtools index " + userBam
			os.system(command)

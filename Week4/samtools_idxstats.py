import sys
import re
import os

#copied and modified from run_RSeQC.py (in Week 5)

alignmentFolder = "hg19_STAR_Alignment"
output_folder = "samtools_idxstats_output"

command = "mkdir " + output_folder
os.system(command)

fileResults = os.listdir(alignmentFolder)

for file in fileResults:
	result = re.search(".bam$",file)
	bam_file = os.path.join(alignmentFolder, file)
	
	if result:
		sample = re.sub(".bam$","",file)
		print sample
		
		statFile = output_folder + "/" + sample + "_alignment_summary-ALL.txt"
		command = "/opt/samtools/samtools idxstats " + bam_file + " > " + statFile
		os.system(command)
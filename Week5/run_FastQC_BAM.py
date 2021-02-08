import sys
import re
import os

#copied and modified from run_RSeQC.py (this week)

alignmentFolder = "../Week4/hg19_STAR_Alignment"
output_folder = "FastQC_Output_BAM"

command = "mkdir " + output_folder
os.system(command)

fileResults = os.listdir(alignmentFolder)

for file in fileResults:
	result = re.search(".bam$",file)
	bam_file = os.path.join(alignmentFolder, file)
	
	if result:
		sample = re.sub(".bam$","",file)
		print sample
		
		command = "/opt/FastQC/fastqc -f bam_mapped -o "+output_folder+ " " + bam_file
		os.system(command)
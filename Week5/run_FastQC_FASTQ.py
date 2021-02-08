import sys
import re
import os
import subprocess

#copied and modified STAR_PE_Alignment.py (Week 4) and run_FastQC_BAM.py (this week)

finishedSamples = ()
readsFolder = "../Week4"
output_folder = "FastQC_Output_FASTQ"

command = "mkdir " + output_folder
os.system(command)

fileResults = os.listdir(readsFolder)

for file in fileResults:
	result = re.search("(.*_\d).fastq.gz$",file)
	fastq_file = os.path.join(readsFolder, file)
	
	if result:
		sample = result.group(1)
		
		if (sample not in finishedSamples):
			print sample
				
			command = "/opt/FastQC/fastqc -o "+output_folder+ " " + fastq_file
			os.system(command)
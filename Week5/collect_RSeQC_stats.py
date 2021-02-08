import sys
import re
import os

#NOTE: This code is based upon scripts in the following locations:
#https://storage.googleapis.com/rna_deg/bucket_file_list.txt (from https://sourceforge.net/projects/rnaseq-deg-methodlimit/files/)
#https://github.com/cwarden45/RNAseq_templates/tree/master/TopHat_Workflow/collect_RSeQC_stats.py

alignmentFolder = "../Week4/hg19_STAR_Alignment"
qcFolder = "RSeQC_Output"
		
fileResults = os.listdir(alignmentFolder)

statsFile = "RSeQC_stats.txt"
statHandle = open(statsFile, 'w')
text = "Sample\tReverse.Strand.Bias\tmedTIN\n"
statHandle.write(text)

for file in fileResults:
	result = re.search("(.*).bam$",file)
	
	if result:
		sample = result.group(1)
		print sample
		text = sample
		
		strandStat = qcFolder + "/" + sample + "_infer_strand.txt"
		inHandle = open(strandStat)
		lines = inHandle.readlines()
					
		for line in lines:
			line = re.sub("\n","",line)
			line = re.sub("\r","",line)
			

			result2 = re.search("Fraction of reads explained by \"1\+-,1-\+,2\+\+,2--\": (.*)",line)
				
			if result2:
				formatFrac =  '{0:.2g}'.format(float(result2.group(1)))
				text = text + "\t" +  formatFrac

			
	
		tinStat = qcFolder + "/" + sample + ".summary.txt"
		lineCount = 0
		inHandle = open(tinStat)
		lines = inHandle.readlines()
					
		for line in lines:
			line = re.sub("\n","",line)
			line = re.sub("\r","",line)
			
			lineCount += 1
			
			if lineCount == 2:
				lineInfo = line.split("\t")
				text = text + "\t" + '{0:.3g}'.format(float(lineInfo[2]))
				
		text = text + "\n"
		statHandle.write(text)

import sys
import re
import os

#NOTE: This code is based upon scripts in the following locations:
#https://storage.googleapis.com/rna_deg/bucket_file_list.txt (from https://sourceforge.net/projects/rnaseq-deg-methodlimit/files/)
#https://github.com/cwarden45/RNAseq_templates/tree/master/TopHat_Workflow/run_RSeQC.py


alignmentFolder = "../Week4/hg19_STAR_Alignment"
bed_file = "/home/cwarden/Ref/RSeQC/hg19.HouseKeepingGenes.bed"
output_folder = "RSeQC_Output"

command = "mkdir " + output_folder
os.system(command)

fileResults = os.listdir(alignmentFolder)

finishedSamples = ()

for file in fileResults:
	result = re.search(".bam$",file)
	fullPath = os.path.join(alignmentFolder, file)
	
	if result:
		sample = re.sub(".bam$","",file)
		sortResult = re.search(".name.sort.bam",file)
		if (sample not in finishedSamples) and (not sortResult):
			print sample
			
			print "Determining Strand for Housekeeping Genes"
			strandStat = output_folder + "/" + sample + "_infer_strand.txt"
			command = "infer_experiment.py -r " + bed_file + " -i " + fullPath + " > " + strandStat
			os.system(command)

			print "Calculating TIN scores"
			command = "tin.py -r " + bed_file + " -i " + fullPath
			os.system(command)
			
			tinOut1 = sample + ".summary.txt"
			command = "mv " + tinOut1 + " " + output_folder + "/" + tinOut1
			os.system(command)
			
			tinOut2 = sample + ".tin.xls"
			command = "mv " + tinOut2 + " " + output_folder + "/" + tinOut2
			os.system(command)
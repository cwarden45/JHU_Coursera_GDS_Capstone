import sys
import re
import os

#portions copied and modified from run_FastQC_FASTQ.py and collect_RSeQC_stats.py (this week)

readsFolder = "../Week4"
qcFolder = "FastQC_Output_FASTQ"
		
fileResults = os.listdir(readsFolder)

statsFile = "FastQC_stats-extra_3var.txt"
statHandle = open(statsFile, 'w')
text = "Sample\tPercent_Average_Q30\tAverage_GC\tPercent_Unique\n"
statHandle.write(text)

for file in fileResults:
	result = re.search("(.*_\d).fastq.gz$",file)
	fastq_file = os.path.join(readsFolder, file)
	
	if result:
		sample = result.group(1)
		print sample
		text = sample
		
		Q30_count = 0
		Q30_total = 0
		GC_total = 0
		GC_count = 0
		percent_unique = 0
		
		extracted_fastqc = qcFolder + "/" + sample + "_fastqc"		
		
		fullInfo = extracted_fastqc + "/fastqc_data.txt"
		
		seq_qual_section_flag=0
		gc_section_flag=0
		dup_section_flag=0
		
		inHandle = open(fullInfo)
		lines = inHandle.readlines()
					
		for line in lines:
			line = re.sub("\n","",line)
			line = re.sub("\r","",line)

			resultQual = re.search(">>Per sequence quality scores",line)
			resultGC = re.search(">>Per sequence GC content",line)
			resultDup = re.search(">>Sequence Duplication Levels",line)
			resultReset = re.search(">>END_MODULE",line)

			if resultQual:
				print "...start reading section to calculate percent average Q30"
				seq_qual_section_flag=1
			elif resultGC:
				print "...start reading section to calculate average GC"
				gc_section_flag=1
			elif resultDup:
				print "...start reading section to calculate duplicated / unique sequence content"
				dup_section_flag=1
			elif resultReset:
				if seq_qual_section_flag==1:
					print "...finish reading section to calculate percent average Q30"
					seq_qual_section_flag=0
				
				if gc_section_flag==1:
					print "...finish reading section to calculate average GC"
					gc_section_flag=0

				if dup_section_flag==1:
					print "...finish reading section to calculate duplicated / unique sequence content"
					dup_section_flag=0
			else:
				if seq_qual_section_flag:
					headerResult = re.search("^#",line)
					if not headerResult:
						lineInfo = line.split("\t")
						temp_qual = int(lineInfo[0])
						temp_count = float(lineInfo[1])
						
						if temp_qual >= 30:
							Q30_count = Q30_count + temp_count
						
						Q30_total = Q30_total + temp_count

				if gc_section_flag:
					headerResult = re.search("^#",line)
					if not headerResult:
						lineInfo = line.split("\t")
						temp_value = float(lineInfo[0])
						temp_count = float(lineInfo[1])
						
						GC_total = GC_total + temp_value*temp_count
						GC_count = GC_count + temp_count

				if dup_section_flag:
					dup1result = re.search("^1\t",line)
					if dup1result:
						lineInfo = line.split("\t")
						percent_unique=float(lineInfo[2])
		
		inHandle.close()
		
		if Q30_total == 0:
			print "Issue with Q30 calculation:"
			print "Total = " + str(Q30_total)
			print "Count = " + str(Q30_count)
			sys.exit()
		percent_average_Q30 = 100 * Q30_count / Q30_total
		
		if GC_count == 0:
			print "Issue with GC calculation:"
			print "Total = " + str(GC_total)
			print "Count = " + str(GC_count)
			sys.exit()
		average_GC = GC_total / GC_count
		
		text = text+"\t"+'{0:.4g}'.format(percent_average_Q30)+"\t"+'{0:.4g}'.format(average_GC)+"\t"+'{0:.4g}'.format(percent_unique)+"\n"
		statHandle.write(text)
		
statHandle.close()

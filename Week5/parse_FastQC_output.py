import sys
import re
import os

#portions copied and modified from run_FastQC_FASTQ.py and collect_RSeQC_stats.py (this week)

readsFolder = "../Week4"
qcFolder = "FastQC_Output_FASTQ"
		
fileResults = os.listdir(readsFolder)

statsFile = "FastQC_stats.txt"
statHandle = open(statsFile, 'w')
text = "Sample\tTotal_Reads\tBasic_Statistics\tBase_Qual\tSeq_Qual\tSeq_Content\tGC\tN_content\tSeq_Length\tSeq_Dup\tEnrich_String\tAdapter_Content\n"
statHandle.write(text)

for file in fileResults:
	result = re.search("(.*_\d).fastq.gz$",file)
	fastq_file = os.path.join(readsFolder, file)
	
	if result:
		sample = result.group(1)
		print sample
		text = sample
		
		extracted_fastqc = qcFolder + "/" + sample + "_fastqc"
		
		if not os.path.isdir(extracted_fastqc):
			fastqc_zip = extracted_fastqc + ".zip"
			command = "unzip "+fastqc_zip+ " -d " + qcFolder
			os.system(command)			
		
		fullInfo = extracted_fastqc + "/fastqc_data.txt"
		fastqcSummary = extracted_fastqc + "/summary.txt"
		
		#extract read counts
		inHandle = open(fullInfo)
		lines = inHandle.readlines()
					
		for line in lines:
			line = re.sub("\n","",line)
			line = re.sub("\r","",line)
			

			result2 = re.search("Total Sequences\t(\d+)",line)
			if result2:
				starting_count = result2.group(1)
				text = text + "\t" +  starting_count
		
		inHandle.close()
		
		#extract summary flags
		Basic_Statistics=""
		Base_Qual=""
		Seq_Qual=""
		Seq_Content=""
		GC=""
		N_content=""
		Seq_Length=""
		Seq_Dup=""
		Enrich_String=""
		Adapter_Content=""


		inHandle = open(fastqcSummary)
		lines = inHandle.readlines()
					
		for line in lines:
			line = re.sub("\n","",line)
			line = re.sub("\r","",line)
			
			result2 = re.search("(\w+)\tBasic Statistics",line)
			if result2:
				Basic_Statistics=result2.group(1)

			result2 = re.search("(\w+)\tPer base sequence quality",line)
			if result2:
				Base_Qual=result2.group(1)

			result2 = re.search("(\w+)\tPer sequence quality scores",line)
			if result2:
				Seq_Qual=result2.group(1)

			result2 = re.search("(\w+)\tPer base sequence content",line)
			if result2:
				Seq_Content=result2.group(1)

			result2 = re.search("(\w+)\tPer sequence GC content",line)
			if result2:
				GC=result2.group(1)

			result2 = re.search("(\w+)\tPer base N content",line)
			if result2:
				N_content=result2.group(1)

			result2 = re.search("(\w+)\tSequence Length Distribution",line)
			if result2:
				Seq_Length=result2.group(1)

			result2 = re.search("(\w+)\tSequence Duplication Levels",line)
			if result2:
				Seq_Dup=result2.group(1)

			result2 = re.search("(\w+)\tOverrepresented sequences",line)
			if result2:
				Enrich_String=result2.group(1)

			result2 = re.search("(\w+)\tAdapter Content",line)
			if result2:
				Adapter_Content=result2.group(1)
		
		inHandle.close()

		text = text+"\t"+Basic_Statistics+"\t"+Base_Qual+"\t"+Seq_Qual+"\t"+Seq_Content+"\t"+GC+"\t"+N_content+"\t"+Seq_Length+"\t"+Seq_Dup+"\t"+Enrich_String+"\t"+Adapter_Content
				
		text = text + "\n"
		statHandle.write(text)
		
statHandle.close()

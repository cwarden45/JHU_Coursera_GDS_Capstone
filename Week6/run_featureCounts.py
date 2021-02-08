import sys
import re
import os

#copied and modified from run_FastQC_BAM.py (Week 5)

alignmentFolder = "../Week4/hg19_STAR_Alignment"
GTF="/mnt/usb8/RNAseq_DEG/Bioconductor_GTF/TxDb_hg19_gene.gtf"

fileResults = os.listdir(alignmentFolder)

for file in fileResults:
	result = re.search(".bam$",file)
	bam_file = os.path.join(alignmentFolder, file)
	
	if result:
		sample = re.sub(".bam$","",file)
		print sample
		
		outputPrefix = sample
		#add "-p" to count fragments instead of reads for paired-end data
		#similarly, add "-B" and "-C" to be more restrictive (require alignment of both ends, and require alignment on same strand on same chromosome, respectively)
		command = "/opt/subread-2.0.1-source/bin/featureCounts -p -B -C -s 0 -a "+GTF+ " -o "+outputPrefix+ " " + bam_file
		os.system(command)
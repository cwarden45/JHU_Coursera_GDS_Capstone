flagstat_folder = "samtools_flagstat_output"
metadata_file = "metadata_full.txt"
fastqc_file = "../Week5/FastQC_stats.txt"

###metadata processing
metadata.table = read.table(metadata_file, head=T, sep="\t")
print(dim(metadata.table))
#extra statistics for report (before filtering for downloaded samples)
metadata.table = metadata.table[-grep("\\dC1$",metadata.table$SampleID.Paper),]#count donors rather than samples - so, exclude cytosol samples for this purpose.
print(dim(metadata.table))
print(table(metadata.table$Race))
print(table(metadata.table$Sex))

pointCol = rep("black",nrow(metadata.table))
pointCol[metadata.table$Group == "Adult"]="blue"
pointCol[metadata.table$Group == "Fetal"]="red"
pointCol[metadata.table$Group == "Infant"]="orange"

metadata.table$Group = as.character(metadata.table$Group)
metadata.table$Group[metadata.table$Group == "Fetal:Cytosol"]="Fetal"
metadata.table$Group[metadata.table$Group == "Fetal:Nucleus"]="Fetal"
metadata.table$Group[metadata.table$Group == "Adult:Cytosol"]="Adult"
metadata.table$Group[metadata.table$Group == "Adult:Nucleus"]="Adult"
metadata.table$Group[metadata.table$Group == "50+"]="Senior"
group_levels = c("Fetal","Infant","Child","Teen","Adult","Senior")
metadata.table$Group = factor(metadata.table$Group,
								levels = group_levels)
print(table(metadata.table$Group))

png("age_check.png")
plot(jitter(as.numeric(metadata.table$Group)),metadata.table$Age,
	pch=16, col=pointCol, xaxt='n',
	xlab="Group", ylab="Age")
mtext(group_levels, side=1, at=c(1:length(group_levels)),
		las=2, line=1)
dev.off()

metadata.table = metadata.table[metadata.table$Download == "YES",]
print(dim(metadata.table))

print(table(metadata.table$Race, metadata.table$Group))
print(table(metadata.table$Sex, metadata.table$Group))

###FastQC processing
fastqc.table = read.table(fastqc_file, head=T, sep="\t")
R1 = paste(metadata.table$SRR,"_1",sep="")

total_pairs = fastqc.table$Total_Reads[match(R1,fastqc.table$Sample)]

###samtools processing
R1_counts = c()
supplemental_counts = c()

for (i in 1:nrow(metadata.table)){
	sampleID = metadata.table$SRR[i]
	flagstat_file = paste(flagstat_folder,"/",sampleID,"_alignment_summary.txt",sep="")
	flagstat.content = readLines(flagstat_file)
	
	reg.obj = gregexpr("^(\\d+)",flagstat.content[7])
	R1_counts[i] = as.numeric(unlist(regmatches(flagstat.content[7], reg.obj)))

	reg.obj = gregexpr("^(\\d+)",flagstat.content[2])
	supplemental_counts[i] = as.numeric(unlist(regmatches(flagstat.content[2], reg.obj)))
}#end for (i in 1:nrow(metadata.table))

output.table = data.frame(Sample=metadata.table$SRR, Group=metadata.table$Group,
							Total=total_pairs,
							Aligned_Count=R1_counts,
							Aligned_Percent=round(100*R1_counts/total_pairs, digits=4),
							Supplemental=supplemental_counts)
write.table(output.table,"alignment_stats-samtools_flagstat.txt", quote=F, sep="\t", row.names=F)

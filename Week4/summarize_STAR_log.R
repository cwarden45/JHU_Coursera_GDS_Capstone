alignment_folder = "hg19_STAR_Alignment"
GitHub_folder = "STAR_Alignment_Log_Files"
metadata_file = "metadata_full.txt"
fastqc_file = "../Week5/FastQC_stats.txt"

command = paste("mkdir ",GitHub_folder,sep="")
system(command)

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
Unique_Aligned_Count = c()
Unique_Aligned_Percent = c()

for (i in 1:nrow(metadata.table)){
	sampleID = metadata.table$SRR[i]
	log_file = paste(alignment_folder,"/",sampleID,"/",sampleID,"_Log.final.out",sep="")
	log.content = readLines(log_file)

	log_copy = paste(GitHub_folder,"/",sampleID,"_Log.final.out",sep="")
	command = paste("cp ",log_file," ",log_copy,sep="")
	system(command)
	
	reg.obj = gregexpr("(\\d+)",log.content[9])
	Unique_Aligned_Count[i] = as.numeric(unlist(regmatches(log.content[9], reg.obj)))

	reg.obj = gregexpr("(\\d+\\.\\d+)",log.content[10])
	Unique_Aligned_Percent[i] = as.numeric(unlist(regmatches(log.content[10], reg.obj)))
}#end for (i in 1:nrow(metadata.table))

output.table = data.frame(Sample=metadata.table$SRR, Group=metadata.table$Group,
							Total=total_pairs,
							Unique_Aligned_Count, Unique_Aligned_Percent)
write.table(output.table,"alignment_stats.txt", quote=F, sep="\t", row.names=F)

pointCol = rep("black",nrow(metadata.table))
pointCol[metadata.table$Group == "Adult"]="blue"
pointCol[metadata.table$Group == "Fetal"]="red"
pointCol[metadata.table$Group == "Infant"]="orange"

png("alignment_comparison-starting.png")
plot(metadata.table$totalMapped.paper,total_pairs, pch=16, col=pointCol,
	xlab="Reported Mapped (Paper)",
	ylab="SRA Downloaded (This Project)")
abline(a=0, b=1,col="gray")
legend("top",legend=c("x=y","linear fit"),
		col=c("gray","black"), lwd=2, ncol=2,
		xpd=T, inset = -0.1)
abline(lm(total_pairs~metadata.table$totalMapped.paper),col="black", lwd=2)
dev.off()

png("alignment_comparison-STAR_log_unique.png")
plot(metadata.table$totalMapped.paper, Unique_Aligned_Count, pch=16, col=pointCol,
	xlab="Reported Mapped (Paper)",
	ylab="STAR Re-Aligned (This Project)")
abline(a=0, b=1,col="gray")
legend("top",legend=c("x=y","linear fit"),
		col=c("gray","black"), lwd=2, ncol=2,
		xpd=T, inset = -0.1)
abline(lm(Unique_Aligned_Count~metadata.table$totalMapped.paper),col="black", lwd=2)
dev.off()

metadata_file = "../Week4/metadata_full.txt"
FastQC_file = "FastQC_stats-extra_3var.txt"

metadata.table = read.table(metadata_file, head=T, sep="\t")
print(dim(metadata.table))
metadata.table = metadata.table[metadata.table$Download == "YES",]
print(dim(metadata.table))

FastQC.table = read.table(FastQC_file, head=T, sep="\t")
R1 = paste(metadata.table$SRR,"_1",sep="")
print(dim(FastQC.table))
FastQC.table = FastQC.table[match(R1,FastQC.table$Sample),]
print(dim(FastQC.table))

#I realize some groups are not represented, but I hope this can help visualize that.
metadata.table$Group = as.character(metadata.table$Group)
metadata.table$Group[metadata.table$Group == "Fetal:Cytosol"]="Fetal"
metadata.table$Group[metadata.table$Group == "Fetal:Nucleus"]="Fetal"
metadata.table$Group[metadata.table$Group == "Adult:Cytosol"]="Adult"
metadata.table$Group[metadata.table$Group == "Adult:Nucleus"]="Adult"
metadata.table$Group[metadata.table$Group == "50+"]="Senior"
group_levels = c("Fetal","Infant","Child","Teen","Adult","Senior")
metadata.table$Group = factor(metadata.table$Group,
								levels = group_levels)

pointCol = rep("black",nrow(metadata.table))
pointCol[metadata.table$Group == "Adult"]="blue"
pointCol[metadata.table$Group == "Fetal"]="red"
pointCol[metadata.table$Group == "Infant"]="orange"

png("FastQC_selected_summarized_values.png", height=400, width=1200)
par(mfcol=c(1,3))

plot(jitter(as.numeric(metadata.table$Group)),FastQC.table$Percent_Average_Q30,
	pch=16, col=pointCol, xaxt='n', xlim=c(0,7), cex=2, cex.main=2,
	xlab="", ylab="Percent Average Q30 (R1)", main="Forward Read Q30 by Group")
mtext(group_levels, side=1, at=c(1:length(group_levels)),
		las=2, line=1)
		
plot(jitter(as.numeric(metadata.table$Group)),FastQC.table$Average_GC,
	pch=16, col=pointCol, xaxt='n', xlim=c(0,7),cex=2, cex.main=2,
	xlab="", ylab="Average GC (R1)", main="Forward Read GC by Group")
mtext(group_levels, side=1, at=c(1:length(group_levels)),
		las=2, line=1)
		
plot(jitter(as.numeric(metadata.table$Group)),FastQC.table$Percent_Unique,
	pch=16, col=pointCol, xaxt='n', xlim=c(0,7),cex=2, cex.main=2,
	xlab="", ylab="Percent Unique / Not Duplicated (R1)", main="Forward Unique Read by Group")
mtext(group_levels, side=1, at=c(1:length(group_levels)),
		las=2, line=1)

dev.off()

#I copied and modified snippets of code from https://github.com/cwarden45/RNAseq_templates/blob/master/TopHat_Workflow/DEG_goseq.R
group.2group = metadata.table$Group[metadata.table$Group != "Infant"]
Q30.2group = FastQC.table$Percent_Average_Q30[metadata.table$Group != "Infant"]
GC.2group = FastQC.table$Average_GC[metadata.table$Group != "Infant"]
unique_R1.2group = FastQC.table$Percent_Unique[metadata.table$Group != "Infant"]

fit = aov(Q30.2group ~ group.2group)
result = summary(fit)
aov.pvalue = result[[1]][['Pr(>F)']][1]
print(aov.pvalue)

fit = aov(GC.2group ~ group.2group)
result = summary(fit)
aov.pvalue = result[[1]][['Pr(>F)']][1]
print(aov.pvalue)

fit = aov(unique_R1.2group ~ group.2group)
result = summary(fit)
aov.pvalue = result[[1]][['Pr(>F)']][1]
print(aov.pvalue)
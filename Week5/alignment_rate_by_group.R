metadata_file = "../Week7/metadata_updated.txt"

metadata.table = read.table(metadata_file, head=T, sep="\t")
print(dim(metadata.table))
metadata.table = metadata.table[metadata.table$Download == "YES",]
print(dim(metadata.table))

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

png("percent_STAR_unique_aligned_reads-by-Group.png")
plot(jitter(as.numeric(metadata.table$Group)),metadata.table$STAR_Unique,
	pch=16, col=pointCol, xaxt='n', xlim=c(0,7), cex=1.5, cex.main=1.2,
	xlab="", ylab="STAR Unique Alignment Rate", main="STAR Unique Alignment Rate by Group")
mtext(group_levels, side=1, at=c(1:length(group_levels)),
		las=2, line=1)
dev.off()

#I copied and modified a snippet of code from https://github.com/cwarden45/RNAseq_templates/blob/master/TopHat_Workflow/DEG_goseq.R
group.2group = metadata.table$Group[metadata.table$Group != "Infant"]
aln_rate.2group = metadata.table$STAR_Unique[metadata.table$Group != "Infant"]
fit = aov(aln_rate.2group ~ group.2group)
result = summary(fit)
aov.pvalue = result[[1]][['Pr(>F)']][1]

print(result)
print(aov.pvalue)

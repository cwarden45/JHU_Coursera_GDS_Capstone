metadata_file = "../Week4/metadata_full.txt"
RSeQC_file = "RSeQC_stats.txt"

RSeQC.table = read.table(RSeQC_file, head=T, sep="\t")

metadata.table = read.table(metadata_file, head=T, sep="\t")
print(dim(metadata.table))
metadata.table = metadata.table[match(RSeQC.table$Sample,metadata.table$SRR),]
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

#correlation
png("RIN_vs_TIN.png")
plot(metadata.table$RIN, RSeQC.table$medTIN,
	xlab="RIN",ylab="Median TIN",
	xlim=c(0,10),ylim=c(0,100),
	pch=16, col=pointCol)
legend("bottom",legend=c("Fetal","Infant","Adult"),
		col=c("red","orange","blue"), pch=16,
		ncol=3, cex=1.3, inset = 0.05)
dev.off()

#separate plots
png("RIN_TIN_by_Group.png", height=400, width=800)
par(mfcol=c(1,2))

plot(jitter(as.numeric(metadata.table$Group)),metadata.table$RIN,
	pch=16, col=pointCol, xaxt='n', xlim=c(0,7), ylim=c(0,10),
	xlab="", ylab="Provided RIN")
mtext(group_levels, side=1, at=c(1:length(group_levels)),
		las=2, line=1)

plot(jitter(as.numeric(metadata.table$Group)),RSeQC.table$medTIN,
	pch=16, col=pointCol, xaxt='n', xlim=c(0,7), ylim=c(0,100),
	xlab="", ylab="RNA-Seq Calculated TIN")
mtext(group_levels, side=1, at=c(1:length(group_levels)),
		las=2, line=1)

dev.off()
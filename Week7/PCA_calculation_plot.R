input.file = "featureCounts_quantified_CPM.txt"
meta.file = "../Week6/alignment_stats_with_quant_total.txt"
rounding_factor = 0.1
pca.text.file = paste("PCA_featureCounts_log2_",rounding_factor,"_quantified_CPM.txt",sep="")
pca.plot.file = paste("PCA_featureCounts_log2_",rounding_factor,"_quantified_CPM.png",sep="")

#input.file = "featureCounts_quantified_CPM.txt"
#meta.file = "../Week6/alignment_stats_with_quant_total.txt"
#rounding_factor = 1
#pca.text.file = paste("PCA_featureCounts_log2_",rounding_factor,"_quantified_CPM.txt",sep="")
#pca.plot.file = paste("PCA_featureCounts_log2_",rounding_factor,"_quantified_CPM.png",sep="")

#input.file = "featureCounts_quantified_CPM.txt"
#meta.file = "../Week6/alignment_stats_with_quant_total.txt"
#rounding_factor = 10
#pca.text.file = paste("PCA_featureCounts_log2_",rounding_factor,"_quantified_CPM.txt",sep="")
#pca.plot.file = paste("PCA_featureCounts_log2_",rounding_factor,"_quantified_CPM.png",sep="")


input.table = read.table(input.file, head=T, sep="\t")
meta.table = read.table(meta.file, head=T, sep="\t")

CPM.mat = input.table[,match(meta.table$Sample,names(input.table))]

CPM.log2.mat = log2(CPM.mat+rounding_factor)

#copy and modify portion of https://github.com/cwarden45/RNAseq_templates/blob/master/TopHat_Workflow/qc.R
pca.values = prcomp(CPM.log2.mat)
pc.values = data.frame(pca.values$rotation)
variance.explained = (pca.values$sdev)^2 / sum(pca.values$sdev^2)
pca.table = data.frame(PC = 1:length(variance.explained), percent.variation = variance.explained, t(pc.values))

write.table(pca.table, pca.text.file, quote=F, row.names=F, sep="\t")
	
plotCol = rep("gray",nrow(meta.table))
plotCol[meta.table$Group == "Adult"]="blue"
plotCol[meta.table$Group == "Fetal"]="red"
plotCol[meta.table$Group == "Infant"]="orange"

png(file=pca.plot.file)
plot(pc.values$PC1, pc.values$PC2, col = plotCol, xlab = paste("PC1 (",round(100* variance.explained[1] , digits = 2),"%)", sep = ""),
		ylab = paste("PC2 (",round(100* variance.explained[2] , digits = 2),"%)", sep = ""), pch=16,
	    main = paste("Plot of PC1 vs PC2 from log2(CPM + ",rounding_factor,") Values",sep=""))
legend("bottomright",legend=c("Adult","Fetal","Infant"),
			col=c("blue","red","orange"),  pch=16)
dev.off()

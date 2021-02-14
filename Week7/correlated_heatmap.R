meta.file = "../Week4/metadata_full.txt"
meta.flagstat = "../Week4/alignment_stats-samtools_flagstat.txt"
meta.GenomicAlignments = "../Week4/GenomicAlignments_aligned_reads_counts.txt"
meta.FastQC = "../Week5/FastQC_stats-extra_3var.txt"
meta.RSeQC = "../Week5/RSeQC_stats.txt"
meta.quant = "../Week6/alignment_stats_with_quant_total.txt"#includes the STAR log unique alignment stats
pca.file1 = "PCA_featureCounts_log2_0.1_quantified_CPM.txt"
pca.file2 = "PCA_featureCounts_log2_1_quantified_CPM.txt"
pca.file3 = "PCA_featureCounts_log2_10_quantified_CPM.txt"
combined_meta = "metadata_updated.txt"
heatmap.file = "correlation_heatmap.png"

### first create larger metadata table ###

meta.in = read.table(meta.file, head=T, sep="\t")
print(dim(meta.in))
meta.in = meta.in[meta.in$Download == "YES",]
print(dim(meta.in))

#provided variables
meta.in$Group = as.character(meta.in$Group)
group_levels = c("Fetal","Infant","Child","Teen","Adult","Senior")
Age_Group = factor(meta.in$Group, levels = group_levels)
print(table(Age_Group))

meta.out = data.frame(meta.in,
					Age_Group, Age_Float=meta.in$Age,
					Mapped_Provided=meta.in$totalMapped.paper)

#samtools flagstat
flagstat.table = read.table(meta.flagstat, head=T, sep="\t")
print(dim(flagstat.table))
flagstat.table = flagstat.table[match(meta.in$SRR,flagstat.table$Sample),]
print(dim(flagstat.table))
meta.out = data.frame(meta.out, Flagstat_Percent=flagstat.table$Aligned_Percent)

#STAR log and quantified percent
quant.table = read.table(meta.quant, head=T, sep="\t")
print(dim(quant.table))
quant.table = quant.table[match(meta.in$SRR,quant.table$Sample),]
print(dim(quant.table))
meta.out = data.frame(meta.out, STAR_Unique=quant.table$Unique_Aligned_Percent,
						Percent_Quant=quant.table$percent_quantified)

#GenomicAlignments stats
GA.table = read.table(meta.GenomicAlignments, head=T, sep="\t")
print(dim(GA.table))
GA.table = GA.table[match(meta.in$SRR,GA.table$Sample),]
print(dim(GA.table))

Percent_GA = 100 * GA.table$aligned.reads /quant.table$Total
meta.out = data.frame(meta.out, Percent_GA)

#FastQC stats
FastQC.table = read.table(meta.FastQC, head=T, sep="\t")
print(dim(FastQC.table))
FastQC.table = FastQC.table[match(paste(meta.in$SRR,"_1",sep=""),FastQC.table$Sample),]
print(dim(FastQC.table))

meta.out = data.frame(meta.out, FastQC.table[,2:4])

#RSeQC stats
RSeQC.table = read.table(meta.RSeQC, head=T, sep="\t")
print(dim(RSeQC.table))
RSeQC.table = RSeQC.table[match(meta.in$SRR,RSeQC.table$Sample),]
print(dim(RSeQC.table))

meta.out = data.frame(meta.out, TIN=RSeQC.table$medTIN)

#PCA stats
PC.table1 = read.table(pca.file1, head=T, sep="\t")
print(dim(PC.table1))
PC.mat1 = PC.table1[,3:ncol(PC.table1)]
rownames(PC.mat1)=paste("PC",PC.table1$PC,sep="")
PC.mat1 = t(PC.mat1)
PC.mat1 = data.frame(Sample=rownames(PC.mat1),PC.mat1)
print(dim(PC.mat1))
PC.mat1 = PC.mat1[match(meta.in$SRR,PC.mat1$Sample),]
print(dim(PC.mat1))
meta.out = data.frame(meta.out, PC1_0.1=PC.mat1$PC1,
						PC2_0.1=PC.mat1$PC2)
						
PC.table2 = read.table(pca.file2, head=T, sep="\t")
print(dim(PC.table2))
PC.mat2 = PC.table2[,3:ncol(PC.table2)]
rownames(PC.mat2)=paste("PC",PC.table2$PC,sep="")
PC.mat2 = t(PC.mat2)
PC.mat2 = data.frame(Sample=rownames(PC.mat2),PC.mat2)
print(dim(PC.mat2))
PC.mat2 = PC.mat2[match(meta.in$SRR,PC.mat2$Sample),]
print(dim(PC.mat2))
meta.out = data.frame(meta.out, PC1_1=PC.mat2$PC1,
						PC2_1=PC.mat2$PC2)

PC.table3 = read.table(pca.file3, head=T, sep="\t")
print(dim(PC.table3))
PC.mat3 = PC.table3[,3:ncol(PC.table3)]
rownames(PC.mat3)=paste("PC",PC.table3$PC,sep="")
PC.mat3 = t(PC.mat3)
PC.mat3 = data.frame(Sample=rownames(PC.mat3),PC.mat3)
print(dim(PC.mat3))
PC.mat3 = PC.mat3[match(meta.in$SRR,PC.mat3$Sample),]
print(dim(PC.mat3))
meta.out = data.frame(meta.out, PC1_10=PC.mat3$PC1,
						PC2_10=PC.mat3$PC2)

write.table(meta.out, combined_meta, quote=F, sep="\t", row.names=F)
### Now, create heatmap of correlation coefficients ###

#leave out age group (only for visualization later on)
plot.vars = c("Age_Float","RIN",
				"Mapped_Provided",
				"Flagstat_Percent","STAR_Unique","Percent_GA",
				"Percent_Average_Q30","Average_GC","Percent_Unique",
				"TIN",
				"Percent_Quant",
				"PC1_0.1","PC2_0.1",
				"PC1_1","PC2_1",
				"PC1_10","PC2_10")
				
plot.mat = meta.out[,match(plot.vars,names(meta.out))]
rownames(plot.mat)=meta.in$SRR

cor.mat = cor(as.matrix(plot.mat))

library(gplots)
#copied and modified from https://github.com/cwarden45/RNAseq_templates/blob/master/TopHat_Workflow/DEG_goseq.R
png(heatmap.file)
		heatmap.2(cor.mat,
			col=colorpanel(33, low="purple", mid="grey", high="red"),
			density.info="none", key=TRUE,
			trace="none", margins = c(10,10))
dev.off()
FDR.cutoff = 0.01
fc.cutoff = 2
expr.cutoff = 5

deg.file = paste("edgeR_DEG_FDR",FDR.cutoff,"_fc",fc.cutoff,"_expr",expr.cutoff,".txt",sep="")
cpm.file = "../Week7/featureCounts_quantified_CPM.txt"
meta.file = "../Week7/metadata_updated.txt"
OUTPUT.FILE = paste("edgeR_DEG_FDR",FDR.cutoff,"_fc",fc.cutoff,"_expr",expr.cutoff,"-heatmap.png",sep="")

source("heatmap.3.R")
library(gplots)
library("RColorBrewer")

cpm.table = read.table(cpm.file, head=T, sep="\t")
meta.table = read.table(meta.file, head=T, sep="\t")

CPM.mat = cpm.table[,match(meta.table$SRR,names(cpm.table))]
CPM.log2.mat = log2(CPM.mat+1)
rownames(CPM.log2.mat)=cpm.table$Gene

deg.table = read.table(deg.file, head=T, sep="\t")
print(dim(deg.table))
deg.table = deg.table[abs(deg.table$Log2.Fold.Change) >= log2(fc.cutoff),]
print(dim(deg.table))
deg.table = deg.table[deg.table$FDR < FDR.cutoff,]
print(dim(deg.table))

deg.genes = as.character(deg.table$Gene)
deg.CPM = CPM.log2.mat[match(deg.genes, rownames(CPM.log2.mat)),]

groupCol = rep("gray",nrow(meta.table))
groupCol[meta.table$Group == "Adult"]="blue"
groupCol[meta.table$Group == "Fetal"]="red"
groupCol[meta.table$Group == "Infant"]="orange"

compCol = groupCol
compCol[meta.table$Group == "Infant"] = "black"

RIN = meta.table$RIN

#portions of code are copied over and modified from https://github.com/cwarden45/RNAseq_templates/blob/master/TopHat_Workflow/DEG_goseq.R
standardize.expr = function(expr)
{
	center.expr = as.numeric(expr) - mean(as.numeric(expr), na.rm=T)
	norm.expr = center.expr / sd(center.expr, na.rm=T)
	return(norm.expr)
}#end def standardize.expr

rinCol = rep("black",times=nrow(meta.table))
continuous.color.breaks = 10

continuous.plot.min = min(RIN, na.rm=T)
continuous.plot.max = max(RIN, na.rm=T)
				
continuous.plot.range = continuous.plot.max - continuous.plot.min
continuous.plot.interval = continuous.plot.range / continuous.color.breaks
				
color.range = colorRampPalette(c("orange","black","cyan"))(n = continuous.color.breaks)
continuous.plot.breaks = continuous.plot.min + continuous.plot.interval*(0:continuous.color.breaks)
for (j in 1:continuous.color.breaks){
	#print(paste(continuous.plot.breaks[j],"to",continuous.plot.breaks[j+1]))
	rinCol[(RIN >= continuous.plot.breaks[j]) &(RIN <= continuous.plot.breaks[j+1])] = color.range[j]
}#end for (j in 1:continuous.color.breaks)

row_annotation = data.frame(label1 = compCol, label2 = groupCol, label3 = rinCol)
row_annotation = as.matrix(t(row_annotation))
rownames(row_annotation) <- c("AvF","Group","RIN")

std.expr = apply(deg.CPM, 1, standardize.expr)
colnames(std.expr) = rep("", length(deg.genes))
rownames(std.expr) = meta.table$SRR

png(OUTPUT.FILE)
heatmap.3(std.expr,
			col=colorpanel(33, low="blue", mid="black", high="red"),
			density.info="none", key=TRUE,
			RowSideColors=row_annotation, trace="none",
			margins = c(10,15),RowSideColorsSize=4, dendrogram="both")
legend("right",legend=c(round(continuous.plot.max,digits=1),rep("",length(color.range)-2),round(continuous.plot.min,digits=1)),
			col=rev(color.range),  pch=15, y.intersp = 0.4, cex=0.8, pt.cex=1.5)
legend("topright", legend=c("Fetal","Infant","Adult"),
				col=c("Red","Orange","Blue"), pch=15, cex=0.7)
dev.off()
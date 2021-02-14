#FDR.cutoff = 0.01
#fc.cutoff = 1
#expr.cutoff = 0

#FDR.cutoff = 0.01
#fc.cutoff = 1
#expr.cutoff = 5

FDR.cutoff = 0.01
fc.cutoff = 2
expr.cutoff = 5


input.file = "../Week6/featureCounts_counts.txt"
meta.file = "../Week6/alignment_stats_with_quant_total.txt"
#don't filter fold-change or FDR, but we do filter the genes being tested
output.file = paste("limma-voom_DEG_FDR",FDR.cutoff,"_fc",fc.cutoff,"_expr",expr.cutoff,".txt",sep="")
output.file2 = paste("limma-voom_DEG_FDR",FDR.cutoff,"_fc",fc.cutoff,"_expr",expr.cutoff,"-voom_plot.png",sep="")


library(edgeR)
library(limma)

input.table = read.table(input.file, head=T, sep="\t")
Gene = input.table$Gene
meta.table = read.table(meta.file, head=T, sep="\t")
print(dim(meta.table))
meta.table = meta.table[meta.table$Group != "Infant",]
print(dim(meta.table))

count.mat = input.table[,match(meta.table$Sample,names(input.table))]
rownames(count.mat)=Gene

print(dim(count.mat))
#calculate CPM
CPM.mat = matrix(ncol=ncol(count.mat), nrow=nrow(count.mat),
				dimnames = list(Genes=Gene, Samples=colnames(count.mat)))

for (i in 1:ncol(count.mat)){
	temp_counts = as.numeric(count.mat[,i])
	quantified_millions = meta.table$quantified_count[i]/1000000
	
	CPM.mat[,i]= temp_counts / quantified_millions
}#end for (i in 1:ncol(count.mat))

#calculate maximum CPM per gene
max_CPM = apply(CPM.mat, 1, max)

#filter genes whose maximum CPM is not above some minimum threshold
count.mat = count.mat[max_CPM >= expr.cutoff,]
print(dim(count.mat))

group = factor(meta.table$Group, levels = c("Fetal","Adult"))
y = DGEList(counts=count.mat, group=group)
design = model.matrix(~group)

png(output.file2)
v = voom(y, design, plot=TRUE)
dev.off()

fit = lmFit(v, design)
fit = eBayes(fit)

result = topTable(fit, coef=ncol(design), number=nrow(count.mat))
result = result[sort(rownames(result)),]

FDR = p.adjust(result$P.Value,"BH")

output.table = data.frame(Gene=rownames(result),
							Log2.Fold.Change = result$logFC,
							p.value = result$P.Value,
							FDR)
write.table(output.table, output.file, quote=F, sep="\t", row.names=F)

print(output.table[output.table$Gene == "SOX11",])
print(dim(output.table))
output.table = output.table[output.table$FDR < FDR.cutoff,]
print(dim(output.table))

print(paste("Adult Up: ", nrow(output.table[output.table$Log2.Fold.Change >= log2(fc.cutoff),]),sep=""))
print(paste("Fetal Up: ", nrow(output.table[output.table$Log2.Fold.Change <= -log2(fc.cutoff),]),sep=""))

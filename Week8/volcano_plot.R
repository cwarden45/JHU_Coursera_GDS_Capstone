FDR.cutoff = 0.01
fc.cutoff = 2
expr.cutoff = 5

library("plotrix")

input.file = paste("edgeR_DEG_FDR",FDR.cutoff,"_fc",fc.cutoff,"_expr",expr.cutoff,".txt",sep="")
output.file = paste("edgeR_DEG_FDR",FDR.cutoff,"_fc",fc.cutoff,"_expr",expr.cutoff,".png",sep="")

input.table = read.table(input.file, head=T, sep="\t")

pointCol = rep("black",nrow(input.table))
pointCol[input.table$FDR < FDR.cutoff]="red"
pointCol[abs(input.table$Log2.Fold.Change) >= log2(fc.cutoff)]="orange"
pointCol[(abs(input.table$Log2.Fold.Change) >= log2(fc.cutoff))&(input.table$FDR < FDR.cutoff)]="darkgreen"

plot.genes = union(c("SOX11"),as.character(input.table$Gene[-log10(input.table$p.value) > 110]))
print(plot.genes)

text.X = input.table$Log2.Fold.Change[match(plot.genes,input.table$Gene)]
text.X[plot.genes == "ST8SIA2"] = text.X[plot.genes == "ST8SIA2"]-1.2
text.Y = -log10(input.table$p.value)[match(plot.genes,input.table$Gene)]+8

png(output.file)
plot(input.table$Log2.Fold.Change,
		-log10(input.table$p.value),
		xlab="log2FoldChange",ylab="-log10(pvalue)",
		pch=16, col=pointCol, ylim=c(0,max(-log10(input.table$p.value))+15),
		main = "Volcano Plot (Positive FC = Adult-Up, Negative FC = Fetal=Up)")
text(x=text.X, y=text.Y, labels=plot.genes, font=2, cex=0.8)
#I found the following function here: https://stackoverflow.com/questions/22265704/drawing-circle-in-r/22266006
draw.circle(x=text.X[plot.genes == "SOX11"],
			y=text.Y[plot.genes == "SOX11"],
			radius=1.2, border="darkgreen", lwd=1.5)
dev.off()
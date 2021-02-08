alignment_stat_file = "alignment_stats_with_quant_total.txt"
GenomicAlignments_count_file = "../Week4/GenomicAlignments_aligned_reads_counts.txt"

quantified.alignment.table = read.table(alignment_stat_file, head=T, sep="\t")
GenomicAlignments.table = read.table(GenomicAlignments_count_file, head=T, sep="\t")

recalculated_aligned_reads = GenomicAlignments.table$aligned.reads[match(quantified.alignment.table$Sample, GenomicAlignments.table$Sample)]

pointCol = rep("black",nrow(quantified.alignment.table))
pointCol[quantified.alignment.table$Group == "Adult"]="blue"
pointCol[quantified.alignment.table$Group == "Fetal"]="red"
pointCol[quantified.alignment.table$Group == "Infant"]="orange"

png("quantified_starting_read_cor.png")
plot(quantified.alignment.table$quantified_count, quantified.alignment.table$Total,
	xlab="Total featureCounts Quantified Reads",
	ylab="Total Starting Reads",
	pch=16, col=pointCol, cex=1.2)
abline(a=0, b=1,col="gray")
legend("top",legend=c("x=y","linear fit"),
		col=c("gray","black"), lwd=2, ncol=2,
		xpd=T, inset = -0.1)
abline(lm(quantified.alignment.table$Total~quantified.alignment.table$quantified_count),col="black", lwd=1)
dev.off()

png("quantified_GenomicAlignments_cor.png")
plot(quantified.alignment.table$quantified_count, recalculated_aligned_reads,
	xlab="Total featureCounts Quantified Reads",
	ylab="Total GenomicAlignments Aligned Reads",
	pch=16, col=pointCol)
abline(a=0, b=1,col="gray")
legend("top",legend=c("x=y","linear fit"),
		col=c("gray","black"), lwd=2, ncol=2,
		xpd=T, inset = -0.1)
abline(lm(recalculated_aligned_reads~quantified.alignment.table$quantified_count),col="black", lwd=1)
dev.off()

percent_GenomicAlignment = 100 * recalculated_aligned_reads /quantified.alignment.table$Total

png("quantified_GenomicAlignments_cor-PERCENT.png")
plot(quantified.alignment.table$percent_quantified, percent_GenomicAlignment,
	xlab="Percent featureCounts Quantified Reads",
	ylab="Percent GenomicAlignments Aligned Reads",
	pch=16, col=pointCol, xlim=c(0,100), ylim=c(0,100))
abline(lm(percent_GenomicAlignment~quantified.alignment.table$percent_quantified),col="black", lwd=1)
dev.off()

png("quantified_STAR_reported_unique_alignment_cor-PERCENT.png")
plot(quantified.alignment.table$percent_quantified, quantified.alignment.table$Unique_Aligned_Percent,
	xlab="Percent featureCounts Quantified Reads",
	ylab="Percent STAR Unique Aligned Reads",
	pch=16, col=pointCol, xlim=c(0,100), ylim=c(0,100))
abline(lm(quantified.alignment.table$Unique_Aligned_Percent~quantified.alignment.table$percent_quantified),col="black", lwd=1)
dev.off()
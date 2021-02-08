alignment_stat_file = "../Week4/alignment_stats.txt"
count_file = "featureCounts_counts.txt"

alignment.stats = read.table(alignment_stat_file, head=T, sep="\t")

sampleIDs  = as.character(alignment.stats$Sample)

for (i in 1:length(sampleIDs)){
	print(sampleIDs[i])
	featureCounts_table = read.table(sampleIDs[i], head=T, sep="\t")
	
	temp_gene = featureCounts_table$Geneid
	temp_length = featureCounts_table$Length
	temp_counts = featureCounts_table[,7]
	
	if(i == 1){
		count.table = data.frame(Gene=temp_gene, Length=temp_length, Count=temp_counts)
	}else{
		temp.table = data.frame(Gene=temp_gene, Length=temp_length, Count=temp_counts)
		
		matched_genes = temp_gene[match(count.table$Gene,temp_gene, nomatch=0)]
		if(length(matched_genes) != length(temp_gene)){
			print("Issue with mismatched sets of gene annotations!")
			stop()
		}#end if(length(matched_genes) != length(temp_gene))
		
		#shouldn't really need to do this
		temp_counts = temp_counts[match(count.table$Gene,temp_gene)]
		
		count.table = data.frame(count.table, Count=temp_counts)
	}#end else
}#end for (i in 1:length(sampleIDs))

colnames(count.table) =  c("Gene","Length",sampleIDs)
count.table = count.table[order(as.character(count.table$Gene)),]
write.table(count.table, count_file, quote=F, sep="\t", row.names=F)

quantified_count = apply(count.table[,3:ncol(count.table)], 2, sum)

percent_quantified = 100 * quantified_count / alignment.stats$Total

extended.stat.file = data.frame(alignment.stats, quantified_count, percent_quantified)
write.table(extended.stat.file, "alignment_stats_with_quant_total.txt", quote=F, sep="\t", row.names=F)
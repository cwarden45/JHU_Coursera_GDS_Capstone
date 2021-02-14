input.file = "../Week6/featureCounts_counts.txt"
output.file = "featureCounts_quantified_CPM.txt"
output.plot = "featureCounts_raw_counts_violin.png"
output.plot2 = "featureCounts_quantified_CPM_violin.png"
meta.file = "../Week6/alignment_stats_with_quant_total.txt"

library(reshape2)
library(ggplot2)

input.table = read.table(input.file, head=T, sep="\t")
Gene = input.table$Gene
meta.table = read.table(meta.file, head=T, sep="\t")

count.mat = input.table[,match(meta.table$Sample,names(input.table))]

### plot raw values ###
plotCol = rep("gray",nrow(meta.table))
plotCol[meta.table$Group == "Adult"]="blue"
plotCol[meta.table$Group == "Fetal"]="red"
plotCol[meta.table$Group == "Infant"]="orange"

count.log2=log2(count.mat+1)
#copied and modified from https://stackoverflow.com/questions/13250872/reshaping-data-to-plot-in-r-using-ggplot2
count.log2.reshape = reshape(data.frame(Gene, count.log2),
					dir="long", idvar="Gene",
					timevar=c("Sample"),times=names(count.log2),
					varying=names(count.log2), v.names=c("Log2Counts"))

Group.reshape = rep(NA,nrow(count.log2.reshape))
for (i in 1:nrow(meta.table)){
	Group.reshape[count.log2.reshape$Sample == meta.table$Sample[i]]=as.character(meta.table$Group[i])
}#end for (i in 1:nrow(meta.table))
count.log2.reshape = data.frame(count.log2.reshape, Group=Group.reshape)

#copied and modified from http://www.sthda.com/english/wiki/ggplot2-violin-plot-quick-start-guide-r-software-and-data-visualization
plot_obj = ggplot(count.log2.reshape, aes(x=Sample, y=Log2Counts, fill=Group))
plot_obj = plot_obj + geom_violin(trim=FALSE)
plot_obj = plot_obj + geom_boxplot(width=0.1, fill="white", col="black")
plot_obj = plot_obj + scale_fill_manual(values=c("blue","red","orange"))
#copied and modified from https://stackoverflow.com/questions/1330989/rotating-and-spacing-axis-labels-in-ggplot2 and 
plot_obj = plot_obj + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
								legend.position = "top",
								panel.background = element_rect(fill = "darkgrey", color = "darkgrey", size=0.1, linetype="solid"),
								panel.grid.major = element_line(size = 0.2, linetype = 'solid',colour = "grey"), 
								panel.grid.minor = element_line(size = 0.1, linetype = 'solid',colour = "grey"))
ggsave(output.plot, width = 14, height = 7)

### calculate and plot CPM values ###
CPM.mat = matrix(ncol=ncol(count.mat), nrow=nrow(count.mat),
				dimnames = list(Genes=Gene, Samples=colnames(count.mat)))

for (i in 1:ncol(count.mat)){
	temp_counts = as.numeric(count.mat[,i])
	quantified_millions = meta.table$quantified_count[i]/1000000
	
	CPM.mat[,i]= temp_counts / quantified_millions
}#end for (i in 1:ncol(count.mat))

gene.info = input.table[,1:2]
output.table = data.frame(gene.info, CPM.mat)
write.table(output.table, output.file, quote=F, sep="\t", row.names=F)

#copy and modify plotting code from above
CPM.log2=log2(CPM.mat+0.1)
CPM.log2.reshape = reshape(data.frame(Gene, CPM.log2),
					dir="long", idvar="Gene",
					timevar=c("Sample"),times=colnames(CPM.log2),
					varying=colnames(CPM.log2), v.names=c("Log2CPM"))

Group.reshape = rep(NA,nrow(CPM.log2.reshape))
for (i in 1:nrow(meta.table)){
	Group.reshape[CPM.log2.reshape$Sample == meta.table$Sample[i]]=as.character(meta.table$Group[i])
}#end for (i in 1:nrow(meta.table))
CPM.log2.reshape = data.frame(CPM.log2.reshape, Group=Group.reshape)

plot_obj = ggplot(CPM.log2.reshape, aes(x=Sample, y=Log2CPM, fill=Group))
plot_obj = plot_obj + geom_violin(trim=FALSE)
plot_obj = plot_obj + geom_boxplot(width=0.1, fill="white", col="black")
plot_obj = plot_obj + scale_fill_manual(values=c("blue","red","orange"))
plot_obj = plot_obj + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
								legend.position = "top",
								panel.background = element_rect(fill = "darkgrey", color = "darkgrey", size=0.1, linetype="solid"),
								panel.grid.major = element_line(size = 0.2, linetype = 'solid',colour = "grey"), 
								panel.grid.minor = element_line(size = 0.1, linetype = 'solid',colour = "grey"))
ggsave(output.plot2, width = 14, height = 7)
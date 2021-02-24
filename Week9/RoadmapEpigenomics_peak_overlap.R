#setwd("E:\\Courses\\Coursera\\JHU_Genomic_Data_Science\\Course8-Capstone\\Week9")

RNAdeg.file = "../Week8/edgeR_DEG_FDR0.01_fc2_expr5.txt"

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(biomaRt)
library(AnnotationHub)

###import and filter DEG stat file (to identify symbols for genes of interest)
RNAdeg.table = read.table(RNAdeg.file, head=T, sep="\t")
print(dim(RNAdeg.table))
RNAdeg.table = RNAdeg.table[RNAdeg.table$FDR < 0.01,]
print(dim(RNAdeg.table))
RNAdeg.table = RNAdeg.table[abs(RNAdeg.table$Log2.Fold.Change) >= 1,]
print(dim(RNAdeg.table))

#copy and modify content to get gene IDs: https://github.com/cwarden45/RNAseq_templates/blob/master/Genome_Ref_Code/get_exon_annotations.R
gene.symbols = keys(org.Hs.eg.db, keytype="SYMBOL")
genecols = c("SYMBOL", "ENTREZID","GENENAME")
genetable = select(org.Hs.eg.db, keys=gene.symbols, columns=genecols, keytype="SYMBOL")

geneIDs = genetable$ENTREZID[match(as.character(RNAdeg.table$Gene), as.character(genetable$SYMBOL))]
RNAdeg.table = RNAdeg.table[!is.na(geneIDs),]
geneIDs = geneIDs[!is.na(geneIDs)]
print(dim(RNAdeg.table))
rownames(RNAdeg.table)=geneIDs
print(RNAdeg.table[RNAdeg.table$Gene == "SOX11",])

RNAseq.UP = rownames(RNAdeg.table)[RNAdeg.table$Log2.Fold.Change > 0]
RNAseq.DOWN = rownames(RNAdeg.table)[RNAdeg.table$Log2.Fold.Change < 0]

###find RNA-Seq promoters
UCSC_promoters = promoters(TxDb.Hsapiens.UCSC.hg19.knownGene,
							columns=c("gene_id"))
print(length(UCSC_promoters))
							
UCSC_promoters.UP = UCSC_promoters[as.character(UCSC_promoters$gene_id) %in% RNAseq.UP,]
print(length(UCSC_promoters.UP))

UCSC_promoters.DOWN = UCSC_promoters[as.character(UCSC_promoters$gene_id) %in% RNAseq.DOWN,]
print(length(UCSC_promoters.DOWN))
print("6664" %in% as.character(UCSC_promoters$gene_id))

###compare H3K4me3 peaks (brain)
ah = AnnotationHub()
ah = subset(ah, species == "Homo sapiens")

ah_results = query(ah, "EpigenomeRoadMap")
ah_results = query(ah_results, "H3K4me3")
ah_resultsFetal = query(ah_results, "E081")
ah_resultsAdult = query(ah_results, "E073")

fetal.narrow_peak = ah_resultsFetal[2]
fetal.narrow_peak.gr = fetal.narrow_peak[[1]]
print(length(fetal.narrow_peak.gr))
fetal.narrow_peak.gr = fetal.narrow_peak.gr[fetal.narrow_peak.gr$qValue > 5,]
print(length(fetal.narrow_peak.gr))

adult.narrow_peak = ah_resultsAdult[2]
adult.narrow_peak.gr = adult.narrow_peak[[1]]
print(length(adult.narrow_peak.gr))
adult.narrow_peak.gr = adult.narrow_peak.gr[adult.narrow_peak.gr$qValue > 5,]
print(length(adult.narrow_peak.gr))

###compare regions
RNAadult.H3K4me3adult = subsetByOverlaps(UCSC_promoters.UP, adult.narrow_peak.gr)
print(length(RNAadult.H3K4me3adult))
RNAadult.H3K4me3fetal = subsetByOverlaps(UCSC_promoters.UP, fetal.narrow_peak.gr)
print(length(RNAadult.H3K4me3fetal))

RNAfetal.H3K4me3adult = subsetByOverlaps(UCSC_promoters.DOWN, adult.narrow_peak.gr)
print(length(RNAfetal.H3K4me3adult))
print("6664" %in% as.character(RNAfetal.H3K4me3adult$gene_id))#this is NOT what we want, at least if we expect a difference
RNAfetal.H3K4me3fetal = subsetByOverlaps(UCSC_promoters.DOWN, fetal.narrow_peak.gr)
print(length(RNAfetal.H3K4me3fetal))
print("6664" %in% as.character(RNAfetal.H3K4me3adult$gene_id))#this IS what we want

fisher.mat = matrix(ncol=2,nrow=2)
colnames(fisher.mat) = c("RNA.adult","RNA.fetal")
rownames(fisher.mat) = c("H3K4me3.adult","H3K4me3.fetal")
fisher.mat[1,1]=length(RNAadult.H3K4me3adult)
fisher.mat[1,2]=length(RNAfetal.H3K4me3adult)
fisher.mat[2,1]=length(RNAadult.H3K4me3fetal)
fisher.mat[2,2]=length(RNAfetal.H3K4me3fetal)
print(fisher.mat)

print(fisher.test(fisher.mat))

###compare H3K4me3 peaks (liver)
ah = AnnotationHub()
ah = subset(ah, species == "Homo sapiens")

ah_results = query(ah, "EpigenomeRoadMap")
ah_results = query(ah_results, "H3K4me3")
ah_resultsLiver = query(ah_results, "E066")

liver.narrow_peak = ah_resultsLiver[2]
liver.narrow_peak.gr = liver.narrow_peak[[1]]
print(length(liver.narrow_peak.gr))
liver.narrow_peak.gr = liver.narrow_peak.gr[liver.narrow_peak.gr$qValue > 5,]
print(length(liver.narrow_peak.gr))

print("6664" %in% as.character(RNAadult.H3K4me3adult$gene_id))#not supposed to be here to begin with
DualAdult.H3K4me3liver = subsetByOverlaps(RNAadult.H3K4me3adult, liver.narrow_peak.gr)
print(length(DualAdult.H3K4me3liver))
print(100*length(DualAdult.H3K4me3liver)/length(RNAadult.H3K4me3adult))
print("6664" %in% as.character(DualAdult.H3K4me3liver$gene_id))

print("6664" %in% as.character(RNAfetal.H3K4me3fetal$gene_id))
DualFetal.H3K4me3liver = subsetByOverlaps(RNAfetal.H3K4me3fetal, liver.narrow_peak.gr)
print(length(DualFetal.H3K4me3liver))
print(100*length(DualFetal.H3K4me3liver)/length(RNAfetal.H3K4me3fetal))
print("6664" %in% as.character(DualAdult.H3K4me3liver$gene_id))

starting.adult.geneID = unique(as.character(RNAadult.H3K4me3adult$gene_id))
starting.fetal.geneID = unique(as.character(RNAfetal.H3K4me3fetal$gene_id))

liver.adult.geneID = unique(as.character(DualAdult.H3K4me3liver$gene_id))
liver.fetal.geneID = unique(as.character(DualFetal.H3K4me3liver$gene_id))

nonoverlapping.adult.geneID = starting.adult.geneID[-match(liver.adult.geneID, starting.adult.geneID)]
nonoverlapping.adult.symbol = as.character(RNAdeg.table$Gene[match(nonoverlapping.adult.geneID,rownames(RNAdeg.table))])

nonoverlapping.fetal.geneID = starting.fetal.geneID[-match(liver.fetal.geneID, starting.fetal.geneID)]
nonoverlapping.fetal.symbol = as.character(RNAdeg.table$Gene[match(nonoverlapping.fetal.geneID,rownames(RNAdeg.table))])

output.table = data.frame(Gene=c(nonoverlapping.adult.symbol,nonoverlapping.fetal.symbol),
						Status=c(rep("Adult",length(nonoverlapping.adult.symbol)),rep("Fetal",length(nonoverlapping.fetal.symbol))))
write.table(output.table,"H3K4me3_candidates_exclude_liver.txt", quote=F, sep="\t", row.names=F)
## Steps for Analysis:

**1)** Carry out tasks directly in assignment using `RoadmapEpigenomics_peak_overlap.R`

**Unfiltered Overlap (*Brain*)**:

<table>
  <tbody>
    <tr>
      <th align="center"></th>
      <th align="center"><b>Adult RNA-Seq Promoter</b></th>
	  <th align="center"><b>Fetal RNA-Seq Promoter</b></th>
    </tr>
     <tr>
      <th align="center"><b>Adult H3K4me3</b></th>
      <th align="center">9,200 promoter overlaps</th>
	  <th align="center">8,365 promoter overlaps</th>
    </tr>
     <tr>
      <th align="center"><b>Fetal H3K4me3</b></th>
      <th align="center">4,616 promoter overlaps</th>
	  <th align="center">6,925 promoter overlaps</th>
    </tr>
</tbody>
</table>

**Overlap with qValue > 5 (*Brain*)**:

<table>
  <tbody>
    <tr>
      <th align="center"></th>
      <th align="center"><b>Adult RNA-Seq Promoter</b></th>
	  <th align="center"><b>Fetal RNA-Seq Promoter</b></th>
    </tr>
     <tr>
      <th align="center"><b>Adult H3K4me3</b></th>
      <th align="center">8,607 promoter overlaps</th>
	  <th align="center">8,109 promoter overlaps</th>
    </tr>
     <tr>
      <th align="center"><b>Fetal H3K4me3</b></th>
      <th align="center">1,864 promoter overlaps</th>
	  <th align="center">4,194 promoter overlaps</th>
    </tr>
</tbody>
</table>

***Non-specific* Overlap (*Liver*, with qValue > 5)**:

<table>
  <tbody>
    <tr>
      <th align="center"></th>
      <th align="center"><b>Adult<br>(RNA-Seq + H3K4me3)</b></th>
	  <th align="center"><b>Fetal<br>(RNA-Seq + H3K4me3)</b></th>
    </tr>
     <tr>
      <th align="center"><b>Adult H3K4me3</b></th>
      <th align="center">8,607 promoter overlaps</th>
	  <th align="center">4,194 promoter overlaps</th>
    </tr>
     <tr>
      <th align="center"><b>Fetal H3K4me3</b></th>
	     <th align="center">6,637 promoter overlaps<br>(77.1% <i>not</i> specific)</th>
	<th align="center"> 3,991 promoter overlaps<br>(95.2% <i>not</i> specific)</th>
    </tr>
</tbody>
</table>

There is a lot genes that do not appear to have liver-specific expression, but I can check the ***filtered*** genes that have overlapping brain peaks but *not* liver peaks (in the next step).

*R*: version 4.0.3 (different verison from other weeks, due to AnnotationHub dependencies)

*biomaRt package*: version 2.46.0

*TxDb.Hsapiens.UCSC.hg19.knownGene package*: version 3.2.2

*org.Hs.eg.db package*: version 3.12.0

*AnnotationHub package*: version 2.22.0

**2)** Carry out enrichment of filtered gene list using [Enrichr](https://maayanlab.cloud/Enrichr/) (on 2/23/2021)

The full results were downloaded and uploaded to this page.  I have summarized the top categories (by p-value) below:

<table>
  <tbody>
    <tr>
      <th align="center"></th>
      <th align="center"><b>Adult Enrichment</b></th>
	<th align="center"><b>Fetal Enrichment</b></th>
    </tr>
     <tr>
      <th align="center"><b>BioPlanet 2019</b></th>
      <th align="center">neuronal system</th>
	<th align="center">axon guidance</th>
    </tr>
     <tr>
      <th align="center"><b>Elseivier Pathway Collection</b></th>
      <th align="center">Proteins Involved in Epilepsy</th>
	<th align="center">Glioblastoma, Proneural Subtype</th>
    </tr>
     <tr>
      <th align="center"><b>Reactome 2016</b></th>
      <th align="center">Neuronal System</th>
	<th align="center">Neuronal System</th>
    </tr>
     <tr>
      <th align="center"><b>GO Biological Process 2018</b></th>
      <th align="center">chemical synaptic transmission</th>
	<th align="center">central nervous system development</th>
    </tr>
     <tr>
      <th align="center"><b>Human Gene Atlas</b></th>
      <th align="center">Amygdala</th>
	<th align="center">Fetalbrain</th>
    </tr>
</tbody>
</table>

So, with somewhat limited expertise in neuroscience, I think that the filtered gene list may be useful for identifying genes with relevant function in the dorsolateral prefrontal cortex of the brain.

## Analysis / Notes Beyond Report:

Additional public data datasets:

[GSE160810](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE160810) - fetal, pediatric, and adult scRNA-Seq of brain tissue

[GSE73721](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE73721) - among other samples, fetal and adult brain tissue (specifically, astrocytes for fetal tissue)

Perhaps as expected by the assignment, you can find fetal tissues processed with various techologies under [University of Washington Human Reference Epigenome Mapping Project / GSE18927](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18927), [Broad Institute Human Reference Epigenome Mapping Project / GSE17312](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE17312), and [	UCSD Human Reference Epigenome Mapping Project / GSE16256](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE16256).

[GSE121723](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE121723) - glioblastoma (brain tumor) samples including purtabation of SOX10 (although not exactly SOX11)

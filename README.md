# Brassicaceae GRN
- Comparative genomics (especially RNA-seq and DAP-seq) analyses of ABA responses among Brassicaceae species (crucifers).  
- Currently, housing scripts and files associated with the first preprint of this project ([Sun et al., 2020](https://doi.org/10.1101/2020.11.18.349449))
- For Python scripts, type the script followed by '-h' for more details on pre-requisite, input, parameters/options, and output.  
- For shell scripts, see annotations inside for more details.
- Go to:  
[Correlation matrices for the Phylogenetically informed Profiling (PiP) analysis](https://github.com/dinnenylab/BrassicaceaeGRN#correlation-matrices-for-the-phylogenetically-informed-profiling-pip-analysis)  
[Merging and annotating DAP-seq peaks](https://github.com/dinnenylab/BrassicaceaeGRN#merging-and-annotating-dap-seq-peaks)
---
## Correlation matrices for the Phylogenetically informed Profiling (PiP) analysis
### 1. Preparing input files
#### 1.1 Gene expression of ortholog pairs input
- RNA-seq reads were processed with `RNA_Seq2020.sh` and [DESeq2](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html) to estimate log<sub>2</sub> fold-change (LFC) in response to ABA treatment.  
- Ortholog pairs were determined between each pair of species, as reciprocal best homologs using [CLfinder](https://github.com/ohdongha/OrthNet#running-clfinder), and organized as follows (i.e. tab-delimited with an ortholog pair ID followed by gene IDs for each species pairs; species other than the pair were filled with placeholders, e.g. "na"):
```
P000001	AT1G01010.1	Si_s2233_00040.a	na	na
P000002	AT1G01020.1	Si_s2233_00060.m	na	na
P000003	AT1G01030.1	Si_s2233_00070.a	na	na
... (a lot of rows) ...
P114315	na	na	Sp5g19900.1	Thhalv10018708m.v1.0
P114316	na	na	Sp5g19260.1	Thhalv10018768m.v1.0
P114317	na	na	Sp5g19300.1	Thhalv10019022m.v1.0
``` 
- DESeq2 LFC columns for each sample, after filtering, were merged to the ortholog pair columns above using `merge_by_NthCol.py` (see final input files in the [PiP_example_ABA_response_4crucifers](https://github.com/dinnenylab/BrassicaceaeGRN/tree/master/PiP_example_ABA_response_4crucifers) folder) 

#### 1.2 GO annotation input
- GO annotation for _Arabidopsis thaliana_ (obtained from [GO consortium](http://geneontology.org/) on July 1st, 2020) was transferred to their homologs in non-model species, if the homolog pair showed protein alignments (e<10<sup>-5</sup>) in total covering over 70% of the length of both query and subject genes.
- GO annotation for all species were organized as follows (i.e. tab-delimited with GO term, followed by gene IDs separated by "|"; GO terms can appear multiple times for different species) 
```
response_to_abscisic_acid	AT1G27730.1|AT2G21230.3|AT3G46750.1|...
response_to_abscisic_acid	SI_S1667_00090.M|SI_S2478_01200.M|SI_S2478_01890.M|...
response_to_abscisic_acid	SP7G00210.1|SP6G19050.1|SP4G25450.1|...
response_to_abscisic_acid	THHALV10019152M.V1.0|THHALV10016322M.V1.0|THHALV10028988M.V1.0|...
...
```
- To achieve this, Arabidopsis-homolog gene ID pairs were converted to a custom GO annotation file using `transfer_GO_annotation_4BinGO.py` following instructions in the [BiNGO manual](https://www.psb.ugent.be/cbd/papers/BiNGO/Customize.html). 
- Using the custom GO annotation, a GO annotation input for a species can be created by running BiNGO with all gene IDs as input and options to print all GO terms without a statistical test. The first and last columns of such BiNGO output were concatanated to create the final GO annotation input as in the [PiP_example_ABA_response_4crucifers](https://github.com/dinnenylab/BrassicaceaeGRN/tree/master/PiP_example_ABA_response_4crucifers) folder.

*Note: GO IDs can be used instead of GO terms, e.g. "GO:0009737" or "9737" instead of "response to abscisic acid."  
*Note: not only GO terms, any gene sets (e.g. ortholog pairs with a certain _cis_-regulatory motif in their promoters) can be added to the GO annotation input.

### 2. Pairwise comparisons of gene expression
- Once the two input files are ready, run `PiP_correlation_matrix_pairwise_orthologs.py` to obtain the matrix of correlations of gene expression LFC between ortholog pairs of all species pairs, for all GO terms and samples.
- For example, the following prints (to STDOUT) both Pearson and Spearman's correlation coefficients and p-values for root 3hr (R03) ABA-responsive LFCs of ortholog pairs, for all six species pairs (among At, Si, Sp, and Es) and all non-redundant GO terms of size 20~3000:
```
PiP_correlation_matrix_pairwise_orthologs.py -N 6 -1 -n \
      -s AtSiSpEs_GOid-GeneIDs_non_redundant_size20-3000.txt \
      AtSiSpEs_rcBHpairs_R03.lfc_filtered.minOneDEG+.txt
```
- The script can also draw "pairs plots" (tested with python 3.8 and seaborn 0.10), comparing LFC values of orthologs pairs of all species pairs, for each GO term and sample.  See the script help ('-h') for details. 

### 3. Statistics and ranking
- Matrices resulted from the step 2 contain columns of p-values estimated for correlation coefficients. These p-values were corrected for multiple testing (g-test), using `stat_multiple_test_correction.py`. See the script help ('-h') for details. 
- Variance of Spearman's correlations across all six species pairs (_V_) were plotted against the median number of ortholog pairs annotated with the GO term (_m_), to detect the presence of outliers (i.e. GO terms with large modifications in gene regulation). 
- In the _V_ * _m_ plot, smaller GO terms tended to have larger _V_, likely due to noises. Empirically, we ranked GO-samples based on higher _V_ * log<sub>2</sub> _m_ values, to identify larger (in effect sizes) modifications among species pairs and, likely, lineages.
---
## Merging and annotating DAP-seq peaks
### 1. Merging peak positions co-occuring among replicates and transcription factors 
- DAP-seq reads were processed with `DAP_Seq2020.sh` to call significant peaks using [GEM](https://groups.csail.mit.edu/cgs/gem/), and the resulting [.narrowPeak](https://genome.ucsc.edu/FAQ/FAQformat.html#format12) files were used for downstream analysis.   
- We used an additional filter based on DAP-seq replicates to detect high-confidence binding events. After concatenating all DAP-seq replicated experiments, we selected GEM-called binding events whose center positions appearing within six nucleotide positions in at least two replicates using [bedtools merge](https://bedtools.readthedocs.io/en/latest/content/tools/merge.html#). Detailed shell scripts for this process are in `mergebed_[At|Si|Sp|Es]_YS.txt`.   
- Overlaps among replicated-merged high-confidence binding events were identified by again `bedtools merge` and the center coordinates of these ABFs-binding peak positions were used for peak annotation and counting overlaps among ABFs, as depicted in the following cartoon: 
<img src="https://user-images.githubusercontent.com/748486/111260241-77969500-85ee-11eb-95e2-0d48e74069dc.png" width="500">

### 2. Annotating DAP-seq peaks 
#### 2.1 Marking ABFs-binding DAP-seq peak positions coinciding with an ACGT or an ABRE 
- First, we marked positions of all ACGT (i.e. ABRE core) in a genome as follows (example for At shown):
```
echo -e 'ACGT\tACGT' > ACGT.txt
find_motifs_in_sequences.py ACGT.txt At_genome.fa ACGT_in_At.txt
```
*Note: the genome sequence (e.g. At_genome.fa) must be formatted as a "one-liner" fasta, i.e. for each sequence, after the sequence header (the line starts with ">"), all nuceoltide sequences should be in the next line without a space.  The easiest way to achieve this is to run ```fasta_formatter``` included in the [FASTX-Toolkit](http://hannonlab.cshl.edu/fastx_toolkit/commandline.html#fasta_formatter_usage). The compiled ```fasta_formatter``` binary from the FASTX-Toolkit is included here (AGP Lincense). If ```fasta_formatter``` does not work, you could also try this [awk trick](https://www.biostars.org/p/9262/#9264).
- Second, we marked the positions of all ABREs in a genome, with redundant positions called by multiple ABRE definitions consolidated, as follows:
```
find_motifs_in_sequences.py ABRE.txt At_genome.fa ABRE_in_At.txt
sort -k1,1 -k2,2n ABRE_in_At.txt > ABRE_in_At.sorted.txt
genomic_regions_collapse_overlaps.py ABRE_in_At.sorted.txt 1 temp.At
grep -vP "\t0$" temp.At > ABRE_in_At.sorted.collapsed.txt
```
- Third, we marked all ACGTs that constitute the core of an ABRE, as follows:
```
genomic_regions_mark_overlaps.py ABRE_in_At.sorted.collapsed.txt 1 \
                  ACGT_in_At.txt ACGT_in_At.ABRE_marked.txt
awk '{print "ACGT_"$1":"$2"-"$3"="$5"\t"$1"\t"$2"\t"$3}' ACGT_in_At.ABRE_marked.txt \
                  | sed "s/=_na_//g" > ACGT_in_At.ABRE_marked.uID.txt
```
- Finally, we marked ABFs-binding peak positions that coincide with an ACGT or an ABRE:
```
genomic_regions_mark_overlaps.py -r ACGT_in_At.ABRE_marked.uID.txt 1 \
                  Peaks_in_At.bed Peaks_in_At.ACGT-ABRE_marked.bed 
```

#### 2.2 Counting ABFs-binding DAP-seq peak positions adjacent to a gene model
- We identified genomic regions adjacent to gene models, as 5' distal (5pD), 5' proximal (5pP), exons, introns, 3' proximal (3pP), and 3' distal (see `genomic_region_annotate.py -h` for details). Exons and introns can be combined to gene body (gCDS):
```
count_chrom_sizes.py At.genome.fa At.genome.chr.list
genomic_regions_annotate.py At.gtf At.genome.chr.list At.genomicRegions.list 
```
- To ABFs-binding peak positions, we added they occured in a genomic region adjacent to a gene model:  
```
genomic_regions_mark_overlaps.py At.genomicRegions.list 1 \
                  Peaks_in_At.ACGT-ABRE_marked.bed Peaks_in_At.annotated.bed
```  
*Note: a peak position may appear in multiple genomic regions of closely located neighboring genes. In this case, all overlaps associated with all neighboring genes are reported separately. The resulting annotation was used to identify all ABFs-binding events occuring adjacent to each gene. 

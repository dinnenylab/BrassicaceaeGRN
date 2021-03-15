# Brassicaceae GRN
- Comparative genomics (especially RNA-seq and DAP-seq) analyses of ABA responses among Brassicaceae species (crucifers).  
- Currently, housing scripts and files associated with the first preprint of this project ([Sun et al., 2020](https://doi.org/10.1101/2020.11.18.349449))

## Custom scripts
For Python scripts, type the script followed by '-h' for more details on pre-requisite, input, parameters/options, and output.
For shell scripts, see annotations inside for more details.

### Correlation matrices for the Phylogenetically informed Profiling (PiP) analysis
#### 1. Preparing input files
1.1 Gene expression of ortholog pairs input
- RNA-seq reads were processed with `RNA_Seq2020.sh` and [DESeq2](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html) to estimate log<sub>2</sub> fold-change (LFC) in response to ABA treatment.  
- Ortholog pairs were determined between each pair of species, as reciprocal best homologs, using [CLfinder](https://github.com/ohdongha/OrthNet#running-clfinder), organized as follows (i.e. tab-delimited with Ortholog pair ID followed by gene IDs for each species pairs; species other than the pair filled with placeholders, e.g. "na"):
```
P000001	AT1G01010.1	Si_s2233_00040.a	na	na
P000002	AT1G01020.1	Si_s2233_00060.m	na	na
P000003	AT1G01030.1	Si_s2233_00070.a	na	na
... (a lot of rows) ...
P114315	na	na	Sp5g19900.1	Thhalv10018708m.v1.0
P114316	na	na	Sp5g19260.1	Thhalv10018768m.v1.0
P114317	na	na	Sp5g19300.1	Thhalv10019022m.v1.0
``` 
- DESeq2 LFC columns, after filtering, for each sample were merged to the ortholog pair columns using `merge_by_NthCol.py` (see resulted final input files in the `PiP_example_ABA_response_4crucifers` folder) 

1.2 GO annotation input
- GO annotation for _Arabidopsis thaliana_ (obtained from [GO consortium](http://geneontology.org/) on July 1st, 2020) was transferred to their homologs in non-model species, if the homolog pair showed protein alignments (e<10<sup>-5</sup>) in total covering over 70% of the length of both query and subject genes.
- GO annotation for all species were organized as follows (i.e. tab-delimited with GO term, followed by gene IDs separated by "|"; GO terms can appear multiple times for different species) 
```
response_to_abscisic_acid	AT1G27730.1|AT2G21230.3|AT3G46750.1|...
response_to_abscisic_acid	SI_S1667_00090.M|SI_S2478_01200.M|SI_S2478_01890.M|...
response_to_abscisic_acid	SP7G00210.1|SP6G19050.1|SP4G25450.1|...
response_to_abscisic_acid	THHALV10019152M.V1.0|THHALV10016322M.V1.0|THHALV10028988M.V1.0|...
...
```
- To achieve this, Arabidopsis-homolog gene ID pairs were converted to a custom GO annotation file using `transfer_GO_annotation_4BinGO.py` following instruction in the [BiNGO manual](https://www.psb.ugent.be/cbd/papers/BiNGO/Customize.html). 
- Using the custom GO annotation, a GO annotation input for a species can be created by running BiNGO with all gene IDs as input and options to print all GO terms without a statistical test. The first and last columns of such BiNGO output were concatanated to create the final GO annotation input as in the `PiP_example_ABA_response_4crucifers` folder. 

#### 2. Pairwise comparison of gene expression
- Once the two input files were ready, run `PiP_correlation_matrix_pairwise_orthologs.py` to obtain the matrix of correlations of gene expression LFC between ortholog pairs of all species pairs, for all GO terms and samples.
- For example, the following will print both Pearson and Spearman's correlation coefficients and p-values for all six pairs among At, Si, Sp, and Es, for all non-redundant GO terms of size 20~3000, to stdout:
```
PiP_correlation_matrix_pairwise_orthologs.py -N 6 -1n -s AtSiSpEs_GOid-GeneIDs_non_redundant_size20-3000.txt AtSiSpEs_rcBHpairs_R03.lfc_filtered.minOneDEG+.txt
```
- The script can also draw "pairs plots" (tested with python 3.8 and seaborn 0.10), plotting LFC values for orthologs pairs of all species pairs.  See the script help ('-h') for details.  

#### 3. Statistical test
- Matrices resulted from the step 2 contain columns of p-values estimated for correlation coefficients. These p-values were corrected for multiple testing (g-test), using `stat_multiple_test_correction.py`. See the script help ('-h') for details. 
- Variance of correlations across all six species pairs (V) were plotted against the median number of genes annotated with the GO term (m), to detect the presence of outliers (i.e. GO terms with large modifications in gene regulation)

### Merging and annotating DAP-seq peaks
#### Merging peak positions co-occuring among replicates and transcription factors 
- DAP-seq reads were processed with `DAP_Seq2020.sh` to call peaks using GEM, and the resulting `.narrowPeak` files are used for merging peaks in different replicates

#### Annotating DAP-seq peaks 
- See the script help ('-h') for details. 
```
find_motifs_in_promoters.py
genomic_regions_annotate.py
genomic_regions_collapse_overlaps.py
genomic_regions_extract_sequences.py
genomic_regions_mark_regions_included_in_others.py
```

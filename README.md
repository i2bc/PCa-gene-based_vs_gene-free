Scripts for reproducing results of "Gene-free transcriptomics identifies new biomarkers of prostate cancer risk"

1. script_infer_signature_gene_level.R: Infers a signature using conventional gene expression counting. 

2. script_infer_signature_kmer_level.R: Infers a signature from a set of RNA contigs obtained by running DE-kupl without Differential expression analysis module.

Both 1. &2. The script starts from a count table in the Discovery Set (TCGA-PRAD cohor), then normalisation is performed as count perbillion for k-mers or count per million for genes. Features are ranked according to their F1-score computed by cross validation using a Naive Bayes classifier (NB) and the top 500 features are retained. Among the top 500, features are selected using lasso logistic regression combined with stability selection. The signature is then used to build a predictor in the Validation Set (ICGC-PRAD cohor).

3. box_plot_signature.R: This script draws the box plot of signature in Discovery set and Validation set.

4. useful_functions.R: collection of functions that are used for infer signature.

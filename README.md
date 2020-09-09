# PCa-gene-based vs gene-free Comparison

Scripts for reproducing results of "A Comparative Analysis of Reference-Free and Conventional Transcriptome Signatures for Prostate Cancer Prognosis".

1. infer_gene-based_risk_signature.R: Infers a risk signature using conventional gene expression counting. 

2. infer_gene-free_risk_signature.R: Infers a risk signature from a set of RNA contigs obtained by, running DE-kupl without Differential expression analysis module.

3. infer_gene-based_relaspe_signature.R: The same as script 1 but in relapse condition.

4. infer_gene-free_relaspe_signature.R: The same as script 2 but in relapse condition.


All four scripts starts from a count table in the Discovery Set (TCGA-PRAD cohor), then normalisation is performed as count perbillion for k-mers or count per million for genes. Features are ranked according to their F1-score computed by cross validation using a Naive Bayes classifier (NB) and the top 500 features are retained. Among the top 500, features are selected using lasso logistic regression combined with stability selection. The signature is then used to build a predictor in the two independent RNA-seq datasets (ICGC-PRAD and Stelloo cohors) for validation.

5. box_plot_risk_signature.R: This script draws the box plot of risk signature in Discovery set and Validation set.

6. box_plot_relapse_signature.R: This script draws the box plot of relapse signature in Discovery set and Validation set.

7. useful_functions.R: collection of functions that are used for infer signature in both tumor risk and relapse models.


## KaMRaT Usage

KaMRaT means "k-mer Matrix Reduction Toolbox", or "k-mer Matrix, Really Tremendous !".

The toolbox contains for now modules below:

- kamratMerge which takes a k-mer count table as input, and merge k-mers into longer contigs according to their overlap
- kamratNorm which takes a k-mer count table as input, and normalize k-mer counts
- kamratReduce which takes a k-mer count table as input, and evaluates the performance of each k-mer according to different metrics

### Prerequest

The kamratReduce module is dependent on MLPack library. It could be installed from [conda cloud](https://anaconda.org/conda-forge/mlpack).

After installation, please add the following line into your ```.bashrc``` file in the ```home/``` directory:

```bash
export LD_LIBRARY_PATH="/path_to_conda_env/mlpack/lib:$LD_LIBRARY_PATH"
```

For compiling the source, you can use ```compile.bash``` in root directory of KaMRaT:

If you installed MLPack library by conda, please do

```bash
bash compile.bash /path_to_conda_environment_for_MLPack
```
Otherwise, if you installed MLPack library without conda, please do

```bash
bash compile.bash
```

And the executable files are in ```bin/``` directory.

### File of Sample Info

The sample-info file is indicated by the option ```-d```. This file aims to indicate which columns in the k-mer count matrix should be considered as sample columns. Please do not add any header line in the file, since the column contents are already defined as below.

If provided, this file may contains one or two columns:

- if the file contains only one column, it indicates sample names, and all samples are considered as in same condition
- if the file contains two columns, the first column correspond to sample names, and the second conrresponds to conditions

If not provided, all columns apart from the first one in the k-mer count matrix are considered as samples.

### kamratMerge

Merge Usage

```text
kamratMerge [-h] [-k k_length] [-m min_overlap] [-nxj] [-d sample_info] [-i interv_method] [-q quant_mode] [-t tmp_dir] kmer_count_path
```

Merge Parameters

```text
-h         Print the helper
-n         If the k-mers are generated from non-stranded RNA-seq data
-k INT     k-mer length (max_value: 32) [31]
-m INT     Min assembly overlap (max_value: k) [15]
-d STRING  Sample-info path, either list or table with sample names as the first column
           if absent, all columns except the first one in k-mer count table are taken as samples
-i STRING  Intervention method (none, mac:0.25) [none]
           the threshold can be precised after a ':' symbol
-q STRING  Quantification mode (rep, mean) [rep]
           the column name for selecting representative k-mer can be precised after ':' symbol, e.g. rep:pvalue
           if no column name precised, the firstly input k-mer of a contig will be taken as representative k-mer
-j         Adjacent k-mer comparison (valid only with intervention) [false]
           if absent, the counts of representative k-mers or mean counts are taken according to quantification mode
-x         Query on disk [false]
-t STRING  Temporary directory [./]
```

Intervention Method

```text
none: without any intervention, only consider overlap for extension
mac: mean absolute contrast, mac(k1, k2) = mean(abs(k1-k2)./(k1+k2)), where k1, k2 are count vectors of two k-mers and './' means division for each components.
```

### kamratNorm

Norm Usage

```text
kamratNorm -b CHAR [-d STRING] [-T STRING] STRING
```

Norm Parameters

```text
-h         Print the helper
-b CHAR    BASE for normalization: count per BASE (MANDATORY, could be B, M, or K)
-d STRING  Sample-info path, if absent, all columns will be printed
-T STRING  Transformation before evaluation (e.g. log)"
```

### kamratReduce

Reduce Usage

```text
kamratReduce [-d samp_info_path -m eval_method:fold_num -s sort_mode -N top_num -T transf_mode -C count_offset] kmer_count_path
```

Reduce Parameters

```text
-h           Print the helper
-d STRING    Path to sample-condition or sample file, without header line
             if absent, all except the first column in k-mer count table will be regarded as samples
-m STRING    Evaluation method name and their sub-evaluation-command, seperated by ':'
-s STRING    Sorting mode, default value depends on evaluation method (c.f [SORT MODE])
-N INT       Number of top features to print
-T STRING    Transformation before evaluation (e.g. log)
-C DOUBLE    Count offset to be added before transformation and evaluation
-r INT       Accepted minimum recurrence number among samples
-a INT       Accepted minimum count abundance for counting sample recurrence
```

Reduce Evaluation Methods

```text
nb           Naive Bayes classification between conditions, the fold number may be precised after ':' (if not precised, fold_num=2)
lr           Logistic regression (slower than Naive Bayes) between conditions, the fold number may be precised after ':' (if not precised, fold_num=2)
sd           Standard deviation
rsd          Relative standard deviation
mc           Contrast of mean between conditions
rsdc         Contrast of relative standard deviation between conditions
ttest        T-test between conditions
es           Effect size between conditions
lfc:mean     Log2 fold change by group mean, 'mean' can be omitted as default value
lfc:median   Log2 fold change by group median
user:name    User-defined method, where name indicates a column in the k-mer count table
```

Sorting Modes

```text
dec          Sorting by decreasing order                              (as default for nb, lr, sd, rsd, user:name)
dec:abs      Sorting by decreasing order but on the absolute value    (as default for mc, rsdc, es, lfc:mean, lfc:median)
inc          Sorting by increasing order                              (as default for ttest)
inc:abs      Sorting by increasing order but on the absolute value
```

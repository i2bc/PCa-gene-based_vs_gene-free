# KaMRaT Usage

KaMRaT means "k-mer Matrix Reduction Toolbox", or "k-mer Matrix, Really Tremendous !".

The toolbox contains for now modules below:

- kamratMask which takes a k-mer count matrix as input, and selects/removes the k-mers in a given fasta file
- kamratMerge which takes a k-mer count matrix as input, and merge k-mers into longer contigs according to their overlap
- kamratNorm which takes a k-mer count matrix as input, and normalize k-mer counts
- kamratReduce which takes a k-mer count matrix as input, and evaluates the performance of each k-mer according to different metrics

## Prerequisite

The kamratReduce module is dependent on MLPack library. It could be installed from [conda cloud](https://anaconda.org/conda-forge/mlpack).

After installation from conda, please add the following line into your ```.bashrc``` file in the ```home/``` directory:

```bash
export LD_LIBRARY_PATH="/path_to_conda_env/mlpack/lib:$LD_LIBRARY_PATH"
```

For compiling the source, you can use ```compile.bash``` in the root directory of KaMRaT:

```bash
bash compile.bash /path_to_conda_environment_for_MLPack
```

And the executable files are in ```bin/``` directory.

## File of Sample Info

The sample-info file is given after the option ```-d```. This file aims to indicate which columns in the k-mer count matrix should be considered as samples. Please do not add any header line in the file, since the columns are defined as below :

- the first column indicates always sample names, if the file contains only one column, all samples are considered as in same condition
- if the file contains a second column, it conrresponds to condition for each sample

This sample-info file option can be omitted. In this case, all columns apart from the first one in the k-mer count matrix are considered as samples.

## kamratMask

Mask usage:

```text
kmerFilter [-h] [-n] [-k k_length] [-l] -f filter_path kmer_count_matrix_path
```

Mask parameter:

```text
-n         if the k-mers are generated from unstranded RNA-seq data
-k INT     the length of k-mers [31]
-l         if to REMOVE the k-mers in filter list (to keep "liquid" in real filtering experiment)
               do not put this parameter if to SELECT the k-mers in filter list (to keep "solid" in real filtering experiment)
-f STRING  the filter list path, could be either a fasta file or a list WITHOUT header
```

## kamratMerge

Merge usage

```text
kamratMerge [-h] [-n] [-k k_length] [-m min_overlap] [-d sample_info] [-i interv_method] [-q quant_mode] [-t tmp_dir] kmer_count_path
```

Merge parameters

```text
-h         Print the helper
-n         If the k-mers are generated from non-stranded RNA-seq data
-k INT     k-mer length (max_value: 32) [31]
-m INT     Min assembly overlap (max_value: k) [15]
-d STRING  Sample-info path, either list or table with sample names as the first column
               if absent, all columns except the first one in k-mer count matrix are taken as samples
-i STRING  Intervention method [none]
               none             without any intervention using count vector, only considering overlap between contigs
               mac:threshold    with intervention by mean absolute contrast between contigs' count vectors, by default, accept merge if mac <= 0.25
-q STRING  Quantification mode [rep]
               rep:colname      use the count vector of k-mer with lowest value in the 'colname' column for representing the corresponding contig's count vector
                                if no colname provided, the input order of k-mers will be taken for deciding representative k-mer (first in, first represent)
               mean             calculate the average of all composite k-mers of the corresponding contig as its count vector
-x         Query on disk [false]
-t STRING  Temporary directory [./]
```

## kamratNorm

Norm usage

```text
kamratNorm -b CHAR [-d STRING] [-T STRING] kmer_count_matrix_path
```

Norm parameters

```text
-h         Print the helper
-b CHAR    BASE for normalization (MANDATORY)
               B: count per billion
               M: count per million
               K: count per thousand
-d STRING  Sample-info path
               if absent, all columns will be printed
-T STRING  Transformation before evaluation (e.g. log)
```

## kamratReduce

Reduce usage

```text
kamratReduce [-d samp_info_path] [-m eval_method:fold_num] [-s sort_mode] [-N top_num] [-T transf_mode] [-C count_offset] kmer_count_path
```

Reduce parameters

```text
-h           Print the helper
-d STRING    Path to sample-condition or sample file, without header line
                 if absent, all except the first column in k-mer count matrix will be regarded as samples
-m STRING    Evaluation method name and their sub-evaluation-command, seperated by ':' (cf. Reduce evaluation methods)
-s STRING    Sorting mode, default value depends on evaluation method (c.f Sorting mode)
-N INT       Number of top features to print
-T STRING    Transformation before evaluation (e.g. -T log)
-C DOUBLE    Count offset to be added before transformation and evaluation
-r INT       Accepted minimum recurrence number among samples
-a INT       Accepted minimum count abundance for counting sample recurrence
```

Reduce Evaluation Methods

```text
nb           Naïve Bayes classification between conditions, fold number for cross-validation may be precised after ':' (if not precised, fold number = 2)
lr           Logistic regression (slower than Naive Bayes) between conditions, the fold number for cross-validation may be precised after ':' (if not precised, fold_num=2)
sd           Standard deviation
rsd          Relative standard deviation
ttest        p-value by t-test between conditions
es           Effect size between conditions
lfc:mean     Log2 fold change by group mean, 'mean' can be omitted as default value
lfc:median   Log2 fold change by group median
user:name    User-defined method, where name indicates a column in the k-mer count matrix
```

Sorting Modes

```text
dec          Sorting by decreasing order                              (as default for nb, lr, sd, rsd, user:name)
dec:abs      Sorting by decreasing order but on the absolute value    (as default for es, lfc:mean, lfc:median)
inc          Sorting by increasing order                              (as default for ttest)
inc:abs      Sorting by increasing order but on the absolute value
```

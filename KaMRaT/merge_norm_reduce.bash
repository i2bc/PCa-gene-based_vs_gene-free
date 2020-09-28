#!/bin/bash

raw_count_path=$1
smp_condition_path=$2
out_dir=$3

/store/EQUIPES/SSFA/haoliang-shared/KaMRaT_for_Paper_Ha/PCa-gene-based_vs_gene-free/KaMRaT/bin/kamratMerge -d $smp_condition_path -i mac -q mean -n -t $out_dir <(zcat $raw_count_path) > $out_dir/raw-counts-merge-TCGA-LRHR-3-10.tsv

/store/EQUIPES/SSFA/haoliang-shared/KaMRaT_for_Paper_Ha/PCa-gene-based_vs_gene-free/KaMRaT/bin/kamratNorm -b B -d $smp_condition_path -T log $out_dir/raw-counts-merge-TCGA-LRHR-3-10.tsv > $out_dir/raw-counts-merge-norm-TCGA-LRHR-3-10.tsv

/store/EQUIPES/SSFA/haoliang-shared/KaMRaT_for_Paper_Ha/PCa-gene-based_vs_gene-free/KaMRaT/bin/kamratReduce -d $smp_condition_path -N 500 -m nb:5 $out_dir/raw-counts-merge-norm-TCGA-LRHR-3-10.tsv > $out_dir/top_contig_merge_norm.nb5.out

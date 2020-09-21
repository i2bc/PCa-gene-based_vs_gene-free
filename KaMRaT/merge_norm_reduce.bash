#!/bin/bash

raw_count_path=/store/EQUIPES/SSFA/MEMBERS/thi-ngoc-ha.nguyen/Prostate_TCGA_LRvsHR/DEkupl_LRvsHR_limma_nonMask/Result_DEkupl_checking/raw-counts-TCGA-LRHR-3-10.tsv.gz         # raw count path
smp_condition_path=/store/EQUIPES/SSFA/MEMBERS/thi-ngoc-ha.nguyen/Prostate_TCGA_LRvsHR/sample_conditions.tsv         # sample condition path
out_dir=/store/EQUIPES/SSFA/haoliang-shared/KaMRaT_for_Paper_Ha/rec_results/          # output directory

source /home/haoliang.xue/.bashrc       # for export LD_LIBRARY_PATH with mlpack conda environment

/store/EQUIPES/SSFA/haoliang-shared/KaMRaT_for_Paper_Ha/PCa-gene-based_vs_gene-free/KaMRaT/bin/kamratMerge -i mac -q mean -n <(zcat $raw_count_path) > $out_dir/raw-counts-merge-TCGA-LRHR-3-10.tsv

/store/EQUIPES/SSFA/haoliang-shared/KaMRaT_for_Paper_Ha/PCa-gene-based_vs_gene-free/KaMRaT/bin/kamratNorm -b B -d $smp_condition_path -T log $out_dir/raw-counts-merge-TCGA-LRHR-3-10.tsv > $out_dir/raw-counts-merge-norm-TCGA-LRHR-3-10.tsv

/store/EQUIPES/SSFA/haoliang-shared/KaMRaT_for_Paper_Ha/PCa-gene-based_vs_gene-free/KaMRaT/bin/kamratReduce -d $smp_condition_path -N 500 -m nb:5 $out_dir/raw-counts-merge-norm-TCGA-LRHR-3-10.tsv > $out_dir/top2k_contig_merge_norm.nb5.out

# kmerFilter Usage

The kmerFilter tool aims to select/remove a list of k-mers in a k-mer count table.

Usage: 

```text
./kmerfilter [-n] [-k k_length] [-l] -f filter_list_path kmer_count_table_path
```

Parameter:

```text
-n         if the k-mers are generated from unstranded RNA-seq data
-k INT     the length of k-mers [31]
-l         if to REMOVE the k-mers in filter list (to keep "liquid" in real filtering experiment)
               do not put this parameter if to SELECT the k-mers in filter list (to keep "solid" in real filtering experiment)
-f STRING  the filter file path, could be either a fasta file or a list WITHOUT header
```

There are 2 principal modes: keep solid and keep liquid.

- keep solid mode (WITHOUT parameter -l) is for the case where we want to KEEP the k-mers from filter\_list. It is in same logic as mask step in DE-kupl
- keep liquid mode (WITH parameter -l) is for the case where we want to REMOVE the k-mers from filter\_list.

For compilation: simply type ```make``` in kmerFilter folder.

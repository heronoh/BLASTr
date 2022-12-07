
# BLASTr

<!-- badges: start -->
<!-- badges: end -->

BLASTr is an R package to facilitate identification of sequences (ASVs, OTUs, ...) derived from metabarcoding analyses. It wraps BLAST+ from NCBI parallelizes the search of each sequence to exploit the computational resources available. It can be used for any sequence and any BLAST+ formated database.

## Installation

You can install the development version of BLASTr from [GitHub](https://github.com/heronoh/BLASTr) with:

``` r
# install.packages("remotes")
remotes::install_github("heronoh/BLASTr")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(BLASTr)
options(BLASTr.num_threads = length(future::availableWorkers()) - 2)
options(BLASTr.db_path = fs::path())

ASVs_test <- c(
  "CTAGCCATAAACTTAAATGAAGCTATACTAAACTCGTTCGCCAGAGTACTACAAGCGAAAGCTTAAAACTCATAGGACTTGGCGGTGTTTCAGACCCAC",
  "CTAGCCATAAACTTAAATGAAGCTATACTAAACTCGTTCGCCAGAGTACTACAAGTGAAAGCTTAAAACTCATAGGACTTGGCGGTGTTTCAGACCCAC",
  "GCCAAATTTGTGTTTTGTCCTTCGTTTTTAGTTAATTGTTACTGGCAAATGACTAACGACAAATGATAAATTACTAATAC",
  "AACATTGTATTTTGTCTTTGGGGCCTGGGCAGGTGCAGTAGGAACTTCACTTAGAATAATTATTCGTACTGAGCTTGGGCATCCAGGAAGACTTATCGGGGATGATCAAATCTATAATGTAATTGTTACAGCACATGCATTTGTGATAATTTTTTTTATAGTAATACCTATTATGATT",
  "ACTATACCTATTATTCGGCGCATGAGCTGGAGTCCTAGGCACAGCTCTAAGCCTCCTTATTCGAGCCGAGCTGGGCCAGCCAGGCAACCTTCTAGGTAACGACCACATCTACAACGTTATCGTCACAGCCCATGCATTTGTAATAATCTTCTTCATAGTAATACCCATCATAATCGGAGGCTTTGGCAACTGACTAGTTCCCCTAATAATCGGTGCCCCCGATATG",
  "TTAGCCATAAACATAAAAGTTCACATAACAAGAACTTTTGCCCGAGAACTACTAGCAACAGCTTAAAACTCAAAGGACTTGGCGGTGCTTTATATCCAC"
)


blast_res <- BLASTr::parallel_blast(
  asvs = ASVs_test[1],
  db_path = "/data/databases/nt/nt",
  out_file = NULL,
  out_RDS = NULL,
  perc_id = 80,
  num_threads = 1,
  perc_qcov_hsp = 80,
  num_alignments = 3,
  blast_type = "blastn"
)

```

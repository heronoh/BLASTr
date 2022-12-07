
# BLASTr

<!-- badges: start -->
<!-- badges: end -->

The BLASTr package is a powerful tool for performing BLAST operations from within the R environment.
Despite being initially developed for metagenomic applications, the BLASTr package is flexible.
It can be used with any installed NCBI BLAST+ search strategy, configuration, and any database on UNIX and Windows platforms through Windows Subsystem for Linux.
This makes it a versatile tool that can be used in various applications and contexts requiring sequence identification on tabular outputs.
Additionally, the package includes documentation and functions for parsing and analyzing the results of BLAST searches, making it easier for users to extract useful information from their BLAST results.
Overall, the BLASTr package is a valuable tool for bioinformaticians and researchers who need to perform BLAST operations from within R.

## Installation

### R-Universe

The main way to install `BLASTr` is through:

``` r
install.packages('BLASTr', repos = "https://heronoh.r-universe.dev")
```

### Development version

You can install the development version of BLASTr from [GitHub](https://github.com/heronoh/BLASTr) with:

``` r
# install.packages("remotes")
remotes::install_github("heronoh/BLASTr")
```

## Documentation

The main place to look for help and documentation is <heronoh.github.io/BLASTr/>

## Usage

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

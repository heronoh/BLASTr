
# BLASTr

<!-- badges: start -->
<!-- badges: end -->

The goal of BLASTr is to ...

## Installation

You can install the development version of BLASTr from [GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("villabioinfo/BLASTr")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(BLASTr)
options(BLASTr.num_threads = length(future::availableWorkers()) - 2)
options(BLASTr.db_path = fs::path())

run_blast()

```


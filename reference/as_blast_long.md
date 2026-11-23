# Convert BLASTr results to the long format

Turns wide BLASTr results (one row per sequence, hits `1_` to `n_`, from
any version, including `k_indentity`) into the long format of
[`parallel_blast()`](https://heronoh.github.io/BLASTr/reference/parallel_blast.md),
ranking the hits with the same rules, so that old results can be read
with
[`blast_summary()`](https://heronoh.github.io/BLASTr/reference/blast_summary.md)
and [`get_tax()`](https://heronoh.github.io/BLASTr/reference/get_tax.md)
without running BLAST again. The BLAST order is taken from
`k_blast_rank` or, in older results, from `k`. Long results are returned
as they are.

## Usage

``` r
as_blast_long(
  x,
  rank_by = c("similarity", "identity", "bitscore"),
  tolerance = 0.5,
  perc_qcov_hsp = 80,
  uninformative = NULL
)
```

## Arguments

- x:

  BLASTr results, wide or long. The sequence is in `Sequence` (or
  `sequence`).

- rank_by, tolerance:

  As in
  [`check_hits()`](https://heronoh.github.io/BLASTr/reference/check_hits.md).

- perc_qcov_hsp:

  Coverage under which a hit is flagged `low_qcov`.

- uninformative:

  `NULL` (default) keeps the flags `informative` of `x` (all hits
  informative when `x` has none); a table of
  [`uninformative_references()`](https://heronoh.github.io/BLASTr/reference/uninformative_references.md)
  flags the hits again with it, and ranks them among the informative
  ones, as
  [`parallel_blast()`](https://heronoh.github.io/BLASTr/reference/parallel_blast.md)
  does.

## Value

Long results (see
[`parallel_blast()`](https://heronoh.github.io/BLASTr/reference/parallel_blast.md)),
with all the hits of `x`.

## Examples

``` r
wide <- data.frame(
  Sequence = "ACGT", `1_subject` = "A1", `1_indentity` = 99, `1_qcovhsp` = 100,
  `1_staxid` = "9606", `2_subject` = "A2", `2_indentity` = 100,
  `2_qcovhsp` = 100, `2_staxid` = "9598", check.names = FALSE
)
as_blast_long(wide)
#> # A tibble: 2 × 25
#>   Sequence  rank tie   at_best low_qcov informative uninformative_category
#>   <chr>    <int> <lgl> <lgl>   <lgl>    <lgl>       <chr>                 
#> 1 ACGT         1 FALSE TRUE    FALSE    TRUE        NA                    
#> 2 ACGT         2 FALSE FALSE   FALSE    TRUE        NA                    
#> # ℹ 18 more variables: blast_rank <int>, subject <chr>, `subject header` <chr>,
#> #   staxid <chr>, identity <dbl>, `Similarity score` <dbl>, qcovhsp <dbl>,
#> #   nident <dbl>, length <dbl>, mismatches <dbl>, gaps <dbl>,
#> #   `query start` <dbl>, `query end` <dbl>, `subject start` <dbl>,
#> #   `subject end` <dbl>, `subject strand` <chr>, `e-value` <dbl>,
#> #   bitscore <dbl>
```

# Summary of the BLAST hits of each sequence

One row per sequence of a long
[`parallel_blast()`](https://heronoh.github.io/BLASTr/reference/parallel_blast.md)
result (or of
[`as_blast_long()`](https://heronoh.github.io/BLASTr/reference/as_blast_long.md)):
`Check ID` and the other summary columns described in
[`check_hits()`](https://heronoh.github.io/BLASTr/reference/check_hits.md),
and `Best low qcov` (the best hit covers less of the query than
`perc_qcov_hsp`: the identification rests on a partial alignment).
Sequences without hits have `N hits = 0` and `NA` elsewhere.

## Usage

``` r
blast_summary(
  x,
  check_n = NULL,
  tolerance = NULL,
  rank_by = NULL,
  search_alignments = NULL
)
```

## Arguments

- x:

  Long BLASTr results.

- check_n, tolerance, rank_by, search_alignments:

  As in
  [`check_hits()`](https://heronoh.github.io/BLASTr/reference/check_hits.md).
  `NULL` (default) uses the values of the search, kept in the attribute
  `BLASTr_metadata`, or else `10`, `0.5`, `"similarity"` and `NA`.

## Value

A tibble with `Sequence` and the summary columns.

## Examples

``` r
if (FALSE) { # \dontrun{
res <- parallel_blast(seqs, db_path)
dplyr::left_join(blast_summary(res), res[res$rank %in% 1L, ], by = "Sequence")
} # }
```

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

The long results keep only the first `num_alignments` hits and the ones
tied with the best, but `N hits`, `N taxids` and `N uninformative` count
all the subjects BLAST returned (up to `search_alignments`), as do
`N hits same taxid as 1st` and `Best uninformative score` when `rank_by`
is the one of the search:
[`parallel_blast()`](https://heronoh.github.io/BLASTr/reference/parallel_blast.md)
keeps these counts in the attribute `BLASTr_metadata` (`hit_counts`).
For sequences without them (results of older versions, other tables)
they are computed from the rows of `x`. So `N hits` equal to
`search_alignments` means BLAST returned as many subjects as it was
asked for: the list may have been cut.

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

## See also

Other result functions:
[`as_blast_long()`](https://heronoh.github.io/BLASTr/reference/as_blast_long.md),
[`blast_lca()`](https://heronoh.github.io/BLASTr/reference/blast_lca.md),
[`check_hits()`](https://heronoh.github.io/BLASTr/reference/check_hits.md),
[`similarity_score()`](https://heronoh.github.io/BLASTr/reference/similarity_score.md),
[`uninformative_references()`](https://heronoh.github.io/BLASTr/reference/uninformative_references.md)

## Examples

``` r
# hits of two sequences: the first one tied between Homo and Pan
wide <- data.frame(
  Sequence = c("ACGT", "TTGA"),
  `1_subject` = c("A1", "B1"), `1_identity` = c(100, 99), `1_qcovhsp` = 100,
  `1_staxid` = "9606",
  `2_subject` = c("A2", "B2"), `2_identity` = c(99.8, 90), `2_qcovhsp` = 100,
  `2_staxid` = "9598",
  check.names = FALSE
)
hits <- as_blast_long(wide)
blast_summary(hits)
#> # A tibble: 2 × 12
#>   Sequence `Check ID` `N hits` `N taxids` `N hits same taxid as 1st`
#>   <chr>    <chr>         <int>      <int>                      <int>
#> 1 ACGT     check ID          2          2                          1
#> 2 TTGA     best ID           2          2                          1
#> # ℹ 7 more variables: `N hits at best` <int>, `N taxids at best` <int>,
#> #   `Taxids at best` <chr>, `At best truncated` <lgl>, `N uninformative` <int>,
#> #   `Best uninformative score` <dbl>, `Best low qcov` <lgl>

# needs BLAST+ (see install_dependencies())
if (has_blast()) {
  fasta <- system.file("extdata", "minimal_db_blast.fasta", package = "BLASTr")
  db <- file.path(tempdir(), "minimal_db_blast")
  make_blast_db(fasta, db_path = db)
  seqs <- readLines(system.file("extdata", "seqs_test.txt", package = "BLASTr"))
  res <- parallel_blast(seqs, db)
  blast_summary(res)
}
```

# Check the agreement of the best BLAST hits

For each query sequence, orders its hits by `rank_by`, decides whether
the best ones agree on the taxon (`Check ID`) and summarises all of
them.

The hits are ordered by `rank_by` (Similarity score, identity or
bitscore) and then by identity, coverage, bitscore and BLAST order. Hits
within `tolerance` of the best value are "at best": they are tied
candidates for the identification.

`Check ID` looks at the first `check_n` hits and compares only the ones
at best: a hit clearly worse than the best one is most likely another
taxon, not a competing identification.

- `"best ID"`: all the hits at best among the first `check_n` have the
  taxid of the first hit.

- `"check ID"`: some of them has another taxid or an unknown one (`0` or
  empty, as in databases built without a taxid map).

- `NA`: no hits.

The other columns use all the hits of the sequence: `N hits`, `N taxids`
(distinct known taxids), `N hits same taxid as 1st`, `N hits at best`,
`N taxids at best`, `Taxids at best` (the known taxids of the hits at
best, in rank order, separated by `";"`, for a lowest common ancestor)
and `At best truncated` (`TRUE` when the last subject BLAST returned is
at best: more tied subjects may exist beyond `search_alignments`; `NA`
when `search_alignments` is unknown).

## Usage

``` r
check_hits(
  hits,
  check_n = 10L,
  tolerance = 0.5,
  rank_by = c("similarity", "identity", "bitscore"),
  search_alignments = NA_integer_
)
```

## Arguments

- hits:

  Data frame with one row per hit and the columns `query`, `identity`,
  `staxid` and, as needed by `rank_by`, `qcovhsp` or `bitscore`. The
  BLAST order is taken from `blast_rank` or `rank`, or else from the row
  order.

- check_n:

  Number of best hits used by `Check ID`.

- tolerance:

  Difference from the best value, in the units of `rank_by` (points of
  Similarity score or identity, bits of bitscore), under which hits are
  tied.

- rank_by:

  Criterion to order the hits: `"similarity"` (Similarity score, see
  [`similarity_score()`](https://heronoh.github.io/BLASTr/reference/similarity_score.md)),
  `"identity"` or `"bitscore"`.

- search_alignments:

  Number of subjects BLAST was asked for, to compute
  `At best truncated`.

## Value

A tibble with `query` and the summary columns.

## Examples

``` r
hits <- data.frame(
  query = "q1", identity = c(100, 99, 98.7), qcovhsp = 100,
  staxid = c("9606", "9598", "9597")
)
check_hits(hits)
#> # A tibble: 1 × 9
#>   query `Check ID` `N hits` `N taxids` N hits same taxid as 1…¹ `N hits at best`
#>   <chr> <chr>         <int>      <int>                    <int>            <int>
#> 1 q1    best ID           3          3                        1                1
#> # ℹ abbreviated name: ¹​`N hits same taxid as 1st`
#> # ℹ 3 more variables: `N taxids at best` <int>, `Taxids at best` <chr>,
#> #   `At best truncated` <lgl>
```

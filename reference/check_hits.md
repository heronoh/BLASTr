# Check the agreement of the best BLAST hits

For each query sequence, decides whether its best hits agree on the
taxon (`Check ID`) and summarises all of its hits.

`Check ID` looks at the first `check_n` hits. Only the hits whose
identity is within `tolerance` percentage points of the best identity
among them are compared: a hit clearly less identical than the best one
is most likely another taxon, not a competing identification. Coverage
is not used here, so a good reference that covers only part of the query
does not hide or create a disagreement.

- `"best ID"`: all the compared hits have the taxid of the first hit.

- `"check ID"`: some compared hit has another taxid or an unknown one
  (`0` or empty, as in databases built without a taxid map), or the
  first hit is not among the most identical ones.

- `NA`: no hits.

The other columns use all the hits of the sequence: `N hits`, `N taxids`
(distinct known taxids), `N hits same taxid as 1st` and
`N taxids at best identity` (distinct taxids among the hits within
`tolerance` of the best identity of all hits).

## Usage

``` r
check_hits(hits, check_n = 3L, tolerance = 0.5)
```

## Arguments

- hits:

  Data frame with one row per hit and the columns `query`, `rank`,
  `identity` and `staxid`.

- check_n:

  Number of best hits used by `Check ID`.

- tolerance:

  Identity difference (percentage points) under which two hits are
  considered tied.

## Value

A tibble with `query` and the summary columns.

## Examples

``` r
hits <- data.frame(
  query = "q1", rank = 1:3,
  identity = c(100, 99, 98.7), staxid = c("9606", "9598", "9597")
)
check_hits(hits)
#> # A tibble: 1 × 6
#>   query `Check ID` `N hits` `N taxids` `N hits same taxid as 1st`
#>   <chr> <chr>         <int>      <int>                      <int>
#> 1 q1    best ID           3          3                          1
#> # ℹ 1 more variable: `N taxids at best identity` <int>
```

# Similarity score of BLAST hits

Combines identity and query coverage in one value, weighting identity
ten times more: `(10 * identity + qcovhsp) / 11`. It was called "BLASTn
pseudo-score" in earlier versions of the metabarcoding pipeline.

## Usage

``` r
similarity_score(identity, qcovhsp)
```

## Arguments

- identity:

  Percent identity of the hit (`pident`).

- qcovhsp:

  Percent of the query covered by the hit (`qcovhsp`).

## Value

Numeric vector, same length as the inputs.

## Examples

``` r
similarity_score(identity = c(100, 98.5), qcovhsp = c(100, 85))
#> [1] 100.00000  97.27273
```

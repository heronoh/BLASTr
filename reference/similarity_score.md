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

## See also

Other result functions:
[`as_blast_long()`](https://heronoh.github.io/BLASTr/reference/as_blast_long.md),
[`blast_lca()`](https://heronoh.github.io/BLASTr/reference/blast_lca.md),
[`blast_summary()`](https://heronoh.github.io/BLASTr/reference/blast_summary.md),
[`check_hits()`](https://heronoh.github.io/BLASTr/reference/check_hits.md),
[`uninformative_references()`](https://heronoh.github.io/BLASTr/reference/uninformative_references.md)

## Examples

``` r
similarity_score(identity = c(100, 98.5), qcovhsp = c(100, 85))
#> [1] 100.00000  97.27273
```

# Consensus lowest common ancestor of the best BLAST hits

For each sequence, the taxon shared by the hits tied with the best one
(`at_best`, see
[`parallel_blast()`](https://heronoh.github.io/BLASTr/reference/parallel_blast.md)):
the deepest node of the NCBI taxonomy that contains at least `consensus`
of them. With `consensus = 1` it is the strict lowest common ancestor;
below 1, a few mislabelled references among many tied hits (e.g. one
"Picramnia" among 44 *Picea*) do not push the identification up to a
high rank.

The lineages come from
[`get_tax()`](https://heronoh.github.io/BLASTr/reference/get_tax.md)
(`lineage_taxids`), so the ancestor is exact even with homonyms and
clades without rank. Hits without a known taxid (`0`, empty, not in the
taxonomy) do not vote, nor do the hits flagged as not `informative` (see
[`uninformative_references()`](https://heronoh.github.io/BLASTr/reference/uninformative_references.md)).
A hit with several taxids votes with the first one.

The default `0.9` was chosen on sequences of known taxonomy (see
`dev/benchmark/R/10_rank_calibration.R`): compared with the strict
ancestor it resolved 3 points more genera with almost no more errors,
while also accepting one dissenting hit among few ("N - 1") resolved
more but about half of the extra ones wrong. With up to 9 tied hits,
`0.9` is the strict ancestor.

## Usage

``` r
blast_lca(
  x,
  consensus = 0.9,
  named_first = TRUE,
  taxonomy = NULL,
  taxonomy_db = NULL,
  taxdump = NULL
)
```

## Arguments

- x:

  BLASTr results: long, as returned by
  [`parallel_blast()`](https://heronoh.github.io/BLASTr/reference/parallel_blast.md),
  or wide (converted with
  [`as_blast_long()`](https://heronoh.github.io/BLASTr/reference/as_blast_long.md)).

- consensus:

  Fraction of the tied hits (with known taxid) that the ancestor must
  contain, above 0.5. Defaults to `0.9`.

- named_first:

  If `TRUE` (default), the tied hits named only above genus (e.g.
  "Bacterium", "alpha proteobacterium XYZ", "Entolomataceae sp.") do not
  vote when a named species is tied with them: their lineage is right
  but shallow, and would pull the ancestor up. Hits in open nomenclature
  within a genus ("Picea sp.", "Carex cf. flava") still vote: they are
  other species of that genus (or of another one), and keep the ancestor
  at the genus (or above) instead of on the named species. They count in
  `N tied`; `LCA support` is among the voters.

- taxonomy:

  The table of
  [`get_tax()`](https://heronoh.github.io/BLASTr/reference/get_tax.md)
  for the taxids of `x`. `NULL` (default) calls
  `get_tax(x, taxonomy_db, taxdump = taxdump)`. Its rank columns define
  the ranks of `LCA ranked ...`: keep only some (e.g. `species`,
  `genus`, `family`, `order`) to place the ranked ancestor on them.

- taxonomy_db, taxdump:

  Passed to
  [`get_tax()`](https://heronoh.github.io/BLASTr/reference/get_tax.md).

## Value

A tibble with one row per sequence: `Sequence`, `LCA taxid`, `LCA name`,
`LCA rank` (`"genus"`, `"family"`, ...; `NA` for a node without one of
the ranks of
[`get_tax()`](https://heronoh.github.io/BLASTr/reference/get_tax.md),
e.g. a clade or a superfamily), `LCA ranked taxid`, `LCA ranked name`,
`LCA ranked rank` (the ancestor itself when it has one of those ranks,
otherwise the deepest node above it that has one; a species in open
nomenclature, e.g. "Picea sp. X", goes to its genus: for tables by
rank), `LCA support` (fraction of the tied hits it contains), `N tied`
(hits at best) and `N tied no taxid`. `NA` for sequences without hits.
Join `get_tax(result$`LCA taxid`)` for its ranks and lineage.

## Examples

``` r
if (FALSE) { # \dontrun{
res <- parallel_blast(seqs, db_path)
lca <- blast_lca(res)
dplyr::left_join(lca, get_tax(lca$`LCA taxid`), by = c(`LCA taxid` = "taxid"))
} # }
```

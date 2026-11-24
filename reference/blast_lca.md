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

## See also

Other result functions:
[`as_blast_long()`](https://heronoh.github.io/BLASTr/reference/as_blast_long.md),
[`blast_summary()`](https://heronoh.github.io/BLASTr/reference/blast_summary.md),
[`check_hits()`](https://heronoh.github.io/BLASTr/reference/check_hits.md),
[`similarity_score()`](https://heronoh.github.io/BLASTr/reference/similarity_score.md),
[`uninformative_references()`](https://heronoh.github.io/BLASTr/reference/uninformative_references.md)

## Examples

``` r
# a simplified NCBI taxdump: Hominidae with the genera Homo and Pan
taxdump <- file.path(tempdir(), "taxdump")
dir.create(taxdump, showWarnings = FALSE)
dmp <- function(...) paste0(paste(..., sep = "\t|\t"), "\t|")
writeLines(c(
  dmp(1, 1, "no rank"), dmp(2759, 1, "superkingdom"), dmp(9604, 2759, "family"),
  dmp(9605, 9604, "genus"), dmp(9606, 9605, "species"),
  dmp(9596, 9604, "genus"), dmp(9598, 9596, "species")
), file.path(taxdump, "nodes.dmp"))
writeLines(c(
  dmp(1, "root", "", "scientific name"), dmp(2759, "Eukaryota", "", "scientific name"),
  dmp(9604, "Hominidae", "", "scientific name"), dmp(9605, "Homo", "", "scientific name"),
  dmp(9606, "Homo sapiens", "", "scientific name"), dmp(9596, "Pan", "", "scientific name"),
  dmp(9598, "Pan troglodytes", "", "scientific name")
), file.path(taxdump, "names.dmp"))
tax_db <- make_taxonomy_db(tempfile(fileext = ".sqlite"), taxdump = taxdump)
#> ✔ Taxonomy database: /tmp/Rtmp39IRAe/file19216b1f6cf1.sqlite (7 taxa).

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
lca <- blast_lca(hits, taxonomy_db = tax_db)
lca
#> # A tibble: 2 × 10
#>   Sequence `LCA taxid` `LCA name`   `LCA rank` `LCA ranked taxid`
#>   <chr>    <chr>       <chr>        <chr>      <chr>             
#> 1 ACGT     9604        Hominidae    family     9604              
#> 2 TTGA     9606        Homo sapiens species    9606              
#> # ℹ 5 more variables: `LCA ranked name` <chr>, `LCA ranked rank` <chr>,
#> #   `LCA support` <dbl>, `N tied` <int>, `N tied no taxid` <int>
# ranks and lineage of each ancestor
dplyr::left_join(lca, get_tax(lca$`LCA taxid`, taxonomy_db = tax_db), by = c(`LCA taxid` = "taxid"))
#> # A tibble: 2 × 28
#>   Sequence `LCA taxid` `LCA name`   `LCA rank` `LCA ranked taxid`
#>   <chr>    <chr>       <chr>        <chr>      <chr>             
#> 1 ACGT     9604        Hominidae    family     9604              
#> 2 TTGA     9606        Homo sapiens species    9606              
#> # ℹ 23 more variables: `LCA ranked name` <chr>, `LCA ranked rank` <chr>,
#> #   `LCA support` <dbl>, `N tied` <int>, `N tied no taxid` <int>, found <lgl>,
#> #   sci_name <chr>, rank <chr>, domain <chr>, kingdom <chr>, phylum <chr>,
#> #   subphylum <chr>, class <chr>, subclass <chr>, order <chr>, suborder <chr>,
#> #   family <chr>, subfamily <chr>, genus <chr>, species <chr>, lineage <chr>,
#> #   lineage_taxids <chr>, merged_into <chr>
```

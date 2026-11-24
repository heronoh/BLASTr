# Taxonomy of NCBI taxids from a local database

Looks up the NCBI taxonomy of taxids in a local SQLite database: no
internet, no `efetch`, thousands of taxids in a few seconds. Returns the
minimal table, with one row per distinct taxid, to be joined to the
BLAST hits by taxid.

`x` can be a vector of taxids or any result of
[`parallel_blast()`](https://heronoh.github.io/BLASTr/reference/parallel_blast.md)
(long or wide) or
[`parallel_blast_hsps()`](https://heronoh.github.io/BLASTr/reference/parallel_blast_hsps.md):
the taxids are taken from the `staxid`, `k_staxid` and `Taxids at best`
columns. Unknown taxids (`0`, empty, `N/A`) are dropped.

The database comes with the BLAST databases downloaded by
[`download_blast_db()`](https://heronoh.github.io/BLASTr/reference/ncbi_blast_dbs.md)
(`taxonomy = TRUE`). Without it, `get_tax()` builds it only when asked
(`download = TRUE`, a `taxdump`, or a yes in an interactive session),
from the NCBI taxdump (about 80 MB, downloaded once, then reused);
otherwise it stops and says how. If `x` is a BLASTr result, the taxdump
is the first one the NCBI published after the date of its BLAST
database, so the taxonomy matches the database; the database is kept in
the BLASTr cache as `taxonomy_<database date>.sqlite`. With only taxids
and no database set, the newest database in the cache is used, or one is
built from the current taxdump.

The database has the tables `nodes` (`id`, `rank`, `parent`) and `names`
(`id`, `name`, `scientific`, `type`), as built by
[`make_taxonomy_db()`](https://heronoh.github.io/BLASTr/reference/make_taxonomy_db.md)
(or by `taxonomizr::prepareDatabase()`, whose file can be used as is).
Build it from the NCBI taxdump of the same date as the BLAST database:
taxids newer than the taxdump are not found.

## Usage

``` r
get_tax(
  x,
  taxonomy_db = NULL,
  ranks = tax_ranks_default,
  taxdump = NULL,
  download = FALSE
)
```

## Arguments

- x:

  Taxids (character or numeric) or a BLASTr result.

- taxonomy_db:

  Path to the SQLite taxonomy database. `NULL` uses the option
  `BLASTr.taxonomy_db` or the environment variable `BLASTR_TAXONOMY_DB`;
  with none set, a database in the BLASTr cache (see below). If the file
  does not exist, it is built.

- ranks:

  Ranks returned as columns. `"domain"` is also filled from
  `"superkingdom"`, its name in taxdumps before 2025.

- taxdump:

  NCBI taxdump used to build a missing database (see
  [`make_taxonomy_db()`](https://heronoh.github.io/BLASTr/reference/make_taxonomy_db.md)).
  `NULL`: the taxdump matching the BLAST database of `x` (the first
  monthly NCBI archive after its date), or the current one.

- download:

  Download the NCBI taxdump and build the database when it does not
  exist. `FALSE` (default) stops with a message instead (or asks, in an
  interactive session).

## Value

A tibble with one row per taxid: `taxid`, `found` (`FALSE` when the
taxid is not in the database), `sci_name`, `rank` (rank of the taxid
itself), one column per rank in `ranks` (`NA` when the lineage has no
node of that rank), `lineage` (all the nodes from the root to the taxid,
without "root" and "cellular organisms", separated by `"; "`),
`lineage_taxids` (their taxids, in the same order, for an exact lowest
common ancestor) and `merged_into` (the current taxid when the NCBI
merged the taxid into another one; then the taxonomy is that of the
current taxid; needs a database with the table `merged`, as built by
[`make_taxonomy_db()`](https://heronoh.github.io/BLASTr/reference/make_taxonomy_db.md)).

## See also

Other taxonomy functions:
[`get_tax_online`](https://heronoh.github.io/BLASTr/reference/get_tax_online.md),
[`make_taxonomy_db()`](https://heronoh.github.io/BLASTr/reference/make_taxonomy_db.md)

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
#> ✔ Taxonomy database: /tmp/Rtmp39IRAe/file19211b790915.sqlite (7 taxa).

# 1234 is not in this taxonomy: found = FALSE, with a warning
get_tax(c("9606", "9598", "1234"), taxonomy_db = tax_db)
#> Warning: ! 1 taxid not found in the taxonomy database: 1234
#> ℹ Use a taxonomy database as recent as the BLAST database (see
#>   `make_taxonomy_db()`).
#> # A tibble: 3 × 19
#>   taxid found sci_name      rank  domain kingdom phylum subphylum class subclass
#>   <chr> <lgl> <chr>         <chr> <chr>  <chr>   <chr>  <chr>     <chr> <chr>   
#> 1 9606  TRUE  Homo sapiens  spec… Eukar… NA      NA     NA        NA    NA      
#> 2 9598  TRUE  Pan troglody… spec… Eukar… NA      NA     NA        NA    NA      
#> 3 1234  FALSE NA            NA    NA     NA      NA     NA        NA    NA      
#> # ℹ 9 more variables: order <chr>, suborder <chr>, family <chr>,
#> #   subfamily <chr>, genus <chr>, species <chr>, lineage <chr>,
#> #   lineage_taxids <chr>, merged_into <chr>

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
tax <- get_tax(hits, taxonomy_db = tax_db)
dplyr::left_join(hits, tax, by = c(staxid = "taxid"))
#> # A tibble: 4 × 43
#>   Sequence rank.x tie   at_best low_qcov informative uninformative_category
#>   <chr>     <int> <lgl> <lgl>   <lgl>    <lgl>       <chr>                 
#> 1 ACGT          1 FALSE TRUE    FALSE    NA          NA                    
#> 2 ACGT          2 FALSE TRUE    FALSE    NA          NA                    
#> 3 TTGA          1 FALSE TRUE    FALSE    NA          NA                    
#> 4 TTGA          2 FALSE FALSE   FALSE    NA          NA                    
#> # ℹ 36 more variables: blast_rank <int>, subject <chr>, `subject header` <chr>,
#> #   staxid <chr>, identity <dbl>, `Similarity score` <dbl>, qcovhsp <dbl>,
#> #   nident <dbl>, length <dbl>, mismatches <dbl>, gaps <dbl>,
#> #   `query start` <dbl>, `query end` <dbl>, `subject start` <dbl>,
#> #   `subject end` <dbl>, `subject strand` <chr>, `e-value` <dbl>,
#> #   bitscore <dbl>, found <lgl>, sci_name <chr>, rank.y <chr>, domain <chr>,
#> #   kingdom <chr>, phylum <chr>, subphylum <chr>, class <chr>, …
```

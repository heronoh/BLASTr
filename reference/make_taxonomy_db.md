# Build the local taxonomy database

Builds the SQLite database used by
[`get_tax()`](https://heronoh.github.io/BLASTr/reference/get_tax.md)
from the NCBI taxdump (`nodes.dmp`, `names.dmp` and `merged.dmp`),
downloading it if not given. The tables `nodes` and `names` follow the
`taxonomizr` layout, so a database built by
`taxonomizr::prepareDatabase()` also works (without resolving merged
taxids); `merged` (`old_id`, `new_id`) maps the taxids the NCBI merged
into others.

Use a taxdump of the same date as the BLAST database (or newer): taxids
created after the taxdump are not found.

## Usage

``` r
make_taxonomy_db(path, taxdump = NULL, overwrite = FALSE)
```

## Arguments

- path:

  Path of the SQLite file to create.

- taxdump:

  NCBI taxdump: a `taxdump.tar.gz` or `taxdmp_*.zip` file, its URL, or a
  directory with `nodes.dmp`, `names.dmp` and (optional) `merged.dmp`.
  `NULL` downloads the current taxdump from the NCBI.

- overwrite:

  Replace `path` if it exists.

## Value

The path of the database, invisibly.

## See also

Other taxonomy functions:
[`get_tax()`](https://heronoh.github.io/BLASTr/reference/get_tax.md),
[`get_tax_online`](https://heronoh.github.io/BLASTr/reference/get_tax_online.md),
[`has_taxonomy_db()`](https://heronoh.github.io/BLASTr/reference/has_taxonomy_db.md)

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
#> ✔ Taxonomy database: /tmp/RtmpYW6x5r/file19ca2b29229f.sqlite (7 taxa).
get_tax("9606", taxonomy_db = tax_db)
#> # A tibble: 1 × 19
#>   taxid found sci_name     rank   domain kingdom phylum subphylum class subclass
#>   <chr> <lgl> <chr>        <chr>  <chr>  <chr>   <chr>  <chr>     <chr> <chr>   
#> 1 9606  TRUE  Homo sapiens speci… Eukar… NA      NA     NA        NA    NA      
#> # ℹ 9 more variables: order <chr>, suborder <chr>, family <chr>,
#> #   subfamily <chr>, genus <chr>, species <chr>, lineage <chr>,
#> #   lineage_taxids <chr>, merged_into <chr>

if (FALSE) { # \dontrun{
# from the current NCBI taxdump (downloaded, about 70 MB)
make_taxonomy_db("~/blast_databases/taxonomy.sqlite")
options(BLASTr.taxonomy_db = "~/blast_databases/taxonomy.sqlite")
} # }
```

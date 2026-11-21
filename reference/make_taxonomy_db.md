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

## Examples

``` r
if (FALSE) { # \dontrun{
make_taxonomy_db("/data/databases/tax/taxonomy.sqlite")
options(BLASTr.taxonomy_db = "/data/databases/tax/taxonomy.sqlite")
} # }
```

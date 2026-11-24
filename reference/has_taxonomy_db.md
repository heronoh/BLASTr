# Is there a taxonomy database for get_tax()?

Checks, without downloading or building anything, whether
[`get_tax()`](https://heronoh.github.io/BLASTr/reference/get_tax.md) has
its taxonomy database: the one it would use with the same `taxonomy_db`,
or, for the results of a search in `db_path`, the one of the date of
those BLAST databases in the BLASTr cache. Useful to stop at the start
of a long analysis instead of after the BLAST search, as
[`has_blast()`](https://heronoh.github.io/BLASTr/reference/has_blast.md)
does for BLAST+.

The database counts when it is an SQLite file with the tables `nodes`
and `names`. When there is none,
[`download_blast_db()`](https://heronoh.github.io/BLASTr/reference/ncbi_blast_dbs.md)
brings it with the BLAST database, and `get_tax(..., download = TRUE)`
builds it from the NCBI taxdump.

## Usage

``` r
has_taxonomy_db(taxonomy_db = NULL, db_path = NULL)
```

## Arguments

- taxonomy_db:

  Path to the SQLite taxonomy database, as in
  [`get_tax()`](https://heronoh.github.io/BLASTr/reference/get_tax.md).
  `NULL` uses the option `BLASTr.taxonomy_db` or the environment
  variable `BLASTR_TAXONOMY_DB`; with none set, a database in the BLASTr
  cache.

- db_path:

  Path(s) of the BLAST database(s) to be searched, as in
  [`parallel_blast()`](https://heronoh.github.io/BLASTr/reference/parallel_blast.md):
  the database in the cache is then the one of their date, the one
  [`get_tax()`](https://heronoh.github.io/BLASTr/reference/get_tax.md)
  uses for the results of that search. `NULL`: the newest database in
  the cache, the one
  [`get_tax()`](https://heronoh.github.io/BLASTr/reference/get_tax.md)
  uses for taxids alone. Ignored when a `taxonomy_db` is set.

## Value

`TRUE` or `FALSE`, with the path of the database looked for in the
attribute `"path"`.

## See also

Other taxonomy functions:
[`get_tax()`](https://heronoh.github.io/BLASTr/reference/get_tax.md),
[`get_tax_online`](https://heronoh.github.io/BLASTr/reference/get_tax_online.md),
[`make_taxonomy_db()`](https://heronoh.github.io/BLASTr/reference/make_taxonomy_db.md)

## Examples

``` r
has_taxonomy_db()
#> [1] FALSE
#> attr(,"path")
#> [1] "/home/runner/.cache/R/BLASTr/taxonomy/taxonomy_2026-11-24.sqlite"

# a database that does not exist: FALSE, and nothing is created
x <- has_taxonomy_db(file.path(tempdir(), "taxonomy.sqlite"))
x
#> [1] FALSE
#> attr(,"path")
#> [1] "/tmp/RtmpYW6x5r/taxonomy.sqlite"
attr(x, "path")
#> [1] "/tmp/RtmpYW6x5r/taxonomy.sqlite"
```

# Taxonomy of NCBI taxids from a local database

Looks up the NCBI taxonomy of taxids in a local SQLite database: no
internet, no `efetch`, thousands of taxids in a few seconds. Returns the
minimal table, with one row per distinct taxid, to be joined to the
BLAST hits by taxid.

`x` can be a vector of taxids or any result of
[`parallel_blast()`](https://heronoh.github.io/BLASTr/reference/parallel_blast.md)
(long or wide) or
[`parallel_blast_hsps()`](https://heronoh.github.io/BLASTr/reference/parallel_blast_hsps.md):
the taxids are taken from the `staxid`, `k_staxid` and
`Taxids at best identity` columns. Unknown taxids (`0`, empty, `N/A`)
are dropped.

When the database does not exist, it is built from the NCBI taxdump
(downloaded once, then reused). If `x` is a BLASTr result, the taxdump
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
get_tax(x, taxonomy_db = NULL, ranks = tax_ranks_default, taxdump = NULL)
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

## Examples

``` r
if (FALSE) { # \dontrun{
options(BLASTr.taxonomy_db = "/data/databases/tax/nameNode.sqlite")
get_tax(c("9606", "10090"))

res <- parallel_blast(seqs, db_path = "/data/databases/core_nt/core_nt")
tax <- get_tax(res)
dplyr::left_join(res, tax, by = c(staxid = "taxid"))
} # }
```

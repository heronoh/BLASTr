# BLAST databases of the NCBI

`list_ncbi_blast_dbs()` lists the preformatted BLAST databases the NCBI
publishes (<https://ftp.ncbi.nlm.nih.gov/blast/db/>), with their size
and date. `download_blast_db()` downloads one of them (e.g. `"core_nt"`,
`"16S_ribosomal_RNA"`, `"ITS_eukaryote_sequences"`), checks the md5 of
each file, unpacks it and returns the path to give to
[`parallel_blast()`](https://heronoh.github.io/BLASTr/reference/parallel_blast.md).

The database goes to `<dir>/<db>_<date>/` (the NCBI date of the
database), so a newer release does not overwrite the one older results
were searched against. An interrupted download is resumed: files already
unpacked (with a matching md5) are not downloaded again. With
`taxonomy = TRUE`, the taxonomy database of
[`get_tax()`](https://heronoh.github.io/BLASTr/reference/get_tax.md)
matching the date of the BLAST database is built too, so the
identifications can be made offline.

`core_nt` has about 300 GB unpacked (and is downloaded in volumes of a
few GB, unpacked one at a time); see `list_ncbi_blast_dbs()` for the
sizes.

## Usage

``` r
list_ncbi_blast_dbs()

download_blast_db(
  db = "core_nt",
  dir,
  taxdb = TRUE,
  taxonomy = TRUE,
  keep_archives = FALSE,
  check_space = TRUE
)
```

## Arguments

- db:

  Name of the database, as in `list_ncbi_blast_dbs()$db`.

- dir:

  Directory where the database directory is created.

- taxdb:

  Also download `taxdb` (scientific names of the taxids, used by BLAST
  in `sscinames`) into the same directory.

- taxonomy:

  Build the local taxonomy database of
  [`get_tax()`](https://heronoh.github.io/BLASTr/reference/get_tax.md)
  for the date of this database (kept in the BLASTr cache).

- keep_archives:

  Keep the `.tar.gz` files after unpacking.

- check_space:

  Stop when the disk does not have room for the unpacked database plus
  one archive.

## Value

`list_ncbi_blast_dbs()`: a tibble with `db`, `type` (`"nucl"`/`"prot"`),
`description`, `n_seqs`, `n_letters`, `n_volumes`, `size_gb` (unpacked),
`download_gb` and `last_updated`. `download_blast_db()`: the path of the
database (directory and name, as `db_path` of
[`parallel_blast()`](https://heronoh.github.io/BLASTr/reference/parallel_blast.md)),
invisibly.

## Examples

``` r
if (FALSE) { # \dontrun{
list_ncbi_blast_dbs()
db <- download_blast_db("16S_ribosomal_RNA", dir = "/data/databases")
parallel_blast(seqs, db_path = db, total_cores = 8)
} # }
```

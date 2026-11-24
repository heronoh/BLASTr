# Build a BLAST database from a FASTA file

Formats a FASTA file as a BLAST database with `makeblastdb`, optionally
with the taxid of each sequence (`taxid_map`), so that the searches
report `staxid`.

## Usage

``` r
make_blast_db(
  fasta_path,
  db_path = NULL,
  db_type = "nucl",
  taxid_map = NULL,
  parse_seqids = TRUE,
  verbose = c("silent", "cmd", "output", "full"),
  env_name = "blastr-blast-env"
)
```

## Arguments

- fasta_path:

  Path to the input FASTA file.

- db_path:

  Path (prefix) of the output BLAST database. Defaults to `fasta_path`,
  i.e. the database files are written next to the FASTA file
  (`<fasta_path>.nsq`, `<fasta_path>.nin`, ...). The function stops with
  an error if the directory of `db_path` does not exist or cannot be
  written to.

- db_type:

  Type of database to create, either "nucl" or "prot".

- taxid_map:

  Optional path to a file mapping sequence IDs to taxonomy IDs (two
  columns, separated by a space or a tab), passed to
  `makeblastdb -taxid_map`.

- parse_seqids:

  Whether to parse sequence IDs.

- verbose:

  Should `[condathis::run()]` internal command be shown? With `"silent"`
  (default) nothing is printed.

- env_name:

  The name of the conda environment with the parameter (i.e.
  "blastr-blast-env").

## Value

The result of the `makeblastdb` command.

## See also

Other database functions:
[`ncbi_blast_dbs`](https://heronoh.github.io/BLASTr/reference/ncbi_blast_dbs.md),
[`reference_dbs`](https://heronoh.github.io/BLASTr/reference/reference_dbs.md)

## Examples

``` r
# needs BLAST+ (see install_dependencies())
if (has_blast()) {
  fasta <- system.file("extdata", "minimal_db_blast.fasta", package = "BLASTr")
  make_blast_db(fasta, db_path = file.path(tempdir(), "minimal_db_blast"))

  # without `db_path`, the database is written next to the FASTA file
  fasta_copy <- file.path(tempdir(), "my_db.fasta")
  file.copy(fasta, fasta_copy, overwrite = TRUE)
  make_blast_db(fasta_copy)
  list.files(tempdir(), pattern = "^my_db")
}
```

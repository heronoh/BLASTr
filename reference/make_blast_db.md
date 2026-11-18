# Make BLAST Database

Create a BLAST database from a FASTA file.

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

  Optional path to a file mapping sequence IDs to taxonomy IDs. Must
  have the same radical name than the .fasta db file

- parse_seqids:

  Whether to parse sequence IDs.

- verbose:

  Should `[condathis::run()]` internal command be shown?

- env_name:

  The name of the conda environment with the parameter (i.e.
  "blastr-blast-env").

## Value

The result of the `makeblastdb` command.

## Examples

``` r
if (FALSE) { # \dontrun{
make_blast_db(
  fasta_path = fs::path_package(
    "BLASTr", "extdata", "minimal_db_blast",
    ext = "fasta"
  ),
  db_path = fs::path_temp("minimal_db_blast"),
  db_type = "nucl"
)

# Without `db_path`, the database is written next to the FASTA file:
fasta_copy <- fs::path_temp("my_db", ext = "fasta")
fs::file_copy(
  fs::path_package("BLASTr", "extdata", "minimal_db_blast", ext = "fasta"),
  fasta_copy
)
make_blast_db(fasta_path = fasta_copy)
} # }
```

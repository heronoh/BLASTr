# Header of a subject in a BLAST database

Retrieve complete sequence header from BLAST DB based on SubjectIDs

## Usage

``` r
get_fasta_header(
  id,
  db_path,
  env_name = "blastr-blast-env",
  verbose = c("silent", "cmd", "output", "full")
)
```

## Arguments

- id:

  SubjectID from BLAST results or any NCBI Nucleotide identifier.

- db_path:

  Complete path do formatted BLAST database.

- env_name:

  The name of the conda environment used to run command-line tools (i.e.
  "blastr-blast-env").

- verbose:

  Should condathis::run() internal command be shown?

## Value

Complete identifier for the SubjectID as is on the database.

## Examples

``` r
# needs BLAST+ (see install_dependencies())
if (has_blast()) {
  fasta <- system.file("extdata", "minimal_db_blast.fasta", package = "BLASTr")
  db <- file.path(tempdir(), "minimal_db_blast")
  make_blast_db(fasta, db_path = db)
  seqs <- readLines(system.file("extdata", "seqs_test.txt", package = "BLASTr"))
  get_fasta_header(id = "AP011979.1", db_path = db)
}
```

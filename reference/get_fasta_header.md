# Get header from SubjectID

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
if (FALSE) { # \dontrun{
dna_fasta_path <- fs::path_package(
  "BLASTr", "extdata", "minimal_db_blast",
  ext = "fasta"
)
temp_db_path <- fs::path_temp("minimal_db_blast")
make_blast_db(fasta_path = dna_fasta_path, db_path = temp_db_path)
get_fasta_header(id = "AP011979.1", db_path = temp_db_path)
} # }
```

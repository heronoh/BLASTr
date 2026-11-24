# Search one sequence with BLAST, as a table

Retrieve BLAST results from a given query sequence

## Usage

``` r
get_blast_results(
  query_seqs,
  db_path,
  num_threads = 1L,
  blast_type = "blastn",
  perc_id = 80L,
  perc_qcov_hsp = 80L,
  num_alignments = 4L,
  verbose = c("silent", "cmd", "output", "full"),
  env_name = "blastr-blast-env"
)
```

## Arguments

- query_seqs:

  Vector of sequences to be BLASTed.

- db_path:

  Complete path do formatted BLAST database.

- num_threads:

  Number of threads to run BLAST on. Passed on to BLAST+ argument
  `-num_threads`.

- blast_type:

  One of the available BLAST+ search engines, \#' one of:
  `c("blastn", "blastp", "blastx", "tblastn", "tblastx")`.

- perc_id:

  Lowest identity percentage cutoff. Passed on to BLAST+
  `-perc_identity`.

- perc_qcov_hsp:

  Lowest query coverage per HSP percentage cutoff. Passed on to BLAST+
  `-qcov_hsp_perc`.

- num_alignments:

  Number of alignments to retrieve from BLAST. Max = 6.

- verbose:

  Should condathis::run() internal command be shown?

- env_name:

  The name of the conda environment used to run command-line tools (i.e.
  "blastr-blast-env").

## Value

A `tibble` with the results of BLASTn for each sequence.

## Examples

``` r
# needs BLAST+ (see install_dependencies())
if (has_blast()) {
  fasta <- system.file("extdata", "minimal_db_blast.fasta", package = "BLASTr")
  db <- file.path(tempdir(), "minimal_db_blast")
  make_blast_db(fasta, db_path = db)
  seqs <- readLines(system.file("extdata", "seqs_test.txt", package = "BLASTr"))
  get_blast_results(query_seqs = seqs[1], db_path = db)
}
```

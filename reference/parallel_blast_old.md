# Run parallel BLAST for a set of sequences.

**\[deprecated\]**

This functions was depreacated in version 0.1.7. The recommended
function to use is
[`parallel_blast()`](https://heronoh.github.io/BLASTr/reference/parallel_blast.md).

The new implementation is more resilient to internal errors in BLAST
execution.

## Usage

``` r
parallel_blast_old(
  asvs,
  db_path,
  ...,
  out_file = NULL,
  out_RDS = NULL,
  total_cores = 1L,
  num_threads = 1L,
  blast_type = "blastn",
  perc_id = 80L,
  perc_qcov_hsp = 80L,
  num_alignments = 4L,
  verbose = "silent",
  env_name = "blastr-blast-env",
  query_seqs = NULL
)
```

## Arguments

- asvs:

  Character vector with sequences

- db_path:

  Complete path do formatted BLAST database.

- ...:

  These dots are for future extensions and must be empty.

- out_file:

  Complete path to output .csv file on an existing folder.

- out_RDS:

  Complete path to output RDS file on an existing folder.

- total_cores:

  Total available cores to run BLAST in parallel. Check your max with
  *future::availableCores()*.

- num_threads:

  Number of threads to run BLAST on. Passed on to BLAST+ argument
  `-num_threads`.

- blast_type:

  BLAST+ executable to be used on search.

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

- query_seqs:

  the same as `asvs` (to provide API compatibility with
  [`parallel_blast()`](https://heronoh.github.io/BLASTr/reference/parallel_blast.md)).

## Value

A tibble with the BLAST tabular output.

## Examples

``` r
if (FALSE) { # \dontrun{
dna_fasta_path <- fs::path_package(
  "BLASTr", "extdata", "minimal_db_blast",
  ext = "fasta"
)
temp_db_path <- fs::path_temp("minimal_db_blast")
make_blast_db(fasta_path = dna_fasta_path, db_path = temp_db_path)
asvs_string <- c(
  "CTAGCCATAAACTTAAATGAAGCTATACTAA",
  "ACTCGTTCGCCAGAGTACTACAAGCGAAAG"
)
blast_res <- parallel_blast_old(
  asvs = asvs_string, # vector of sequences to be searched
  db_path = temp_db_path, # path to a formatted blast database
  out_file = NULL, # path to a .csv file to be created with results (on an existing folder)
  out_RDS = NULL, # path to a .RDS file to be created with results (on an existing folder)
  perc_id = 80, # minimum identity percentage cutoff
  total_cores = 1, # Number of BLAST process to start in parallel
  perc_qcov_hsp = 80, # minimum percentage coverage of query sequence by subject sequence cutoff
  num_threads = 1, # number of threads/cores to run each blast on
  # maximum number of alignments/matches to retrieve results for each query sequence
  num_alignments = 3,
  blast_type = "blastn" # blast search engine to use
)
} # }
```

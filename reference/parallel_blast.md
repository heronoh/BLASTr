# Run Parallelized BLAST

Run parallel BLAST for a set of sequences

## Usage

``` r
parallel_blast(
  query_seqs,
  db_path,
  ...,
  total_cores = 1L,
  num_threads = 1L,
  blast_type = "blastn",
  perc_id = 80L,
  perc_qcov_hsp = 80L,
  num_alignments = 4L,
  retry_times = 3L,
  mt_mode = c("2", "1", "0"),
  verbose = c("progress", "silent", "cmd", "output", "full"),
  env_name = "blastr-blast-env",
  asvs = deprecated(),
  out_file = deprecated(),
  out_RDS = deprecated()
)
```

## Arguments

- query_seqs:

  Character vector with sequences.

- db_path:

  Path to the formatted BLAST database.

- ...:

  These dots are for future extensions and must be empty.

- total_cores:

  Number of parallel BLAST processes to run.

- num_threads:

  Number of threads/cores to run each BLAST process on.

- blast_type:

  BLAST+ executable to be used on search.

- perc_id:

  Lowest identity percentage cutoff.

- perc_qcov_hsp:

  Lowest query coverage per HSP percentage cutoff.

- num_alignments:

  Number of alignments to retrieve results for each query sequence.
  Defaults to `4`.

- retry_times:

  Number of times to retry failed BLAST jobs. Defaults to `3` attempts.

- mt_mode:

  Multithreading mode to be used by BLAST+. One of: `c("2", "1", "0")`.
  See BLAST+ manual for details.

- verbose:

  Strategy used for showing outputting internal commands. Defaults to
  "progress".

- env_name:

  The name of the conda environment used to run command-line tools.
  Defaults to `"blastr-blast-env"`.

- asvs:

  Deprecated. Same as `query_seqs`. Use `query_seqs` instead.

- out_file:

  Deprecated. Path to output `.csv` file on an existing directory.

- out_RDS:

  Deprecated. Path to output `RDS` file on an existing directory.

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
blast_res <- parallel_blast(
  query_seqs = asvs_string,
  db_path = temp_db_path
)
blast_res
} # }
```

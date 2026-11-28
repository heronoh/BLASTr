# Log of BLAST searches

When a log directory is set, each search made by
[`parallel_blast()`](https://heronoh.github.io/BLASTr/reference/parallel_blast.md)
adds one line to `BLASTr_logs/blastr_runs.tsv` inside it (the folder and
file are created when needed). The directory is taken, in this order,
from the `log_dir` argument of
[`parallel_blast()`](https://heronoh.github.io/BLASTr/reference/parallel_blast.md),
the option `BLASTr.log_dir` and the environment variable
`BLASTR_LOG_DIR`; with none of them, nothing is logged. In background
jobs, the line is written when the results are collected with
[`blast_job_result()`](https://heronoh.github.io/BLASTr/reference/blast_job.md).

Each line has:

- the run: `run_id`, `start`, `end`, `host`, `user`, `mode`
  (foreground/background), `function`
  ([`parallel_blast()`](https://heronoh.github.io/BLASTr/reference/parallel_blast.md)
  or
  [`parallel_blast_hsps()`](https://heronoh.github.io/BLASTr/reference/parallel_blast_hsps.md)),
  `blastr_version`, `blast_version`, `job_dir`;

- the sequences: `n_input` (as given), `n_seqs` (unique ones searched),
  `n_invalid` and their length distribution (`len_min`, `len_q25`,
  `len_median`, `len_mean`, `len_q75`, `len_max`);

- the database(s): `db_path`, `n_dbs`, `db_gb` (sum), `db_volumes`,
  `ram_gb`, `fits_ram`;

- the strategy and settings: `total_cores`, `processes`, `threads`,
  `mt_mode`, `batch_size`, `n_batches`, `task`, `num_alignments`,
  `search_alignments`, `perc_id`, `perc_qcov_hsp`, `check_n`,
  `check_tolerance`;

- the cost: `wall_s`, `cpu_s` (CPU time of the BLAST processes),
  `cpu_h`, `cpu_efficiency` (`cpu_s` over `wall_s` times the threads in
  use), `n_batches_failed`, `n_seqs_failed`;

- the identifications: `n_with_hits`, `pct_with_hits`, `pct_best_id`,
  `pct_check_id` (both over all sequences), `id1_median` (median
  identity of the first hit) and `pct_id1_ge_99`, `pct_id1_ge_97`
  (sequences whose first hit has identity of at least 99% / 97%). For
  [`parallel_blast_hsps()`](https://heronoh.github.io/BLASTr/reference/parallel_blast_hsps.md),
  `pct_best_id` and `pct_check_id` are `NA` and the first hit is the
  first HSP of the best subject.

When a newer BLASTr adds columns, the file is rewritten with all of them
(missing values as `NA`).

## Usage

``` r
blast_log_path(log_dir = NULL)

read_blast_log(log_dir = NULL)
```

## Arguments

- log_dir:

  Directory where the `BLASTr_logs` folder is kept. `NULL` uses the
  option and the environment variable described above.

## Value

`blast_log_path()` returns the path of the log file, or `NULL` when no
log directory is set. `read_blast_log()` returns the log as a tibble.

## See also

Other search functions:
[`blast_job()`](https://heronoh.github.io/BLASTr/reference/blast_job.md),
[`blast_strategy()`](https://heronoh.github.io/BLASTr/reference/blast_strategy.md),
[`exit_codes()`](https://heronoh.github.io/BLASTr/reference/exit_codes.md),
[`parallel_blast()`](https://heronoh.github.io/BLASTr/reference/parallel_blast.md),
[`parallel_blast_hsps()`](https://heronoh.github.io/BLASTr/reference/parallel_blast_hsps.md)

## Examples

``` r
log_dir <- tempfile("blastr_logs_")
blast_log_path(log_dir)
#> /tmp/RtmpROkON7/blastr_logs_1970333a1633/BLASTr_logs/blastr_runs.tsv

# needs BLAST+ (see install_dependencies())
if (has_blast()) {
  fasta <- system.file("extdata", "minimal_db_blast.fasta", package = "BLASTr")
  db <- file.path(tempdir(), "minimal_db_blast")
  make_blast_db(fasta, db_path = db)
  seqs <- readLines(system.file("extdata", "seqs_test.txt", package = "BLASTr"))
  res <- parallel_blast(seqs, db, log_dir = log_dir)
  read_blast_log(log_dir)
}
```

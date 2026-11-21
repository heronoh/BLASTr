# Run Parallelized BLAST

Search a set of sequences against one or more BLAST databases, using the
cores and RAM given. The search strategy is chosen automatically (see
[`blast_strategy()`](https://heronoh.github.io/BLASTr/reference/blast_strategy.md)):
the sequences are grouped in batches, and each `blastn` process searches
a whole batch, so the database is read once per batch instead of once
per sequence.

The progress bar shows the estimated time left (ETA), computed from the
batches already finished. With `background = TRUE` the search runs in a
separate process and the function returns at once; see
[`blast_job()`](https://heronoh.github.io/BLASTr/reference/blast_job.md).

## Usage

``` r
parallel_blast(
  query_seqs,
  db_path,
  ...,
  total_cores = 1L,
  ram_gb = NULL,
  num_alignments = 10L,
  search_alignments = 50L,
  perc_id = 80L,
  perc_qcov_hsp = 80L,
  search_qcov_hsp = 50L,
  check_n = 10L,
  check_tolerance = 0.5,
  rank_by = c("similarity", "identity", "bitscore"),
  output = c("long", "wide"),
  num_threads = NULL,
  batch_size = NULL,
  mt_mode = NULL,
  task = NULL,
  background = FALSE,
  job_dir = NULL,
  log_dir = NULL,
  retry_times = 3L,
  blast_type = "blastn",
  verbose = c("progress", "silent", "cmd", "output", "full"),
  env_name = "blastr-blast-env",
  out_file = deprecated(),
  out_RDS = deprecated()
)
```

## Arguments

- query_seqs:

  Character vector with sequences.

- db_path:

  Path(s) to the formatted BLAST database(s). Several databases (a
  vector, or one string separated by spaces) are searched together, the
  smallest first.

- ...:

  These dots are for future extensions and must be empty.

- total_cores:

  Number of CPU cores (threads) to use in total.

- ram_gb:

  RAM available for the search, in GB. Used to decide whether the
  database fits in memory. `NULL` uses the total RAM of the machine.

- num_alignments:

  Number of hits to keep for each sequence. Defaults to `10`.

- search_alignments:

  Number of subjects BLAST keeps during the search (`-max_target_seqs`).
  BLAST cuts its list of subjects before the `perc_qcov_hsp` filter, so
  asking only for `num_alignments` can lose valid hits. BLASTr searches
  `max(num_alignments, search_alignments)` subjects and keeps the
  `num_alignments` best ones that pass the filters. In the benchmark, 50
  instead of 10 added about 5% to the time.

- perc_id:

  Lowest identity percentage cutoff.

- perc_qcov_hsp:

  Query coverage per HSP (percent) expected of a good hit. Hits below it
  are kept and flagged (`low_qcov`), not removed: removing them inside
  BLAST hides a closer reference that covers a bit less of the query
  (common with partial reference sequences) and promotes a worse one.
  The coverage still weighs in the rank through the Similarity score.
  Defaults to `80`.

- search_qcov_hsp:

  Lowest query coverage per HSP (percent) used by BLAST
  (`-qcov_hsp_perc`); hits below it are not returned. Defaults to `50`.
  Set it to `perc_qcov_hsp` to remove the hits below `perc_qcov_hsp`, as
  in earlier versions.

- check_n:

  Number of best hits compared by `Check ID` (see
  [`check_hits()`](https://heronoh.github.io/BLASTr/reference/check_hits.md)).
  Defaults to `10`.

- check_tolerance:

  Difference from the best hit, in the units of `rank_by`, under which
  hits are tied ("at best"). Defaults to `0.5`.

- rank_by:

  Criterion to rank the hits of each sequence: `"similarity"` (default,
  the Similarity score, see
  [`similarity_score()`](https://heronoh.github.io/BLASTr/reference/similarity_score.md)),
  `"identity"` or `"bitscore"`. Ties are broken by identity, coverage,
  bitscore and BLAST order. BLAST orders its hits by e-value/bitscore,
  which favours longer alignments, so its first hit is not always the
  most similar one.

- output:

  `"long"` (default): one row per hit. `"wide"`: one row per sequence,
  with the summary of
  [`blast_summary()`](https://heronoh.github.io/BLASTr/reference/blast_summary.md)
  and the first `num_alignments` hits side by side (`1_` is the best
  one).

- num_threads, batch_size, mt_mode:

  Override the automatic strategy: threads per `blastn` process,
  sequences per batch and BLAST+ multithreading mode (`0`, `1` or `2`).
  `NULL` (default) chooses automatically.

- task:

  BLAST task (`-task`), e.g. `"megablast"` (the `blastn` default),
  `"dc-megablast"` or `"blastn"`. `NULL` uses the BLAST default.

- background:

  If `TRUE`, start the search in the background and return a job (see
  [`blast_job()`](https://heronoh.github.io/BLASTr/reference/blast_job.md))
  instead of the results.

- job_dir:

  Directory for the job files. Defaults to a temporary directory
  (removed at the end) in the foreground, and to a new directory in the
  BLASTr cache in the background.

- log_dir:

  Directory where a `BLASTr_logs` folder keeps a log of the searches
  (number and lengths of the sequences, database size, strategy, time,
  CPU time and percentage of identifications). `NULL` uses the option
  `BLASTr.log_dir` or the environment variable `BLASTR_LOG_DIR`; with
  none set, nothing is logged. See
  [`blast_log_path()`](https://heronoh.github.io/BLASTr/reference/blast_log_path.md).

- retry_times:

  Number of times to retry failed BLAST jobs. Defaults to `3` attempts.

- blast_type:

  BLAST+ executable to be used on search.

- verbose:

  `"progress"` (default) shows the progress and ETA; `"silent"` shows
  nothing. `"cmd"`, `"output"` and `"full"` also show the chosen
  strategy.

- env_name:

  The name of the conda environment used to run command-line tools.
  Defaults to `"blastr-blast-env"`.

- out_file:

  Deprecated. Path to output `.csv` file on an existing directory.

- out_RDS:

  Deprecated. Path to output `RDS` file on an existing directory.

## Value

With `output = "long"`, a tibble with one row per hit, best first:
`Sequence`, `rank` (1 = best; tied hits share the rank), `tie` (the rank
is shared), `at_best` (within `check_tolerance` of the best hit),
`low_qcov` (`qcovhsp` below `perc_qcov_hsp`), `blast_rank` (order given
by BLAST), `subject`, `subject header`, `staxid`, `identity`,
`Similarity score`, `qcovhsp`, `nident`, `length`, `mismatches`, `gaps`,
`query start`, `query end`, `subject start`, `subject end`,
`subject strand`, `e-value` and `bitscore`. The first `num_alignments`
hits are kept, plus all the hits at best. A sequence without hits has
one row with `NA` hit fields. The summary per sequence (`Check ID`, ...)
comes from
[`blast_summary()`](https://heronoh.github.io/BLASTr/reference/blast_summary.md).

With `output = "wide"`, one row per sequence: `Sequence`, the summary of
[`blast_summary()`](https://heronoh.github.io/BLASTr/reference/blast_summary.md)
and, for each hit `k` from 1 to `num_alignments` in rank order, the hit
fields prefixed by `k_`.

The attribute `BLASTr_metadata` keeps the parameters, the strategy, the
BLAST command (`blast_cmd`, without the sequences), the BLAST version,
the date of each database and the exit codes (see
[`exit_codes()`](https://heronoh.github.io/BLASTr/reference/exit_codes.md)).
With `background = TRUE`, a `blastr_job`.

## Examples

``` r
if (FALSE) { # \dontrun{
dna_fasta_path <- fs::path_package(
  "BLASTr", "extdata", "minimal_db_blast",
  ext = "fasta"
)
temp_db_path <- fs::path_temp("minimal_db_blast")
make_blast_db(fasta_path = dna_fasta_path, db_path = temp_db_path)
query_seqs_string <- c(
  "CTAGCCATAAACTTAAATGAAGCTATACTAA",
  "ACTCGTTCGCCAGAGTACTACAAGCGAAAG"
)
blast_res <- parallel_blast(
  query_seqs = query_seqs_string,
  db_path = temp_db_path
)
blast_res

# in the background
job <- parallel_blast(query_seqs_string, temp_db_path, background = TRUE)
blast_job_status(job)
blast_res <- blast_job_result(job)
} # }
```

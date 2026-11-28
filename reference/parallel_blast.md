# Search sequences with BLAST in parallel

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
  search_alignments = 100L,
  perc_id = 80L,
  perc_qcov_hsp = 80L,
  search_qcov_hsp = 50L,
  check_n = 10L,
  check_tolerance = 0.5,
  rank_by = c("similarity", "identity", "bitscore"),
  uninformative = uninformative_references(),
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

  Number of subjects BLAST keeps during the search (`-num_alignments`).
  BLAST cuts its list of subjects before the `perc_qcov_hsp` filter, so
  asking only for `num_alignments` can lose valid hits. BLASTr searches
  `max(num_alignments, search_alignments)` subjects and keeps the
  `num_alignments` best ones that pass the filters. Subjects tied with
  the best one beyond this limit are lost (`At best truncated`, see
  [`check_hits()`](https://heronoh.github.io/BLASTr/reference/check_hits.md)).
  Defaults to `100`: in the benchmark on `core_nt`, 50 instead of 10
  added about 5% to the time, and 100 instead of 50 from 3% to 18%,
  while halving the sequences with truncated ties (from 170 to 94 in
  1260 GenBank references) and the ones whose ancestor differs from a
  search with 500 subjects (from 37 to 18).

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
  most similar one. In a calibration on plant barcodes of known genus
  (ITS2, rbcL) and on 1257 GenBank references of six markers searched
  without themselves, ranking by Similarity score or identity picked the
  right taxon significantly more often than bitscore (genus: 72% and 76%
  against 62%), with no significant difference between them.

- uninformative:

  Table of the references that do not inform the identification
  (uncultured, environmental, metagenomes, vectors, ...; see
  [`uninformative_references()`](https://heronoh.github.io/BLASTr/reference/uninformative_references.md)):
  the table shipped with BLASTr (default), a path to another one, a data
  frame, or `NULL` to treat all hits alike (`informative` NA, not
  evaluated). The hits are flagged (`informative`,
  `uninformative_category`) and ranked among the informative ones, so
  `rank`, `at_best`, `Check ID` and
  [`blast_lca()`](https://heronoh.github.io/BLASTr/reference/blast_lca.md)
  use only those; a sequence whose hits are all uninformative keeps its
  best ones, ranked among themselves. Uninformative hits take subjects
  of the search, so a sequence of a poorly known group needs a larger
  `search_alignments`.

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
`low_qcov` (`qcovhsp` below `perc_qcov_hsp`), `informative` and
`uninformative_category` (see `uninformative`), `blast_rank` (order
given by BLAST), `subject`, `subject header`, `staxid`, `identity`,
`Similarity score`, `qcovhsp`, `nident`, `length`, `mismatches`, `gaps`,
`query start`, `query end`, `subject start`, `subject end`,
`subject strand`, `e-value`, `bitscore`, `at_best_refs` and `btop`. The
first `num_alignments` hits are kept (informative ones, when there are),
plus all the hits at best (by either rule), the hits with the taxid of a
closest reference (see below) and the uninformative hits among the first
`num_alignments` of all (after the ranked ones, with `rank` `NA`). A
sequence without hits has one row with `NA` hit fields. The summary per
sequence (`Check ID`, ...) comes from
[`blast_summary()`](https://heronoh.github.io/BLASTr/reference/blast_summary.md).

`at_best_refs` marks the hits tied by a second rule, by references, kept
beside `at_best` to compare the two (it does not change `rank`,
`at_best` or the summary;
[`blast_lca()`](https://heronoh.github.io/BLASTr/reference/blast_lca.md)
uses it with `ties = "at_best_refs"`). After the rule of OBITools'
`obitag`: the closest references to the query are found by distance
(differences in the alignment plus query bases outside it, over the
query length), and the references as close to one of them as the query
is are tied with it, comparing references with each other on the query
region both cover (from `btop`, BLAST's summary of the alignment). The
window is at least 0.5% of the query, even when the query is identical
to its closest reference: in markers whose species differ by few bases
(e.g. long rbcL amplicons), an identical reference of a neighbour
species is not taken for the query's own. The price: a query identical
to a reference of its own species, with another species within 0.5%,
goes up to the genus. A reference that shares less than half of the
query with the closest one is tied. Unlike `at_best`, the window does
not grow with the length of the sequence: 0.5 point of Similarity score
is under one base in 100 bp and four in 800 bp. `NA` for results without
`btop` (searches of earlier versions).

With `output = "wide"`, one row per sequence: `Sequence`, the summary of
[`blast_summary()`](https://heronoh.github.io/BLASTr/reference/blast_summary.md)
and, for each hit `k` from 1 to `num_alignments` in rank order, the hit
fields prefixed by `k_`.

The attribute `BLASTr_metadata` keeps the parameters, the strategy, the
BLAST command (`blast_cmd`, without the sequences), the BLAST version,
the date of each database, the exit codes (see
[`exit_codes()`](https://heronoh.github.io/BLASTr/reference/exit_codes.md))
and, per sequence, the counts over all the subjects BLAST returned
(`hit_counts`, used by
[`blast_summary()`](https://heronoh.github.io/BLASTr/reference/blast_summary.md)).
With `background = TRUE`, a `blastr_job`.

## See also

Other search functions:
[`blast_job()`](https://heronoh.github.io/BLASTr/reference/blast_job.md),
[`blast_log_path()`](https://heronoh.github.io/BLASTr/reference/blast_log_path.md),
[`blast_strategy()`](https://heronoh.github.io/BLASTr/reference/blast_strategy.md),
[`exit_codes()`](https://heronoh.github.io/BLASTr/reference/exit_codes.md),
[`parallel_blast_hsps()`](https://heronoh.github.io/BLASTr/reference/parallel_blast_hsps.md)

## Examples

``` r
# needs BLAST+ (see install_dependencies())
if (has_blast()) {
  fasta <- system.file("extdata", "minimal_db_blast.fasta", package = "BLASTr")
  db <- file.path(tempdir(), "minimal_db_blast")
  make_blast_db(fasta, db_path = db)
  seqs <- readLines(system.file("extdata", "seqs_test.txt", package = "BLASTr"))
  res <- parallel_blast(seqs, db)
  # one row per hit, best first; one row with NA hit fields for a sequence without hits
  res[, c("Sequence", "rank", "at_best", "subject", "identity", "qcovhsp", "Similarity score")]
  # one row per sequence
  blast_summary(res)
}
```

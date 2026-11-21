# Changelog

## BLASTr (development version)

### Breaking changes

- The term “ASV” was dropped from the package API and documentation,
  since it denotes a sequence type specific to some pipelines. Sequences
  are now referred to only as `seqs`.
- [`parallel_blast_old()`](https://heronoh.github.io/BLASTr/reference/parallel_blast_old.md)
  lost its `asvs` argument; use `query_seqs`. This also fixes
  [`parallel_blast_old()`](https://heronoh.github.io/BLASTr/reference/parallel_blast_old.md),
  which still called
  [`get_blast_results()`](https://heronoh.github.io/BLASTr/reference/get_blast_results.md)
  with the removed `asv` argument.
- Example data `inst/extdata/asvs_test.txt` was renamed to
  `seqs_test.txt`.
- [`parallel_blast()`](https://heronoh.github.io/BLASTr/reference/parallel_blast.md)
  was rewritten. Instead of one `blastn` per sequence, the sequences are
  grouped in batches and each `blastn` searches a whole batch, so the
  database is read once per batch (in the benchmark on `core_nt`, 100
  sequences took 79 s instead of 944 s, with the same hits):
  - `total_cores` is now the total number of threads to use. The number
    of processes, threads per process, batch size and `mt_mode` are
    chosen automatically from `total_cores`, the new `ram_gb`, the size
    of the database(s) and the number of sequences (see the new
    [`blast_strategy()`](https://heronoh.github.io/BLASTr/reference/blast_strategy.md)).
    `num_threads`, `batch_size` and `mt_mode` (now `NULL` by default)
    only override the automatic choice.
  - `num_alignments` now defaults to `10` (was `4`). With few
    alignments, BLAST cuts the hit list before the `perc_qcov_hsp`
    filter and valid hits can be lost.
  - Hit column `k_indentity` was renamed to `k_identity` (typo).
  - Subjects are reported by accession (`saccver`), as before.

### New features

- [`parallel_blast()`](https://heronoh.github.io/BLASTr/reference/parallel_blast.md)
  results have a hit summary after `Sequence`: `Check ID` (“best ID” /
  “check ID”, agreement of the taxids of the best `check_n` hits, by
  default 10, comparing only hits within `check_tolerance` identity
  points of the best one), `N hits`, `N taxids`,
  `N hits same taxid as 1st`, `N hits at best identity` (1 = a single
  reference at the best identity) and `N taxids at best identity` (see
  [`check_hits()`](https://heronoh.github.io/BLASTr/reference/check_hits.md)).
- New `k_Similarity score` column for each hit,
  `(10 * identity + qcovhsp) / 11` (formerly “BLASTn pseudo-score”; see
  [`similarity_score()`](https://heronoh.github.io/BLASTr/reference/similarity_score.md)).
- Progress bar with the estimated time left (ETA), computed from the
  finished batches. In non-interactive sessions, one line per batch.
- `parallel_blast(background = TRUE)` starts the search in a separate
  process and returns a job at once;
  [`blast_job_status()`](https://heronoh.github.io/BLASTr/reference/blast_job.md),
  [`blast_job_result()`](https://heronoh.github.io/BLASTr/reference/blast_job.md),
  [`blast_job_cancel()`](https://heronoh.github.io/BLASTr/reference/blast_job.md)
  and
  [`blast_job_resume()`](https://heronoh.github.io/BLASTr/reference/blast_job.md)
  follow it. The job lives in a directory (`job_dir`), so it survives
  the end of the R session and can be opened again with
  [`blast_job()`](https://heronoh.github.io/BLASTr/reference/blast_job.md).
- Several databases (a vector or a space-separated string in `db_path`)
  are searched together, smallest first; their sizes are added when
  choosing the strategy.
- A sequence that makes BLAST fail is retried alone, so the others in
  its batch are not affected.
- [`parallel_blast()`](https://heronoh.github.io/BLASTr/reference/parallel_blast.md)
  searches `search_alignments` subjects (default 50) and keeps the
  `num_alignments` best ones that pass the identity and coverage
  filters: BLAST cuts its subject list before the coverage filter, so
  asking only for `num_alignments` could lose valid hits.
- New `task` argument in
  [`parallel_blast()`](https://heronoh.github.io/BLASTr/reference/parallel_blast.md)
  (BLAST `-task`, e.g. `"dc-megablast"`).
- New
  [`parallel_blast_hsps()`](https://heronoh.github.io/BLASTr/reference/parallel_blast_hsps.md):
  the raw hits in long format, one row per HSP, with the BLAST fields
  chosen in `fields` and no HSP limit or coverage filter unless asked.
  Searches a database or, with `subject`, the sequences of a FASTA file
  (BLAST `-subject`). Same batches, strategy, ETA, background jobs and
  log as
  [`parallel_blast()`](https://heronoh.github.io/BLASTr/reference/parallel_blast.md).
- Optional log of the searches: with `log_dir` (or the option
  `BLASTr.log_dir` / the environment variable `BLASTR_LOG_DIR`), each
  search adds a line to `BLASTr_logs/blastr_runs.tsv` with the number
  and lengths of the sequences, database size, strategy, time, CPU time
  and percentage of identifications. See
  [`blast_log_path()`](https://heronoh.github.io/BLASTr/reference/blast_log_path.md)
  and
  [`read_blast_log()`](https://heronoh.github.io/BLASTr/reference/blast_log_path.md).

## BLASTr 0.1.7

Release Date: 2025-11-04

Development Changelog:
[dev](https://github.com/heronoh/BLASTr/compare/v0.1.6...v0.1.7)

### Added

- New `retry_times` argument in
  [`parallel_blast()`](https://heronoh.github.io/BLASTr/reference/parallel_blast.md)
  to specify the number of retries for failed BLAST jobs. Defaults to 10
  times.

- New `verbose = "progress"` in
  [`parallel_blast()`](https://heronoh.github.io/BLASTr/reference/parallel_blast.md)
  to show a progress bar for the BLAST jobs.

- New
  [`get_blastr_cache()`](https://heronoh.github.io/BLASTr/reference/get_blastr_cache.md)
  for printing the directory used as cache for `BLASTr`.

- The previous version of
  [`parallel_blast()`](https://heronoh.github.io/BLASTr/reference/parallel_blast.md)
  is now available as
  [`parallel_blast_old()`](https://heronoh.github.io/BLASTr/reference/parallel_blast_old.md).
  It is kept for backward compatibility, but users are encouraged to
  switch to the new
  [`parallel_blast()`](https://heronoh.github.io/BLASTr/reference/parallel_blast.md)
  implementation.

- New `exit_codes` function to extract exit codes and error message from
  internal BLAST runs.

### Changed

- [`parallel_blast()`](https://heronoh.github.io/BLASTr/reference/parallel_blast.md)
  will not throw errors immediately when a BLAST job fails. Instead, it
  will log the error and continue with the remaining jobs.

- [`parallel_blast()`](https://heronoh.github.io/BLASTr/reference/parallel_blast.md)
  output has a new attributes data frame with columns:

  - `query_seq`: The query sequence that was processed.
  - `exit_code`: The exit code of the BLAST command.
  - `stderr`: The standard error output from the BLAST command.

- [`parallel_blast()`](https://heronoh.github.io/BLASTr/reference/parallel_blast.md)
  will automatically retry failed BLAST jobs up to 3 times (or the value
  of `retry_times`) before giving up.

- Argument `asvs` is renamed to `query_seqs` in
  [`parallel_blast()`](https://heronoh.github.io/BLASTr/reference/parallel_blast.md).

- Arguments `asvs`, `out_file` and `out_RDS` in
  [`parallel_blast()`](https://heronoh.github.io/BLASTr/reference/parallel_blast.md)
  are now deprecated.

## BLASTr 0.1.6

Release Date: 2025-10-28

Development Changelog:
[0.1.6](https://github.com/heronoh/BLASTr/compare/v0.1.5...v0.1.6)

### Fixed

- Fix hard coded `num_threads = 1L` in
  [`parallel_blast()`](https://heronoh.github.io/BLASTr/reference/parallel_blast.md).
  [`parallel_blast()`](https://heronoh.github.io/BLASTr/reference/parallel_blast.md)
  now correctly uses the user provided `num_threads` argument.

- Fix parallel functions not using the the custom compute profile.
  [`purrr::map`](https://purrr.tidyverse.org/reference/map.html) do not
  have an argument for specifying the `mirai` `.compute` argument
  correctly. Now all parallel functions use the “default” mirai compute
  profile `.compute = NULL`.

- Fix
  [`get_tax_by_taxID()`](https://heronoh.github.io/BLASTr/reference/get_tax_by_taxID.md)
  when multiple tax IDs are provided.

- Fix
  [`parse_fasta()`](https://heronoh.github.io/BLASTr/reference/parse_fasta.md)
  concatenate lines for multiline FASTA.

## BLASTr 0.1.5

Release Date: 2025-08-27

Development Changelog:
[0.1.5](https://github.com/heronoh/BLASTr/compare/v0.1.4...v0.1.5)

### Added

- New parallel framework based on native `purrr` functions, using
  `mirai` and `carrier` in the background.

- New
  [`make_blast_db()`](https://heronoh.github.io/BLASTr/reference/make_blast_db.md)
  function to create custom BLAST databases.

### Changed

- Internal `entrez-direct` version bump to `24.0`.

- The internal environments used by BLASTr have been renamed to
  `blastr-blast-env` and `blastr-entrez-env`.

- Remove of the `furrr` based functions in favor of the new `purrr`
  native parallel.

- All [`options()`](https://rdrr.io/r/base/options.html) based
  parameters have been removed in favor of explicit arguments in the
  functions.

## BLASTr 0.1.4

### Changed

- All functions have better defaults and less mandatory arguments.

### Fixed

- Fix progress bars when parallel processing is used.

## BLASTr 0.1.3

### Added

- Refactor handling of internal errors.

- `shell_exec` is not available anymore.

## BLASTr 0.1.2

### Added

- New
  [`install_dependencies()`](https://heronoh.github.io/BLASTr/reference/install_dependencies.md)
  automatically create environments.

- Initial support for `condathis` for managing dependency installation.

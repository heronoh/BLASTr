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
- [`parallel_blast()`](https://heronoh.github.io/BLASTr/reference/parallel_blast.md)
  returns the long format by default: one row per hit, ranked by BLASTr
  (`rank`, `tie`, `at_best`, `blast_rank`), and one row with `NA` hit
  fields for each sequence without hits. The wide format (one row per
  sequence, hits `1_` to `n_`) is still available with
  `output = "wide"`; its `k_` hits now follow BLASTr’s rank (BLAST’s
  order is in `k_blast_rank`). The summary per sequence (`Check ID`, …)
  comes from the new
  [`blast_summary()`](https://heronoh.github.io/BLASTr/reference/blast_summary.md).
- Hits are ranked by `rank_by` (default `"similarity"`, the Similarity
  score; or `"identity"`, `"bitscore"`), no longer by BLAST’s order,
  which favours long alignments. `check_tolerance` is in the units of
  `rank_by`. The default was calibrated on sequences of known taxonomy:
  Similarity score and identity picked the right taxon significantly
  more often than bitscore (scripts in `dev/benchmark/R/10_*` to
  `12_*`).
- `perc_qcov_hsp` (80) no longer removes hits: BLAST searches with the
  new `search_qcov_hsp` (default 50) and the hits below `perc_qcov_hsp`
  are kept and flagged (`low_qcov`; summary `Best low qcov`). Removing
  them inside BLAST hid closer references that cover a bit less of the
  query (common with partial references) and promoted worse ones: in 66
  ITS2 and rbcL sequences, 65 got a hit of another genus. With
  `search_qcov_hsp = perc_qcov_hsp` the search is as before.
- Summary columns `N hits at best identity` and
  `N taxids at best identity` were renamed to `N hits at best` and
  `N taxids at best` (they follow `rank_by`).
- [`get_tax_by_taxID()`](https://heronoh.github.io/BLASTr/reference/get_tax_online.md)
  and
  [`get_tax_by_name()`](https://heronoh.github.io/BLASTr/reference/get_tax_online.md)
  were rewritten: they ask the NCBI online through the E-utilities over
  HTTP (no Entrez Direct), in batches of 200 taxids and within the NCBI
  rate limit (optional API key), and return the same columns as
  [`get_tax()`](https://heronoh.github.io/BLASTr/reference/get_tax.md),
  for when the local taxonomy database is not there yet.
  [`get_tax_by_name()`](https://heronoh.github.io/BLASTr/reference/get_tax_online.md)
  matches scientific names first, then synonyms and common names
  (`matched_by`).
  [`parallel_get_tax()`](https://heronoh.github.io/BLASTr/reference/get_tax_online.md)
  is deprecated (use
  [`get_tax_by_taxID()`](https://heronoh.github.io/BLASTr/reference/get_tax_online.md)).
- [`install_dependencies()`](https://heronoh.github.io/BLASTr/reference/install_dependencies.md)
  no longer installs Entrez Direct (`efetch`): nothing uses it now.
- The internal, unexported `search_primers_on_fq()` was removed; the
  check of which primers are found in the reads is now
  `detect_primers()` in mtbrcdR.

### New features

- A second tie rule, by references, kept beside the current one to
  compare the two on real data (the default does not change).
  [`parallel_blast()`](https://heronoh.github.io/BLASTr/reference/parallel_blast.md)
  asks BLAST for the `btop` of each hit (a short summary of the
  alignment) and its long results gain two columns at the end: `btop`
  and `at_best_refs`, the hits tied by the new rule. After the rule of
  OBITools’ `obitag`: the references as close to the closest reference
  as the query is are tied with it, comparing references with each other
  on the region of the query both cover; the window is at least 0.5% of
  the query, also when the query is identical to its closest reference.
  The rows tied by the new rule, and the hits with the taxid of a
  closest reference, are kept too, after the first `num_alignments`.
  `blast_lca(ties = "at_best_refs")` gives the ancestor of those hits
  and two more columns, `Species flag` and `Flagged species`:
  `"single reference"` (the species of the closest reference has only
  references identical to it) or `"beyond spread"` (the query is farther
  from it than the references of the species are from each other), when
  the query is not identical to it; the flag does not change the
  ancestor. On a benchmark of 47 primer pairs of vertebrates,
  invertebrates and plants, with the species cut of the mtbrcdR
  (Similarity score 98), the new rule called a species wrongly less
  often (by 0.3 to 1.4 points) and rightly less often (by 1.2 to 2.1
  points) than the current one; flagged species calls were wrong 3 to 8
  times as often as the others (species in the database, without the
  cut). Searches of earlier versions have no `btop`: `at_best_refs` is
  `NA` and `blast_lca(ties = "at_best_refs")` stops with an error.

- Reference databases:
  [`list_ncbi_blast_dbs()`](https://heronoh.github.io/BLASTr/reference/ncbi_blast_dbs.md)
  lists the BLAST databases of the NCBI with their sizes and dates;
  [`download_blast_db()`](https://heronoh.github.io/BLASTr/reference/ncbi_blast_dbs.md)
  downloads one (e.g. `"core_nt"`), checking the md5 of each file, into
  a directory named after its date, resumes interrupted downloads, and
  by default also gets `taxdb` and builds the taxonomy database of
  [`get_tax()`](https://heronoh.github.io/BLASTr/reference/get_tax.md)
  for that date.
  [`make_unite_db()`](https://heronoh.github.io/BLASTr/reference/reference_dbs.md)
  (UNITE general FASTA release, all eukaryotes or fungi) and
  [`make_silva_db()`](https://heronoh.github.io/BLASTr/reference/reference_dbs.md)
  (SILVA SSU Ref NR99) download and format them as BLAST databases with
  NCBI taxids, so
  [`get_tax()`](https://heronoh.github.io/BLASTr/reference/get_tax.md)
  and
  [`blast_lca()`](https://heronoh.github.io/BLASTr/reference/blast_lca.md)
  work as with the NCBI databases; UNITE names are matched to the NCBI
  taxonomy (the deepest name found; homonyms resolved by the higher
  ranks), SILVA gives the NCBI taxid of each record. A table next to the
  database records the taxid of each sequence and how it was found.

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

- All the hits tied with the best one (`at_best`) are kept, even beyond
  `num_alignments`. New summary columns `Taxids at best` (the tied
  taxids, for a lowest common ancestor) and `At best truncated` (the tie
  may go beyond the `search_alignments` subjects BLAST returned).

- New hit fields `nident` and `subject strand`.

- New
  [`as_blast_long()`](https://heronoh.github.io/BLASTr/reference/as_blast_long.md):
  converts wide results of any version (including `k_indentity`) to the
  long format, ranked with the same rules.

- New
  [`get_tax()`](https://heronoh.github.io/BLASTr/reference/get_tax.md):
  taxonomy of taxids from a local SQLite database (NCBI taxdump,
  `taxonomizr` layout), in a fraction of a second for thousands of
  taxids. Takes taxids or any BLASTr result and returns one row per
  taxid: ranks, `lineage`, `lineage_taxids` and `merged_into` (merged
  taxids). When the database does not exist, it is built from the NCBI
  taxdump matching the date of the BLAST database of the search (or the
  one given), once, and kept in the BLASTr cache.
  [`make_taxonomy_db()`](https://heronoh.github.io/BLASTr/reference/make_taxonomy_db.md)
  builds it explicitly.

- Uninformative references:
  [`parallel_blast()`](https://heronoh.github.io/BLASTr/reference/parallel_blast.md)
  gains `uninformative`, the table of
  [`uninformative_references()`](https://heronoh.github.io/BLASTr/reference/uninformative_references.md)
  (shipped in `inst/extdata/uninformative_references.tsv`: uncultured,
  metagenomes and MAGs, environmental samples, synthetic constructs and
  vectors, unidentified, chimeras, and wrong references by accession).
  The hits are flagged (`informative`, `uninformative_category`) and
  ranked among the informative ones, so `rank`, `at_best`, `Check ID`
  and
  [`blast_lca()`](https://heronoh.github.io/BLASTr/reference/blast_lca.md)
  use only those; a sequence with only uninformative hits keeps its best
  ones.
  [`blast_summary()`](https://heronoh.github.io/BLASTr/reference/blast_summary.md)
  gains `N uninformative` and `Best uninformative score` (a sign of a
  taxon missing from the database). The md5 of the table goes to the
  metadata and to the log; `as_blast_long(uninformative = )` flags older
  results.

- [`blast_lca()`](https://heronoh.github.io/BLASTr/reference/blast_lca.md)
  gains `named_first`: hits named only above genus (“Bacterium”,
  “Entolomataceae sp.”) do not vote when a named species is tied with
  them; and a species in open nomenclature (“Picea sp. X”) is ranked at
  its genus in `LCA ranked ...`.

- New
  [`blast_lca()`](https://heronoh.github.io/BLASTr/reference/blast_lca.md):
  consensus lowest common ancestor of the hits tied with the best one
  (by default the deepest taxon holding 90% of them, so a few
  mislabelled references do not push the identification up), from the
  exact lineages of
  [`get_tax()`](https://heronoh.github.io/BLASTr/reference/get_tax.md).
  When the ancestor has none of the ranks of
  [`get_tax()`](https://heronoh.github.io/BLASTr/reference/get_tax.md)
  (a clade, a superfamily), `LCA ranked taxid/name/rank` give the
  deepest ranked node above it, for tables by rank.

- The BLAST command of each search (without the sequences) is kept in
  `BLASTr_metadata$blast_cmd`, with the BLAST version and the date of
  each database; it is shown in interactive sessions and written to the
  log.

- Optional log of the searches: with `log_dir` (or the option
  `BLASTr.log_dir` / the environment variable `BLASTR_LOG_DIR`), each
  search adds a line to `BLASTr_logs/blastr_runs.tsv` with the number
  and lengths of the sequences, database size, strategy, time, CPU time
  and percentage of identifications. See
  [`blast_log_path()`](https://heronoh.github.io/BLASTr/reference/blast_log_path.md)
  and
  [`read_blast_log()`](https://heronoh.github.io/BLASTr/reference/blast_log_path.md).

### Fixed

- [`parse_fasta()`](https://heronoh.github.io/BLASTr/reference/parse_fasta.md)
  failed on R 4.3 (`could not find function "declare"`): it called
  [`declare()`](https://rdrr.io/r/base/declare.html), which only exists
  from R 4.4, to keep the `carrier` import. `carrier` was not used
  otherwise and left `Imports`.
- [`blast_summary()`](https://heronoh.github.io/BLASTr/reference/blast_summary.md):
  `N hits`, `N taxids`, `N uninformative`, `N hits same taxid as 1st`
  and `Best uninformative score` were computed only on the rows kept in
  the long results (the first `num_alignments` and the hits at best),
  not on all the subjects BLAST returned, as documented.
  [`parallel_blast()`](https://heronoh.github.io/BLASTr/reference/parallel_blast.md)
  now computes them on all the subjects and keeps them in
  `BLASTr_metadata$hit_counts`, which
  [`blast_summary()`](https://heronoh.github.io/BLASTr/reference/blast_summary.md)
  uses; `N hits` equal to `search_alignments` now shows that BLAST
  returned a full list. The summary of the wide results and of the log
  is also computed on all the subjects. Results of older versions keep
  the old counts.
- [`get_tax()`](https://heronoh.github.io/BLASTr/reference/get_tax.md)
  ignored the `Taxids at best` column of wide results (it looked for its
  old name, `Taxids at best identity`, still accepted).
- [`parse_fasta()`](https://heronoh.github.io/BLASTr/reference/parse_fasta.md)
  names the missing file in its error message.
- [`make_blast_db()`](https://heronoh.github.io/BLASTr/reference/make_blast_db.md)
  accepts a `taxid_map` with any name (it required the base name of the
  FASTA file), checks that it exists, and prints the `makeblastdb`
  command only when `verbose` is not `"silent"`.
- The documentation of `search_alignments` in
  [`parallel_blast()`](https://heronoh.github.io/BLASTr/reference/parallel_blast.md)
  named `-max_target_seqs`; the option passed to BLAST is
  `-num_alignments`.

### Changed

- [`parallel_blast()`](https://heronoh.github.io/BLASTr/reference/parallel_blast.md):
  `search_alignments` now defaults to `100` (was `50`). In 1260 GenBank
  references searched in `core_nt`, it halved the sequences whose ties
  with the best hit were cut at the limit (170 to 94) and those whose
  ancestor differs from a search with 500 subjects (37 to 18), for 3% to
  18% more time.
- BLAST+ 2.17.0 (was 2.16.0), installed by
  [`install_dependencies()`](https://heronoh.github.io/BLASTr/reference/install_dependencies.md)
  from bioconda in the conda environment `blastr-blast-env`; an
  environment created by an older BLASTr is recreated. On 300 GenBank
  references searched in `core_nt`, 2.17.0 gave the same hits, summary
  and ancestors as 2.16.0, in about the same time.
- Which BLAST+ is used: the one of the conda environment; without it,
  the one of the system, only if it is 2.16.0 or newer (BLAST+ 2.12, for
  example, does not accept `-mt_mode 2`, used on large databases);
  otherwise an error that says to run
  [`install_dependencies()`](https://heronoh.github.io/BLASTr/reference/install_dependencies.md).
  Searches no longer install or update BLAST+. The source is kept in
  `BLASTr_metadata$blast_source` (`"conda"` or `"system"`).
- New
  [`has_blast()`](https://heronoh.github.io/BLASTr/reference/has_blast.md):
  whether BLAST+ 2.16.0 or newer is ready (conda environment or system),
  without installing anything. The examples and the vignette run their
  BLAST code only when it is `TRUE`.
- Tests: the ones that run BLAST+ use the BLAST+ already installed and
  are skipped without it (before, they were skipped on CRAN, and the
  test setup created a new conda environment at each run). A real BLAST+
  output (`inst/extdata/blast_output_example.tsv`, made by
  `data-raw/blast_output_example.R`) tests reading, ranking and summary
  where BLAST+ is not installed.
- [`get_tax()`](https://heronoh.github.io/BLASTr/reference/get_tax.md)
  no longer downloads the NCBI taxdump on its own: without the taxonomy
  database it stops and says how to get it, or asks in an interactive
  session. The new argument `download = TRUE` downloads the taxdump
  (about 80 MB) and builds the database, for those who use only
  [`get_tax()`](https://heronoh.github.io/BLASTr/reference/get_tax.md);
  [`download_blast_db()`](https://heronoh.github.io/BLASTr/reference/ncbi_blast_dbs.md)
  (`taxonomy = TRUE`),
  [`make_unite_db()`](https://heronoh.github.io/BLASTr/reference/reference_dbs.md)
  and
  [`make_silva_db()`](https://heronoh.github.io/BLASTr/reference/reference_dbs.md)
  still bring it with the BLAST databases.
- New
  [`has_taxonomy_db()`](https://heronoh.github.io/BLASTr/reference/has_taxonomy_db.md):
  whether
  [`get_tax()`](https://heronoh.github.io/BLASTr/reference/get_tax.md)
  has its taxonomy database, without downloading or building anything;
  with `db_path`, the database of the date of the BLAST database(s) to
  be searched. A pipeline can check it at the start instead of stopping
  after the BLAST search.
- [`get_tax()`](https://heronoh.github.io/BLASTr/reference/get_tax.md)
  takes `taxonomy_db = NA` (or `""`) as not set, like `NULL`: before,
  `NA` was used as the path of the database.

### Documentation

- Examples that run in 22 of the 25 help pages of exported functions (5
  before), with the small database in `inst/extdata` and a small
  taxonomy written in the example; they run only where `blastn` is found
  and write only to [`tempdir()`](https://rdrr.io/r/base/tempfile.html).
  Titles, families of functions and the package page
  ([`?BLASTr`](https://heronoh.github.io/BLASTr/reference/BLASTr-package.md))
  were rewritten.
- New vignette “Get started with BLASTr”
  ([`vignette("BLASTr")`](https://heronoh.github.io/BLASTr/articles/BLASTr.md));
  the vignette on searches in the background is now
  [`vignette("background-searches")`](https://heronoh.github.io/BLASTr/articles/background-searches.md).
- `DESCRIPTION`: `SystemRequirements` (BLAST+, Quarto) and `biocViews`.

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
  [`get_tax_by_taxID()`](https://heronoh.github.io/BLASTr/reference/get_tax_online.md)
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

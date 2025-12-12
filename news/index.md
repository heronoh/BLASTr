# Changelog

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

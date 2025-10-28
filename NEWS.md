## BLASTr 0.1.6

Release Date: 2025-10-28

Development Changelog: [0.1.6](https://github.com/heronoh/BLASTr/compare/v0.1.5...v0.1.6)

### Fixed

* Fix hard coded `num_threads = 1L` in `parallel_blast()`.
  `parallel_blast()` now correctly uses the user provided `num_threads`
  argument.

* Fix parallel functions not using the the custom compute profile.
  `purrr::map` do not have an argument for specifying the `mirai` `.compute`
  argument correctly.
  Now all parallel functions use the "default" mirai compute profile
  `.compute = NULL`.

* Fix `get_tax_by_taxID()` when multiple tax IDs are provided.

* Fix `parse_fasta()` concatenate lines for multiline FASTA.

## BLASTr 0.1.5

Release Date: 2025-08-27

Development Changelog: [0.1.5](https://github.com/heronoh/BLASTr/compare/v0.1.4...v0.1.5)

### Added

* New parallel framework based on native `purrr` functions,
  using `mirai` and `carrier` in the background.

* New `make_blast_db()` function to create custom BLAST databases.

### Changed

* Internal `entrez-direct` version bump to `24.0`.

* The internal environments used by BLASTr have been renamed to `blastr-blast-env` and `blastr-entrez-env`.

* Remove of the `furrr` based functions in favor of the new `purrr` native parallel.

* All `options()` based parameters have been removed in favor of explicit arguments in the functions.

## BLASTr 0.1.4

### Changed

* All functions have better defaults and less mandatory arguments.

### Fixed

* Fix progress bars when parallel processing is used.

## BLASTr 0.1.3

### Added

* Refactor handling of internal errors.

* `shell_exec` is not available anymore.

## BLASTr 0.1.2

### Added

* New `install_dependencies()` automatically create environments.

* Initial support for `condathis` for managing dependency installation.


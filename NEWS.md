# BLASTr 0.1.5 (dev)

## Breaking Changes

* The internal environments used by BLASTr have been renamed to `blastr-blast-env` and `blastr-entrez-env`.

* Remove of the `furrr` based functions in favor of the new `purrr` native parallel.

* All `options()` based parameters have been removed in favor of explicit arguments in the functions.

### New features

* New parallel framework based on native `purrr` functions,
  using `mirai` and `carrier` in the background.

### Minor improvements and fixes

* Internal `entrez-direct` version bump to `24.0`.

# BLASTr 0.1.3

* Refactor handling of internal errors.

* `shell_exec` is not available anymore.

# BLASTr 0.1.2

* New `install_dependencies()` automatically create environments.

* Initial support for `condathis` for managing dependency installation.

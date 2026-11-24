# Path in the BLASTr cache directory

Determines the appropriate user cache directory for the `BLASTr` package
based on the operating system. On macOS, it avoids using paths with
spaces and follows the XDG base directory specification.

## Usage

``` r
get_blastr_cache(..., ext = NULL)
```

## Arguments

- ...:

  Additional path components to append to the cache directory.

- ext:

  An optional file extension (e.g., "rds", "csv") to append to the final
  path.

## Value

A path character string representing the path to the user cache
directory for the `BLASTr` package.

## Details

This function uses the `[tools::R_user_dir()]` function to determine the
user cache directory.

## See also

Other setup functions:
[`has_blast()`](https://heronoh.github.io/BLASTr/reference/has_blast.md),
[`install_dependencies()`](https://heronoh.github.io/BLASTr/reference/install_dependencies.md)

## Examples

``` r
get_blastr_cache()
#> /home/runner/.cache/R/BLASTr
# where get_tax() keeps the taxonomy of a BLAST database of a given date
get_blastr_cache("taxonomy", "taxonomy_2026-01-01", ext = "sqlite")
#> /home/runner/.cache/R/BLASTr/taxonomy/taxonomy_2026-01-01.sqlite
```

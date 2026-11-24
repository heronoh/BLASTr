# Is BLAST+ ready for BLASTr?

Checks, without installing or downloading anything, whether BLASTr has a
`blastn` it runs with (BLAST+ 2.16.0 or newer): the one of its conda
environment (installed by
[`install_dependencies()`](https://heronoh.github.io/BLASTr/reference/install_dependencies.md),
BLAST+ 2.17.0) or, when that environment has none, the one of the
system. Useful to run code only where BLAST+ is ready, as in the
examples.

## Usage

``` r
has_blast(env_name = "blastr-blast-env")
```

## Arguments

- env_name:

  Name of the conda environment.

## Value

`TRUE` or `FALSE`.

## See also

Other setup functions:
[`get_blastr_cache()`](https://heronoh.github.io/BLASTr/reference/get_blastr_cache.md),
[`install_dependencies()`](https://heronoh.github.io/BLASTr/reference/install_dependencies.md)

## Examples

``` r
has_blast()
#> [1] FALSE
```

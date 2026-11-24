# Install BLAST+ in a conda environment

Installs BLAST+ 2.17.0 (`blastn`, `makeblastdb`, ...) from bioconda into
the conda environment `blastr-blast-env`, through the condathis package.
BLASTr runs BLAST+ from this environment; without it, from the system,
if it is BLAST+ 2.16.0 or newer (see
[`has_blast()`](https://heronoh.github.io/BLASTr/reference/has_blast.md)).
When the environment exists with another version of BLAST+ (created by
an older BLASTr), it is recreated with BLAST+ 2.17.0. Searches never
install or update BLAST+: only this function does. The online taxonomy
functions
([`get_tax_by_taxID()`](https://heronoh.github.io/BLASTr/reference/get_tax_online.md),
[`get_tax_by_name()`](https://heronoh.github.io/BLASTr/reference/get_tax_online.md))
query the NCBI over HTTP and need no external tool.

## Usage

``` r
install_dependencies(
  verbose = c("silent", "cmd", "output", "full"),
  force = FALSE
)
```

## Arguments

- verbose:

  A character string specifying the verbosity level during environment
  creation. Options include `"silent"`, `"cmd"`, `"output"`, or
  `"full"`. Default is `"silent"`.

- force:

  Recreate the conda environment even if it exists. Default is `FALSE`.

## Value

Invisibly returns `TRUE` after attempting to install the dependencies.

## See also

Other setup functions:
[`get_blastr_cache()`](https://heronoh.github.io/BLASTr/reference/get_blastr_cache.md),
[`has_blast()`](https://heronoh.github.io/BLASTr/reference/has_blast.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Install dependencies with default settings
install_dependencies()

# Install dependencies with verbose output
install_dependencies(verbose = "output")

# Force re-installation of dependencies
install_dependencies(force = TRUE)
} # }
```

# Install Required Command-line Tool Dependencies

Installs BLAST+ (`blastn`, `makeblastdb`, ...) into the conda
environment `blastr-blast-env` if it is not already there. The online
taxonomy functions
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

  A logical value indicating whether to force the re-creation of the
  Conda environments even if they already exist. Default is `FALSE`.

## Value

Invisibly returns `TRUE` after attempting to install the dependencies.

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

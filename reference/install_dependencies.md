# Install Required Command-line Tool Dependencies

Installs the necessary command-line tools ('blast' and 'efetch') into
conda environments if they are not already available. This ensures that
all dependencies required by the package are present and up-to-date.

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

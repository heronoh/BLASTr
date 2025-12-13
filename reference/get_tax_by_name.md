# Get taxonomy ranks using *NCBI Taxonomy Tax ID*

Recover complete taxonomy for a given *NCBI Taxonomy Tax ID*

## Usage

``` r
get_tax_by_name(
  organisms_names,
  parse_result = TRUE,
  verbose = FALSE,
  env_name = "entrez-env"
)
```

## Arguments

- organisms_names:

  Vector of *NCBI Taxonomy Tax ID* to retrieve taxonomy for.

- parse_result:

  Should the taxonomy be returned as the `efetch` returns it or should
  it be parsed into a tibble.

- verbose:

  Character string indicating whether to print verbose messages during
  the process. Default is `silent`.

- env_name:

  Character string specifying the name of the conda environment where
  `efetch` is installed. Default is `"blastr-entrez-env"`.

- total_cores:

  Number of threads to use. Defaults to 1.

- show_command:

  Should the `efetch` command on command line be printed for debugging.
  Defaults to TRUE.

## Value

A tibble with all the taxonomic ranks for the corresponding taxID.

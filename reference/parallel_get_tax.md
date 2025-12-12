# Retrieve Taxonomic Ranks for a List of NCBI Taxonomy Tax IDs in Parallel

Retrieves taxonomy ranks for a list of NCBI Taxonomy Tax IDs using
parallel processing. The function queries the NCBI database and can
retry fetching taxonomy information multiple times if initial attempts
fail.

## Usage

``` r
parallel_get_tax(
  organisms_taxIDs,
  parse_result = TRUE,
  total_cores = 1L,
  retry_times = 10L,
  verbose = c("silent", "cmd", "output", "full"),
  env_name = "blastr-entrez-env"
)
```

## Arguments

- organisms_taxIDs:

  A character vector of NCBI Taxonomy Tax IDs to retrieve taxonomy
  information for.

- parse_result:

  Logical indicating whether to parse the taxonomy information into a
  tibble (`TRUE`, default) or return the raw output as returned by
  `efetch` (`FALSE`).

- total_cores:

  Integer specifying the number of cores to use for parallel processing.
  Defaults to `1`.

- retry_times:

  Integer specifying the number of times to retry fetching taxonomy
  information if it fails. Defaults to `10`.

- verbose:

  Character string indicating whether to print verbose messages during
  the process. Default is `silent`.

- env_name:

  Character string specifying the name of the conda environment where
  `efetch` is installed. Default is `"blastr-entrez-env"`.

## Value

A tibble containing the taxonomic ranks for the given Tax IDs.

## Examples

``` r
if (FALSE) { # \dontrun{
# Retrieve taxonomy for multiple Tax IDs using 2 processes
tax_ids <- c("9606", "10090", "10116") # Human, Mouse, Rat
tax_info <- parallel_get_tax(tax_ids, total_cores = 2)

# Retrieve unparsed taxonomy result
tax_info_unparsed <- parallel_get_tax(tax_ids, parse_result = FALSE)

# Change the number of retry attempts
tax_info <- parallel_get_tax(tax_ids, retry_times = 2)

# Enable verbose output
tax_info <- parallel_get_tax(tax_ids, verbose = "output")
} # }
```

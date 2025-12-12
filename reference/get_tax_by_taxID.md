# Retrieve Taxonomic Ranks Using NCBI Taxonomy Tax IDs

Retrieves complete taxonomy information for given NCBI Taxonomy Tax IDs
by querying the NCBI database using the `efetch` command.

## Usage

``` r
get_tax_by_taxID(
  organisms_taxIDs,
  parse_result = TRUE,
  verbose = c("silent", "cmd", "output", "full"),
  env_name = "blastr-entrez-env"
)
```

## Arguments

- organisms_taxIDs:

  A character vector of NCBI Taxonomy Tax IDs for which to retrieve
  taxonomy information.

- parse_result:

  Logical indicating whether to parse the taxonomy information into a
  tibble (`TRUE`, default) or return the raw output as returned by
  `efetch` (`FALSE`).

- verbose:

  Character indicating whether to print verbose messages during the
  process. Default is `"silent"`.

- env_name:

  Character string specifying the name of the conda environment where
  `efetch` is installed. Default is `"blastr-entrez-env"`.

## Value

A tibble containing the taxonomic ranks for the given Tax IDs.

## Examples

``` r
if (FALSE) { # \dontrun{
# Retrieve taxonomy for a single Tax ID
tax_info <- get_tax_by_taxID("9606") # Human

# Retrieve taxonomy for multiple Tax IDs
tax_ids <- c("9606", "10090", "10116") # Human, Mouse, Rat
tax_info <- get_tax_by_taxID(tax_ids)

# Get unparsed taxonomy result
raw_tax_info <- get_tax_by_taxID("9606", parse_result = FALSE)

# Enable verbose output
tax_info <- get_tax_by_taxID("9606", verbose = "output")
} # }
```

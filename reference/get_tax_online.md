# Taxonomy from the NCBI online, without a local database

Look up taxonomy in the NCBI Taxonomy database through the E-utilities
(`efetch`, `esearch`), over the internet. For when the local database of
[`get_tax()`](https://heronoh.github.io/BLASTr/reference/get_tax.md) is
not there yet, or for a few taxa;
[`get_tax()`](https://heronoh.github.io/BLASTr/reference/get_tax.md) is
much faster for many taxids and works offline once built.

- `get_tax_by_taxID()`: taxonomy of taxids (or of the taxids of a BLASTr
  result, as in
  [`get_tax()`](https://heronoh.github.io/BLASTr/reference/get_tax.md)).

- `get_tax_by_name()`: taxids and taxonomy of names, matched first as
  scientific names, then as any name the NCBI knows (synonyms, common
  names).

Both return the same columns as
[`get_tax()`](https://heronoh.github.io/BLASTr/reference/get_tax.md), so
the tables can be used in place of each other. The taxonomy is the
current one at the NCBI, which may be newer than the BLAST database.

The NCBI accepts 3 requests per second, or 10 with an API key (see
<https://support.nlm.nih.gov/kbArticle/?pn=KA-05317>): set it in
`api_key` or in the environment variable `ENTREZ_KEY` or `NCBI_API_KEY`.
Taxids are sent in batches of 200.

## Usage

``` r
get_tax_by_taxID(
  organisms_taxIDs,
  ranks = tax_ranks_default,
  api_key = NULL,
  ...
)

get_tax_by_name(
  organisms_names,
  ranks = tax_ranks_default,
  api_key = NULL,
  ...
)

parallel_get_tax(organisms_taxIDs, ...)
```

## Arguments

- organisms_taxIDs:

  Taxids (character or numeric), or a BLASTr result.

- ranks:

  Ranks returned as columns, as in
  [`get_tax()`](https://heronoh.github.io/BLASTr/reference/get_tax.md).

- api_key:

  NCBI API key (optional).

- ...:

  Ignored (arguments of older versions).

- organisms_names:

  Names of taxa (e.g. `"Homo sapiens"`, `"Salminus"`).

## Value

A tibble as in
[`get_tax()`](https://heronoh.github.io/BLASTr/reference/get_tax.md):
`taxid`, `found`, `sci_name`, `rank`, one column per rank, `lineage`,
`lineage_taxids` and `merged_into`. `get_tax_by_name()` adds, before
them, `name` (as given) and `matched_by` (`"scientific name"`,
`"other name"` or `NA` when not found); a name with more than one taxid
(a homonym) gets one row per taxid.

## See also

Other taxonomy functions:
[`get_tax()`](https://heronoh.github.io/BLASTr/reference/get_tax.md),
[`has_taxonomy_db()`](https://heronoh.github.io/BLASTr/reference/has_taxonomy_db.md),
[`make_taxonomy_db()`](https://heronoh.github.io/BLASTr/reference/make_taxonomy_db.md)

## Examples

``` r
# \donttest{
# needs internet
if (requireNamespace("curl", quietly = TRUE) && curl::has_internet()) {
  get_tax_by_taxID(c("9606", "10090"))
  get_tax_by_name(c("Homo sapiens", "Salminus brasiliensis"))
}
#> # A tibble: 2 × 21
#>   name     matched_by taxid found sci_name rank  domain kingdom phylum subphylum
#>   <chr>    <chr>      <chr> <lgl> <chr>    <chr> <chr>  <chr>   <chr>  <chr>    
#> 1 Homo sa… scientifi… 9606  TRUE  Homo sa… spec… Eukar… Metazoa Chord… Craniata 
#> 2 Salminu… scientifi… 9302… TRUE  Salminu… spec… Eukar… Metazoa Chord… Craniata 
#> # ℹ 11 more variables: class <chr>, subclass <chr>, order <chr>,
#> #   suborder <chr>, family <chr>, subfamily <chr>, genus <chr>, species <chr>,
#> #   lineage <chr>, lineage_taxids <chr>, merged_into <chr>
# }
```

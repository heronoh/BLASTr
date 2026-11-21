# Deprecated taxonomy functions

**\[deprecated\]**

`get_tax_by_taxID()`, `parallel_get_tax()` and `get_tax_by_name()`
queried the NCBI online (`efetch`), one taxid at a time. Use
[`get_tax()`](https://heronoh.github.io/BLASTr/reference/get_tax.md),
which reads a local taxonomy database and is much faster. These
functions now call it and return its table.

## Usage

``` r
get_tax_by_taxID(organisms_taxIDs, ..., taxonomy_db = NULL)

parallel_get_tax(organisms_taxIDs, ..., taxonomy_db = NULL)

get_tax_by_name(organisms_names, ..., taxonomy_db = NULL)
```

## Arguments

- organisms_taxIDs:

  Taxids.

- ...:

  Ignored (arguments of the old versions).

- taxonomy_db:

  See
  [`get_tax()`](https://heronoh.github.io/BLASTr/reference/get_tax.md).

- organisms_names:

  Scientific (or other) names of taxa.

## Value

See
[`get_tax()`](https://heronoh.github.io/BLASTr/reference/get_tax.md).

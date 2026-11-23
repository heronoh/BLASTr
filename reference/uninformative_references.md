# References that do not inform the identification

The table of reference sequences that
[`parallel_blast()`](https://heronoh.github.io/BLASTr/reference/parallel_blast.md)
sets aside when ranking the hits: sequences of uncultured or
unidentified organisms, metagenomes and MAGs, environmental samples,
synthetic constructs and vectors, chimeras, and specific references
known to be wrong. Their taxid says little about the query (often only
"Bacteria" or "environmental samples"), so a hit to one of them should
not decide the identification when a hit to a named organism is also
there.

Each row is a rule:

- `field`: `"header"` (`pattern` is a Perl regular expression, matched
  against the subject header; use `(?i:...)` groups for case, as the
  rules are joined with `"|"`) or `"accession"` (`pattern` is a subject
  accession, with or without version).

- `category`: the reason, reported in `uninformative_category`.

- `note`: free text.

Organisms named only at a high rank or in open nomenclature (e.g.
"Bacterium", "Picea sp.") are not here: their lineage is right, only
shallow.
[`blast_lca()`](https://heronoh.github.io/BLASTr/reference/blast_lca.md)
sets them aside when a named species is tied with them.

## Usage

``` r
uninformative_references(x = NULL)
```

## Arguments

- x:

  `NULL` (default) for the table shipped with BLASTr, a path to a
  tab-separated file with the same columns, or a data frame.

## Value

A tibble with `field`, `pattern`, `category` and `note`, with the
attribute `md5` (of its content, to record which table was used).

## Examples

``` r
uninformative_references()
#> # A tibble: 14 × 4
#>    field     pattern                                              category note 
#>    <chr>     <chr>                                                <chr>    <chr>
#>  1 header    "(?i:uncultured)"                                    uncultu… uncu…
#>  2 header    "^MAG"                                               metagen… meta…
#>  3 header    "(?i:metagenome)"                                    metagen… meta…
#>  4 header    "(?i:\\benvironmental (sample|DNA)\\b)"              environ… envi…
#>  5 header    "^Leaf litter"                                       environ… envi…
#>  6 header    "16S rRNA amplicon fragment"                         environ… envi…
#>  7 header    "Complete LGCome-Assembled"                          environ… envi…
#>  8 header    "root endophyte"                                     environ… envi…
#>  9 header    "root associated fungus"                             environ… envi…
#> 10 header    "(?i:synthetic construct)"                           synthet… synt…
#> 11 header    "(?i:\\b(cloning|expression|shuttle|transformation|… synthet… vect…
#> 12 header    "(?i:\\bunidentified\\b)"                            unident… unid…
#> 13 header    "(?i:\\bchimeric\\b)"                                chimeric chim…
#> 14 accession "XR_004027927.1"                                     bad ref… PRED…
```

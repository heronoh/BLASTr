# Retrieve Exit Codes and Standard Error from BLASTr Results

This function extracts the exit codes and standard error messages from
the results of a BLAST search performed using the
[`parallel_blast()`](https://heronoh.github.io/BLASTr/reference/parallel_blast.md)
function.

## Usage

``` r
exit_codes(blast_res)
```

## Arguments

- blast_res:

  A data frame containing the results of a BLAST search.

## Value

A data frame with the exit codes and standard error messages.

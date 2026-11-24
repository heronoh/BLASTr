# Exit codes and standard error of the BLAST runs

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

  Results of
  [`parallel_blast()`](https://heronoh.github.io/BLASTr/reference/parallel_blast.md).

## Value

A tibble with one row per sequence: `query_seq`, `exit_code` of the
`blastn` run of its batch (`NA` when cancelled) and `stderr`.

## See also

Other search functions:
[`blast_job()`](https://heronoh.github.io/BLASTr/reference/blast_job.md),
[`blast_log_path()`](https://heronoh.github.io/BLASTr/reference/blast_log_path.md),
[`blast_strategy()`](https://heronoh.github.io/BLASTr/reference/blast_strategy.md),
[`parallel_blast()`](https://heronoh.github.io/BLASTr/reference/parallel_blast.md),
[`parallel_blast_hsps()`](https://heronoh.github.io/BLASTr/reference/parallel_blast_hsps.md)

## Examples

``` r
# needs BLAST+ (see install_dependencies())
if (has_blast()) {
  fasta <- system.file("extdata", "minimal_db_blast.fasta", package = "BLASTr")
  db <- file.path(tempdir(), "minimal_db_blast")
  make_blast_db(fasta, db_path = db)
  seqs <- readLines(system.file("extdata", "seqs_test.txt", package = "BLASTr"))
  res <- parallel_blast(seqs, db)
  exit_codes(res)
}
```

# Read the sequences of a FASTA file

Reads a FASTA file and returns its sequences, without the headers, to
pass to
[`parallel_blast()`](https://heronoh.github.io/BLASTr/reference/parallel_blast.md)

## Usage

``` r
parse_fasta(file_path)
```

## Arguments

- file_path:

  Path to FASTA file.

## Value

vector containing the sequences.

## Examples

``` r
fasta_path <- fs::path_package("BLASTr", "extdata", "minimal_db_blast", ext = "fasta")
seqs <- parse_fasta(file_path = fasta_path)
```

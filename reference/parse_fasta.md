# Parsing fasta file to sequences only

Extract sequences From FASTA file

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

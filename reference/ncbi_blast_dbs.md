# BLAST databases of the NCBI

`list_ncbi_blast_dbs()` lists the preformatted BLAST databases the NCBI
publishes (<https://ftp.ncbi.nlm.nih.gov/blast/db/>), with their size
and date. `download_blast_db()` downloads one of them (e.g. `"core_nt"`,
`"16S_ribosomal_RNA"`, `"ITS_eukaryote_sequences"`), checks the md5 of
each file, unpacks it and returns the path to give to
[`parallel_blast()`](https://heronoh.github.io/BLASTr/reference/parallel_blast.md).

The database goes to `<dir>/<db>_<date>/` (the NCBI date of the
database), so a newer release does not overwrite the one older results
were searched against. An interrupted download is resumed: files already
unpacked (with a matching md5) are not downloaded again. With
`taxonomy = TRUE`, the taxonomy database of
[`get_tax()`](https://heronoh.github.io/BLASTr/reference/get_tax.md)
matching the date of the BLAST database is built too, so the
identifications can be made offline.

`core_nt` has about 300 GB unpacked (and is downloaded in volumes of a
few GB, unpacked one at a time); see `list_ncbi_blast_dbs()` for the
sizes.

## Usage

``` r
list_ncbi_blast_dbs()

download_blast_db(
  db = "core_nt",
  dir,
  taxdb = TRUE,
  taxonomy = TRUE,
  keep_archives = FALSE,
  check_space = TRUE
)
```

## Arguments

- db:

  Name of the database, as in `list_ncbi_blast_dbs()$db`.

- dir:

  Directory where the database directory is created.

- taxdb:

  Also download `taxdb` (scientific names of the taxids, used by BLAST
  in `sscinames`) into the same directory.

- taxonomy:

  Build the local taxonomy database of
  [`get_tax()`](https://heronoh.github.io/BLASTr/reference/get_tax.md)
  for the date of this database (kept in the BLASTr cache).

- keep_archives:

  Keep the `.tar.gz` files after unpacking.

- check_space:

  Stop when the disk does not have room for the unpacked database plus
  one archive.

## Value

`list_ncbi_blast_dbs()`: a tibble with `db`, `type` (`"nucl"`/`"prot"`),
`description`, `n_seqs`, `n_letters`, `n_volumes`, `size_gb` (unpacked),
`download_gb` and `last_updated`. `download_blast_db()`: the path of the
database (directory and name, as `db_path` of
[`parallel_blast()`](https://heronoh.github.io/BLASTr/reference/parallel_blast.md)),
invisibly.

## See also

Other database functions:
[`make_blast_db()`](https://heronoh.github.io/BLASTr/reference/make_blast_db.md),
[`reference_dbs`](https://heronoh.github.io/BLASTr/reference/reference_dbs.md)

## Examples

``` r
# \donttest{
# needs internet
if (requireNamespace("curl", quietly = TRUE) && curl::has_internet()) {
  dbs <- list_ncbi_blast_dbs()
  dbs[dbs$download_gb < 1, c("db", "type", "description", "download_gb")]
}
#> # A tibble: 19 × 4
#>    db                      type  description                         download_gb
#>    <chr>                   <chr> <chr>                                     <dbl>
#>  1 16S_ribosomal_RNA       nucl  16S ribosomal RNA (Bacteria and Ar…        0.07
#>  2 18S_fungal_sequences    nucl  18S ribosomal RNA sequences (SSU) …        0.07
#>  3 28S_fungal_sequences    nucl  28S ribosomal RNA sequences (LSU) …        0.07
#>  4 ITS_RefSeq_Fungi        nucl  Internal transcribed spacer region…        0.07
#>  5 ITS_eukaryote_sequences nucl  ITS eukaryote BLAST                        0.08
#>  6 LSU_eukaryote_rRNA      nucl  Large subunit ribosomal nucleic ac…        0.06
#>  7 LSU_prokaryote_rRNA     nucl  Large subunit ribosomal nucleic ac…        0.06
#>  8 SSU_eukaryote_rRNA      nucl  Small subunit ribosomal nucleic ac…        0.06
#>  9 env_nt                  nucl  environmental samples                      0.42
#> 10 mito                    nucl  NCBI Genomic Mitochondrial Referen…        0.2 
#> 11 pdbnt                   nucl  PDB nucleotide database                    0.07
#> 12 ref_viroids_rep_genomes nucl  Refseq viroids reference genomes           0.07
#> 13 ref_viruses_rep_genomes nucl  Refseq viruses reference genomes           0.21
#> 14 refseq_select_rna       nucl  RefSeq Select RNA sequences                0.12
#> 15 taxdb                   nucl  Taxonomy database                          0.07
#> 16 landmark                prot  Landmark database for SmartBLAST           0.21
#> 17 pataa                   prot  Protein sequences derived from the…        0.84
#> 18 pdbaa                   prot  PDB protein database                       0.13
#> 19 swissprot               prot  Non-redundant UniProtKB/SwissProt …        0.23
# }
if (FALSE) { # \dontrun{
# downloads the database and the taxonomy of the same date
db <- download_blast_db("16S_ribosomal_RNA", dir = "~/blast_databases")
parallel_blast(seqs, db_path = db, total_cores = 8)
} # }
```

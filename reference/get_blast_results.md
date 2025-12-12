# Run BLASTn

Retrieve BLAST results from a given ASV

## Usage

``` r
get_blast_results(
  asv,
  db_path,
  num_threads = 1L,
  blast_type = "blastn",
  perc_id = 80L,
  perc_qcov_hsp = 80L,
  num_alignments = 4L,
  verbose = c("silent", "cmd", "output", "full"),
  env_name = "blastr-blast-env"
)
```

## Arguments

- asv:

  Vector of sequences to be BLASTed.

- db_path:

  Complete path do formatted BLAST database.

- num_threads:

  Number of threads to run BLAST on. Passed on to BLAST+ argument
  `-num_threads`.

- blast_type:

  One of the available BLAST+ search engines, \#' one of:
  `c("blastn", "blastp", "blastx", "tblastn", "tblastx")`.

- perc_id:

  Lowest identity percentage cutoff. Passed on to BLAST+
  `-perc_identity`.

- perc_qcov_hsp:

  Lowest query coverage per HSP percentage cutoff. Passed on to BLAST+
  `-qcov_hsp_perc`.

- num_alignments:

  Number of alignments to retrieve from BLAST. Max = 6.

- verbose:

  Should condathis::run() internal command be shown?

- env_name:

  The name of the conda environment used to run command-line tools (i.e.
  "blastr-blast-env").

## Value

A `tibble` with the results of BLASTn for each sequence.

## Examples

``` r
if (FALSE) { # \dontrun{
dna_fasta_path <- fs::path_package(
  "BLASTr", "extdata", "minimal_db_blast",
  ext = "fasta"
)
temp_db_path <- fs::path_temp("minimal_db_blast")
make_blast_db(fasta_path = dna_fasta_path, db_path = temp_db_path)
asvs_string <- "CTAGCCATAAACTTAAATGAAGCTATACTAAACTCGTTCGCCAG
AGTACTACAAGCGAAAGCTTAAAACTCATAGGACTTGGCGGTGTTTCAGACCCAC"

get_blast_results(asv = asvs_string, db_path = temp_db_path)
} # }
```

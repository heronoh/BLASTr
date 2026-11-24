# Get started with BLASTr

BLASTr searches many sequences with BLAST+ from R and returns the hits
ranked, a summary per sequence and the consensus lowest common ancestor
of the best hits, with the NCBI taxonomy. This vignette goes from a
FASTA file to a summary per sequence, with the small database that comes
with the package. It runs in a few seconds and needs no internet.

``` r

library(BLASTr)
```

## Install

BLASTr is installed from GitHub; BLAST+ 2.17.0 is installed by
[`install_dependencies()`](https://heronoh.github.io/BLASTr/reference/install_dependencies.md)
in a conda environment of its own (`blastr-blast-env`), through the
condathis package. Without it, BLASTr uses the BLAST+ of the system, if
it is 2.16.0 or newer.

``` r

# install.packages("pak")
pak::pak("heronoh/BLASTr")
BLASTr::install_dependencies()
```

The code below runs only where BLAST+ is ready
([`has_blast()`](https://heronoh.github.io/BLASTr/reference/has_blast.md)).

``` r

blast_ready <- has_blast()
blast_ready
#> [1] FALSE
```

## A BLAST database from a FASTA file

[`make_blast_db()`](https://heronoh.github.io/BLASTr/reference/make_blast_db.md)
formats a FASTA file as a BLAST database. The package has a small one,
`minimal_db_blast.fasta`, with 47 GenBank sequences. For real analyses,
[`download_blast_db()`](https://heronoh.github.io/BLASTr/reference/ncbi_blast_dbs.md)
downloads any database of the NCBI (such as `core_nt`) with the taxonomy
of the same date, and
[`make_unite_db()`](https://heronoh.github.io/BLASTr/reference/reference_dbs.md)
and
[`make_silva_db()`](https://heronoh.github.io/BLASTr/reference/reference_dbs.md)
build UNITE and SILVA as BLAST databases (see their help pages).

``` r

fasta <- system.file("extdata", "minimal_db_blast.fasta", package = "BLASTr")
db <- file.path(tempdir(), "minimal_db_blast")
make_blast_db(fasta, db_path = db)
```

## Searching

[`parallel_blast()`](https://heronoh.github.io/BLASTr/reference/parallel_blast.md)
takes a character vector of sequences and the path of the database. It
chooses how to split the search (how many `blastn` processes, threads
per process and sequences per batch) from the size of the database and
the cores (`total_cores`) and RAM available;
[`blast_strategy()`](https://heronoh.github.io/BLASTr/reference/blast_strategy.md)
shows the choice before the search.

``` r

seqs <- readLines(system.file("extdata", "seqs_test.txt", package = "BLASTr"))
length(seqs)
blast_strategy(n_seqs = length(seqs), db_path = db, total_cores = 2)
res <- parallel_blast(seqs, db, total_cores = 2)
```

By default, BLAST keeps up to 100 subjects per sequence
(`search_alignments`) with at least 80% identity (`perc_id`) and 50%
query coverage (`search_qcov_hsp`).

## Reading the results

The result has one row per hit, best first, and one row with empty hit
fields for each sequence without hits:

``` r

res[, c("Sequence", "rank", "tie", "at_best", "subject", "identity", "qcovhsp", "Similarity score", "low_qcov")]
```

- `Similarity score` is `(10 * identity + qcovhsp) / 11` (see
  [`similarity_score()`](https://heronoh.github.io/BLASTr/reference/similarity_score.md)).
  The hits are ranked by it (`rank_by = "similarity"`), then by
  identity, coverage, bitscore and BLAST order. `blast_rank` keeps the
  order of BLAST, which favours longer alignments.
- `at_best` marks the hits within `check_tolerance` (0.5 points) of the
  best one: the tied candidates for the identification. `tie` marks hits
  with the same value as another one.
- `low_qcov` marks hits that cover less than `perc_qcov_hsp` (80%) of
  the sequence; they are kept and ranked, but the identification of a
  sequence whose best hit has `low_qcov` rests on a partial alignment.
- `informative` is `FALSE` for references that do not inform the
  identification (uncultured, environmental, metagenomes, vectors…; see
  [`uninformative_references()`](https://heronoh.github.io/BLASTr/reference/uninformative_references.md));
  they get no `rank` and come after the informative hits.
- The first 10 ranked hits (`num_alignments`) are kept, plus all the
  hits at best.

## One row per sequence

[`blast_summary()`](https://heronoh.github.io/BLASTr/reference/blast_summary.md)
summarises the hits of each sequence:

``` r

blast_summary(res)
```

`Check ID` is `"best ID"` when all the hits at best among the first 10
have the taxid of the first hit, and `"check ID"` otherwise (another
taxid, or an unknown one, as here: the example database has no taxids).
`At best truncated` is `TRUE` when the last subject BLAST returned is
tied with the best one, so more tied subjects may exist beyond
`search_alignments`.

## Taxonomy and the lowest common ancestor

With a database that has taxids
([`download_blast_db()`](https://heronoh.github.io/BLASTr/reference/ncbi_blast_dbs.md),
[`make_unite_db()`](https://heronoh.github.io/BLASTr/reference/reference_dbs.md),
[`make_silva_db()`](https://heronoh.github.io/BLASTr/reference/reference_dbs.md),
or `make_blast_db(taxid_map = )`),
[`get_tax()`](https://heronoh.github.io/BLASTr/reference/get_tax.md)
adds the NCBI taxonomy of each taxid from a local SQLite database (it
comes with
[`download_blast_db()`](https://heronoh.github.io/BLASTr/reference/ncbi_blast_dbs.md);
otherwise `get_tax(res, download = TRUE)` downloads the NCBI taxdump of
the date of the BLAST database once), and
[`blast_lca()`](https://heronoh.github.io/BLASTr/reference/blast_lca.md)
gives the deepest taxon that contains at least 90% of the hits tied with
the best one:

``` r

tax <- get_tax(res)
dplyr::left_join(res, tax, by = c(staxid = "taxid"))
lca <- blast_lca(res)
```

The help pages of
[`get_tax()`](https://heronoh.github.io/BLASTr/reference/get_tax.md) and
[`blast_lca()`](https://heronoh.github.io/BLASTr/reference/blast_lca.md)
have examples that run with a small taxonomy written on the spot.

## Next

- [`vignette("background-searches", package = "BLASTr")`](https://heronoh.github.io/BLASTr/articles/background-searches.md):
  long searches in the background, resumed after a restart.
- [`?parallel_blast_hsps`](https://heronoh.github.io/BLASTr/reference/parallel_blast_hsps.md):
  every HSP, with the BLAST fields of your choice.

## Session information

``` r

sessionInfo()
#> R version 4.6.2 (2026-10-29)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 26.04.1 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.32.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] BLASTr_0.1.7
#> 
#> loaded via a namespace (and not attached):
#>  [1] digest_0.6.39   fastmap_1.2.0   xfun_0.62       glue_1.8.1     
#>  [5] knitr_1.52      htmltools_0.5.9 rmarkdown_2.32  lifecycle_1.0.5
#>  [9] cli_3.6.6       vctrs_0.7.3     compiler_4.6.2  tools_4.6.2    
#> [13] evaluate_1.0.5  pillar_1.11.1   yaml_2.3.12     otel_0.3.0     
#> [17] rlang_1.3.0     jsonlite_2.0.0  fs_2.1.0
```

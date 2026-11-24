# BLASTr: Parallel Taxonomic Classification of Metabarcoding Sequences

BLASTr runs BLAST+ from R on many sequences at once and turns the hits
into an identification that can be checked. The search is split in
batches and processes chosen from the size of the database and the cores
and RAM available
([`blast_strategy()`](https://heronoh.github.io/BLASTr/reference/blast_strategy.md));
the hits of each sequence are ranked by the Similarity score, and the
hits tied with the best one give the summary per sequence (`Check ID`)
and the consensus lowest common ancestor, with the NCBI taxonomy from a
local database.

## Main functions

- [`make_blast_db()`](https://heronoh.github.io/BLASTr/reference/make_blast_db.md),
  [`download_blast_db()`](https://heronoh.github.io/BLASTr/reference/ncbi_blast_dbs.md),
  [`make_unite_db()`](https://heronoh.github.io/BLASTr/reference/reference_dbs.md),
  [`make_silva_db()`](https://heronoh.github.io/BLASTr/reference/reference_dbs.md):
  BLAST databases, from a FASTA file or downloaded.

- [`parallel_blast()`](https://heronoh.github.io/BLASTr/reference/parallel_blast.md):
  the search; long results, one row per hit.

- [`blast_summary()`](https://heronoh.github.io/BLASTr/reference/blast_summary.md):
  one row per sequence, with `Check ID`.

- [`blast_lca()`](https://heronoh.github.io/BLASTr/reference/blast_lca.md):
  consensus lowest common ancestor of the tied hits.

- [`get_tax()`](https://heronoh.github.io/BLASTr/reference/get_tax.md)
  and
  [`make_taxonomy_db()`](https://heronoh.github.io/BLASTr/reference/make_taxonomy_db.md):
  NCBI taxonomy of the taxids, from a local SQLite database.

- [`parallel_blast_hsps()`](https://heronoh.github.io/BLASTr/reference/parallel_blast_hsps.md):
  every HSP, with the BLAST fields chosen.

- [`blast_job()`](https://heronoh.github.io/BLASTr/reference/blast_job.md):
  long searches in the background.

## External programs

BLAST+ (`blastn`, `makeblastdb`, `blastdbcmd`), installed in a conda
environment by
[`install_dependencies()`](https://heronoh.github.io/BLASTr/reference/install_dependencies.md)
(through the condathis package).

## Vignettes

[`vignette("BLASTr", package = "BLASTr")`](https://heronoh.github.io/BLASTr/articles/BLASTr.md)
(get started) and
[`vignette("background-searches", package = "BLASTr")`](https://heronoh.github.io/BLASTr/articles/background-searches.md)
(long searches in the background).

## See also

Useful links:

- <https://github.com/heronoh/BLASTr>

- <https://heronoh.github.io/BLASTr/>

- Report bugs at <https://github.com/heronoh/BLASTr/issues>

## Author

**Maintainer**: Heron Hilário <heronoh@gmail.com>
([ORCID](https://orcid.org/0000-0002-2548-1053))

Authors:

- Heron Hilário <heronoh@gmail.com>
  ([ORCID](https://orcid.org/0000-0002-2548-1053))

- Lucio Queiroz <luciorqueiroz@gmail.com>
  ([ORCID](https://orcid.org/0000-0002-6090-1834))

#' @description
#' BLASTr runs BLAST+ from R on many sequences at once and turns the hits into
#' an identification that can be checked. The search is split in batches and
#' processes chosen from the size of the database and the cores and RAM
#' available ([blast_strategy()]); the hits of each sequence are ranked by the
#' Similarity score, and the hits tied with the best one give the summary per
#' sequence (`Check ID`) and the consensus lowest common ancestor, with the
#' NCBI taxonomy from a local database.
#'
#' @section Main functions:
#' * [make_blast_db()], [download_blast_db()], [make_unite_db()],
#'   [make_silva_db()]: BLAST databases, from a FASTA file or downloaded.
#' * [parallel_blast()]: the search; long results, one row per hit.
#' * [blast_summary()]: one row per sequence, with `Check ID`.
#' * [blast_lca()]: consensus lowest common ancestor of the tied hits.
#' * [get_tax()] and [make_taxonomy_db()]: NCBI taxonomy of the taxids, from
#'   a local SQLite database.
#' * [parallel_blast_hsps()]: every HSP, with the BLAST fields chosen.
#' * [blast_job()]: long searches in the background.
#'
#' @section External programs:
#' BLAST+ (`blastn`, `makeblastdb`, `blastdbcmd`), installed in a conda
#' environment by [install_dependencies()] (through the condathis package).
#'
#' @section Vignettes:
#' `vignette("BLASTr", package = "BLASTr")` (get started) and
#' `vignette("background-searches", package = "BLASTr")` (long searches in
#' the background).
#'
#' @keywords internal
#' @importFrom lifecycle deprecated
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL

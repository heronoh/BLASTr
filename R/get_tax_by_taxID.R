#' @title Deprecated taxonomy functions
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' `get_tax_by_taxID()`, `parallel_get_tax()` and `get_tax_by_name()`
#' queried the NCBI online (`efetch`), one taxid at a time. Use [get_tax()],
#' which reads a local taxonomy database and is much faster. These functions
#' now call it and return its table.
#'
#' @param organisms_taxIDs Taxids.
#' @param organisms_names Scientific (or other) names of taxa.
#' @param taxonomy_db See [get_tax()].
#' @param ... Ignored (arguments of the old versions).
#'
#' @returns See [get_tax()].
#' @keywords internal
#' @name deprecated-tax
NULL

#' @rdname deprecated-tax
#' @export
get_tax_by_taxID <- function(organisms_taxIDs, ..., taxonomy_db = NULL) {
  lifecycle::deprecate_warn("0.1.8", "get_tax_by_taxID()", "get_tax()")
  get_tax(organisms_taxIDs, taxonomy_db = taxonomy_db)
}

#' @rdname deprecated-tax
#' @export
parallel_get_tax <- function(organisms_taxIDs, ..., taxonomy_db = NULL) {
  lifecycle::deprecate_warn("0.1.8", "parallel_get_tax()", "get_tax()")
  get_tax(organisms_taxIDs, taxonomy_db = taxonomy_db)
}

#' @rdname deprecated-tax
#' @export
get_tax_by_name <- function(organisms_names, ..., taxonomy_db = NULL) {
  lifecycle::deprecate_warn("0.1.8", "get_tax_by_name()", "get_tax()")
  db <- resolve_taxonomy_db(taxonomy_db)
  con <- DBI::dbConnect(RSQLite::SQLite(), db, flags = RSQLite::SQLITE_RO)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  nm <- unique(as.character(organisms_names))
  ids <- DBI::dbGetQuery(
    con,
    sprintf("SELECT name, id FROM names WHERE name IN (%s)", paste(DBI::dbQuoteString(con, nm), collapse = ","))
  )
  out <- get_tax(unique(as.character(ids$id)), taxonomy_db = db)
  dplyr::left_join(tibble::tibble(name = ids$name, taxid = as.character(ids$id)), out, by = "taxid")
}

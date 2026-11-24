#' @title Taxonomy from the NCBI online, without a local database
#'
#' @description
#' Look up taxonomy in the NCBI Taxonomy database through the E-utilities
#' (`efetch`, `esearch`), over the internet. For when the local database of
#' [get_tax()] is not there yet, or for a few taxa; [get_tax()] is much faster
#' for many taxids and works offline once built.
#'
#' * `get_tax_by_taxID()`: taxonomy of taxids (or of the taxids of a BLASTr
#'   result, as in [get_tax()]).
#' * `get_tax_by_name()`: taxids and taxonomy of names, matched first as
#'   scientific names, then as any name the NCBI knows (synonyms, common
#'   names).
#'
#' Both return the same columns as [get_tax()], so the tables can be used in
#' place of each other. The taxonomy is the current one at the NCBI, which
#' may be newer than the BLAST database.
#'
#' The NCBI accepts 3 requests per second, or 10 with an API key (see
#' <https://support.nlm.nih.gov/kbArticle/?pn=KA-05317>): set it in `api_key`
#' or in the environment variable `ENTREZ_KEY` or `NCBI_API_KEY`. Taxids are
#' sent in batches of 200.
#'
#' @param organisms_taxIDs Taxids (character or numeric), or a BLASTr result.
#' @param organisms_names Names of taxa (e.g. `"Homo sapiens"`, `"Salminus"`).
#' @param ranks Ranks returned as columns, as in [get_tax()].
#' @param api_key NCBI API key (optional).
#' @param ... Ignored (arguments of older versions).
#'
#' @returns A tibble as in [get_tax()]: `taxid`, `found`, `sci_name`, `rank`,
#'   one column per rank, `lineage`, `lineage_taxids` and `merged_into`.
#'   `get_tax_by_name()` adds, before them, `name` (as given) and
#'   `matched_by` (`"scientific name"`, `"other name"` or `NA` when not
#'   found); a name with more than one taxid (a homonym) gets one row per
#'   taxid.
#'
#' @examples
#' \donttest{
#' # needs internet
#' if (requireNamespace("curl", quietly = TRUE) && curl::has_internet()) {
#'   get_tax_by_taxID(c("9606", "10090"))
#'   get_tax_by_name(c("Homo sapiens", "Salminus brasiliensis"))
#' }
#' }
#'
#' @family taxonomy functions
#' @name get_tax_online
NULL

#' @rdname get_tax_online
#' @export
get_tax_by_taxID <- function(organisms_taxIDs, ranks = tax_ranks_default, api_key = NULL, ...) {
  taxids <- extract_taxids(organisms_taxIDs)
  ranks <- as.character(ranks)
  if (length(taxids) == 0L) {
    return(online_tax_table(character(0L), list(), ranks))
  }
  api_key <- ncbi_api_key(api_key)
  parts <- split(taxids, ceiling(seq_along(taxids) / 200L))
  taxa <- unlist(lapply(parts, function(p) {
    x <- eutils_get("efetch", list(db = "taxonomy", id = paste(p, collapse = ","), retmode = "xml"), api_key)
    parse_taxa_xml(x)
  }), recursive = FALSE)
  out <- online_tax_table(taxids, taxa, ranks)
  miss <- out$taxid[!out$found]
  if (length(miss) > 0L) {
    cli::cli_warn(c(`!` = "{length(miss)} taxid{?s} not found at the NCBI: {utils::head(miss, 10)}{if (length(miss) > 10) ', ...' else ''}"))
  }
  out
}

#' @rdname get_tax_online
#' @export
get_tax_by_name <- function(organisms_names, ranks = tax_ranks_default, api_key = NULL, ...) {
  nm <- unique(trimws(as.character(organisms_names)))
  nm <- nm[!is.na(nm) & nzchar(nm)]
  api_key <- ncbi_api_key(api_key)
  search <- function(term) {
    x <- eutils_get("esearch", list(db = "taxonomy", term = term, retmax = "20"), api_key)
    xml2::xml_text(xml2::xml_find_all(x, "/eSearchResult/IdList/Id"))
  }
  hits <- lapply(nm, function(n) {
    ids <- search(paste0("\"", n, "\"[Scientific Name]"))
    by <- "scientific name"
    if (length(ids) == 0L) {
      ids <- search(paste0("\"", n, "\"[All Names]"))
      by <- "other name"
    }
    if (length(ids) == 0L) {
      return(tibble::tibble(name = n, matched_by = NA_character_, taxid = NA_character_))
    }
    tibble::tibble(name = n, matched_by = by, taxid = ids)
  })
  hits <- dplyr::bind_rows(tibble::tibble(name = character(0L), matched_by = character(0L), taxid = character(0L)), hits)
  ids <- unique(hits$taxid[!is.na(hits$taxid)])
  tax <- if (length(ids) > 0L) {
    suppressWarnings(get_tax_by_taxID(ids, ranks = ranks, api_key = api_key))
  } else {
    online_tax_table(character(0L), list(), ranks)
  }
  out <- dplyr::left_join(hits, tax, by = "taxid")
  out$found[is.na(out$found)] <- FALSE
  miss <- unique(out$name[is.na(out$matched_by)])
  if (length(miss) > 0L) {
    cli::cli_warn(c(`!` = "{length(miss)} name{?s} not found at the NCBI: {utils::head(miss, 10)}{if (length(miss) > 10) ', ...' else ''}"))
  }
  out
}

#' @rdname get_tax_online
#' @export
parallel_get_tax <- function(organisms_taxIDs, ...) {
  lifecycle::deprecate_warn("0.1.8", "parallel_get_tax()", "get_tax_by_taxID()")
  get_tax_by_taxID(organisms_taxIDs, ...)
}

ncbi_api_key <- function(api_key = NULL) {
  if (!is.null(api_key) && nzchar(api_key)) return(api_key)
  k <- Sys.getenv("ENTREZ_KEY")
  if (!nzchar(k)) k <- Sys.getenv("NCBI_API_KEY")
  if (nzchar(k)) k else NULL
}

# One E-utilities request (GET), with the NCBI rate limit and retries;
# returns the parsed XML.
eutils_get <- function(tool, query, api_key = NULL, tries = 4L) {
  if (!is.null(api_key)) query$api_key <- api_key
  query$tool <- "BLASTr"
  url <- paste0(
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/", tool, ".fcgi?",
    paste(names(query), vapply(query, utils::URLencode, character(1L), reserved = TRUE), sep = "=", collapse = "&")
  )
  wait <- if (is.null(api_key)) 0.34 else 0.1
  last <- getOption("BLASTr.eutils_last", 0)
  dt <- as.numeric(Sys.time()) - last
  if (dt < wait) Sys.sleep(wait - dt)
  tmp <- tempfile(fileext = ".xml")
  on.exit(unlink(tmp), add = TRUE)
  for (i in seq_len(tries)) {
    options(BLASTr.eutils_last = as.numeric(Sys.time()))
    ok <- tryCatch(
      {
        utils::download.file(url, tmp, quiet = TRUE, mode = "wb")
        TRUE
      },
      error = function(e) FALSE, warning = function(w) FALSE
    )
    if (ok) {
      x <- tryCatch(xml2::read_xml(tmp), error = function(e) NULL)
      if (!is.null(x)) return(x)
    }
    Sys.sleep(2^i)
  }
  cli::cli_abort(
    c(`x` = "The NCBI E-utilities ({.field {tool}}) did not answer after {tries} tries.", `i` = "Check the internet connection."),
    class = "blastr_error_eutils"
  )
}

# Taxa of an efetch taxonomy XML: one list per <Taxon> at the top level.
parse_taxa_xml <- function(x) {
  lapply(xml2::xml_find_all(x, "/TaxaSet/Taxon"), function(t) {
    txt <- function(path) xml2::xml_text(xml2::xml_find_all(t, path))
    lin <- xml2::xml_find_all(t, "./LineageEx/Taxon")
    list(
      taxid = txt("./TaxId"),
      name = txt("./ScientificName"),
      rank = txt("./Rank"),
      aka = txt("./AkaTaxIds/TaxId"),
      lin_id = xml2::xml_text(xml2::xml_find_first(lin, "./TaxId")),
      lin_name = xml2::xml_text(xml2::xml_find_first(lin, "./ScientificName")),
      lin_rank = xml2::xml_text(xml2::xml_find_first(lin, "./Rank"))
    )
  })
}

# The table of get_tax() from the parsed taxa, one row per requested taxid.
online_tax_table <- function(taxids, taxa, ranks) {
  key <- unlist(lapply(seq_along(taxa), function(i) {
    stats::setNames(rep(i, 1L + length(taxa[[i]]$aka)), c(taxa[[i]]$taxid, taxa[[i]]$aka))
  }))
  idx <- unname(key[as.character(taxids)])
  skip <- c("root", "cellular organisms")
  rows <- lapply(seq_along(taxids), function(j) {
    i <- idx[j]
    if (is.na(i)) {
      return(c(list(found = FALSE, sci_name = NA_character_, rank = NA_character_, merged_into = NA_character_),
        stats::setNames(as.list(rep(NA_character_, length(ranks))), ranks),
        list(lineage = NA_character_, lineage_taxids = NA_character_)))
    }
    t <- taxa[[i]]
    ids <- c(t$lin_id, t$taxid)
    nms <- c(t$lin_name, t$name)
    rks <- c(t$lin_rank, t$rank)
    pick <- function(k) {
      w <- which(rks == k)
      if (length(w)) nms[max(w)] else NA_character_
    }
    rk <- lapply(ranks, function(k) {
      v <- pick(k)
      if (k == "domain" && is.na(v)) v <- pick("superkingdom")
      v
    })
    keep <- !nms %in% skip
    list2 <- list(
      found = TRUE, sci_name = t$name, rank = t$rank,
      merged_into = if (taxids[j] != t$taxid) t$taxid else NA_character_
    )
    c(list2, stats::setNames(rk, ranks), list(
      lineage = if (any(keep)) paste(nms[keep], collapse = "; ") else NA_character_,
      lineage_taxids = if (any(keep)) paste(ids[keep], collapse = "; ") else NA_character_
    ))
  })
  col <- function(k, type) vapply(rows, function(r) r[[k]], type)
  tibble::tibble(
    taxid = as.character(taxids),
    found = col("found", logical(1L)),
    sci_name = col("sci_name", character(1L)),
    rank = col("rank", character(1L)),
    !!!stats::setNames(lapply(ranks, function(k) col(k, character(1L))), ranks),
    lineage = col("lineage", character(1L)),
    lineage_taxids = col("lineage_taxids", character(1L)),
    merged_into = col("merged_into", character(1L))
  )
}

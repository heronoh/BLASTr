#' @title Consensus lowest common ancestor of the best BLAST hits
#'
#' @description
#' For each sequence, the taxon shared by the hits tied with the best one
#' (`at_best`, see [parallel_blast()]): the deepest node of the NCBI
#' taxonomy that contains at least `consensus` of them. With
#' `consensus = 1` it is the strict lowest common ancestor; below 1, a few
#' mislabelled references among many tied hits (e.g. one "Picramnia" among
#' 44 *Picea*) do not push the identification up to a high rank.
#'
#' The lineages come from [get_tax()] (`lineage_taxids`), so the ancestor is
#' exact even with homonyms and clades without rank. Hits without a known
#' taxid (`0`, empty, not in the taxonomy) do not vote. A hit with several
#' taxids votes with the first one.
#'
#' The default `0.9` was chosen on sequences of known taxonomy (see
#' `dev/benchmark/R/10_rank_calibration.R`): compared with the strict
#' ancestor it resolved 3 points more genera with almost no more errors,
#' while also accepting one dissenting hit among few ("N - 1") resolved more
#' but about half of the extra ones wrong. With up to 9 tied hits, `0.9` is
#' the strict ancestor.
#'
#' @param x BLASTr results: long, as returned by [parallel_blast()], or wide
#'   (converted with [as_blast_long()]).
#' @param consensus Fraction of the tied hits (with known taxid) that the
#'   ancestor must contain, above 0.5. Defaults to `0.9`.
#' @param taxonomy The table of [get_tax()] for the taxids of `x`. `NULL`
#'   (default) calls `get_tax(x, taxonomy_db, taxdump = taxdump)`.
#' @param taxonomy_db,taxdump Passed to [get_tax()].
#'
#' @returns A tibble with one row per sequence: `Sequence`, `LCA taxid`,
#'   `LCA name`, `LCA rank` (`"genus"`, `"family"`, ...; `NA` for a node
#'   without one of the ranks of [get_tax()], e.g. a clade or a superfamily),
#'   `LCA ranked taxid`, `LCA ranked name`, `LCA ranked rank` (the ancestor
#'   itself when it has one of those ranks, otherwise the deepest node above
#'   it that has one: for tables by rank), `LCA support` (fraction of the
#'   tied hits it contains), `N tied` (hits at best) and `N tied no taxid`.
#'   `NA` for sequences without hits. Join `get_tax(result$`LCA taxid`)` for
#'   its ranks and lineage.
#'
#' @examples
#' \dontrun{
#' res <- parallel_blast(seqs, db_path)
#' lca <- blast_lca(res)
#' dplyr::left_join(lca, get_tax(lca$`LCA taxid`), by = c(`LCA taxid` = "taxid"))
#' }
#' @export
blast_lca <- function(x, consensus = 0.9, taxonomy = NULL, taxonomy_db = NULL, taxdump = NULL) {
  if (!is.numeric(consensus) || length(consensus) != 1L || consensus <= 0.5 || consensus > 1) {
    cli::cli_abort(c(`x` = "{.arg consensus} must be a number above 0.5 and at most 1."))
  }
  meta <- attr(x, "BLASTr_metadata")
  x <- as_blast_long(x)
  seqs <- unique(x$Sequence)
  tied <- x[!is.na(x$at_best) & x$at_best, c("Sequence", "staxid"), drop = FALSE]
  tied$taxid <- trimws(sub(";.*", "", as.character(tied$staxid)))
  tied$taxid[is.na(tied$taxid) | tied$taxid %in% unknown_taxids] <- NA_character_

  if (is.null(taxonomy)) {
    q <- tied[!is.na(tied$taxid), "taxid", drop = FALSE]
    names(q) <- "staxid"
    attr(q, "BLASTr_metadata") <- meta
    taxonomy <- if (nrow(q) > 0L) {
      get_tax(q, taxonomy_db = taxonomy_db, taxdump = taxdump)
    } else {
      tibble::tibble(taxid = character(0L), lineage_taxids = character(0L))
    }
  }
  k <- match(tied$taxid, taxonomy$taxid)
  lin <- strsplit(ifelse(is.na(k), NA_character_, taxonomy$lineage_taxids[k]), "; ", fixed = TRUE)
  known <- !vapply(lin, function(l) length(l) == 0L || all(is.na(l)), logical(1L))

  # node names and ranks, for the ancestors found
  by_seq <- split(seq_len(nrow(tied)), factor(tied$Sequence, levels = seqs))
  res <- lapply(seqs, function(s) {
    i <- by_seq[[s]]
    n_tied <- length(i)
    i <- i[known[i]]
    n <- length(i)
    if (n == 0L) {
      return(list(taxid = NA_character_, support = NA_real_, n_tied = n_tied, n_unknown = n_tied))
    }
    l <- lin[i]
    counts <- table(unlist(lapply(l, unique)))
    ok <- names(counts)[counts >= ceiling(consensus * n - 1e-9)]
    if (length(ok) == 0L) {
      return(list(taxid = NA_character_, support = NA_real_, n_tied = n_tied, n_unknown = n_tied - n))
    }
    # nodes held by more than half of the hits lie on one path: the deepest
    depth <- vapply(ok, function(t) max(vapply(l, function(p) {
      m <- match(t, p)
      if (is.na(m)) 0L else m
    }, integer(1L))), integer(1L))
    best <- ok[which.max(depth)]
    list(taxid = best, support = as.numeric(counts[[best]]) / n, n_tied = n_tied, n_unknown = n_tied - n)
  })
  out <- tibble::tibble(
    Sequence = seqs,
    `LCA taxid` = vapply(res, `[[`, character(1L), "taxid"),
    `LCA support` = round(vapply(res, `[[`, numeric(1L), "support"), 4),
    `N tied` = vapply(res, function(r) as.integer(r$n_tied), integer(1L)),
    `N tied no taxid` = vapply(res, function(r) as.integer(r$n_unknown), integer(1L))
  )
  # name and rank of each ancestor, from the lineage of one of its hits
  rank_cols <- intersect(c(tax_ranks_default, "superkingdom"), names(taxonomy))
  # and the deepest node at or above it with one of those ranks
  none <- rep(NA_character_, 5L)
  info <- lapply(seq_along(res), function(r) {
    t <- res[[r]]$taxid
    if (is.na(t)) return(none)
    i <- by_seq[[seqs[r]]]
    i <- i[known[i]]
    for (h in i) {
      p <- match(t, lin[[h]])
      if (!is.na(p)) {
        row <- taxonomy[k[h], , drop = FALSE]
        nms <- strsplit(row$lineage, "; ", fixed = TRUE)[[1L]]
        rank_of <- function(q) {
          rk <- rank_cols[vapply(rank_cols, function(cc) identical(row[[cc]], nms[q]), logical(1L))]
          if (q == length(nms) && "rank" %in% names(row)) rk <- row$rank
          if (length(rk)) rk[1L] else NA_character_
        }
        rk <- rank_of(p)
        q <- p
        while (q > 0L && !rank_of(q) %in% rank_cols) q <- q - 1L
        if (q == 0L) return(c(nms[p], rk, NA_character_, NA_character_, NA_character_))
        return(c(nms[p], rk, lin[[h]][q], nms[q], rank_of(q)))
      }
    }
    none
  })
  out$`LCA name` <- vapply(info, `[`, character(1L), 1L)
  out$`LCA rank` <- vapply(info, `[`, character(1L), 2L)
  out$`LCA ranked taxid` <- vapply(info, `[`, character(1L), 3L)
  out$`LCA ranked name` <- vapply(info, `[`, character(1L), 4L)
  out$`LCA ranked rank` <- vapply(info, `[`, character(1L), 5L)
  out$`N tied`[out$`N tied` == 0L] <- NA_integer_
  out[, c(
    "Sequence", "LCA taxid", "LCA name", "LCA rank", "LCA ranked taxid", "LCA ranked name", "LCA ranked rank",
    "LCA support", "N tied", "N tied no taxid"
  )]
}

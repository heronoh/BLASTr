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
#' taxid (`0`, empty, not in the taxonomy) do not vote, nor do the hits
#' flagged as not `informative` (see [uninformative_references()]). A hit
#' with several taxids votes with the first one.
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
#' @param named_first If `TRUE` (default), the tied hits named only above
#'   genus (e.g. "Bacterium", "alpha proteobacterium XYZ", "Entolomataceae
#'   sp.") do not vote when a named species is tied with them: their lineage
#'   is right but shallow, and would pull the ancestor up. Hits in open
#'   nomenclature within a genus ("Picea sp.", "Carex cf. flava") still
#'   vote: they are other species of that genus (or of another one), and
#'   keep the ancestor at the genus (or above) instead of on the named
#'   species. They count in `N tied`; `LCA support` is among the voters.
#' @param taxonomy The table of [get_tax()] for the taxids of `x`. `NULL`
#'   (default) calls `get_tax(x, taxonomy_db, taxdump = taxdump)`. Its rank
#'   columns define the ranks of `LCA ranked ...`: keep only some (e.g.
#'   `species`, `genus`, `family`, `order`) to place the ranked ancestor on
#'   them.
#' @param taxonomy_db,taxdump Passed to [get_tax()].
#'
#' @returns A tibble with one row per sequence: `Sequence`, `LCA taxid`,
#'   `LCA name`, `LCA rank` (`"genus"`, `"family"`, ...; `NA` for a node
#'   without one of the ranks of [get_tax()], e.g. a clade or a superfamily),
#'   `LCA ranked taxid`, `LCA ranked name`, `LCA ranked rank` (the ancestor
#'   itself when it has one of those ranks, otherwise the deepest node above
#'   it that has one; a species in open nomenclature, e.g. "Picea sp. X",
#'   goes to its genus: for tables by rank), `LCA support` (fraction of the
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
blast_lca <- function(x, consensus = 0.9, named_first = TRUE, taxonomy = NULL, taxonomy_db = NULL, taxdump = NULL) {
  if (!is.numeric(consensus) || length(consensus) != 1L || consensus <= 0.5 || consensus > 1) {
    cli::cli_abort(c(`x` = "{.arg consensus} must be a number above 0.5 and at most 1."))
  }
  meta <- attr(x, "BLASTr_metadata")
  x <- as_blast_long(x)
  seqs <- unique(x$Sequence)
  voter <- !is.na(x$at_best) & x$at_best
  if (!is.null(x[["informative"]])) voter <- voter & !(x[["informative"]] %in% FALSE)
  tied <- x[voter, c("Sequence", "staxid"), drop = FALSE]
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
  lin <- strsplit(as.character(ifelse(is.na(k), NA_character_, taxonomy$lineage_taxids[k])), "; ", fixed = TRUE)
  known <- !vapply(lin, function(l) length(l) == 0L || all(is.na(l)), logical(1L))
  has_tax <- known
  # hits named only above genus ("Bacterium", "Entolomataceae sp.") do not
  # vote when a named species is tied with them; "Picea sp." still votes
  # (another species of Picea)
  if (isTRUE(named_first) && nrow(tied) > 0L) {
    sp <- if (!is.null(taxonomy[["species"]])) {
      taxonomy$species[k]
    } else if (!is.null(taxonomy[["sci_name"]]) && !is.null(taxonomy[["rank"]])) {
      ifelse(taxonomy$rank[k] %in% "species", taxonomy$sci_name[k], NA_character_)
    } else {
      rep("named", length(k))
    }
    named <- known & !is.na(sp) & !grepl(open_nomenclature_pattern, sp, perl = TRUE)
    in_genus <- if (!is.null(taxonomy[["genus"]])) !is.na(taxonomy$genus[k]) else rep(TRUE, length(k))
    has_named <- stats::ave(named, factor(tied$Sequence), FUN = any)
    known <- known & (named | in_genus | !has_named)
  }

  # node names and ranks, for the ancestors found
  by_seq <- split(seq_len(nrow(tied)), factor(tied$Sequence, levels = seqs))
  res <- lapply(seqs, function(s) {
    i <- by_seq[[s]]
    n_tied <- length(i)
    n_unknown <- sum(!has_tax[i])
    i <- i[known[i]]
    n <- length(i)
    if (n == 0L) {
      return(list(taxid = NA_character_, support = NA_real_, n_tied = n_tied, n_unknown = n_unknown))
    }
    l <- lin[i]
    counts <- table(unlist(lapply(l, unique)))
    ok <- names(counts)[counts >= ceiling(consensus * n - 1e-9)]
    if (length(ok) == 0L) {
      return(list(taxid = NA_character_, support = NA_real_, n_tied = n_tied, n_unknown = n_unknown))
    }
    # nodes held by more than half of the hits lie on one path: the deepest
    depth <- vapply(ok, function(t) max(vapply(l, function(p) {
      m <- match(t, p)
      if (is.na(m)) 0L else m
    }, integer(1L))), integer(1L))
    best <- ok[which.max(depth)]
    list(taxid = best, support = as.numeric(counts[[best]]) / n, n_tied = n_tied, n_unknown = n_unknown)
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
        # a species in open nomenclature ("Picea sp. X") is placed on its genus
        while (q > 0L && (!rank_of(q) %in% rank_cols ||
          (rank_of(q) %in% "species" && grepl(open_nomenclature_pattern, nms[q], perl = TRUE)))) {
          q <- q - 1L
        }
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

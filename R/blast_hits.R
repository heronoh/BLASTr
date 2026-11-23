# Columns of the BLAST tabular output (`-outfmt`) and their names in the
# BLASTr results.
blast_out_cols <- c(
  query = "qaccver",
  subject = "saccver",
  identity = "pident",
  length = "length",
  mismatches = "mismatch",
  gaps = "gapopen",
  `query start` = "qstart",
  `query end` = "qend",
  `subject start` = "sstart",
  `subject end` = "send",
  `e-value` = "evalue",
  bitscore = "bitscore",
  qcovhsp = "qcovhsp",
  nident = "nident",
  `subject strand` = "sstrand",
  staxid = "staxid",
  `subject header` = "stitle"
)

blast_outfmt <- function() {
  paste("6", paste(blast_out_cols, collapse = " "))
}

# Fields of each hit in the wide results, in display order.
hit_fields <- c(
  "subject header", "Similarity score", "subject", "identity", "length",
  "mismatches", "gaps", "query start", "query end", "subject start",
  "subject end", "e-value", "bitscore", "qcovhsp", "staxid", "nident",
  "subject strand", "blast_rank", "low_qcov", "informative",
  "uninformative_category"
)

# Columns of the long results, in order.
long_cols <- c(
  "Sequence", "rank", "tie", "at_best", "low_qcov", "informative",
  "uninformative_category", "blast_rank", "subject",
  "subject header", "staxid", "identity", "Similarity score", "qcovhsp",
  "nident", "length", "mismatches", "gaps", "query start", "query end",
  "subject start", "subject end", "subject strand", "e-value", "bitscore"
)

hit_num_cols <- c(
  "identity", "length", "mismatches", "gaps", "query start", "query end",
  "subject start", "subject end", "e-value", "bitscore", "qcovhsp", "nident",
  "Similarity score"
)

# Columns that summarise the hits of each sequence.
hit_summary_cols <- c(
  "Check ID", "N hits", "N taxids", "N hits same taxid as 1st",
  "N hits at best", "N taxids at best", "Taxids at best", "At best truncated",
  "Best low qcov", "N uninformative", "Best uninformative score"
)

#' @title Similarity score of BLAST hits
#'
#' @description
#' Combines identity and query coverage in one value, weighting identity
#' ten times more: `(10 * identity + qcovhsp) / 11`. It was called
#' "BLASTn pseudo-score" in earlier versions of the metabarcoding pipeline.
#'
#' @param identity Percent identity of the hit (`pident`).
#' @param qcovhsp Percent of the query covered by the hit (`qcovhsp`).
#'
#' @returns Numeric vector, same length as the inputs.
#'
#' @examples
#' similarity_score(identity = c(100, 98.5), qcovhsp = c(100, 85))
#' @export
similarity_score <- function(identity, qcovhsp) {
  (10 * as.numeric(identity) + as.numeric(qcovhsp)) / 11
}


# Value used to rank the hits.
rank_key <- function(h, rank_by) {
  switch(rank_by,
    similarity = similarity_score(h[["identity"]], h[["qcovhsp"]]),
    identity = as.numeric(h[["identity"]]),
    bitscore = as.numeric(h[["bitscore"]])
  )
}

unknown_taxids <- c("", "0", "N/A", "NA")

# Order the hits of each query by BLASTr's criterion and add `rank`, `tie`,
# `at_best`, `.pos` (position in that order) and `.all_pos` (position among
# all the hits, informative or not). With an `informative` column, the
# hits are ranked among the informative ones (or, when a query has none,
# among all); the others follow, with `rank` and `tie` NA and `at_best`
# FALSE.
rank_blast_hits <- function(hits, rank_by = "similarity", tolerance = 0.5) {
  if (nrow(hits) == 0L) {
    hits$rank <- integer(0L)
    hits$tie <- logical(0L)
    hits$at_best <- logical(0L)
    hits$.pos <- integer(0L)
    hits$.all_pos <- integer(0L)
    return(hits)
  }
  num <- function(x) if (is.null(x)) rep(NA_real_, nrow(hits)) else as.numeric(x)
  key <- round(rank_key(hits, rank_by), 6)
  inf <- if (is.null(hits[["informative"]])) rep(TRUE, nrow(hits)) else !(hits[["informative"]] %in% FALSE)
  pool <- inf | !stats::ave(inf, hits$query, FUN = any)
  qo <- match(hits$query, unique(hits$query))
  tiebreak <- list(-num(hits[["identity"]]), -num(hits[["qcovhsp"]]), -num(hits[["bitscore"]]), num(hits[["blast_rank"]]))
  all_o <- do.call(order, c(list(qo, -key), tiebreak, list(na.last = TRUE)))
  all_pos <- integer(nrow(hits))
  all_pos[all_o] <- as.integer(stats::ave(seq_along(all_o), qo[all_o], FUN = seq_along))
  o <- do.call(order, c(list(qo, !pool, -key), tiebreak, list(na.last = TRUE)))
  hits <- hits[o, , drop = FALSE]
  key <- key[o]
  pool <- pool[o]
  pos <- stats::ave(seq_len(nrow(hits)), hits$query, FUN = seq_along)
  grp <- paste(hits$query, pool, key, sep = "\r")
  hits$rank <- ifelse(pool, as.integer(stats::ave(pos, grp, FUN = min)), NA_integer_)
  hits$tie <- ifelse(pool, stats::ave(pos, grp, FUN = length) > 1L, NA)
  best <- stats::ave(ifelse(pool, key, -Inf), hits$query, FUN = function(k) max(k, na.rm = TRUE))
  hits$at_best <- pool & !is.na(key) & key >= best - tolerance - 1e-9
  hits$.pos <- as.integer(pos)
  hits$.all_pos <- all_pos[o]
  hits
}

#' @title Check the agreement of the best BLAST hits
#'
#' @description
#' For each query sequence, orders its hits by `rank_by`, decides whether
#' the best ones agree on the taxon (`Check ID`) and summarises all of them.
#'
#' The hits are ordered by `rank_by` (Similarity score, identity or
#' bitscore) and then by identity, coverage, bitscore and BLAST order. Hits
#' within `tolerance` of the best value are "at best": they are tied
#' candidates for the identification.
#'
#' `Check ID` looks at the first `check_n` hits and compares only the ones at
#' best: a hit clearly worse than the best one is most likely another taxon,
#' not a competing identification.
#'
#' * `"best ID"`: all the hits at best among the first `check_n` have the
#'   taxid of the first hit.
#' * `"check ID"`: some of them has another taxid or an unknown one (`0` or
#'   empty, as in databases built without a taxid map).
#' * `NA`: no hits.
#'
#' The other columns use all the hits of the sequence: `N hits`, `N taxids`
#' (distinct known taxids), `N hits same taxid as 1st`, `N hits at best`,
#' `N taxids at best`, `Taxids at best` (the known taxids of the hits at
#' best, in rank order, separated by `";"`, for a lowest common ancestor) and
#' `At best truncated` (`TRUE` when the last subject BLAST returned is at
#' best: more tied subjects may exist beyond `search_alignments`; `NA` when
#' `search_alignments` is unknown).
#'
#' With an `informative` column (see [uninformative_references()]), all of
#' this is computed among the informative hits, or among all the hits when
#' a query has no informative one (and then `Check ID` is `"check ID"`).
#' `N uninformative` counts the other hits and `Best uninformative score`
#' is the best value of `rank_by` among them: above the best informative
#' hit, the query is closer to an uncultured or environmental sequence than
#' to any named organism (a sign of a taxon missing from the database).
#'
#' @param hits Data frame with one row per hit and the columns `query`,
#'   `identity`, `staxid` and, as needed by `rank_by`, `qcovhsp` or
#'   `bitscore`. The BLAST order is taken from `blast_rank` or `rank`, or
#'   else from the row order. Optional: `informative` (logical).
#' @param check_n Number of best hits used by `Check ID`.
#' @param tolerance Difference from the best value, in the units of
#'   `rank_by` (points of Similarity score or identity, bits of bitscore),
#'   under which hits are tied.
#' @param rank_by Criterion to order the hits: `"similarity"` (Similarity
#'   score, see [similarity_score()]), `"identity"` or `"bitscore"`.
#' @param search_alignments Number of subjects BLAST was asked for, to
#'   compute `At best truncated`.
#'
#' @returns A tibble with `query` and the summary columns.
#'
#' @examples
#' hits <- data.frame(
#'   query = "q1", identity = c(100, 99, 98.7), qcovhsp = 100,
#'   staxid = c("9606", "9598", "9597")
#' )
#' check_hits(hits)
#' @export
check_hits <- function(
  hits,
  check_n = 10L,
  tolerance = 0.5,
  rank_by = c("similarity", "identity", "bitscore"),
  search_alignments = NA_integer_
) {
  rank_by <- rlang::arg_match(rank_by)
  hits <- as.data.frame(hits, check.names = FALSE)
  if (is.null(hits$blast_rank)) {
    hits$blast_rank <- if (!is.null(hits$rank)) hits$rank else stats::ave(seq_len(nrow(hits)), hits$query, FUN = seq_along)
  }
  hits$rank <- NULL
  hits <- rank_blast_hits(hits, rank_by = rank_by, tolerance = tolerance)
  key <- round(rank_key(hits, rank_by), 6)
  inf <- if (is.null(hits[["informative"]])) rep(TRUE, nrow(hits)) else !(hits[["informative"]] %in% FALSE)
  pool <- !is.na(hits$rank)
  tx <- as.character(hits$staxid)
  known <- !is.na(tx) & !tx %in% unknown_taxids
  first_tx <- stats::ave(tx, hits$query, FUN = function(x) x[1L])
  first_known <- stats::ave(known, hits$query, FUN = function(x) x[1L])
  top <- pool & hits$.pos <= check_n
  best_top <- stats::ave(ifelse(top, key, -Inf), hits$query, FUN = max)
  cmp <- top & key >= best_top - tolerance - 1e-9
  d <- data.frame(
    query = hits$query, tx = tx, known = known & pool, cmp = cmp, pool = pool, inf = inf, key = key,
    same1 = known & pool & first_known & tx == first_tx,
    at_best = hits$at_best, blast_rank = hits$blast_rank
  )
  q <- unique(hits$query)
  by_q <- split(d, factor(d$query, levels = q))
  out <- tibble::tibble(
    query = q,
    `Check ID` = vapply(by_q, function(x) {
      c1 <- x[x$cmp, , drop = FALSE]
      if (any(c1$inf) && all(c1$known) && all(c1$tx == x$tx[1L])) "best ID" else "check ID"
    }, character(1L), USE.NAMES = FALSE),
    `N hits` = vapply(by_q, nrow, integer(1L), USE.NAMES = FALSE),
    `N taxids` = vapply(by_q, function(x) length(unique(x$tx[x$known])), integer(1L), USE.NAMES = FALSE),
    `N hits same taxid as 1st` = vapply(by_q, function(x) {
      if (x$known[1L]) sum(x$same1) else NA_integer_
    }, integer(1L), USE.NAMES = FALSE),
    `N hits at best` = vapply(by_q, function(x) sum(x$at_best), integer(1L), USE.NAMES = FALSE),
    `N taxids at best` = vapply(by_q, function(x) length(unique(x$tx[x$known & x$at_best])), integer(1L), USE.NAMES = FALSE),
    `Taxids at best` = vapply(by_q, function(x) {
      t <- unique(trimws(unlist(strsplit(x$tx[x$known & x$at_best], ";", fixed = TRUE))))
      t <- t[!t %in% unknown_taxids]
      if (length(t)) paste(t, collapse = ";") else NA_character_
    }, character(1L), USE.NAMES = FALSE),
    `At best truncated` = vapply(by_q, function(x) {
      if (is.na(search_alignments)) NA else any(x$at_best & x$blast_rank >= search_alignments)
    }, logical(1L), USE.NAMES = FALSE),
    `N uninformative` = vapply(by_q, function(x) sum(!x$inf), integer(1L), USE.NAMES = FALSE),
    `Best uninformative score` = vapply(by_q, function(x) {
      k <- x$key[!x$inf & !is.na(x$key)]
      if (length(k)) max(k) else NA_real_
    }, numeric(1L), USE.NAMES = FALSE)
  )
  out
}

# Read the tabular outputs of a job into one row per hit (`query`,
# `blast_rank`, the hit fields and `Similarity score`).
read_blast_hits <- function(out_files) {
  out_files <- out_files[file.exists(out_files) & file.size(out_files) > 0]
  empty <- tibble::as_tibble(stats::setNames(
    rep(list(character(0L)), length(blast_out_cols)),
    names(blast_out_cols)
  ))
  hits <- lapply(out_files, function(f) {
    readr::read_tsv(
      f,
      col_names = names(blast_out_cols),
      col_types = readr::cols(.default = readr::col_character()),
      quote = "",
      comment = "",
      na = character(0L),
      trim_ws = TRUE,
      progress = FALSE
    )
  })
  hits <- dplyr::bind_rows(c(list(empty), hits))
  for (k in intersect(hit_num_cols, names(hits))) {
    hits[[k]] <- suppressWarnings(as.numeric(hits[[k]]))
  }
  # BLAST lists the hits of each query from best to worst (by its score)
  hits$blast_rank <- as.integer(stats::ave(seq_len(nrow(hits)), hits$query, FUN = seq_along))
  hits$`Similarity score` <- similarity_score(hits$identity, hits$qcovhsp)
  hits
}

# Long results: one row per kept hit, ranked; one row with NA hit fields for
# each sequence without hits. With a table of `uninformative` references,
# the hits are flagged (`informative`, `uninformative_category`) and ranked
# among the informative ones. Kept: the first `num_alignments` ranked hits,
# all the hits at best, and the uninformative hits among the first
# `num_alignments` of all.
blast_long <- function(hits, seq_ids, seqs, num_alignments = Inf, rank_by = "similarity", tolerance = 0.5,
                       perc_qcov_hsp = NA, uninformative = NULL) {
  if (!is.null(uninformative)) {
    hits$uninformative_category <- uninformative_category(hits$subject, hits[["subject header"]], uninformative)
    hits$informative <- is.na(hits$uninformative_category)
  } else if (is.null(hits[["informative"]])) {
    # not evaluated: ranked as informative
    hits$informative <- rep(NA, nrow(hits))
  }
  hits <- rank_blast_hits(hits, rank_by = rank_by, tolerance = tolerance)
  hits$low_qcov <- if (is.na(perc_qcov_hsp)) NA else as.numeric(hits[["qcovhsp"]]) < perc_qcov_hsp
  pool <- !is.na(hits$rank)
  hits <- hits[(pool & (hits$.pos <= num_alignments | hits$at_best)) | (!pool & hits$.all_pos <= num_alignments), , drop = FALSE]
  hits$Sequence <- seqs[match(hits$query, seq_ids)]
  none <- seq_ids[!seq_ids %in% hits$query]
  if (length(none) > 0L) {
    hits <- dplyr::bind_rows(hits, tibble::tibble(query = none, Sequence = seqs[match(none, seq_ids)]))
  }
  hits <- hits[order(match(hits$query, seq_ids), hits$.pos, na.last = TRUE), , drop = FALSE]
  for (k in setdiff(long_cols, names(hits))) {
    hits[[k]] <- if (k %in% hit_num_cols) NA_real_ else NA_character_
  }
  hits$rank <- as.integer(hits$rank)
  hits$blast_rank <- as.integer(hits$blast_rank)
  hits$tie <- as.logical(hits$tie)
  hits$at_best <- as.logical(hits$at_best)
  hits$low_qcov <- as.logical(hits$low_qcov)
  hits$informative <- as.logical(hits$informative)
  hits$uninformative_category <- as.character(hits$uninformative_category)
  tibble::as_tibble(hits[, long_cols])
}

#' @title Summary of the BLAST hits of each sequence
#'
#' @description
#' One row per sequence of a long [parallel_blast()] result (or of
#' [as_blast_long()]): `Check ID` and the other summary columns described in
#' [check_hits()], and `Best low qcov` (the best hit covers less of the
#' query than `perc_qcov_hsp`: the identification rests on a partial
#' alignment). Sequences without hits have `N hits = 0` and `NA` elsewhere.
#'
#' @param x Long BLASTr results.
#' @param check_n,tolerance,rank_by,search_alignments As in [check_hits()].
#'   `NULL` (default) uses the values of the search, kept in the attribute
#'   `BLASTr_metadata`, or else `10`, `0.5`, `"similarity"` and `NA`.
#'
#' @returns A tibble with `Sequence` and the summary columns.
#'
#' @examples
#' \dontrun{
#' res <- parallel_blast(seqs, db_path)
#' dplyr::left_join(blast_summary(res), res[res$rank %in% 1L, ], by = "Sequence")
#' }
#' @export
blast_summary <- function(x, check_n = NULL, tolerance = NULL, rank_by = NULL, search_alignments = NULL) {
  meta <- attr(x, "BLASTr_metadata")
  pick <- function(v, m, default) {
    if (!is.null(v)) return(v)
    if (!is.null(m) && length(m) == 1L && !is.na(m)) return(m)
    default
  }
  check_n <- pick(check_n, meta$check_n, 10L)
  tolerance <- pick(tolerance, meta$check_tolerance, 0.5)
  rank_by <- pick(rank_by, meta$rank_by, "similarity")
  search_alignments <- pick(search_alignments, meta$search_alignments, NA_integer_)
  x <- as_blast_long(x, rank_by = rank_by, tolerance = tolerance)
  h <- x[!is.na(x$blast_rank), , drop = FALSE]
  h$query <- h$Sequence
  s <- check_hits(h, check_n = check_n, tolerance = tolerance, rank_by = rank_by, search_alignments = search_alignments)
  names(s)[names(s) == "query"] <- "Sequence"
  out <- dplyr::left_join(tibble::tibble(Sequence = unique(x$Sequence)), s, by = "Sequence")
  out$`N hits`[is.na(out$`N hits`)] <- 0L
  out$`N uninformative`[is.na(out$`N uninformative`)] <- 0L
  # coverage of the best hit, as ranked in x (the first row of each sequence)
  first <- x[!duplicated(x$Sequence), , drop = FALSE]
  out$`Best low qcov` <- if (is.null(first[["low_qcov"]])) NA else as.logical(first$low_qcov[match(out$Sequence, first$Sequence)])
  out
}

# Wide results: summary, then the first `num_alignments` hits side by side,
# in rank order.
widen_blast_hits <- function(long, summary, num_alignments) {
  out <- summary
  h <- long[!is.na(long$blast_rank), , drop = FALSE]
  pos <- stats::ave(seq_len(nrow(h)), h$Sequence, FUN = seq_along)
  for (r in seq_len(num_alignments)) {
    hr <- h[pos == r, , drop = FALSE]
    i <- match(out$Sequence, hr$Sequence)
    for (f in hit_fields) {
      out[[paste0(r, "_", f)]] <- hr[[f]][i]
    }
  }
  out
}

#' @title Convert BLASTr results to the long format
#'
#' @description
#' Turns wide BLASTr results (one row per sequence, hits `1_` to `n_`, from
#' any version, including `k_indentity`) into the long format of
#' [parallel_blast()], ranking the hits with the same rules, so that old
#' results can be read with [blast_summary()] and [get_tax()] without
#' running BLAST again. The BLAST order is taken from `k_blast_rank` or, in
#' older results, from `k`. Long results are returned as they are.
#'
#' @param x BLASTr results, wide or long. The sequence is in `Sequence` (or
#'   `sequence`).
#' @param rank_by,tolerance As in [check_hits()].
#' @param perc_qcov_hsp Coverage under which a hit is flagged `low_qcov`.
#' @param uninformative `NULL` (default) keeps the flags `informative` of
#'   `x` (`NA`, not evaluated, when `x` has none; ranked as informative); a
#'   table of
#'   [uninformative_references()] flags the hits again with it, and ranks
#'   them among the informative ones, as [parallel_blast()] does.
#'
#' @returns Long results (see [parallel_blast()]), with all the hits of `x`.
#'
#' @examples
#' wide <- data.frame(
#'   Sequence = "ACGT", `1_subject` = "A1", `1_indentity` = 99, `1_qcovhsp` = 100,
#'   `1_staxid` = "9606", `2_subject` = "A2", `2_indentity` = 100,
#'   `2_qcovhsp` = 100, `2_staxid` = "9598", check.names = FALSE
#' )
#' as_blast_long(wide)
#' @export
as_blast_long <- function(x, rank_by = c("similarity", "identity", "bitscore"), tolerance = 0.5, perc_qcov_hsp = 80,
                          uninformative = NULL) {
  rank_by <- rlang::arg_match(rank_by)
  meta <- attr(x, "BLASTr_metadata")
  nm <- names(x)
  if (!"Sequence" %in% nm && "sequence" %in% nm) names(x)[nm == "sequence"] <- "Sequence"
  if (!"Sequence" %in% names(x)) {
    cli::cli_abort(c(`x` = "No {.field Sequence} column in {.arg x}."), class = "blastr_not_blast_result_error")
  }
  if (all(c("rank", "blast_rank") %in% names(x)) && !any(grepl("^[0-9]+_", names(x)))) {
    if (is.null(uninformative)) {
      return(x)
    }
    # flag the hits of a long result again, with this table
    seqs <- unique(as.character(x$Sequence))
    hits <- x[!is.na(x$blast_rank), setdiff(names(x), c("rank", "tie", "at_best", "informative", "uninformative_category")), drop = FALSE]
    hits$query <- as.character(hits$Sequence)
    out <- blast_long(
      hits, seq_ids = seqs, seqs = seqs, num_alignments = Inf, rank_by = rank_by, tolerance = tolerance,
      perc_qcov_hsp = perc_qcov_hsp, uninformative = uninformative
    )
    attr(out, "BLASTr_metadata") <- meta
    return(out)
  }
  names(x) <- sub("^([0-9]+)_indentity$", "\\1_identity", names(x))
  ks <- sort(unique(as.integer(sub("_.*", "", grep("^[0-9]+_", names(x), value = TRUE)))))
  if (length(ks) == 0L) {
    cli::cli_abort(c(`x` = "{.arg x} has no hit columns ({.field 1_subject}, ...)."), class = "blastr_not_blast_result_error")
  }
  seqs <- unique(as.character(x$Sequence))
  hits <- dplyr::bind_rows(lapply(ks, function(k) {
    cols <- grep(paste0("^", k, "_"), names(x), value = TRUE)
    h <- as.data.frame(x[, cols, drop = FALSE], check.names = FALSE)
    names(h) <- sub("^[0-9]+_", "", cols)
    h[] <- lapply(names(h), function(f) {
      if (f %in% hit_num_cols || f == "blast_rank") {
        suppressWarnings(as.numeric(h[[f]]))
      } else if (f %in% c("informative", "low_qcov")) {
        as.logical(h[[f]])
      } else {
        as.character(h[[f]])
      }
    })
    if (is.null(h$blast_rank)) h$blast_rank <- k
    tibble::as_tibble(cbind(query = as.character(x$Sequence), h))
  }))
  has <- function(f) if (is.null(hits[[f]])) rep(FALSE, nrow(hits)) else !is.na(hits[[f]]) & hits[[f]] != ""
  hits <- hits[has("subject") | has("subject header"), , drop = FALSE]
  if (is.null(hits[["Similarity score"]]) || all(is.na(hits[["Similarity score"]]))) {
    hits$`Similarity score` <- similarity_score(hits$identity, hits$qcovhsp)
  }
  out <- blast_long(
    hits, seq_ids = seqs, seqs = seqs, num_alignments = Inf, rank_by = rank_by, tolerance = tolerance,
    perc_qcov_hsp = perc_qcov_hsp, uninformative = uninformative
  )
  attr(out, "BLASTr_metadata") <- meta
  out
}

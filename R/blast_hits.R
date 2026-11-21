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
  "subject end", "e-value", "bitscore", "qcovhsp", "staxid"
)

# Columns that summarise the hits of each sequence.
hit_summary_cols <- c(
  "Check ID", "N hits", "N taxids", "N hits same taxid as 1st",
  "N hits at best identity", "N taxids at best identity"
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

#' @title Check the agreement of the best BLAST hits
#'
#' @description
#' For each query sequence, decides whether its best hits agree on the
#' taxon (`Check ID`) and summarises all of its hits.
#'
#' `Check ID` looks at the first `check_n` hits. Only the hits whose identity
#' is within `tolerance` percentage points of the best identity among them
#' are compared: a hit clearly less identical than the best one is most
#' likely another taxon, not a competing identification. Coverage is not
#' used here, so a good reference that covers only part of the query does
#' not hide or create a disagreement.
#'
#' * `"best ID"`: all the compared hits have the taxid of the first hit.
#' * `"check ID"`: some compared hit has another taxid or an unknown one
#'   (`0` or empty, as in databases built without a taxid map), or the first
#'   hit is not among the most identical ones.
#' * `NA`: no hits.
#'
#' The other columns use all the hits of the sequence: `N hits`, `N taxids`
#' (distinct known taxids), `N hits same taxid as 1st`, `N hits at best
#' identity` (hits within `tolerance` of the best identity of all hits) and
#' `N taxids at best identity` (distinct known taxids among those hits).
#'
#' Together they tell how well a `"best ID"` is supported: `N hits at best
#' identity = 1` means a single reference at the best identity (the others,
#' if any, are clearly less identical); several hits at best identity with
#' `N taxids at best identity = 1` mean several references agreeing; more
#' than one taxid at best identity means a tie between taxa. Sequences with
#' fewer hits than `check_n` are checked with the hits they have.
#'
#' @param hits Data frame with one row per hit and the columns `query`,
#'   `rank`, `identity` and `staxid`.
#' @param check_n Number of best hits used by `Check ID`.
#' @param tolerance Identity difference (percentage points) under which two
#'   hits are considered tied.
#'
#' @returns A tibble with `query` and the summary columns.
#'
#' @examples
#' hits <- data.frame(
#'   query = "q1", rank = 1:3,
#'   identity = c(100, 99, 98.7), staxid = c("9606", "9598", "9597")
#' )
#' check_hits(hits)
#' @export
check_hits <- function(hits, check_n = 10L, tolerance = 0.5) {
  hits <- hits[order(hits$query, hits$rank), , drop = FALSE]
  eps <- 1e-9
  by_query <- split(hits, factor(hits$query, levels = unique(hits$query)))
  out <- lapply(by_query, function(h) {
    id <- as.numeric(h$identity)
    tx <- as.character(h$staxid)
    known <- !is.na(tx) & !tx %in% c("", "0", "N/A")
    top <- seq_len(min(check_n, nrow(h)))
    tied_top <- top[id[top] >= max(id[top]) - tolerance - eps]
    check <- if (1L %in% tied_top && all(known[tied_top]) &&
      all(tx[tied_top] == tx[1L])) "best ID" else "check ID"
    tibble::tibble(
      query = h$query[1L],
      `Check ID` = check,
      `N hits` = nrow(h),
      `N taxids` = length(unique(tx[known])),
      `N hits same taxid as 1st` = if (known[1L]) sum(tx == tx[1L] & known) else NA_integer_,
      `N hits at best identity` = sum(id >= max(id) - tolerance - eps),
      `N taxids at best identity` = length(unique(tx[known & id >= max(id) - tolerance - eps]))
    )
  })
  if (length(out) == 0L) {
    return(tibble::tibble(
      query = character(0L),
      `Check ID` = character(0L),
      `N hits` = integer(0L),
      `N taxids` = integer(0L),
      `N hits same taxid as 1st` = integer(0L),
      `N hits at best identity` = integer(0L),
      `N taxids at best identity` = integer(0L)
    ))
  }
  dplyr::bind_rows(out)
}

# Read the tabular outputs of a job into one row per hit (`query`, `rank`,
# the hit fields and `Similarity score`).
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
  num <- c(
    "identity", "length", "mismatches", "gaps", "query start", "query end",
    "subject start", "subject end", "e-value", "bitscore", "qcovhsp"
  )
  for (k in num) {
    hits[[k]] <- suppressWarnings(as.numeric(hits[[k]]))
  }
  # BLAST lists the hits of each query from best to worst
  hits$rank <- stats::ave(seq_len(nrow(hits)), hits$query, FUN = seq_along)
  hits$`Similarity score` <- similarity_score(hits$identity, hits$qcovhsp)
  hits
}

# One row per query sequence: summary columns, then hits 1..n side by side.
widen_blast_hits <- function(hits, seq_ids, seqs, num_alignments, check_n, tolerance) {
  hits <- hits[hits$rank <= num_alignments, , drop = FALSE]
  summ <- check_hits(hits, check_n = check_n, tolerance = tolerance)
  out <- tibble::tibble(query = seq_ids, Sequence = seqs)
  out <- dplyr::left_join(out, summ, by = "query")
  out$`N hits`[is.na(out$`N hits`)] <- 0L
  for (r in seq_len(num_alignments)) {
    h <- hits[hits$rank == r, , drop = FALSE]
    i <- match(out$query, h$query)
    for (f in hit_fields) {
      out[[paste0(r, "_", f)]] <- h[[f]][i]
    }
  }
  out$query <- NULL
  out
}

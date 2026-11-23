#' @title References that do not inform the identification
#'
#' @description
#' The table of reference sequences that [parallel_blast()] sets aside when
#' ranking the hits: sequences of uncultured or unidentified organisms,
#' metagenomes and MAGs, environmental samples, synthetic constructs and
#' vectors, chimeras, and specific references known to be wrong. Their taxid
#' says little about the query (often only "Bacteria" or "environmental
#' samples"), so a hit to one of them should not decide the identification
#' when a hit to a named organism is also there.
#'
#' Each row is a rule:
#' * `field`: `"header"` (`pattern` is a Perl regular expression, matched
#'   against the subject header; use `(?i:...)` groups for case, as the
#'   rules are joined with `"|"`) or `"accession"` (`pattern` is a subject
#'   accession, with or without version).
#' * `category`: the reason, reported in `uninformative_category`.
#' * `note`: free text.
#'
#' Organisms named only at a high rank or in open nomenclature (e.g.
#' "Bacterium", "Picea sp.") are not here: their lineage is right, only
#' shallow. [blast_lca()] sets them aside when a named species is tied with
#' them.
#'
#' @param x `NULL` (default) for the table shipped with BLASTr, a path to a
#'   tab-separated file with the same columns, or a data frame.
#'
#' @returns A tibble with `field`, `pattern`, `category` and `note`, with
#'   the attribute `md5` (of its content, to record which table was used).
#'
#' @examples
#' uninformative_references()
#' @export
uninformative_references <- function(x = NULL) {
  if (is.null(x)) {
    x <- system.file("extdata", "uninformative_references.tsv", package = "BLASTr", mustWork = TRUE)
  }
  if (is.character(x)) {
    if (length(x) != 1L || !file.exists(x)) {
      cli::cli_abort(c(`x` = "Table of uninformative references not found: {.path {x}}."))
    }
    x <- readr::read_tsv(
      x,
      col_types = readr::cols(.default = readr::col_character()),
      quote = "", comment = "", na = character(0L), progress = FALSE
    )
  }
  x <- tibble::as_tibble(as.data.frame(x, stringsAsFactors = FALSE))
  miss <- setdiff(c("field", "pattern", "category"), names(x))
  if (length(miss) > 0L) {
    cli::cli_abort(c(`x` = "The table of uninformative references lacks {.field {miss}}."))
  }
  if (!"note" %in% names(x)) x$note <- NA_character_
  x <- x[, c("field", "pattern", "category", "note")]
  x[] <- lapply(x, as.character)
  bad <- setdiff(unique(x$field), c("header", "accession"))
  if (length(bad) > 0L) {
    cli::cli_abort(c(`x` = "{.field field} must be {.val header} or {.val accession}, not {.val {bad}}."))
  }
  for (p in x$pattern[x$field == "header"]) {
    ok <- tryCatch({
      grepl(p, "", perl = TRUE)
      TRUE
    }, error = function(e) FALSE, warning = function(w) FALSE)
    if (!ok) cli::cli_abort(c(`x` = "Invalid regular expression in the uninformative references: {.val {p}}."))
  }
  txt <- paste(do.call(paste, c(x, sep = "\t")), collapse = "\n")
  f <- tempfile()
  on.exit(unlink(f), add = TRUE)
  writeLines(txt, f, useBytes = TRUE)
  attr(x, "md5") <- unname(tools::md5sum(f))
  x
}

# Category of each hit (NA: informative), by the first rule that matches;
# accession rules first.
uninformative_category <- function(subject, header, table) {
  out <- rep(NA_character_, length(subject))
  if (is.null(table) || nrow(table) == 0L || length(subject) == 0L) {
    return(out)
  }
  acc <- table[table$field == "accession", , drop = FALSE]
  if (nrow(acc) > 0L) {
    s <- as.character(subject)
    # gi|...|acc| style ids: the accession field
    s <- ifelse(grepl("|", s, fixed = TRUE), sub("^.*\\|([^|]+)\\|?$", "\\1", s), s)
    k <- match(s, acc$pattern)
    nov <- is.na(k)
    k[nov] <- match(sub("\\.[0-9]+$", "", s[nov]), sub("\\.[0-9]+$", "", acc$pattern))
    out[!is.na(k)] <- acc$category[k[!is.na(k)]]
  }
  hd <- table[table$field == "header", , drop = FALSE]
  h <- as.character(header)
  for (r in seq_len(nrow(hd))) {
    todo <- is.na(out) & !is.na(h)
    if (!any(todo)) break
    hit <- grepl(hd$pattern[r], h[todo], perl = TRUE)
    out[which(todo)[hit]] <- hd$category[r]
  }
  out
}

# Names in open nomenclature or only at a high rank: "Picea sp.",
# "Carex cf. flava", "Bacterium", "alpha proteobacterium XYZ",
# "Chlamydiota bacterium hydrae", "Candidate division ...", "uncultured ...".
open_nomenclature_pattern <- paste0(
  "(?i:\\b(sp|cf|aff|nr)\\.(\\s|$))|",
  "(?i:^(uncultured|unidentified|unclassified|environmental)\\b)|",
  "(?i:^(alpha|beta|gamma|delta|epsilon)?\\s?proteobacterium\\b)|",
  "(?i:^(bacterium|archaeon|eukaryote|fungal|fungus)\\b)|",
  "(?i:\\b(bacterium|proteobacterium|archaeon|endophyte)\\b)|",
  "(?i:candidate division)"
)

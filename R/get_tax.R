# Ranks returned by default by get_tax().
tax_ranks_default <- c(
  "domain", "kingdom", "phylum", "subphylum", "class", "subclass", "order",
  "suborder", "family", "subfamily", "genus", "species"
)

#' @title Taxonomy of NCBI taxids from a local database
#'
#' @description
#' Looks up the NCBI taxonomy of taxids in a local SQLite database: no
#' internet, no `efetch`, thousands of taxids in a few seconds. Returns the
#' minimal table, with one row per distinct taxid, to be joined to the BLAST
#' hits by taxid.
#'
#' `x` can be a vector of taxids or any result of [parallel_blast()] (long
#' or wide) or [parallel_blast_hsps()]: the taxids are taken from the
#' `staxid`, `k_staxid` and `Taxids at best` columns. Unknown
#' taxids (`0`, empty, `N/A`) are dropped.
#'
#' When the database does not exist, it is built from the NCBI taxdump
#' (downloaded once, then reused). If `x` is a BLASTr result, the taxdump is
#' the first one the NCBI published after the date of its BLAST database, so
#' the taxonomy matches the database; the database is kept in the BLASTr
#' cache as `taxonomy_<database date>.sqlite`. With only taxids and no
#' database set, the newest database in the cache is used, or one is built
#' from the current taxdump.
#'
#' The database has the tables `nodes` (`id`, `rank`, `parent`) and
#' `names` (`id`, `name`, `scientific`, `type`), as built by
#' [make_taxonomy_db()] (or by `taxonomizr::prepareDatabase()`, whose file
#' can be used as is). Build it from the NCBI taxdump of the same date as the
#' BLAST database: taxids newer than the taxdump are not found.
#'
#' @param x Taxids (character or numeric) or a BLASTr result.
#' @param taxonomy_db Path to the SQLite taxonomy database. `NULL` uses the
#'   option `BLASTr.taxonomy_db` or the environment variable
#'   `BLASTR_TAXONOMY_DB`; with none set, a database in the BLASTr cache
#'   (see below). If the file does not exist, it is built.
#' @param taxdump NCBI taxdump used to build a missing database (see
#'   [make_taxonomy_db()]). `NULL`: the taxdump matching the BLAST database
#'   of `x` (the first monthly NCBI archive after its date), or the current
#'   one.
#' @param ranks Ranks returned as columns. `"domain"` is also filled from
#'   `"superkingdom"`, its name in taxdumps before 2025.
#'
#' @returns A tibble with one row per taxid: `taxid`, `found` (`FALSE` when
#'   the taxid is not in the database), `sci_name`, `rank` (rank of the taxid
#'   itself), one column per rank in `ranks` (`NA` when the lineage has no
#'   node of that rank), `lineage` (all the nodes from the root to the
#'   taxid, without "root" and "cellular organisms", separated by `"; "`),
#'   `lineage_taxids` (their taxids, in the same order, for an exact lowest
#'   common ancestor) and `merged_into` (the current taxid when the NCBI
#'   merged the taxid into another one; then the taxonomy is that of the
#'   current taxid; needs a database with the table `merged`, as built by
#'   [make_taxonomy_db()]).
#'
#' @examples
#' # a simplified NCBI taxdump: Hominidae with the genera Homo and Pan
#' taxdump <- file.path(tempdir(), "taxdump")
#' dir.create(taxdump, showWarnings = FALSE)
#' dmp <- function(...) paste0(paste(..., sep = "\t|\t"), "\t|")
#' writeLines(c(
#'   dmp(1, 1, "no rank"), dmp(2759, 1, "superkingdom"), dmp(9604, 2759, "family"),
#'   dmp(9605, 9604, "genus"), dmp(9606, 9605, "species"),
#'   dmp(9596, 9604, "genus"), dmp(9598, 9596, "species")
#' ), file.path(taxdump, "nodes.dmp"))
#' writeLines(c(
#'   dmp(1, "root", "", "scientific name"), dmp(2759, "Eukaryota", "", "scientific name"),
#'   dmp(9604, "Hominidae", "", "scientific name"), dmp(9605, "Homo", "", "scientific name"),
#'   dmp(9606, "Homo sapiens", "", "scientific name"), dmp(9596, "Pan", "", "scientific name"),
#'   dmp(9598, "Pan troglodytes", "", "scientific name")
#' ), file.path(taxdump, "names.dmp"))
#' tax_db <- make_taxonomy_db(tempfile(fileext = ".sqlite"), taxdump = taxdump)
#'
#' # 1234 is not in this taxonomy: found = FALSE, with a warning
#' get_tax(c("9606", "9598", "1234"), taxonomy_db = tax_db)
#'
#' # hits of two sequences: the first one tied between Homo and Pan
#' wide <- data.frame(
#'   Sequence = c("ACGT", "TTGA"),
#'   `1_subject` = c("A1", "B1"), `1_identity` = c(100, 99), `1_qcovhsp` = 100,
#'   `1_staxid` = "9606",
#'   `2_subject` = c("A2", "B2"), `2_identity` = c(99.8, 90), `2_qcovhsp` = 100,
#'   `2_staxid` = "9598",
#'   check.names = FALSE
#' )
#' hits <- as_blast_long(wide)
#' tax <- get_tax(hits, taxonomy_db = tax_db)
#' dplyr::left_join(hits, tax, by = c(staxid = "taxid"))
#'
#' @family taxonomy functions
#' @export
get_tax <- function(x, taxonomy_db = NULL, ranks = tax_ranks_default, taxdump = NULL) {
  taxids <- extract_taxids(x)
  db <- resolve_taxonomy_db(taxonomy_db, taxdump = taxdump, db_dates = attr(x, "BLASTr_metadata")$db_dates)
  ranks <- as.character(ranks)

  con <- DBI::dbConnect(RSQLite::SQLite(), db, flags = RSQLite::SQLITE_RO)
  on.exit(DBI::dbDisconnect(con), add = TRUE)

  # taxids merged into others (table `merged`, from merged.dmp)
  ids <- suppressWarnings(as.integer(taxids))
  merged_into <- rep(NA_integer_, length(ids))
  if (DBI::dbExistsTable(con, "merged")) {
    present <- sql_in(con, "SELECT id FROM nodes WHERE id IN (%s)", ids)$id
    old <- unique(ids[!is.na(ids) & !ids %in% present])
    if (length(old) > 0L) {
      mg <- sql_in(con, "SELECT old_id, new_id FROM merged WHERE old_id IN (%s)", old)
      merged_into <- as.integer(mg$new_id[match(ids, mg$old_id)])
    }
  }
  ids <- ifelse(is.na(merged_into), ids, merged_into)

  # climb from the taxids to the root, one level of the tree per query
  nodes <- data.frame(id = integer(0L), rank = character(0L), parent = integer(0L))
  todo <- unique(ids)
  todo <- todo[!is.na(todo)]
  while (length(todo) > 0L) {
    got <- sql_in(con, "SELECT id, rank, parent FROM nodes WHERE id IN (%s)", todo)
    got$id <- as.integer(got$id)
    got$parent <- as.integer(got$parent)
    got <- got[!duplicated(got$id), , drop = FALSE]
    nodes <- rbind(nodes, got[, c("id", "rank", "parent")])
    todo <- setdiff(got$parent[got$parent != got$id], nodes$id)
  }
  nm <- sql_in(con, "SELECT id, name FROM names WHERE scientific = 1 AND id IN (%s)", nodes$id)
  nm <- nm[!duplicated(nm$id), , drop = FALSE]

  # lineage matrix: one row per taxid, columns from the taxid up to the root
  found <- !is.na(ids) & ids %in% nodes$id
  levels <- list(ifelse(found, ids, NA_integer_))
  repeat {
    cur <- levels[[length(levels)]]
    p <- nodes$parent[match(cur, nodes$id)]
    p[!is.na(p) & p == cur] <- NA_integer_
    if (all(is.na(p))) break
    levels[[length(levels) + 1L]] <- p
  }
  lin_id <- do.call(cbind, levels)
  lin_name <- matrix(nm$name[match(lin_id, nm$id)], nrow = nrow(lin_id), ncol = ncol(lin_id))
  lin_rank <- matrix(nodes$rank[match(lin_id, nodes$id)], nrow = nrow(lin_id), ncol = ncol(lin_id))
  # nearest node of each rank (the first one going up)
  pick_rank <- function(k) {
    hit <- lin_rank == k & !is.na(lin_rank)
    j <- max.col(hit, ties.method = "first")
    v <- lin_name[cbind(seq_len(nrow(hit)), j)]
    v[rowSums(hit) == 0L] <- NA_character_
    v
  }
  rank_cols <- lapply(ranks, function(k) {
    v <- pick_rank(k)
    if (k == "domain") v <- dplyr::coalesce(v, pick_rank("superkingdom"))
    v
  })
  skip <- c("root", "cellular organisms")
  lin_paths <- lapply(seq_len(nrow(lin_name)), function(i) {
    keep <- rev(!is.na(lin_id[i, ]) & !lin_name[i, ] %in% skip)
    list(name = rev(lin_name[i, ])[keep], id = rev(lin_id[i, ])[keep])
  })
  collapse <- function(v) if (length(v)) paste(v, collapse = "; ") else NA_character_
  lineage <- vapply(lin_paths, function(x) collapse(x$name), character(1L))
  lineage_taxids <- vapply(lin_paths, function(x) collapse(x$id), character(1L))
  out <- tibble::tibble(
    taxid = as.character(taxids),
    found = found,
    sci_name = lin_name[, 1L],
    rank = lin_rank[, 1L],
    !!!stats::setNames(rank_cols, ranks),
    lineage = lineage,
    lineage_taxids = lineage_taxids,
    merged_into = as.character(merged_into)
  )
  miss <- out$taxid[!out$found]
  if (length(miss) > 0L) {
    cli::cli_warn(c(
      `!` = "{length(miss)} taxid{?s} not found in the taxonomy database: {utils::head(miss, 10)}{if (length(miss) > 10) ', ...' else ''}",
      `i` = "Use a taxonomy database as recent as the BLAST database (see {.fun make_taxonomy_db})."
    ))
  }
  out
}

# Run "... IN (%s)" over integer ids, in chunks.
sql_in <- function(con, sql, ids, chunk = 5000L) {
  ids <- unique(as.integer(ids))
  ids <- ids[!is.na(ids)]
  if (length(ids) == 0L) {
    return(DBI::dbGetQuery(con, sprintf(sql, "NULL")))
  }
  parts <- split(ids, ceiling(seq_along(ids) / chunk))
  do.call(rbind, lapply(parts, function(p) {
    DBI::dbGetQuery(con, sprintf(sql, paste(p, collapse = ",")))
  }))
}

# Distinct known taxids of a vector or of a BLASTr result.
extract_taxids <- function(x) {
  if (is.data.frame(x)) {
    # "Taxids at best identity": its name before 0.2.0, in older wide results
    cols <- grep("^([0-9]+_)?staxid$|^Taxids at best( identity)?$", names(x), value = TRUE)
    if (length(cols) == 0L) {
      cli::cli_abort(
        c(`x` = "No taxid column ({.field staxid}, {.field k_staxid} or {.field Taxids at best}) in {.arg x}."),
        class = "blastr_no_taxids_error"
      )
    }
    x <- unlist(lapply(cols, function(k) as.character(x[[k]])), use.names = FALSE)
  }
  x <- trimws(unlist(strsplit(as.character(x), "[;,]"), use.names = FALSE))
  x <- x[!is.na(x) & !x %in% c("", "0", "N/A", "NA")]
  unique(x)
}

# The taxonomy database to use, building it when it does not exist:
# `taxonomy_db` (argument, option BLASTr.taxonomy_db or variable
# BLASTR_TAXONOMY_DB) or else one in the BLASTr cache, named after the date
# of the BLAST database(s) of the search (`db_dates`), or the newest one
# there.
resolve_taxonomy_db <- function(taxonomy_db = NULL, taxdump = NULL, db_dates = NULL) {
  db <- if (is.null(taxonomy_db)) getOption("BLASTr.taxonomy_db") else taxonomy_db
  if (is.null(db) || !nzchar(db)) db <- Sys.getenv("BLASTR_TAXONOMY_DB")
  db_date <- latest_date(db_dates)
  if (!nzchar(db)) {
    cache <- get_blastr_cache("taxonomy")
    if (!is.na(db_date)) {
      db <- fs::path(cache, paste0("taxonomy_", db_date, ".sqlite"))
    } else if (!is.null(taxdump)) {
      db <- fs::path(cache, paste0("taxonomy_", taxdump_label(taxdump), ".sqlite"))
    } else {
      old <- if (fs::dir_exists(cache)) fs::dir_ls(cache, glob = "*/taxonomy_*.sqlite") else character(0L)
      db <- if (length(old) > 0L) sort(old, decreasing = TRUE)[1L] else fs::path(cache, paste0("taxonomy_", Sys.Date(), ".sqlite"))
    }
  }
  if (!fs::file_exists(db)) {
    if (is.null(taxdump)) taxdump <- taxdump_url(db_date)
    cli::cli_inform(c(`i` = "Building the taxonomy database {.path {db}} from {.path {taxdump}}."))
    make_taxonomy_db(db, taxdump = taxdump)
  }
  as.character(db)
}

# Newest date (YYYY-MM-DD) among BLAST database dates.
latest_date <- function(db_dates) {
  d <- suppressWarnings(as.Date(substr(as.character(db_dates), 1L, 10L)))
  d <- d[!is.na(d)]
  if (length(d) == 0L) NA_character_ else format(max(d))
}

taxdump_label <- function(taxdump) {
  if (grepl("^(https?|ftp)://", taxdump)) {
    return(gsub("[^A-Za-z0-9_-]", "_", fs::path_ext_remove(fs::path_ext_remove(basename(taxdump)))))
  }
  format(as.Date(fs::file_info(taxdump)$modification_time))
}

# The NCBI taxdump for a BLAST database of a given date: the first monthly
# archive published after it (so it has every taxid of the database), or the
# current taxdump when there is none yet.
taxdump_url <- function(db_date = NA_character_) {
  current <- "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
  if (is.na(db_date)) return(current)
  d <- as.Date(db_date)
  first <- as.Date(format(d, "%Y-%m-01"))
  month <- if (d == first) first else seq(first, by = "month", length.out = 2L)[2L]
  if (month > Sys.Date()) return(current)
  sprintf("https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_archive/taxdmp_%s.zip", format(month))
}

#' @title Build the local taxonomy database
#'
#' @description
#' Builds the SQLite database used by [get_tax()] from the NCBI taxdump
#' (`nodes.dmp`, `names.dmp` and `merged.dmp`), downloading it if not given.
#' The tables `nodes` and `names` follow the `taxonomizr` layout, so a
#' database built by `taxonomizr::prepareDatabase()` also works (without
#' resolving merged taxids); `merged` (`old_id`, `new_id`) maps the taxids
#' the NCBI merged into others.
#'
#' Use a taxdump of the same date as the BLAST database (or newer): taxids
#' created after the taxdump are not found.
#'
#' @param path Path of the SQLite file to create.
#' @param taxdump NCBI taxdump: a `taxdump.tar.gz` or `taxdmp_*.zip` file,
#'   its URL, or a directory with `nodes.dmp`, `names.dmp` and (optional)
#'   `merged.dmp`. `NULL` downloads the current taxdump from the NCBI.
#' @param overwrite Replace `path` if it exists.
#'
#' @returns The path of the database, invisibly.
#'
#' @examples
#' # a simplified NCBI taxdump: Hominidae with the genera Homo and Pan
#' taxdump <- file.path(tempdir(), "taxdump")
#' dir.create(taxdump, showWarnings = FALSE)
#' dmp <- function(...) paste0(paste(..., sep = "\t|\t"), "\t|")
#' writeLines(c(
#'   dmp(1, 1, "no rank"), dmp(2759, 1, "superkingdom"), dmp(9604, 2759, "family"),
#'   dmp(9605, 9604, "genus"), dmp(9606, 9605, "species"),
#'   dmp(9596, 9604, "genus"), dmp(9598, 9596, "species")
#' ), file.path(taxdump, "nodes.dmp"))
#' writeLines(c(
#'   dmp(1, "root", "", "scientific name"), dmp(2759, "Eukaryota", "", "scientific name"),
#'   dmp(9604, "Hominidae", "", "scientific name"), dmp(9605, "Homo", "", "scientific name"),
#'   dmp(9606, "Homo sapiens", "", "scientific name"), dmp(9596, "Pan", "", "scientific name"),
#'   dmp(9598, "Pan troglodytes", "", "scientific name")
#' ), file.path(taxdump, "names.dmp"))
#' tax_db <- make_taxonomy_db(tempfile(fileext = ".sqlite"), taxdump = taxdump)
#' get_tax("9606", taxonomy_db = tax_db)
#'
#' \dontrun{
#' # from the current NCBI taxdump (downloaded, about 70 MB)
#' make_taxonomy_db("~/blast_databases/taxonomy.sqlite")
#' options(BLASTr.taxonomy_db = "~/blast_databases/taxonomy.sqlite")
#' }
#'
#' @family taxonomy functions
#' @export
make_taxonomy_db <- function(path, taxdump = NULL, overwrite = FALSE) {
  if (fs::file_exists(path) && !isTRUE(overwrite)) {
    cli::cli_abort(c(`x` = "{.path {path}} exists; use {.code overwrite = TRUE}."))
  }
  tmp <- fs::file_temp(pattern = "taxdump_")
  fs::dir_create(tmp)
  on.exit(fs::dir_delete(tmp), add = TRUE)
  if (is.null(taxdump)) taxdump <- taxdump_url()
  if (grepl("^(https?|ftp)://", taxdump)) {
    url <- taxdump
    taxdump <- fs::path(tmp, basename(url))
    cli::cli_inform(c(`i` = "Downloading {.url {url}}."))
    withr::local_options(timeout = max(3600, getOption("timeout")))
    utils::download.file(url, taxdump, mode = "wb", quiet = TRUE)
  }
  files <- c("nodes.dmp", "names.dmp", "merged.dmp")
  dir <- if (fs::is_dir(taxdump)) {
    taxdump
  } else if (grepl("\\.zip$", taxdump)) {
    utils::unzip(taxdump, files = files, exdir = tmp)
    tmp
  } else {
    utils::untar(taxdump, files = files, exdir = tmp)
    tmp
  }
  read_dmp <- function(f, cols, keep) {
    x <- readr::read_delim(
      fs::path(dir, f),
      delim = "\t|\t", col_names = FALSE, col_types = readr::cols(.default = readr::col_character()),
      quote = "", trim_ws = FALSE, progress = FALSE, na = character(0L)
    )
    x <- x[, keep]
    names(x) <- cols
    # the last field of each line ends with "\t|"
    x[[ncol(x)]] <- sub("\t\\|$", "", x[[ncol(x)]])
    x
  }
  nodes <- read_dmp("nodes.dmp", c("id", "parent", "rank"), 1:3)
  names_ <- read_dmp("names.dmp", c("id", "name", "type"), c(1L, 2L, 4L))
  nodes <- data.frame(id = as.integer(nodes$id), rank = nodes$rank, parent = as.integer(nodes$parent))
  names_ <- data.frame(
    id = as.integer(names_$id), name = names_$name,
    scientific = as.integer(names_$type == "scientific name"), type = names_$type
  )
  if (fs::file_exists(path)) fs::file_delete(path)
  fs::dir_create(fs::path_dir(path))
  con <- DBI::dbConnect(RSQLite::SQLite(), path)
  on.exit(DBI::dbDisconnect(con), add = TRUE, after = FALSE)
  DBI::dbWriteTable(con, "nodes", nodes)
  DBI::dbWriteTable(con, "names", names_)
  if (fs::file_exists(fs::path(dir, "merged.dmp"))) {
    mg <- read_dmp("merged.dmp", c("old_id", "new_id"), 1:2)
    DBI::dbWriteTable(con, "merged", data.frame(old_id = as.integer(mg$old_id), new_id = as.integer(mg$new_id)))
    DBI::dbExecute(con, "CREATE INDEX index_merged_old ON merged (old_id)")
  }
  DBI::dbExecute(con, "CREATE INDEX index_nodes_id ON nodes (id)")
  DBI::dbExecute(con, "CREATE INDEX index_nodes_parent ON nodes (parent)")
  DBI::dbExecute(con, "CREATE INDEX index_names_id ON names (id)")
  DBI::dbExecute(con, "CREATE INDEX index_names_name ON names (name, scientific)")
  cli::cli_inform(c(`v` = "Taxonomy database: {.path {path}} ({nrow(nodes)} taxa)."))
  invisible(as.character(path))
}

#' @title BLAST databases of the NCBI
#'
#' @description
#' `list_ncbi_blast_dbs()` lists the preformatted BLAST databases the NCBI
#' publishes (<https://ftp.ncbi.nlm.nih.gov/blast/db/>), with their size and
#' date. `download_blast_db()` downloads one of them (e.g. `"core_nt"`,
#' `"16S_ribosomal_RNA"`, `"ITS_eukaryote_sequences"`), checks the md5 of each
#' file, unpacks it and returns the path to give to [parallel_blast()].
#'
#' The database goes to `<dir>/<db>_<date>/` (the NCBI date of the
#' database), so a newer release does not overwrite the one older results
#' were searched against. An interrupted download is resumed: files already
#' unpacked (with a matching md5) are not downloaded again. With
#' `taxonomy = TRUE`, the taxonomy database of [get_tax()] matching the date
#' of the BLAST database is built too, so the identifications can be made
#' offline.
#'
#' `core_nt` has about 300 GB unpacked (and is downloaded in volumes of a few
#' GB, unpacked one at a time); see `list_ncbi_blast_dbs()` for the sizes.
#'
#' @param db Name of the database, as in `list_ncbi_blast_dbs()$db`.
#' @param dir Directory where the database directory is created.
#' @param taxdb Also download `taxdb` (scientific names of the taxids, used by
#'   BLAST in `sscinames`) into the same directory.
#' @param taxonomy Build the local taxonomy database of [get_tax()] for the
#'   date of this database (kept in the BLASTr cache).
#' @param keep_archives Keep the `.tar.gz` files after unpacking.
#' @param check_space Stop when the disk does not have room for the
#'   unpacked database plus one archive.
#'
#' @returns `list_ncbi_blast_dbs()`: a tibble with `db`, `type`
#'   (`"nucl"`/`"prot"`), `description`, `n_seqs`, `n_letters`, `n_volumes`,
#'   `size_gb` (unpacked), `download_gb` and `last_updated`.
#'   `download_blast_db()`: the path of the database (directory and name, as
#'   `db_path` of [parallel_blast()]), invisibly.
#'
#' @examples
#' \dontrun{
#' list_ncbi_blast_dbs()
#' db <- download_blast_db("16S_ribosomal_RNA", dir = "/data/databases")
#' parallel_blast(seqs, db_path = db, total_cores = 8)
#' }
#' @name ncbi_blast_dbs
NULL

ncbi_blast_url <- "https://ftp.ncbi.nlm.nih.gov/blast/db/"

#' @rdname ncbi_blast_dbs
#' @export
list_ncbi_blast_dbs <- function() {
  js <- read_json_url(paste0(ncbi_blast_url, "blastdb-metadata-1-1.json"))
  num <- function(x, k) vapply(x, function(e) as.numeric(null_to(e[[k]], NA_real_)), numeric(1L))
  chr <- function(x, k) vapply(x, function(e) as.character(null_to(e[[k]], NA_character_)), character(1L))
  out <- tibble::tibble(
    db = chr(js, "dbname"),
    type = ifelse(chr(js, "dbtype") == "Protein", "prot", "nucl"),
    description = chr(js, "description"),
    n_seqs = num(js, "number-of-sequences"),
    n_letters = num(js, "number-of-letters"),
    n_volumes = as.integer(num(js, "number-of-volumes")),
    size_gb = round(num(js, "bytes-total") / 1e9, 2),
    download_gb = round(num(js, "bytes-total-compressed") / 1e9, 2),
    last_updated = as.Date(substr(chr(js, "last-updated"), 1L, 10L))
  )
  out[order(out$type, out$db), ]
}

#' @rdname ncbi_blast_dbs
#' @export
download_blast_db <- function(
  db = "core_nt",
  dir,
  taxdb = TRUE,
  taxonomy = TRUE,
  keep_archives = FALSE,
  check_space = TRUE
) {
  rlang::check_required(dir)
  dbs <- list_ncbi_blast_dbs()
  if (!db %in% dbs$db) {
    cli::cli_abort(c(`x` = "{.val {db}} is not a BLAST database of the NCBI.", `i` = "See {.fun list_ncbi_blast_dbs}."))
  }
  meta_url <- sprintf("%s%s-%s-metadata.json", ncbi_blast_url, db, dbs$type[dbs$db == db])
  meta <- read_json_url(meta_url)
  date <- substr(meta[["last-updated"]], 1L, 10L)
  out_dir <- fs::path(dir, paste0(db, "_", date))
  fs::dir_create(out_dir)
  db_path <- fs::path(out_dir, db)
  files <- sub("^ftp://", "https://", unlist(meta$files))
  if (isTRUE(taxdb)) files <- c(files, paste0(ncbi_blast_url, "taxdb.tar.gz"))

  done_file <- fs::path(out_dir, "BLASTr_download.tsv")
  done <- if (fs::file_exists(done_file)) utils::read.delim(done_file, colClasses = "character") else data.frame(file = character(0L), md5 = character(0L))
  todo <- files[!basename(files) %in% done$file]

  if (isTRUE(check_space) && length(todo) > 0L) {
    need <- as.numeric(meta[["bytes-total"]]) * length(todo) / length(files) +
      as.numeric(meta[["bytes-total-compressed"]]) / max(1L, length(unlist(meta$files)))
    free <- disk_free(out_dir)
    if (!is.na(free) && free < need) {
      cli::cli_abort(c(
        `x` = "Not enough disk in {.path {out_dir}}: {round(free / 1e9, 1)} GB free, {round(need / 1e9, 1)} GB needed.",
        `i` = "Free space or choose another {.arg dir} (or {.code check_space = FALSE})."
      ))
    }
  }
  cli::cli_inform(c(`i` = "{.val {db}} of {date}: {length(files)} file{?s}, {length(todo)} to download, into {.path {out_dir}}."))
  for (u in todo) {
    md5 <- sub("\\s.*$", "", trimws(read_text_url(paste0(u, ".md5"))[1L]))
    tar <- fs::path(out_dir, basename(u))
    if (!fs::file_exists(tar) || unname(tools::md5sum(tar)) != md5) {
      download_checked(u, tar, md5)
    }
    utils::untar(tar, exdir = out_dir)
    if (!isTRUE(keep_archives)) fs::file_delete(tar)
    done <- rbind(done, data.frame(file = basename(u), md5 = md5))
    utils::write.table(done, done_file, sep = "\t", quote = FALSE, row.names = FALSE)
    cli::cli_inform(c(`v` = "{basename(u)} ({nrow(done)}/{length(files)})"))
  }
  writeLines(read_text_url(meta_url), fs::path(out_dir, paste0(db, "-metadata.json")))
  if (isTRUE(taxonomy)) resolve_taxonomy_db(db_dates = date)
  cli::cli_inform(c(`v` = "BLAST database ready: {.path {db_path}}"))
  invisible(as.character(db_path))
}

#' @title UNITE and SILVA as BLAST databases
#'
#' @description
#' Download the reference databases UNITE (ITS of fungi and other
#' eukaryotes) and SILVA (SSU rRNA, 16S/18S) and format them as BLAST
#' databases that work like the ones of the NCBI in BLASTr: each sequence has
#' an NCBI taxid (so [get_tax()], [blast_lca()] and the summaries work), its
#' accession as the subject id, and the classification of UNITE or SILVA in
#' the subject header.
#'
#' * `make_unite_db()`: the UNITE "general FASTA release" (representative and
#'   reference sequences of the species hypotheses, SH), for all eukaryotes or
#'   only fungi. UNITE has no NCBI taxids: each sequence gets the taxid of the
#'   deepest name of its UNITE classification found in the NCBI taxonomy
#'   (the species, else the genus, the family...). Names in open nomenclature
#'   (`"Fungi_sp"`, `"..._gen_Incertae_sedis"`) are skipped, and homonyms are
#'   resolved by the higher ranks of the classification.
#' * `make_silva_db()`: SILVA SSU Ref NR99, with RNA converted to DNA. The
#'   taxid is the NCBI taxid SILVA gives each sequence (from the organism of
#'   the INSDC record), or, when missing, the one found from the names as for
#'   UNITE. The SILVA classification, which may differ from the organism of
#'   the record, is in the header.
#'
#' The table `<name>_taxids.tsv`, next to the database, records for each
#' sequence the taxid, how it was found and the name it was found by.
#'
#' Both databases are licensed under CC BY 4.0: cite the release (the DOI is
#' printed and kept in `<name>_source.txt`).
#'
#' @param dir Directory where the database directory is created.
#' @param set UNITE: `"eukaryotes"` (all eukaryotes) or `"fungi"`.
#' @param doi UNITE: the DOI of a general FASTA release (e.g.
#'   `"10.15156/BIO/3671843"`); `NULL` takes the newest one of `set` listed at
#'   <https://unite.ut.ee/repository.php>.
#' @param singletons UNITE: use the file that also has the singletons (the
#'   "_dev" file of the release).
#' @param release SILVA: release number (e.g. `"138.2"`), or `"current"`.
#' @param taxonomy_db Taxonomy database to find the taxids from names (see
#'   [get_tax()]); `NULL` uses (or builds) the current one.
#'
#' @returns The path of the database (as `db_path` of [parallel_blast()]),
#'   invisibly.
#'
#' @examples
#' \dontrun{
#' its <- make_unite_db("/data/databases")
#' ssu <- make_silva_db("/data/databases")
#' parallel_blast(seqs, db_path = its, total_cores = 8)
#' }
#' @name reference_dbs
NULL

#' @rdname reference_dbs
#' @export
make_unite_db <- function(dir, set = c("eukaryotes", "fungi"), doi = NULL, singletons = FALSE, taxonomy_db = NULL) {
  rlang::check_required(dir)
  set <- rlang::arg_match(set)
  if (is.null(doi)) doi <- unite_latest_doi(set)
  api <- read_json_url(paste0(
    "https://api.plutof.ut.ee/v1/public/dois/?format=vnd.api%2Bjson&identifier=",
    utils::URLencode(doi, reserved = TRUE)
  ))
  att <- api$data[[1L]]$attributes
  media <- att$media[[1L]]
  version <- regmatches(media$name, regexpr("[0-9]{2}\\.[0-9]{2}\\.[0-9]{4}", media$name))
  date <- if (length(version)) format(as.Date(version, "%d.%m.%Y")) else format(Sys.Date())
  name <- paste0("unite_", set, "_", date)
  out_dir <- fs::path(dir, name)
  fs::dir_create(out_dir)

  tmp <- fs::file_temp(pattern = "unite_")
  fs::dir_create(tmp)
  on.exit(fs::dir_delete(tmp), add = TRUE)
  tgz <- fs::path(tmp, media$name)
  download_checked(media$url, tgz)
  utils::untar(tgz, exdir = tmp)
  fa <- fs::dir_ls(tmp, recurse = TRUE, glob = "*.fasta")
  fa <- if (isTRUE(singletons)) fa[grepl("_dev\\.fasta$", fa)] else fa[!grepl("_dev\\.fasta$", fa)]
  if (length(fa) != 1L) cli::cli_abort(c(`x` = "Unexpected files in the UNITE release: {.file {basename(fa)}}."))

  x <- read_fasta_records(fa)
  f <- strsplit(x$header, "|", fixed = TRUE)
  acc <- vapply(f, `[`, character(1L), 2L)
  sh <- vapply(f, `[`, character(1L), 3L)
  kind <- vapply(f, `[`, character(1L), 4L)
  cls <- vapply(f, function(v) v[length(v)], character(1L))
  levels <- unite_levels(cls)
  tax <- names_to_taxids(levels, taxonomy_db)
  ids <- make_unique_ids(acc)
  title <- paste(sh, kind, cls)
  db <- write_ref_db(out_dir, name, ids, title, x$seq, tax, extra = data.frame(sh = sh, unite = cls))
  src <- c(
    paste("UNITE general FASTA release,", set), paste("DOI:", doi),
    paste("File:", basename(fa)), paste("Citation:", att$citation),
    paste("License:", att$licence$rights, att$licence$rights_uri)
  )
  writeLines(src, fs::path(out_dir, paste0(name, "_source.txt")))
  cli::cli_inform(c(`v` = "UNITE database ready: {.path {db}}", `i` = "Cite: {att$citation}"))
  invisible(db)
}

#' @rdname reference_dbs
#' @export
make_silva_db <- function(dir, release = "current", taxonomy_db = NULL) {
  rlang::check_required(dir)
  base <- "https://ftp.arb-silva.de/"
  exports <- if (identical(release, "current")) paste0(base, "current/Exports/") else paste0(base, "release_", release, "/Exports/")
  listing <- paste(read_text_url(exports), collapse = "\n")
  hit <- regmatches(listing, regexpr("SILVA_[0-9.]+_SSURef_NR99_tax_silva_trunc\\.fasta\\.gz", listing))
  if (length(hit) == 0L) cli::cli_abort(c(`x` = "No SSU Ref NR99 FASTA in {.url {exports}}."))
  rel <- sub("^SILVA_([0-9.]+)_.*$", "\\1", hit)
  name <- paste0("silva_ssu_nr99_", rel)
  out_dir <- fs::path(dir, name)
  fs::dir_create(out_dir)

  tmp <- fs::file_temp(pattern = "silva_")
  fs::dir_create(tmp)
  on.exit(fs::dir_delete(tmp), add = TRUE)
  get <- function(url) {
    md5 <- sub("\\s.*$", "", trimws(read_text_url(paste0(url, ".md5"))[1L]))
    f <- fs::path(tmp, basename(url))
    download_checked(url, f, md5)
    f
  }
  fa <- get(paste0(exports, hit))
  map <- get(sprintf("%staxonomy/ncbi/taxmap_embl-ebi_ena_ssu_ref_nr99_%s.txt.gz", exports, rel))

  x <- read_fasta_records(fa)
  id <- sub("\\s.*$", "", x$header)
  cls <- sub("^\\S+\\s*", "", x$header)
  tm <- suppressWarnings(readr::read_tsv(map, col_types = readr::cols(.default = readr::col_character()), progress = FALSE, quote = ""))
  key <- paste(tm$primaryAccession, tm$start, tm$stop, sep = ".")
  taxid <- tm$ncbi_taxonid[match(id, key)]
  taxid[!grepl("^[0-9]+$", taxid)] <- NA_character_ # a few rows of the file are malformed
  levels <- silva_levels(cls)
  tax <- data.frame(taxid = taxid, how = ifelse(is.na(taxid), NA_character_, "SILVA (NCBI taxid of the record)"), name = NA_character_)
  miss <- is.na(tax$taxid) | tax$taxid %in% c("", "0")
  if (any(miss)) tax[miss, ] <- names_to_taxids(levels[miss], taxonomy_db)
  seq <- chartr("Uu", "Tt", x$seq)
  db <- write_ref_db(out_dir, name, gsub(".", "_", id, fixed = TRUE), cls, seq, tax, extra = data.frame(silva_id = id, silva = cls))
  doi <- tryCatch(trimws(read_text_url(paste0(exports, hit, ".doi"))[1L]), error = function(e) NA_character_)
  src <- c(
    paste("SILVA SSU Ref NR99 release", rel), paste("File:", hit), paste("DOI:", doi),
    "License: CC BY 4.0 https://creativecommons.org/licenses/by/4.0/",
    "Cite SILVA: https://www.arb-silva.de/contact/citation"
  )
  writeLines(src, fs::path(out_dir, paste0(name, "_source.txt")))
  cli::cli_inform(c(`v` = "SILVA database ready: {.path {db}}", `i` = "SILVA is CC BY 4.0; cite the release ({doi})."))
  invisible(db)
}

# --- helpers ---------------------------------------------------------------

null_to <- function(x, y) if (is.null(x)) y else x

# DOI of the newest UNITE general FASTA release of a set, from the list at
# unite.ut.ee (newest first).
unite_latest_doi <- function(set) {
  page <- paste(read_text_url("https://unite.ut.ee/repository.php"), collapse = " ")
  page <- gsub("\\s+", " ", gsub("<[^>]*>", " ", page))
  lab <- if (set == "fungi") "Fungi" else "eukaryotes"
  re <- sprintf("UNITE general FASTA release for %s\\. Version [0-9.]+\\. UNITE Community\\. https://doi\\.org/(10\\.15156/BIO/[0-9]+)", lab)
  m <- regmatches(page, regexpr(re, page))
  if (length(m) == 0L) {
    cli::cli_abort(c(`x` = "Could not find the UNITE general FASTA release at {.url https://unite.ut.ee/repository.php}.", `i` = "Give the release with {.arg doi}."))
  }
  sub("^.*doi\\.org/", "", m)
}

read_text_url <- function(url) {
  tmp <- tempfile()
  on.exit(unlink(tmp), add = TRUE)
  download_checked(url, tmp)
  readLines(tmp, warn = FALSE)
}

read_json_url <- function(url) {
  jsonlite::fromJSON(paste(read_text_url(url), collapse = "\n"), simplifyVector = FALSE)
}

# Download with retries; when `md5` is given, the file must match it.
download_checked <- function(url, dest, md5 = NULL, tries = 3L) {
  withr::local_options(timeout = max(24 * 3600, getOption("timeout")))
  for (i in seq_len(tries)) {
    ok <- tryCatch(
      {
        utils::download.file(url, dest, mode = "wb", quiet = TRUE)
        TRUE
      },
      error = function(e) FALSE, warning = function(w) FALSE
    )
    if (ok && (is.null(md5) || unname(tools::md5sum(dest)) == md5)) return(invisible(dest))
    Sys.sleep(5 * i)
  }
  cli::cli_abort(
    c(`x` = "Download failed{if (!is.null(md5)) ' or md5 did not match'}: {.url {url}}"),
    class = "blastr_error_download"
  )
}

disk_free <- function(path) {
  tryCatch(as.numeric(ps::ps_disk_usage(as.character(path))$available), error = function(e) NA_real_)
}

# Records of a FASTA file (any line length): header (without ">") and seq.
read_fasta_records <- function(path) {
  x <- readr::read_lines(path, progress = FALSE)
  h <- startsWith(x, ">")
  grp <- cumsum(h)
  seq <- vapply(split(x[!h], factor(grp[!h], levels = seq_len(sum(h)))), paste, character(1L), collapse = "")
  data.frame(header = substring(x[h], 2L), seq = toupper(gsub("[^A-Za-z]", "", unname(seq))))
}

make_unique_ids <- function(x) {
  x <- gsub("[^A-Za-z0-9_.-]", "_", x)
  d <- stats::ave(seq_along(x), x, FUN = seq_along)
  ifelse(d > 1L, paste0(x, "_", d), x)
}

# Names of a UNITE classification ("k__Fungi;...;s__Entoloma_sp"), top rank
# first, NA where the name says nothing ("Entoloma sp", "Incertae sedis").
unite_levels <- function(cls) {
  lapply(strsplit(cls, ";", fixed = TRUE), function(v) {
    # "Thelephoraceae_gen_Incertae_sedis": a genus unknown within Thelephoraceae
    v <- sub("_(kgd|phy|cls|ord|fam|gen)_Incertae_sedis$", "", sub("^[a-z]+__", "", v))
    v <- gsub("_", " ", v)
    v[grepl("(?i)\\bsp\\.?$|incertae sedis|^unidentified|^unclassified|^unknown", v, perl = TRUE)] <- NA_character_
    v
  })
}

# Names of a SILVA classification ("Bacteria;...;Pseudomonas;uncultured
# bacterium"), with the same rules.
silva_levels <- function(cls) {
  lapply(strsplit(cls, ";", fixed = TRUE), function(v) {
    v <- trimws(v)
    v[grepl("(?i)^uncultured|^unidentified|^unclassified|^unknown|metagenome|incertae sedis|\\bsp\\.?$|^environmental", v, perl = TRUE)] <- NA_character_
    v
  })
}

# NCBI taxid of each classification (a vector of names from the top rank
# down; NA for levels to skip): the deepest name found in the taxonomy
# database; homonyms resolved by the other names of the classification.
names_to_taxids <- function(levels, taxonomy_db = NULL) {
  out <- data.frame(taxid = rep(NA_character_, length(levels)), how = NA_character_, name = NA_character_)
  if (length(levels) == 0L) return(out)
  key <- vapply(levels, paste, character(1L), collapse = ";")
  u <- !duplicated(key)
  ul <- levels[u]
  all_names <- unique(stats::na.omit(unlist(ul)))
  db <- resolve_taxonomy_db(taxonomy_db)
  con <- DBI::dbConnect(RSQLite::SQLite(), db, flags = RSQLite::SQLITE_RO)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  cand <- do.call(rbind, lapply(split(all_names, ceiling(seq_along(all_names) / 5000L)), function(p) {
    DBI::dbGetQuery(con, sprintf(
      "SELECT name, id, scientific, type FROM names WHERE name IN (%s)",
      paste(DBI::dbQuoteString(con, p), collapse = ",")
    ))
  }))
  cand <- cand[!grepl("common name|acronym", cand$type), , drop = FALSE]
  cand$id <- as.character(cand$id)
  lin <- if (nrow(cand) > 0L) suppressWarnings(get_tax(unique(cand$id), taxonomy_db = db)) else NULL
  lin_names <- if (is.null(lin)) list() else stats::setNames(strsplit(lin$lineage, "; ", fixed = TRUE), lin$taxid)
  by_name <- split(cand, cand$name)
  pick <- vapply(ul, function(v) {
    for (k in rev(seq_along(v))) {
      n <- v[k]
      if (is.na(n) || is.null(by_name[[n]])) next
      c1 <- by_name[[n]]
      if (nrow(c1) > 1L) {
        score <- vapply(c1$id, function(i) sum(v[-k] %in% lin_names[[i]]) + 0.5 * c1$scientific[c1$id == i][1L], numeric(1L))
        c1 <- c1[order(-score), , drop = FALSE]
      }
      return(paste(c1$id[1L], n, sep = "\t"))
    }
    NA_character_
  }, character(1L))
  res <- strsplit(pick, "\t", fixed = TRUE)
  ut <- vapply(res, `[`, character(1L), 1L)
  un <- vapply(res, `[`, character(1L), 2L)
  m <- match(key, key[u])
  out$taxid <- ut[m]
  out$name <- un[m]
  out$how <- ifelse(is.na(out$taxid), NA_character_, "name in the NCBI taxonomy")
  out
}

# Write FASTA + taxid map, run makeblastdb, write the taxid table.
write_ref_db <- function(out_dir, name, ids, title, seq, tax, extra = NULL) {
  fasta <- fs::path(out_dir, paste0(name, ".fasta"))
  map <- fs::path(out_dir, paste0(name, ".txt"))
  ok <- nzchar(seq)
  writeLines(paste0(">", ids[ok], " ", title[ok], "\n", seq[ok]), fasta)
  tid <- ifelse(is.na(tax$taxid) | tax$taxid == "", "0", tax$taxid)
  writeLines(paste(ids[ok], tid[ok]), map)
  report <- data.frame(id = ids, taxid = tid, how = tax$how, found_by = tax$name)
  if (!is.null(extra)) report <- cbind(report, extra)
  readr::write_tsv(report, fs::path(out_dir, paste0(name, "_taxids.tsv")), na = "")
  n0 <- sum(tid[ok] == "0")
  if (n0 > 0L) cli::cli_warn(c(`!` = "{n0} of {sum(ok)} sequences without an NCBI taxid (see {.file {name}_taxids.tsv})."))
  db <- fs::path(out_dir, name)
  make_blast_db(fasta_path = fasta, db_path = db, db_type = "nucl", taxid_map = map, parse_seqids = TRUE)
  as.character(db)
}

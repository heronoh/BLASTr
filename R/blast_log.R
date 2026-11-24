# Optional run log: one line per parallel_blast() search, appended to
# <log_dir>/BLASTr_logs/blastr_runs.tsv. Used to follow the cost of the
# searches and as a complement to the BLASTr benchmark.

#' @title Log of BLAST searches
#'
#' @description
#' When a log directory is set, each search made by [parallel_blast()] adds
#' one line to `BLASTr_logs/blastr_runs.tsv` inside it (the folder and file
#' are created when needed). The directory is taken, in this order, from the
#' `log_dir` argument of [parallel_blast()], the option `BLASTr.log_dir` and
#' the environment variable `BLASTR_LOG_DIR`; with none of them, nothing is
#' logged. In background jobs, the line is written when the results are
#' collected with [blast_job_result()].
#'
#' Each line has:
#' * the run: `run_id`, `start`, `end`, `host`, `user`, `mode`
#'   (foreground/background), `function` ([parallel_blast()] or
#'   [parallel_blast_hsps()]), `blastr_version`, `blast_version`, `job_dir`;
#' * the sequences: `n_input` (as given), `n_seqs` (unique ones searched),
#'   `n_invalid` and their length distribution (`len_min`, `len_q25`,
#'   `len_median`, `len_mean`, `len_q75`, `len_max`);
#' * the database(s): `db_path`, `n_dbs`, `db_gb` (sum), `db_volumes`,
#'   `ram_gb`, `fits_ram`;
#' * the strategy and settings: `total_cores`, `processes`, `threads`,
#'   `mt_mode`, `batch_size`, `n_batches`, `task`, `num_alignments`,
#'   `search_alignments`, `perc_id`,
#'   `perc_qcov_hsp`, `check_n`, `check_tolerance`;
#' * the cost: `wall_s`, `cpu_s` (CPU time of the BLAST processes),
#'   `cpu_h`, `cpu_efficiency` (`cpu_s` over `wall_s` times the threads in
#'   use), `n_batches_failed`, `n_seqs_failed`;
#' * the identifications: `n_with_hits`, `pct_with_hits`, `pct_best_id`,
#'   `pct_check_id` (both over all sequences), `id1_median` (median identity
#'   of the first hit) and `pct_id1_ge_99`, `pct_id1_ge_97` (sequences whose
#'   first hit has identity of at least 99% / 97%). For
#'   [parallel_blast_hsps()], `pct_best_id` and `pct_check_id` are `NA` and
#'   the first hit is the first HSP of the best subject.
#'
#' When a newer BLASTr adds columns, the file is rewritten with all of them
#' (missing values as `NA`).
#'
#' @param log_dir Directory where the `BLASTr_logs` folder is kept. `NULL`
#'   uses the option and the environment variable described above.
#'
#' @returns `blast_log_path()` returns the path of the log file, or `NULL`
#'   when no log directory is set. `read_blast_log()` returns the log as a
#'   tibble.
#'
#' @examples
#' log_dir <- tempfile("blastr_logs_")
#' blast_log_path(log_dir)
#'
#' # needs BLAST+ (see install_dependencies())
#' if (has_blast()) {
#'   fasta <- system.file("extdata", "minimal_db_blast.fasta", package = "BLASTr")
#'   db <- file.path(tempdir(), "minimal_db_blast")
#'   make_blast_db(fasta, db_path = db)
#'   seqs <- readLines(system.file("extdata", "seqs_test.txt", package = "BLASTr"))
#'   res <- parallel_blast(seqs, db, log_dir = log_dir)
#'   read_blast_log(log_dir)
#' }
#'
#' @family search functions
#' @export
blast_log_path <- function(log_dir = NULL) {
  log_dir <- resolve_log_dir(log_dir)
  if (is.null(log_dir)) {
    return(NULL)
  }
  fs::path(log_dir, "BLASTr_logs", "blastr_runs.tsv")
}

#' @rdname blast_log_path
#' @export
read_blast_log <- function(log_dir = NULL) {
  path <- blast_log_path(log_dir)
  if (is.null(path) || !fs::file_exists(path)) {
    cli::cli_abort(
      c(`x` = "No BLASTr log found{if (!is.null(path)) paste0(' at ', path)}."),
      class = "blastr_no_log_error"
    )
  }
  readr::read_tsv(path, show_col_types = FALSE, progress = FALSE, quote = "")
}

resolve_log_dir <- function(log_dir = NULL) {
  if (is.null(log_dir)) {
    log_dir <- getOption("BLASTr.log_dir", default = NULL)
  }
  if (is.null(log_dir)) {
    log_dir <- Sys.getenv("BLASTR_LOG_DIR", unset = "")
  }
  if (is.null(log_dir) || is.na(log_dir) || !nzchar(log_dir)) {
    return(NULL)
  }
  as.character(fs::path_abs(fs::path_expand(log_dir)))
}

# CPU time (s) of the blastn processes of a job, from the usage lines the
# scheduler writes at each run (several when the job was resumed).
job_cpu_s <- function(job_dir) {
  f <- fs::path(job_dir, "usage.tsv")
  if (!fs::file_exists(f)) {
    return(NA_real_)
  }
  u <- utils::read.delim(f)
  sum(u$child_user_s + u$child_sys_s)
}

blast_version <- function(blast_bin) {
  v <- tryCatch(
    processx::run(blast_bin, "-version", error_on_status = FALSE, timeout = 30)$stdout,
    error = function(e) ""
  )
  v <- sub("^[^:]*:\\s*", "", strsplit(v, "\n", fixed = TRUE)[[1L]][1L])
  if (length(v) == 0L || is.na(v)) NA_character_ else trimws(v)
}

pct <- function(x, n) if (n > 0) round(100 * x / n, 2) else NA_real_

# Add the line of a finished job to the log. Called once per job.
write_blast_log <- function(job, res, st, b, seqs, summ = NULL) {
  path <- blast_log_path(st$meta$log_dir)
  if (is.null(path)) {
    return(invisible(NULL))
  }
  marker <- fs::path(job$dir, "logged")
  if (fs::file_exists(marker)) {
    return(invisible(path))
  }
  s <- st$strategy
  hsps <- identical(st$meta$output, "hsps")
  len <- nchar(seqs$sequence)
  final <- b[b$status != "split", , drop = FALSE]
  n <- nrow(seqs)
  if (hsps) {
    with_hits <- seqs$sequence %in% res$Sequence
    first <- res[!is.na(res$hit) & res$hit == 1L & res$hsp == 1L, , drop = FALSE]
    id1 <- if ("pident" %in% names(res)) first$pident[match(seqs$sequence, first$Sequence)] else rep(NA_real_, n)
    check <- rep(NA_character_, n)
  } else {
    k <- match(seqs$sequence, summ$Sequence)
    with_hits <- !is.na(summ$`N hits`[k]) & summ$`N hits`[k] > 0L
    check <- summ$`Check ID`[k]
    id1 <- if ("1_identity" %in% names(res)) {
      res$`1_identity`[match(seqs$sequence, res$Sequence)]
    } else {
      # long: the first row of each sequence is its best hit
      first <- res[!duplicated(res$Sequence), , drop = FALSE]
      first$identity[match(seqs$sequence, first$Sequence)]
    }
  }
  threads_used <- s$processes * (if (is.na(s$threads)) 1L else s$threads)
  wall <- max(b$end, na.rm = TRUE) - min(b$start, na.rm = TRUE)
  if (!is.finite(wall)) wall <- NA_real_
  cpu <- job_cpu_s(job$dir)
  q <- if (n > 0L) stats::quantile(len, c(0.25, 0.5, 0.75), names = FALSE) else rep(NA_real_, 3L)

  row <- data.frame(
    check.names = FALSE,
    run_id = paste0(
      format(as.POSIXct(min(b$start, na.rm = TRUE), origin = "1970-01-01"), "%Y%m%d-%H%M%S"),
      "_", basename(job$dir)
    ),
    start = format(as.POSIXct(min(b$start, na.rm = TRUE), origin = "1970-01-01"), "%Y-%m-%d %H:%M:%S"),
    end = format(as.POSIXct(max(b$end, na.rm = TRUE), origin = "1970-01-01"), "%Y-%m-%d %H:%M:%S"),
    host = Sys.info()[["nodename"]],
    user = Sys.info()[["user"]],
    mode = if (isTRUE(st$meta$background)) "background" else "foreground",
    "function" = if (hsps) "parallel_blast_hsps" else "parallel_blast",
    blastr_version = as.character(utils::packageVersion("BLASTr")),
    blast_version = if (is.null(st$meta$blast_version)) blast_version(st$blast_bin) else st$meta$blast_version,
    job_dir = as.character(job$dir),
    n_input = if (is.null(st$meta$n_input)) NA_integer_ else st$meta$n_input,
    n_seqs = n,
    n_invalid = sum(!is_valid_seq(seqs$sequence)),
    len_min = if (n > 0L) min(len) else NA_integer_,
    len_q25 = q[1L],
    len_median = q[2L],
    len_mean = if (n > 0L) round(mean(len), 1) else NA_real_,
    len_q75 = q[3L],
    len_max = if (n > 0L) max(len) else NA_integer_,
    db_path = paste(s$db_path, collapse = " "),
    db_date = paste(if (is.null(st$meta$db_dates)) db_dates(s$db_path) else st$meta$db_dates, collapse = " "),
    n_dbs = length(s$db_path),
    db_gb = round(s$db_gb, 3),
    db_volumes = s$db_volumes,
    ram_gb = round(s$ram_gb, 1),
    fits_ram = s$fits_ram,
    total_cores = if (is.null(st$meta$total_cores)) NA_integer_ else st$meta$total_cores,
    processes = s$processes,
    threads = s$threads,
    mt_mode = s$mt_mode,
    batch_size = s$batch_size,
    n_batches = s$n_batches,
    task = if (is.null(st$meta$task)) NA_character_ else st$meta$task,
    num_alignments = st$meta$num_alignments,
    search_alignments = if (is.null(st$meta$search_alignments)) NA_integer_ else st$meta$search_alignments,
    perc_id = st$meta$perc_id,
    perc_qcov_hsp = st$meta$perc_qcov_hsp,
    search_qcov_hsp = if (is.null(st$meta$search_qcov_hsp)) NA_real_ else st$meta$search_qcov_hsp,
    check_n = st$meta$check_n,
    check_tolerance = st$meta$check_tolerance,
    rank_by = if (is.null(st$meta$rank_by)) NA_character_ else st$meta$rank_by,
    output = if (is.null(st$meta$output)) NA_character_ else st$meta$output,
    uninformative_md5 = if (is.null(st$meta$uninformative_md5)) NA_character_ else st$meta$uninformative_md5,
    wall_s = round(wall, 1),
    cpu_s = round(cpu, 1),
    cpu_h = round(cpu / 3600, 4),
    cpu_efficiency = round(cpu / (wall * threads_used), 3),
    n_batches_failed = sum(final$status == "failed"),
    n_seqs_failed = sum(final$n[final$status == "failed"]),
    n_with_hits = sum(with_hits),
    pct_with_hits = pct(sum(with_hits), n),
    pct_best_id = if (hsps) NA_real_ else pct(sum(check %in% "best ID"), n),
    pct_check_id = if (hsps) NA_real_ else pct(sum(check %in% "check ID"), n),
    id1_median = if (any(!is.na(id1))) stats::median(id1, na.rm = TRUE) else NA_real_,
    pct_id1_ge_99 = pct(sum(id1 >= 99, na.rm = TRUE), n),
    pct_id1_ge_97 = pct(sum(id1 >= 97, na.rm = TRUE), n),
    blast_cmd = if (is.null(st$meta$blast_cmd)) NA_character_ else st$meta$blast_cmd
  )

  fs::dir_create(fs::path_dir(path))
  if (fs::file_exists(path)) {
    header <- strsplit(readLines(path, n = 1L), "\t", fixed = TRUE)[[1L]]
    if (!identical(header, names(row))) {
      # columns changed between BLASTr versions: rewrite with all of them
      old <- utils::read.delim(path, colClasses = "character", check.names = FALSE, quote = "", na.strings = "NA")
      row[] <- lapply(row, as.character)
      row <- dplyr::bind_rows(old, row)
      fs::file_delete(path)
    }
  }
  new_file <- !fs::file_exists(path)
  utils::write.table(
    row, path,
    sep = "\t", quote = FALSE, row.names = FALSE, na = "NA",
    col.names = new_file, append = !new_file
  )
  fs::file_create(marker)
  invisible(path)
}

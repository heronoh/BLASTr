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
#'   (foreground/background), `blastr_version`, `blast_version`,
#'   `job_dir`;
#' * the sequences: `n_input` (as given), `n_seqs` (unique ones searched),
#'   `n_invalid` and their length distribution (`len_min`, `len_q25`,
#'   `len_median`, `len_mean`, `len_q75`, `len_max`);
#' * the database(s): `db_path`, `n_dbs`, `db_gb` (sum), `db_volumes`,
#'   `ram_gb`, `fits_ram`;
#' * the strategy and settings: `total_cores`, `processes`, `threads`,
#'   `mt_mode`, `batch_size`, `n_batches`, `num_alignments`, `perc_id`,
#'   `perc_qcov_hsp`, `check_n`, `check_tolerance`;
#' * the cost: `wall_s`, `cpu_s` (CPU time of the BLAST processes),
#'   `cpu_h`, `cpu_efficiency` (`cpu_s` over `wall_s` times the threads in
#'   use), `n_batches_failed`, `n_seqs_failed`;
#' * the identifications: `n_with_hits`, `pct_with_hits`, `pct_best_id`,
#'   `pct_check_id` (both over all sequences), `id1_median` (median identity
#'   of the first hit) and `pct_id1_ge_99`, `pct_id1_ge_97` (sequences whose
#'   first hit has identity of at least 99% / 97%).
#'
#' @param log_dir Directory where the `BLASTr_logs` folder is kept. `NULL`
#'   uses the option and the environment variable described above.
#'
#' @returns `blast_log_path()` returns the path of the log file, or `NULL`
#'   when no log directory is set. `read_blast_log()` returns the log as a
#'   tibble.
#'
#' @examples
#' \dontrun{
#' options(BLASTr.log_dir = "~/analyses")
#' res <- parallel_blast(query_seqs = seqs, db_path = db, total_cores = 8)
#' read_blast_log()
#' }
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
  readr::read_tsv(path, show_col_types = FALSE, progress = FALSE)
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
write_blast_log <- function(job, res, st, b) {
  path <- blast_log_path(st$meta$log_dir)
  if (is.null(path)) {
    return(invisible(NULL))
  }
  marker <- fs::path(job$dir, "logged")
  if (fs::file_exists(marker)) {
    return(invisible(path))
  }
  s <- st$strategy
  len <- nchar(res$Sequence)
  final <- b[b$status != "split", , drop = FALSE]
  n <- nrow(res)
  with_hits <- !is.na(res$`N hits`) & res$`N hits` > 0L
  id1 <- res$`1_identity`
  wall <- max(b$end, na.rm = TRUE) - min(b$start, na.rm = TRUE)
  if (!is.finite(wall)) wall <- NA_real_
  cpu <- job_cpu_s(job$dir)
  q <- if (n > 0L) stats::quantile(len, c(0.25, 0.5, 0.75), names = FALSE) else rep(NA_real_, 3L)

  row <- data.frame(
    run_id = paste0(
      format(as.POSIXct(min(b$start, na.rm = TRUE), origin = "1970-01-01"), "%Y%m%d-%H%M%S"),
      "_", basename(job$dir)
    ),
    start = format(as.POSIXct(min(b$start, na.rm = TRUE), origin = "1970-01-01"), "%Y-%m-%d %H:%M:%S"),
    end = format(as.POSIXct(max(b$end, na.rm = TRUE), origin = "1970-01-01"), "%Y-%m-%d %H:%M:%S"),
    host = Sys.info()[["nodename"]],
    user = Sys.info()[["user"]],
    mode = if (isTRUE(st$meta$background)) "background" else "foreground",
    blastr_version = as.character(utils::packageVersion("BLASTr")),
    blast_version = blast_version(st$blast_bin),
    job_dir = as.character(job$dir),
    n_input = if (is.null(st$meta$n_input)) NA_integer_ else st$meta$n_input,
    n_seqs = n,
    n_invalid = sum(!is_valid_seq(res$Sequence)),
    len_min = if (n > 0L) min(len) else NA_integer_,
    len_q25 = q[1L],
    len_median = q[2L],
    len_mean = if (n > 0L) round(mean(len), 1) else NA_real_,
    len_q75 = q[3L],
    len_max = if (n > 0L) max(len) else NA_integer_,
    db_path = paste(s$db_path, collapse = " "),
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
    num_alignments = st$meta$num_alignments,
    perc_id = st$meta$perc_id,
    perc_qcov_hsp = st$meta$perc_qcov_hsp,
    check_n = st$meta$check_n,
    check_tolerance = st$meta$check_tolerance,
    wall_s = round(wall, 1),
    cpu_s = round(cpu, 1),
    cpu_h = round(cpu / 3600, 4),
    cpu_efficiency = round(cpu / (wall * s$processes * s$threads), 3),
    n_batches_failed = sum(final$status == "failed"),
    n_seqs_failed = sum(final$n[final$status == "failed"]),
    n_with_hits = sum(with_hits),
    pct_with_hits = pct(sum(with_hits), n),
    pct_best_id = pct(sum(res$`Check ID` %in% "best ID"), n),
    pct_check_id = pct(sum(res$`Check ID` %in% "check ID"), n),
    id1_median = if (any(!is.na(id1))) stats::median(id1, na.rm = TRUE) else NA_real_,
    pct_id1_ge_99 = pct(sum(id1 >= 99, na.rm = TRUE), n),
    pct_id1_ge_97 = pct(sum(id1 >= 97, na.rm = TRUE), n)
  )

  fs::dir_create(fs::path_dir(path))
  new_file <- !fs::file_exists(path)
  utils::write.table(
    row, path,
    sep = "\t", quote = FALSE, row.names = FALSE,
    col.names = new_file, append = !new_file
  )
  fs::file_create(marker)
  invisible(path)
}

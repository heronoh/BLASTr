# Default fields of parallel_blast_hsps(): the 12 "std" ones plus those
# needed to rebuild alignments of split queries (nident, sstrand, qlen) and
# the taxonomy of the subject.
hsp_fields_default <- c(
  "qaccver", "saccver", "pident", "length", "mismatch", "gapopen",
  "qstart", "qend", "sstart", "send", "evalue", "bitscore",
  "nident", "qlen", "slen", "qcovhsp", "sstrand", "staxid", "stitle"
)

# BLAST tabular fields that are numbers.
hsp_numeric_fields <- c(
  "pident", "length", "mismatch", "gapopen", "gaps", "qstart", "qend",
  "sstart", "send", "evalue", "bitscore", "score", "nident", "positive",
  "ppos", "qlen", "slen", "qcovs", "qcovhsp", "qcovus", "frames", "qframe",
  "sframe", "btop_len"
)

#' @title Search sequences with BLAST and return every HSP
#'
#' @description
#' Like [parallel_blast()], but returns the raw BLAST hits in long format:
#' one row per HSP (high-scoring pair), with the BLAST fields chosen in
#' `fields`, instead of one row per sequence. Nothing is filtered or cut
#' beyond what is asked: by default there is no limit of HSPs per subject
#' (`max_hsps`) and no coverage filter (`perc_qcov_hsp`). Use it when the
#' alignments themselves are needed, e.g. to rebuild a hit from the two
#' ends of a concatenated sequence.
#'
#' The search uses the same batches, automatic strategy, progress/ETA,
#' background jobs and log as [parallel_blast()]. Instead of a database,
#' `subject` searches against the sequences of a FASTA file (BLAST
#' `-subject`, no database needed); then each process uses one thread.
#'
#' @inheritParams parallel_blast
#' @param db_path Path(s) to the formatted BLAST database(s). Give either
#'   `db_path` or `subject`.
#' @param subject Path to a FASTA file to search against instead of a
#'   database (BLAST `-subject`).
#' @param fields BLAST tabular output fields (see `blastn -help`, `-outfmt`).
#'   Numeric fields are returned as numbers.
#' @param max_target_seqs Number of subjects kept per sequence. Defaults to
#'   `50`.
#' @param max_hsps Number of HSPs kept per subject. `NULL` (default): no
#'   limit.
#' @param perc_id,perc_qcov_hsp,evalue Lowest identity and query coverage per
#'   HSP (percent) and highest e-value. `NULL` (default): no filter (BLAST
#'   default for `evalue`).
#' @param extra_args Further BLAST arguments, as a character vector, e.g.
#'   `c("-dust", "no")`.
#'
#' @returns A tibble with `Sequence`, `hit` (rank of the subject for that
#'   sequence), `hsp` (rank of the HSP within the subject) and one column per
#'   field. Sequences without hits do not appear. The attribute
#'   `BLASTr_metadata` is as in [parallel_blast()]. With
#'   `background = TRUE`, a `blastr_job` (see [blast_job()]).
#'
#' @examples
#' # needs BLAST+ (see install_dependencies())
#' if (nzchar(Sys.which("blastn"))) {
#'   fasta <- system.file("extdata", "minimal_db_blast.fasta", package = "BLASTr")
#'   db <- file.path(tempdir(), "minimal_db_blast")
#'   make_blast_db(fasta, db_path = db)
#'   seqs <- readLines(system.file("extdata", "seqs_test.txt", package = "BLASTr"))
#'   hsps <- parallel_blast_hsps(seqs, db, perc_id = 80)
#'   hsps[, c("Sequence", "saccver", "pident", "qcovhsp", "evalue")]
#'
#'   # against the sequences of a FASTA file, as BLAST -subject
#'   parallel_blast_hsps(seqs[1], subject = fasta, task = "blastn")
#' }
#'
#' @family search functions
#' @export
parallel_blast_hsps <- function(
  query_seqs,
  db_path = NULL,
  ...,
  subject = NULL,
  total_cores = 1L,
  ram_gb = NULL,
  fields = hsp_fields_default,
  task = NULL,
  max_target_seqs = 50L,
  max_hsps = NULL,
  perc_id = NULL,
  perc_qcov_hsp = NULL,
  evalue = NULL,
  extra_args = character(0L),
  num_threads = NULL,
  batch_size = NULL,
  mt_mode = NULL,
  background = FALSE,
  job_dir = NULL,
  log_dir = NULL,
  retry_times = 3L,
  blast_type = "blastn",
  verbose = c("progress", "silent", "cmd", "output", "full"),
  env_name = "blastr-blast-env"
) {
  rlang::check_dots_empty()
  rlang::check_required(query_seqs)
  verbose <- rlang::arg_match(verbose)
  if (
    rlang::is_null(query_seqs) ||
      isTRUE(length(query_seqs) == 0L) ||
      all(is.na(query_seqs))
  ) {
    cli::cli_abort(
      c(`x` = "No query sequences were provided to {.fun parallel_blast_hsps}."),
      class = "blastr_no_query_seqs_error"
    )
  }
  if (is.null(db_path) == is.null(subject)) {
    cli::cli_abort(
      c(`x` = "Give either {.arg db_path} or {.arg subject}, not both or none."),
      class = "blastr_no_db_path_error"
    )
  }
  if (!is.null(subject) && !fs::file_exists(subject)) {
    cli::cli_abort(
      c(`x` = "The subject file {.path {subject}} does not exist."),
      class = "blastr_no_db_path_error"
    )
  }
  if (!is.null(mt_mode)) {
    mt_mode <- rlang::arg_match(as.character(mt_mode), c("0", "1", "2"))
  }
  fields <- unique(as.character(fields))

  check_cmd(blast_type, env_name = env_name, verbose = "silent")
  seqs <- unique(stringr::str_replace_all(query_seqs[!is.na(query_seqs)], "\\s", ""))

  if (is.null(subject)) {
    strategy <- blast_strategy(
      n_seqs = length(seqs),
      db_path = db_path,
      total_cores = total_cores,
      ram_gb = ram_gb,
      num_threads = num_threads,
      batch_size = batch_size,
      mt_mode = mt_mode
    )
    target <- c(
      "-db", db_arg(strategy$db_path),
      "-mt_mode", as.character(strategy$mt_mode)
    )
  } else {
    strategy <- subject_strategy(
      n_seqs = length(seqs),
      subject = subject,
      total_cores = total_cores,
      batch_size = batch_size
    )
    target <- c("-subject", as.character(fs::path_abs(subject)))
  }
  if (!identical(verbose, "silent")) {
    cli::cli_inform(c(`i` = "Running {.fun parallel_blast_hsps} for {length(seqs)} unique query sequences."))
    print(strategy)
  }

  # the first field identifies the query (BLASTr's own ids, s000001, ...)
  blast_args <- c(
    target,
    "-outfmt", paste("6", paste(c("qaccver", fields), collapse = " ")),
    "-max_target_seqs", as.character(max_target_seqs),
    if (!is.null(max_hsps)) c("-max_hsps", as.character(max_hsps)),
    if (!is.null(perc_id)) c("-perc_identity", as.character(perc_id)),
    if (!is.null(perc_qcov_hsp)) c("-qcov_hsp_perc", as.character(perc_qcov_hsp)),
    if (!is.null(evalue)) c("-evalue", as.character(evalue)),
    if (!is.null(task)) c("-task", task),
    as.character(extra_args)
  )
  meta <- list(
    db_path = strategy$db_path,
    output = "hsps",
    fields = fields,
    num_alignments = as.integer(max_target_seqs),
    search_alignments = as.integer(max_target_seqs),
    task = if (is.null(task)) NA_character_ else task,
    perc_id = if (is.null(perc_id)) NA_real_ else perc_id,
    perc_qcov_hsp = if (is.null(perc_qcov_hsp)) NA_real_ else perc_qcov_hsp,
    blast_type = blast_type,
    check_n = NA_integer_,
    check_tolerance = NA_real_,
    n_input = length(query_seqs),
    total_cores = as.integer(total_cores),
    background = isTRUE(background),
    log_dir = resolve_log_dir(log_dir)
  )

  run_blast_job(
    seqs = seqs,
    strategy = strategy,
    blast_type = blast_type,
    blast_args = blast_args,
    retry_times = retry_times,
    env_name = env_name,
    job_dir = job_dir,
    meta = meta,
    background = background,
    verbose = verbose
  )
}

# Strategy for -subject searches: no database to split, BLAST ignores
# -num_threads, so one single-threaded process per core.
subject_strategy <- function(n_seqs, subject, total_cores = 1L, batch_size = NULL, max_batch_size = 1000L) {
  n_seqs <- max(1L, as.integer(n_seqs))
  processes <- min(max(1L, as.integer(total_cores)), n_seqs)
  plan <- batch_plan(n_seqs, processes, batch_size, max_batch_size)
  structure(
    list(
      processes = as.integer(min(processes, plan$n_batches)),
      threads = NA_integer_,
      mt_mode = NA_integer_,
      batch_size = plan$batch_size,
      n_batches = plan$n_batches,
      n_seqs = n_seqs,
      db_path = as.character(fs::path_abs(subject)),
      db_gb = as.numeric(fs::file_size(subject)) / 1e9,
      db_volumes = 1L,
      ram_gb = NA_real_,
      fits_ram = TRUE,
      reason = "subject FASTA (-subject): one single-threaded process per core"
    ),
    class = "blastr_strategy"
  )
}

# Read the tabular outputs of a parallel_blast_hsps() job.
read_blast_hsps <- function(out_files, seqs, fields) {
  cols <- c(".query", fields)
  out_files <- out_files[file.exists(out_files) & file.size(out_files) > 0]
  empty <- tibble::as_tibble(stats::setNames(rep(list(character(0L)), length(cols)), cols))
  h <- lapply(out_files, function(f) {
    readr::read_tsv(
      f,
      col_names = cols,
      col_types = readr::cols(.default = readr::col_character()),
      quote = "", comment = "", na = character(0L), trim_ws = TRUE,
      progress = FALSE
    )
  })
  h <- dplyr::bind_rows(c(list(empty), h))
  for (k in intersect(fields, hsp_numeric_fields)) {
    h[[k]] <- suppressWarnings(as.numeric(h[[k]]))
  }
  # BLAST lists subjects from best to worst, and the HSPs of each subject
  # together
  subj <- if ("saccver" %in% fields) h$saccver else if ("sseqid" %in% fields) h$sseqid else NULL
  if (is.null(subj)) {
    hit <- rep(NA_integer_, nrow(h))
    hsp <- stats::ave(seq_len(nrow(h)), h$.query, FUN = seq_along)
  } else {
    key <- paste(h$.query, subj, sep = "\r")
    first <- !duplicated(key)
    rank_first <- stats::ave(as.integer(first), h$.query, FUN = cumsum)
    hit <- rank_first[first][match(key, key[first])]
    hsp <- stats::ave(seq_len(nrow(h)), key, FUN = seq_along)
  }
  out <- tibble::tibble(
    Sequence = seqs$sequence[match(h$.query, seqs$query)],
    hit = as.integer(hit),
    hsp = as.integer(hsp)
  )
  dplyr::bind_cols(out, h[fields])
}

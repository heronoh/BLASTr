#' @title Run Parallelized BLAST
#'
#' @description
#' Search a set of sequences against one or more BLAST databases, using the
#' cores and RAM given. The search strategy is chosen automatically (see
#' [blast_strategy()]): the sequences are grouped in batches, and each
#' `blastn` process searches a whole batch, so the database is read once per
#' batch instead of once per sequence.
#'
#' The progress bar shows the estimated time left (ETA), computed from the
#' batches already finished. With `background = TRUE` the search runs in a
#' separate process and the function returns at once; see [blast_job()].
#'
#' @param query_seqs Character vector with sequences.
#' @param db_path Path(s) to the formatted BLAST database(s). Several
#'   databases (a vector, or one string separated by spaces) are searched
#'   together, the smallest first.
#' @param total_cores Number of CPU cores (threads) to use in total.
#' @param ram_gb RAM available for the search, in GB. Used to decide whether
#'   the database fits in memory. `NULL` uses the total RAM of the machine.
#' @param num_alignments Number of hits to keep for each sequence.
#'   Defaults to `10`.
#' @param search_alignments Number of subjects BLAST keeps during the
#'   search (`-max_target_seqs`). BLAST cuts its list of subjects before the
#'   `perc_qcov_hsp` filter, so asking only for `num_alignments` can lose
#'   valid hits. BLASTr searches `max(num_alignments, search_alignments)`
#'   subjects and keeps the `num_alignments` best ones that pass the
#'   filters. In the benchmark, 50 instead of 10 added about 5% to the time.
#' @param perc_id Lowest identity percentage cutoff.
#' @param perc_qcov_hsp Query coverage per HSP (percent) expected of a
#'   good hit. Hits below it are kept and flagged (`low_qcov`), not
#'   removed: removing them inside BLAST hides a closer reference that covers
#'   a bit less of the query (common with partial reference sequences) and
#'   promotes a worse one. The coverage still weighs in the rank through the
#'   Similarity score. Defaults to `80`.
#' @param search_qcov_hsp Lowest query coverage per HSP (percent) used by
#'   BLAST (`-qcov_hsp_perc`); hits below it are not returned. Defaults to
#'   `50`. Set it to `perc_qcov_hsp` to remove the hits below
#'   `perc_qcov_hsp`, as in earlier versions.
#' @param check_n Number of best hits compared by `Check ID`
#'   (see [check_hits()]). Defaults to `10`.
#' @param check_tolerance Difference from the best hit, in the units of
#'   `rank_by`, under which hits are tied ("at best"). Defaults to `0.5`.
#' @param rank_by Criterion to rank the hits of each sequence:
#'   `"similarity"` (default, the Similarity score, see
#'   [similarity_score()]), `"identity"` or `"bitscore"`. Ties are broken by
#'   identity, coverage, bitscore and BLAST order. BLAST orders its hits by
#'   e-value/bitscore, which favours longer alignments, so its first hit is
#'   not always the most similar one. In a calibration on plant barcodes of
#'   known genus (ITS2, rbcL) and on 1257 GenBank references of six markers
#'   searched without themselves, ranking by Similarity score or identity
#'   picked the right taxon significantly more often than bitscore (genus:
#'   72% and 76% against 62%), with no significant difference between them.
#' @param uninformative Table of the references that do not inform the
#'   identification (uncultured, environmental, metagenomes, vectors, ...;
#'   see [uninformative_references()]): the table shipped with BLASTr
#'   (default), a path to another one, a data frame, or `NULL` to treat all
#'   hits alike (`informative` NA, not evaluated). The hits are flagged (`informative`,
#'   `uninformative_category`) and ranked among the informative ones, so
#'   `rank`, `at_best`, `Check ID` and [blast_lca()] use only those; a
#'   sequence whose hits are all uninformative keeps its best ones, ranked
#'   among themselves. Uninformative hits take subjects of the search, so a
#'   sequence of a poorly known group needs a larger `search_alignments`.
#' @param output `"long"` (default): one row per hit. `"wide"`: one row per
#'   sequence, with the summary of [blast_summary()] and the first
#'   `num_alignments` hits side by side (`1_` is the best one).
#' @param num_threads,batch_size,mt_mode Override the automatic strategy:
#'   threads per `blastn` process, sequences per batch and BLAST+
#'   multithreading mode (`0`, `1` or `2`). `NULL` (default) chooses
#'   automatically.
#' @param task BLAST task (`-task`), e.g. `"megablast"` (the `blastn`
#'   default), `"dc-megablast"` or `"blastn"`. `NULL` uses the BLAST default.
#' @param background If `TRUE`, start the search in the background and return
#'   a job (see [blast_job()]) instead of the results.
#' @param job_dir Directory for the job files. Defaults to a temporary
#'   directory (removed at the end) in the foreground, and to a new directory
#'   in the BLASTr cache in the background.
#' @param log_dir Directory where a `BLASTr_logs` folder keeps a log of the
#'   searches (number and lengths of the sequences, database size, strategy,
#'   time, CPU time and percentage of identifications). `NULL` uses the
#'   option `BLASTr.log_dir` or the environment variable `BLASTR_LOG_DIR`;
#'   with none set, nothing is logged. See [blast_log_path()].
#' @param retry_times Number of times to retry failed BLAST jobs.
#'   Defaults to `3` attempts.
#' @param blast_type BLAST+ executable to be used on search.
#' @param verbose `"progress"` (default) shows the progress and ETA;
#'   `"silent"` shows nothing. `"cmd"`, `"output"` and `"full"` also show the
#'   chosen strategy.
#' @param env_name The name of the conda environment used to run
#'  command-line tools. Defaults to `"blastr-blast-env"`.
#'
#' @param out_file Deprecated. Path to output `.csv` file on an existing
#' directory.
#' @param out_RDS Deprecated. Path to output `RDS` file on an existing
#' directory.
#'
#' @inheritParams rlang::args_dots_empty
#'
#' @returns With `output = "long"`, a tibble with one row per hit, best
#'   first: `Sequence`, `rank` (1 = best; tied hits share the rank), `tie`
#'   (the rank is shared), `at_best` (within `check_tolerance` of the best
#'   hit), `low_qcov` (`qcovhsp` below `perc_qcov_hsp`), `informative`
#'   and `uninformative_category` (see `uninformative`), `blast_rank`
#'   (order given by BLAST), `subject`, `subject header`,
#'   `staxid`, `identity`, `Similarity score`, `qcovhsp`, `nident`,
#'   `length`, `mismatches`, `gaps`, `query start`, `query end`,
#'   `subject start`, `subject end`, `subject strand`, `e-value` and
#'   `bitscore`. The first `num_alignments` hits are kept (informative
#'   ones, when there are), plus all the hits at best and the uninformative
#'   hits among the first `num_alignments` of all (after the ranked ones,
#'   with `rank` `NA`). A sequence without hits has one row with `NA` hit fields. The
#'   summary per sequence (`Check ID`, ...) comes from [blast_summary()].
#'
#'   With `output = "wide"`, one row per sequence: `Sequence`, the summary
#'   of [blast_summary()] and, for each hit `k` from 1 to `num_alignments`
#'   in rank order, the hit fields prefixed by `k_`.
#'
#'   The attribute `BLASTr_metadata` keeps the parameters, the strategy, the
#'   BLAST command (`blast_cmd`, without the sequences), the BLAST version,
#'   the date of each database and the exit codes (see [exit_codes()]). With
#'   `background = TRUE`, a `blastr_job`.
#'
#' @examples
#' \dontrun{
#' dna_fasta_path <- fs::path_package(
#'   "BLASTr", "extdata", "minimal_db_blast",
#'   ext = "fasta"
#' )
#' temp_db_path <- fs::path_temp("minimal_db_blast")
#' make_blast_db(fasta_path = dna_fasta_path, db_path = temp_db_path)
#' query_seqs_string <- c(
#'   "CTAGCCATAAACTTAAATGAAGCTATACTAA",
#'   "ACTCGTTCGCCAGAGTACTACAAGCGAAAG"
#' )
#' blast_res <- parallel_blast(
#'   query_seqs = query_seqs_string,
#'   db_path = temp_db_path
#' )
#' blast_res
#'
#' # in the background
#' job <- parallel_blast(query_seqs_string, temp_db_path, background = TRUE)
#' blast_job_status(job)
#' blast_res <- blast_job_result(job)
#' }
#'
#' @importFrom lifecycle deprecated
#'
#' @export
parallel_blast <- function(
  query_seqs,
  db_path,
  ...,
  total_cores = 1L,
  ram_gb = NULL,
  num_alignments = 10L,
  search_alignments = 50L,
  perc_id = 80L,
  perc_qcov_hsp = 80L,
  search_qcov_hsp = 50L,
  check_n = 10L,
  check_tolerance = 0.5,
  rank_by = c("similarity", "identity", "bitscore"),
  uninformative = uninformative_references(),
  output = c("long", "wide"),
  num_threads = NULL,
  batch_size = NULL,
  mt_mode = NULL,
  task = NULL,
  background = FALSE,
  job_dir = NULL,
  log_dir = NULL,
  retry_times = 3L,
  blast_type = "blastn",
  verbose = c("progress", "silent", "cmd", "output", "full"),
  env_name = "blastr-blast-env",
  out_file = deprecated(), # nolint: object_name_linter
  out_RDS = deprecated() # nolint: object_name_linter
) {
  rlang::check_dots_empty()

  if (lifecycle::is_present(out_file)) {
    lifecycle::deprecate_soft("0.1.7", "parallel_blast(out_file)")
  }
  if (lifecycle::is_present(out_RDS)) {
    lifecycle::deprecate_soft("0.1.7", "parallel_blast(out_RDS)")
  }

  rlang::check_required(query_seqs)
  rlang::check_required(db_path)
  verbose <- rlang::arg_match(verbose)
  rank_by <- rlang::arg_match(rank_by)
  output <- rlang::arg_match(output)
  if (!is.null(uninformative)) uninformative <- uninformative_references(uninformative)
  if (!is.null(mt_mode)) {
    mt_mode <- rlang::arg_match(as.character(mt_mode), c("0", "1", "2"))
  }

  if (
    rlang::is_null(query_seqs) ||
      isTRUE(length(query_seqs) == 0L) ||
      all(is.na(query_seqs))
  ) {
    cli::cli_abort(
      c(`x` = "No query sequences were provided to {.fun parallel_blast}."),
      class = "blastr_no_query_seqs_error"
    )
  }
  if (
    rlang::is_null(db_path) ||
      isTRUE(length(db_path) == 0L) ||
      all(is.na(db_path))
  ) {
    cli::cli_abort(
      c(`x` = "{.pkg BLASTr}: No database path was provided to {.fun parallel_blast}."), # nolint: line_length_linter
      class = "blastr_no_db_path_error"
    )
  }

  check_cmd(blast_type, env_name = env_name, verbose = "silent")

  seqs <- unique(stringr::str_replace_all(query_seqs[!is.na(query_seqs)], "\\s", ""))

  strategy <- blast_strategy(
    n_seqs = length(seqs),
    db_path = db_path,
    total_cores = total_cores,
    ram_gb = ram_gb,
    num_threads = num_threads,
    batch_size = batch_size,
    mt_mode = mt_mode
  )
  if (!identical(verbose, "silent")) {
    cli::cli_inform(c(`i` = "Running {.fun parallel_blast} for {length(seqs)} unique query sequences."))
    print(strategy)
    if (!isTRUE(strategy$fits_ram)) {
      cli::cli_warn(c(
        `!` = "The database ({round(strategy$db_gb)} GB) is larger than the RAM ({round(strategy$ram_gb)} GB).",
        `i` = "Each batch will read it from the disk; the search can be slow."
      ))
    }
  }

  blast_args <- c(
    "-db", db_arg(strategy$db_path),
    "-outfmt", blast_outfmt(),
    "-max_hsps", "1",
    "-perc_identity", as.character(perc_id),
    "-qcov_hsp_perc", as.character(min(perc_qcov_hsp, search_qcov_hsp)),
    "-num_alignments", as.character(max(num_alignments, search_alignments)),
    "-mt_mode", as.character(strategy$mt_mode),
    if (!is.null(task)) c("-task", task)
  )

  meta <- list(
    db_path = strategy$db_path,
    output = output,
    rank_by = rank_by,
    num_alignments = as.integer(num_alignments),
    search_alignments = as.integer(max(num_alignments, search_alignments)),
    task = if (is.null(task)) NA_character_ else task,
    perc_id = perc_id,
    perc_qcov_hsp = perc_qcov_hsp,
    search_qcov_hsp = min(perc_qcov_hsp, search_qcov_hsp),
    blast_type = blast_type,
    check_n = as.integer(check_n),
    check_tolerance = check_tolerance,
    uninformative = uninformative,
    uninformative_md5 = if (is.null(uninformative)) NA_character_ else attr(uninformative, "md5"),
    n_input = length(query_seqs),
    total_cores = as.integer(total_cores),
    background = isTRUE(background),
    log_dir = resolve_log_dir(log_dir)
  )

  blast_res <- run_blast_job(
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
  if (inherits(blast_res, "blastr_job")) {
    return(blast_res)
  }

  if (
    lifecycle::is_present(out_file) &&
      !rlang::is_na(out_file) &&
      !rlang::is_null(out_file)
  ) {
    readr::write_csv(x = blast_res, file = out_file, append = FALSE)
  }
  if (
    lifecycle::is_present(out_RDS) &&
      !rlang::is_na(out_RDS) &&
      !rlang::is_null(out_RDS)
  ) {
    readr::write_rds(x = blast_res, file = out_RDS)
  }

  blast_res
}

# Create the job and run it: in the background (returns the job) or in the
# R session with a progress bar (returns the collected results).
run_blast_job <- function(
  seqs,
  strategy,
  blast_type,
  blast_args,
  retry_times,
  env_name,
  job_dir,
  meta,
  background,
  verbose
) {
  if (is.null(job_dir)) {
    job_dir <- if (isTRUE(background)) {
      get_blastr_cache("jobs", format(Sys.time(), "%Y%m%d-%H%M%S"))
    } else {
      fs::file_temp(pattern = "blastr_job_")
    }
    if (!isTRUE(background)) {
      withr::defer(fs::dir_delete(job_dir))
    }
  }
  job <- new_blast_job(
    seqs = seqs,
    strategy = strategy,
    blast_type = blast_type,
    blast_args = blast_args,
    retry_times = retry_times,
    env_name = env_name,
    job_dir = job_dir,
    meta = meta
  )

  if (verbose %in% c("cmd", "output", "full") || (verbose == "progress" && rlang::is_interactive())) {
    cmd <- readRDS(fs::path(job$dir, "job.rds"))$meta$blast_cmd
    cli::cli_inform(c(`i` = "BLAST command: {.code {cmd}}"))
  }

  if (isTRUE(background)) {
    launch_scheduler_bg(job)
    if (!identical(verbose, "silent")) {
      cli::cli_inform(c(
        `v` = "BLAST started in the background: {.path {job$dir}}",
        `i` = "Follow it with {.code blast_job_status(job)}; get the results with {.code blast_job_result(job)}."
      ))
    }
    return(job)
  }

  bar <- new_progress(
    job_progress(read_batches(job)),
    show = !identical(verbose, "silent")
  )
  blastr_scheduler(job$dir, on_progress = bar$update)
  bar$done()
  collect_blast_job(job)
}

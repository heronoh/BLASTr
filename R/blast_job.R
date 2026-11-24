# BLAST jobs: the query sequences split into batches on disk, run by a
# scheduler that keeps at most `processes` blastn running. The scheduler runs
# in the R session (foreground) or in a separate Rscript (background), so a
# long search does not block R and survives the end of the session.
#
# Job directory:
#   job.rds        settings (binary, arguments, processes, threads, ...)
#   seqs.tsv       query id -> sequence
#   batches.tsv    one line per batch: status, exit code, attempts, times
#   queries/<b>.fasta, out/<b>.tsv, out/<b>.err
#   scheduler.R, scheduler.log, scheduler.pid   (background only)
#   cancel         created by blast_job_cancel()

new_blast_job <- function(
  seqs,
  strategy,
  blast_type,
  blast_args,
  retry_times,
  env_name,
  job_dir,
  meta
) {
  job_dir <- fs::path_abs(job_dir)
  if (fs::dir_exists(job_dir) && length(fs::dir_ls(job_dir)) > 0L) {
    cli::cli_abort(
      c(
        `x` = "The job directory {.path {job_dir}} already exists and is not empty.",
        `i` = "Use {.fun blast_job} to open an existing job, or choose another {.arg job_dir}."
      ),
      class = "blastr_job_dir_exists_error"
    )
  }
  fs::dir_create(fs::path(job_dir, c("queries", "out")))

  ids <- sprintf("s%06d", seq_along(seqs))
  utils::write.table(
    data.frame(query = ids, sequence = seqs),
    fs::path(job_dir, "seqs.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )

  # Invalid sequences go alone, so BLAST reports their error and they do
  # not affect the others.
  valid <- is_valid_seq(seqs)
  groups <- split(ids[valid], ceiling(seq_len(sum(valid)) / strategy$batch_size))
  groups <- c(unname(groups), as.list(ids[!valid]))
  batches <- sprintf("b%05d", seq_along(groups))
  for (i in seq_along(groups)) {
    k <- match(groups[[i]], ids)
    writeLines(
      paste0(">", ids[k], "\n", seqs[k]),
      fs::path(job_dir, "queries", batches[i], ext = "fasta")
    )
  }
  utils::write.table(
    data.frame(
      batch = batches,
      n = lengths(groups),
      status = "pending",
      exit_code = NA_integer_,
      attempt = 0L,
      start = NA_real_,
      end = NA_real_
    ),
    fs::path(job_dir, "batches.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )

  bin <- resolve_blast_bin(blast_type, env_name)
  blast_bin <- bin$path
  meta$blast_cmd <- blast_cmd_line(blast_type, blast_args, strategy$threads)
  meta$blast_version <- blast_version(blast_bin)
  # "conda" (the environment of BLASTr) or "system"
  meta$blast_source <- bin$source
  meta$db_dates <- db_dates(strategy$db_path)
  saveRDS(
    list(
      blast_bin = blast_bin,
      blast_args = blast_args,
      processes = strategy$processes,
      threads = strategy$threads,
      retry_times = as.integer(retry_times),
      env = c(
        BLAST_USAGE_REPORT = "false",
        NCBI_DONT_USE_NCBIRC = "true",
        NCBI_DONT_USE_LOCAL_CONFIG = "true"
      ),
      strategy = strategy,
      meta = meta,
      created = Sys.time()
    ),
    fs::path(job_dir, "job.rds")
  )
  structure(list(dir = job_dir), class = "blastr_job")
}

is_valid_seq <- function(seqs) {
  nzchar(seqs) & grepl("^[A-Za-z*-]+$", seqs)
}

# The scheduler. Self-contained (base R + processx): in background mode its
# source is written to the job directory and run by Rscript, which does not
# need BLASTr to be installed.
blastr_scheduler <- function(job_dir, on_progress = NULL) {
  job <- readRDS(file.path(job_dir, "job.rds"))
  bpath <- file.path(job_dir, "batches.tsv")
  read_b <- function() {
    utils::read.delim(
      bpath,
      colClasses = c(
        batch = "character", n = "integer", status = "character",
        exit_code = "integer", attempt = "integer", start = "numeric", end = "numeric"
      )
    )
  }
  write_b <- function(b) {
    tmp <- paste0(bpath, ".tmp")
    utils::write.table(b, tmp, sep = "\t", quote = FALSE, row.names = FALSE)
    file.rename(tmp, bpath)
  }
  qfile <- function(id) file.path(job_dir, "queries", paste0(id, ".fasta"))
  ofile <- function(id) file.path(job_dir, "out", paste0(id, ".tsv"))
  efile <- function(id) file.path(job_dir, "out", paste0(id, ".err"))

  b <- read_b()
  b$status[b$status == "running"] <- "pending" # resumed after a crash
  write_b(b)
  writeLines(as.character(Sys.getpid()), file.path(job_dir, "scheduler.pid"))
  running <- list()
  # CPU time of the finished child processes (blastn), for the run log
  pt0 <- proc.time()
  t0 <- as.numeric(Sys.time())
  on.exit(
    {
      for (p in running) try(p$kill(), silent = TRUE)
      unlink(file.path(job_dir, "scheduler.pid"))
      pt <- proc.time() - pt0
      ufile <- file.path(job_dir, "usage.tsv")
      utils::write.table(
        data.frame(
          start = t0, end = as.numeric(Sys.time()),
          child_user_s = unname(pt[4L]), child_sys_s = unname(pt[5L])
        ),
        ufile,
        sep = "\t", quote = FALSE, row.names = FALSE,
        col.names = !file.exists(ufile), append = file.exists(ufile)
      )
    },
    add = TRUE
  )

  repeat {
    if (file.exists(file.path(job_dir, "cancel"))) {
      for (p in running) try(p$kill(), silent = TRUE)
      running <- list()
      b$status[b$status %in% c("running", "pending")] <- "cancelled"
      write_b(b)
      break
    }
    while (length(running) < job$processes && any(b$status == "pending")) {
      i <- which(b$status == "pending")[1L]
      id <- b$batch[i]
      running[[id]] <- processx::process$new(
        job$blast_bin,
        c(
          job$blast_args,
          "-query", qfile(id),
          "-out", ofile(id),
          # no -num_threads with -subject (BLAST ignores it and warns)
          if (!is.na(job$threads)) c("-num_threads", as.character(job$threads))
        ),
        env = c("current", job$env),
        stdout = NULL,
        stderr = efile(id),
        cleanup = TRUE
      )
      b$status[i] <- "running"
      b$attempt[i] <- b$attempt[i] + 1L
      b$start[i] <- as.numeric(Sys.time())
    }
    write_b(b)
    if (!is.null(on_progress)) on_progress(b)
    if (length(running) == 0L) break

    processx::poll(unname(running), 2000L)
    for (id in names(running)) {
      p <- running[[id]]
      if (p$is_alive()) next
      i <- match(id, b$batch)
      st <- p$get_exit_status()
      b$exit_code[i] <- st
      b$end[i] <- as.numeric(Sys.time())
      running[[id]] <- NULL
      if (identical(as.integer(st), 0L)) {
        b$status[i] <- "done"
      } else if (b$n[i] > 1L) {
        # find the sequence that failed: one batch per sequence
        l <- readLines(qfile(id))
        h <- grep("^>", l)
        new <- sprintf("%s_%03d", id, seq_along(h))
        for (j in seq_along(h)) {
          end <- if (j < length(h)) h[j + 1L] - 1L else length(l)
          writeLines(l[h[j]:end], qfile(new[j]))
        }
        b$status[i] <- "split"
        b <- rbind(b, data.frame(
          batch = new, n = 1L, status = "pending", exit_code = NA_integer_,
          attempt = 0L, start = NA_real_, end = NA_real_
        ))
      } else if (b$attempt[i] <= job$retry_times) {
        b$status[i] <- "pending"
      } else {
        b$status[i] <- "failed"
      }
    }
  }
  invisible(b)
}

launch_scheduler_bg <- function(job) {
  script <- fs::path(job$dir, "scheduler.R")
  writeLines(
    c(
      paste0(
        "blastr_scheduler <- ",
        paste(deparse(blastr_scheduler, width.cutoff = 500L), collapse = "\n")
      ),
      paste0("blastr_scheduler(", deparse(as.character(job$dir)), ")")
    ),
    script
  )
  rscript <- file.path(
    R.home("bin"),
    if (identical(.Platform$OS.type, "windows")) "Rscript.exe" else "Rscript"
  )
  processx::process$new(
    rscript,
    c("--vanilla", script),
    stdout = as.character(fs::path(job$dir, "scheduler.log")),
    stderr = "2>&1",
    cleanup = FALSE,
    supervise = FALSE
  )
  invisible(job)
}

#' @title Background BLAST jobs
#'
#' @description
#' [parallel_blast()] with `background = TRUE` returns a job and the search
#' goes on in a separate process, so R stays free (and the search goes on
#' even if the R session ends). Use these functions to follow the job and
#' get its results:
#'
#' * `blast_job()` opens a job from its directory, e.g. in a new R session.
#' * `blast_job_status()` shows the progress and the estimated time left.
#' * `blast_job_result()` returns the results like [parallel_blast()],
#'   waiting for the job to finish if `wait = TRUE`.
#' * `blast_job_cancel()` stops the job.
#' * `blast_job_resume()` restarts a job whose scheduler was stopped (e.g.
#'   the machine rebooted); finished batches are kept.
#'
#' @param job A job returned by [parallel_blast()] or `blast_job()`, or the
#'   path to its directory.
#' @param job_dir Directory of an existing job.
#' @param wait Wait for the job to finish? If `FALSE` and the job is not
#'   finished, `blast_job_result()` returns `NULL`.
#' @param progress Show a progress bar while waiting.
#'
#' @returns `blast_job()` and `blast_job_resume()` return the job;
#'   `blast_job_status()` a list of class `blastr_job_status`;
#'   `blast_job_result()` a tibble (see [parallel_blast()]);
#'   `blast_job_cancel()` the job, invisibly.
#'
#' @examples
#' # needs BLAST+ (see install_dependencies())
#' if (has_blast()) {
#'   fasta <- system.file("extdata", "minimal_db_blast.fasta", package = "BLASTr")
#'   db <- file.path(tempdir(), "minimal_db_blast")
#'   make_blast_db(fasta, db_path = db)
#'   seqs <- readLines(system.file("extdata", "seqs_test.txt", package = "BLASTr"))
#'   job_dir <- tempfile("blastr_job_")
#'   job <- parallel_blast(seqs, db, background = TRUE, job_dir = job_dir)
#'   blast_job_status(job)
#'   # later, or in another R session:
#'   job <- blast_job(job_dir)
#'   res <- blast_job_result(job)
#'   res
#' }
#'
#' @family search functions
#' @name blast_job
NULL

#' @rdname blast_job
#' @export
blast_job <- function(job_dir) {
  job_dir <- fs::path_abs(job_dir)
  if (!fs::file_exists(fs::path(job_dir, "job.rds"))) {
    cli::cli_abort(
      c(`x` = "{.path {job_dir}} is not a BLASTr job directory."),
      class = "blastr_not_a_job_error"
    )
  }
  structure(list(dir = job_dir), class = "blastr_job")
}

as_blast_job <- function(job) {
  if (inherits(job, "blastr_job")) job else blast_job(job)
}

read_batches <- function(job) {
  utils::read.delim(
    fs::path(job$dir, "batches.tsv"),
    colClasses = c(
      batch = "character", n = "integer", status = "character",
      exit_code = "integer", attempt = "integer", start = "numeric", end = "numeric"
    )
  )
}

scheduler_alive <- function(job) {
  pid_file <- fs::path(job$dir, "scheduler.pid")
  if (!fs::file_exists(pid_file)) {
    return(FALSE)
  }
  pid <- suppressWarnings(as.integer(readLines(pid_file, warn = FALSE)[1L]))
  isTRUE(!is.na(pid) && ps::ps_is_running(ps::ps_handle(pid)))
}

job_progress <- function(b, start_time = NULL) {
  final <- b$status %in% c("done", "failed", "cancelled")
  live <- b$status != "split"
  total <- sum(b$n[live])
  done <- sum(b$n[final])
  if (is.null(start_time)) {
    start_time <- suppressWarnings(min(b$start, na.rm = TRUE))
  }
  now <- as.numeric(Sys.time())
  finished <- all(final[live])
  end_time <- if (finished) suppressWarnings(max(b$end, na.rm = TRUE)) else now
  elapsed <- if (is.finite(start_time)) end_time - start_time else 0
  eta <- if (done > 0 && !finished) elapsed / done * (total - done) else NA_real_
  list(
    total = total,
    done = done,
    failed = sum(b$n[b$status == "failed"]),
    running = sum(b$n[b$status == "running"]),
    batches_done = sum(final[live]),
    batches = sum(live),
    finished = finished,
    cancelled = any(b$status == "cancelled"),
    elapsed = elapsed,
    eta = eta
  )
}

#' @rdname blast_job
#' @export
blast_job_status <- function(job) {
  job <- as_blast_job(job)
  pr <- job_progress(read_batches(job))
  pr$dir <- job$dir
  pr$scheduler_running <- scheduler_alive(job)
  pr$age <- as.numeric(Sys.time()) - as.numeric(file.mtime(fs::path(job$dir, "job.rds")))
  structure(pr, class = "blastr_job_status")
}

fmt_secs <- function(s) {
  if (is.na(s)) {
    return("unknown (after the first batch)")
  }
  s <- round(s)
  if (s < 60) {
    return(paste0(s, "s"))
  }
  if (s < 3600) {
    return(sprintf("%dm %02ds", s %/% 60, s %% 60))
  }
  sprintf("%dh %02dm", s %/% 3600, (s %% 3600) %/% 60)
}

#' @export
print.blastr_job_status <- function(x, ...) {
  state <- if (x$cancelled) {
    "cancelled"
  } else if (x$finished) {
    "finished"
  } else if (x$scheduler_running) {
    "running"
  } else if (x$done == 0 && isTRUE(x$age < 30)) {
    "starting"
  } else {
    "stopped (use blast_job_resume())"
  }
  cli::cli_text("BLASTr job {.path {x$dir}}: {state}")
  cli::cli_text(
    "{x$done}/{x$total} sequences ({round(100 * x$done / max(1, x$total))}%), ",
    "{x$batches_done}/{x$batches} batches; elapsed {fmt_secs(x$elapsed)}",
    "{if (!x$finished) paste0('; ETA ', fmt_secs(x$eta)) else ''}",
    "{if (x$failed > 0) paste0('; ', x$failed, ' failed') else ''}."
  )
  invisible(x)
}

#' @export
print.blastr_job <- function(x, ...) {
  print(blast_job_status(x))
  invisible(x)
}

#' @rdname blast_job
#' @export
blast_job_cancel <- function(job) {
  job <- as_blast_job(job)
  file.create(fs::path(job$dir, "cancel"))
  invisible(job)
}

#' @rdname blast_job
#' @export
blast_job_resume <- function(job) {
  job <- as_blast_job(job)
  if (scheduler_alive(job)) {
    cli::cli_inform(c(`i` = "The job is already running."))
    return(invisible(job))
  }
  if (fs::file_exists(fs::path(job$dir, "cancel"))) {
    fs::file_delete(fs::path(job$dir, "cancel"))
    b <- read_batches(job)
    b$status[b$status == "cancelled"] <- "pending"
    utils::write.table(b, fs::path(job$dir, "batches.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  }
  launch_scheduler_bg(job)
}

#' @rdname blast_job
#' @export
blast_job_result <- function(job, wait = TRUE, progress = rlang::is_interactive()) {
  job <- as_blast_job(job)
  pr <- job_progress(read_batches(job))
  if (!pr$finished) {
    if (!isTRUE(wait)) {
      return(NULL)
    }
    wait_for_job(job, progress = progress)
  }
  collect_blast_job(job)
}

# Show progress until the background scheduler finishes.
wait_for_job <- function(job, progress = TRUE) {
  bar <- new_progress(job_progress(read_batches(job)), progress)
  repeat {
    b <- read_batches(job)
    bar$update(b)
    pr <- job_progress(b)
    if (pr$finished) break
    if (!scheduler_alive(job)) {
      Sys.sleep(3) # the scheduler may be starting
      if (!scheduler_alive(job) && !job_progress(read_batches(job))$finished) {
        cli::cli_abort(
          c(
            `x` = "The BLAST job stopped before finishing.",
            `i` = "Restart it with {.code blast_job_resume(\"{job$dir}\")}."
          ),
          class = "blastr_job_stopped_error"
        )
      }
    }
    Sys.sleep(2)
  }
  bar$done()
}

# Progress display: a cli bar with ETA when interactive; otherwise one line
# each time a batch finishes (for logs).
new_progress <- function(pr0, show) {
  if (!isTRUE(show)) {
    return(list(update = function(b) NULL, done = function() NULL))
  }
  env <- new.env()
  env$last <- -1L
  env$eta <- fmt_secs(NA)
  interactive_bar <- rlang::is_interactive()
  if (interactive_bar) {
    env$id <- cli::cli_progress_bar(
      total = pr0$total,
      format = paste(
        "BLAST {cli::pb_bar} {cli::pb_percent} | {cli::pb_current}/{cli::pb_total} seqs",
        "| {cli::pb_elapsed} | ETA {eta}"
      ),
      clear = FALSE,
      .envir = env
    )
  }
  list(
    update = function(b) {
      pr <- job_progress(b)
      env$eta <- fmt_secs(pr$eta)
      if (interactive_bar) {
        cli::cli_progress_update(id = env$id, set = pr$done, .envir = env)
      } else if (pr$batches_done != env$last) {
        env$last <- pr$batches_done
        cli::cli_inform(c(`i` = paste0(
          "BLAST: {pr$done}/{pr$total} sequences ({round(100 * pr$done / max(1, pr$total))}%), ",
          "{pr$batches_done}/{pr$batches} batches; elapsed {fmt_secs(pr$elapsed)}",
          "{if (!pr$finished) paste0('; ETA ', fmt_secs(pr$eta)) else ''}"
        )))
      }
    },
    done = function() {
      if (interactive_bar) cli::cli_progress_done(id = env$id, .envir = env)
    }
  )
}

# The BLAST command of each batch, without the sequences: the query and
# output files are shown as <batch.fasta> and <batch.tsv>.
blast_cmd_line <- function(blast_type, blast_args, threads = NA) {
  args <- c(
    blast_args, "-query", "<batch.fasta>", "-out", "<batch.tsv>",
    if (!is.na(threads)) c("-num_threads", as.character(threads))
  )
  q <- ifelse(grepl("[^A-Za-z0-9_./<>=:,+-]", args), shQuote(args), args)
  paste(c(blast_type, q), collapse = " ")
}

# Date of each BLAST database: "last-updated" of its .njs file (BLAST
# databases v5), or else the date of its files.
db_dates <- function(db_paths) {
  vapply(db_paths, function(p) {
    njs <- paste0(p, c(".njs", ".nal", ".nin", ".pjs", ".pal", ".pin"))
    njs <- njs[fs::file_exists(njs)]
    j <- njs[grepl("js$", njs)]
    if (length(j) > 0L) {
      l <- grep("\"last-updated\"", readLines(j[1L], n = 50L, warn = FALSE), value = TRUE)
      if (length(l) > 0L) {
        return(sub('.*"last-updated"\\s*:\\s*"([^"]+)".*', "\\1", l[1L]))
      }
    }
    f <- if (length(njs) > 0L) njs[1L] else if (fs::file_exists(p)) p else NA_character_
    if (is.na(f)) NA_character_ else format(fs::file_info(f)$modification_time, "%Y-%m-%dT%H:%M:%S")
  }, character(1L), USE.NAMES = FALSE)
}

# Results of a finished job as returned by parallel_blast().
collect_blast_job <- function(job) {
  st <- readRDS(fs::path(job$dir, "job.rds"))
  seqs <- utils::read.delim(
    fs::path(job$dir, "seqs.tsv"),
    colClasses = "character"
  )
  b <- read_batches(job)
  final <- b[b$status != "split", , drop = FALSE]
  out_files <- fs::path(job$dir, "out", final$batch[final$status == "done"], ext = "tsv")

  # exit code and BLAST messages of each sequence
  q_batch <- lapply(final$batch, function(id) {
    l <- readLines(fs::path(job$dir, "queries", id, ext = "fasta"))
    data.frame(query = sub("^>", "", l[startsWith(l, ">")]), batch = id)
  })
  q_batch <- do.call(rbind, q_batch)
  # a header line inside an invalid sequence is not a query id
  q_batch <- q_batch[q_batch$query %in% seqs$query, , drop = FALSE]
  q_batch <- q_batch[!duplicated(q_batch$query), , drop = FALSE]
  k <- match(seqs$query, q_batch$query)
  bb <- final[match(q_batch$batch[k], final$batch), , drop = FALSE]
  err <- vapply(seq_len(nrow(seqs)), function(i) {
    f <- fs::path(job$dir, "out", q_batch$batch[k[i]], ext = "err")
    if (is.na(k[i]) || !fs::file_exists(f)) {
      return("")
    }
    l <- readLines(f, warn = FALSE)
    if (bb$n[i] > 1L) {
      l <- l[grepl(paste0("\\b", seqs$query[i], "\\b"), l)]
    }
    paste(l, collapse = "\n")
  }, character(1L))

  summ <- NULL
  res <- if (identical(st$meta$output, "hsps")) {
    read_blast_hsps(out_files, seqs, st$meta$fields)
  } else {
    m <- st$meta
    rank_by <- if (is.null(m$rank_by)) "similarity" else m$rank_by
    hits <- read_blast_hits(out_files)
    as_long <- function(n) {
      blast_long(
        hits,
        seq_ids = seqs$query,
        seqs = seqs$sequence,
        num_alignments = n,
        rank_by = rank_by,
        tolerance = m$check_tolerance,
        perc_qcov_hsp = if (is.null(m$perc_qcov_hsp)) NA else m$perc_qcov_hsp,
        uninformative = m$uninformative
      )
    }
    # the summary counts all the subjects BLAST returned, not only the kept rows
    summ <- blast_summary(
      as_long(Inf),
      check_n = m$check_n, tolerance = m$check_tolerance, rank_by = rank_by,
      search_alignments = m$search_alignments
    )
    long <- as_long(m$num_alignments)
    if (identical(m$output, "wide")) widen_blast_hits(long, summ, m$num_alignments) else long
  }
  pr <- job_progress(b)
  log_path <- tryCatch(
    write_blast_log(job, res, st, b, seqs, summ),
    error = function(e) {
      cli::cli_warn(c(`!` = "Could not write the BLASTr log: {conditionMessage(e)}"))
      NULL
    }
  )
  attr(res, "BLASTr_metadata") <- c(
    st$meta,
    list(
      strategy = st$strategy,
      elapsed_s = pr$elapsed,
      hit_counts = if (is.null(summ)) NULL else summ[, c("Sequence", hit_count_cols)],
      job_dir = as.character(job$dir),
      cpu_s = job_cpu_s(job$dir),
      log_file = if (is.null(log_path)) NULL else as.character(log_path),
      exit_codes = tibble::tibble(
        query_seq = seqs$sequence,
        exit_code = ifelse(bb$status == "cancelled", NA_integer_, bb$exit_code),
        stderr = err
      )
    )
  )
  res
}

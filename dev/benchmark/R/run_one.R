# Execute ONE benchmark configuration. Called by bench_lib.R::run_config()
# under /usr/bin/time, so the CPU of every child process is accounted.
#
# Usage: Rscript run_one.R <config.json>
#
# Strategies:
#   A_blastr  BLASTr 0.1.7 parallel_blast(), 1 query per blastn, P processes x T threads
#             (loaded from cfg$pkg_path_old, a worktree at 60ce902)
#   E_blastr  BLASTr parallel_blast() of this repo, automatic strategy with
#             total_cores = P * T (batched; the version being benchmarked)
#   D_raw     raw blastn, 1 query per process, P processes x T threads (no R/mirai/condathis)
#   B_single  raw blastn, all queries in one process, T threads, -mt_mode m
#   C_batch   raw blastn, batches of k queries, P processes x T threads
#
# Writes <out_dir>/hits.tsv (seq_id, rank, subject, staxid, identity, qcovhsp, evalue, bitscore)
# and <out_dir>/jobs.tsv (one line per blastn job: start, end, status).

args <- commandArgs(trailingOnly = TRUE)
cfg <- jsonlite::read_json(args[[1]], simplifyVector = TRUE)
dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)

blast_bin <- path.expand("~/.local/share/R/condathis/envs/blastr-blast-env/bin/blastn")
outfmt <- "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp staxid stitle"

read_fasta <- function(path) {
  l <- readLines(path)
  h <- grep("^>", l)
  ids <- sub("^>", "", l[h])
  seqs <- vapply(seq_along(h), function(i) {
    end <- if (i < length(h)) h[i + 1] - 1 else length(l)
    paste(l[(h[i] + 1):end], collapse = "")
  }, character(1))
  stats::setNames(seqs, ids)
}

queries <- read_fasta(cfg$query_fasta)
if (!is.null(cfg$n_queries) && cfg$n_queries < length(queries)) {
  queries <- queries[seq_len(cfg$n_queries)]
}

blast_args <- function(query_path, out_path, threads) {
  a <- c(
    "-db", cfg$db, "-query", query_path, "-out", out_path,
    "-outfmt", outfmt,
    "-max_hsps", "1",
    "-perc_identity", as.character(cfg$perc_id),
    "-qcov_hsp_perc", as.character(cfg$qcov),
    "-num_alignments", as.character(cfg$num_alignments),
    "-num_threads", as.character(threads)
  )
  if (!is.null(cfg$mt_mode) && !is.na(cfg$mt_mode)) a <- c(a, "-mt_mode", as.character(cfg$mt_mode))
  a
}

blast_env <- c(
  "current",
  BLAST_USAGE_REPORT = "false",
  NCBI_DONT_USE_NCBIRC = "true",
  NCBI_DONT_USE_LOCAL_CONFIG = "true"
)

# Run a list of query chunks with at most P concurrent blastn processes.
run_chunks <- function(chunks, P, threads) {
  tmp <- file.path(cfg$out_dir, "chunks")
  dir.create(tmp, showWarnings = FALSE)
  jobs <- vector("list", length(chunks))
  running <- list()
  next_i <- 1L
  log <- data.frame()
  while (next_i <= length(chunks) || length(running) > 0) {
    while (length(running) < P && next_i <= length(chunks)) {
      q <- file.path(tmp, sprintf("c%05d.fasta", next_i))
      o <- file.path(tmp, sprintf("c%05d.tsv", next_i))
      ch <- chunks[[next_i]]
      writeLines(paste0(">", names(ch), "\n", ch), q)
      p <- processx::process$new(blast_bin, blast_args(q, o, threads),
        env = blast_env, stdout = NULL, stderr = "|", cleanup = TRUE
      )
      running[[as.character(next_i)]] <- list(p = p, start = Sys.time(), n = length(ch), out = o)
      next_i <- next_i + 1L
    }
    procs <- lapply(running, `[[`, "p")
    processx::poll(procs, 1000)
    for (key in names(running)) {
      r <- running[[key]]
      if (!r$p$is_alive()) {
        st <- r$p$get_exit_status()
        err <- tryCatch(r$p$read_all_error(), error = function(e) "")
        log <- rbind(log, data.frame(
          chunk = as.integer(key), n_queries = r$n, start = as.numeric(r$start),
          end = as.numeric(Sys.time()), status = st, stderr = substr(gsub("\\s+", " ", err), 1, 300)
        ))
        running[[key]] <- NULL
      }
    }
  }
  outs <- list.files(tmp, pattern = "\\.tsv$", full.names = TRUE)
  hits <- do.call(rbind, lapply(outs, function(f) {
    if (file.size(f) == 0) return(NULL)
    utils::read.delim(f, header = FALSE, quote = "", comment.char = "", colClasses = "character")
  }))
  unlink(tmp, recursive = TRUE)
  list(hits = hits, log = log)
}

split_k <- function(x, k) split(x, ceiling(seq_along(x) / k))

t0 <- Sys.time()
# wide BLASTr results -> long hits (identity column: k_identity, or
# k_indentity before the rename)
blastr_long <- function(res, n) {
  long <- do.call(rbind, lapply(seq_len(n), function(r) {
    sc <- paste0(r, "_subject")
    if (!sc %in% names(res)) return(NULL)
    idc <- intersect(paste0(r, c("_identity", "_indentity")), names(res))[1]
    data.frame(
      seq_id = ids[res$Sequence], rank = r,
      subject = res[[sc]], staxid = res[[paste0(r, "_staxid")]],
      identity = res[[idc]], qcovhsp = res[[paste0(r, "_qcovhsp")]],
      evalue = res[[paste0(r, "_e-value")]], bitscore = res[[paste0(r, "_bitscore")]]
    )
  }))
  long[!is.na(long$subject), ]
}
ids <- stats::setNames(names(queries), unname(queries))

if (cfg$strategy == "E_blastr") {
  suppressMessages(pkgload::load_all(cfg$pkg_path, quiet = TRUE))
  res <- parallel_blast(
    query_seqs = unname(queries),
    db_path = strsplit(cfg$db, " ", fixed = TRUE)[[1]],
    total_cores = cfg$P * cfg$T,
    perc_id = cfg$perc_id,
    perc_qcov_hsp = cfg$qcov,
    num_alignments = cfg$num_alignments,
    verbose = "silent"
  )
  meta <- attr(res, "BLASTr_metadata")
  st <- meta$strategy
  writeLines(
    sprintf("processes=%d threads=%d mt_mode=%d batch_size=%d n_batches=%d search_alignments=%d",
      st$processes, st$threads, st$mt_mode, st$batch_size, st$n_batches, meta$search_alignments),
    file.path(cfg$out_dir, "strategy.txt")
  )
  long <- blastr_long(res, cfg$num_alignments)
  ex <- meta$exit_codes
  log <- data.frame(chunk = NA, n_queries = length(queries), start = as.numeric(t0),
    end = as.numeric(Sys.time()), status = max(ex$exit_code, na.rm = TRUE), stderr = "")
} else if (cfg$strategy == "A_blastr") {
  pkg_old <- if (!is.null(cfg$pkg_path_old) && dir.exists(cfg$pkg_path_old)) cfg$pkg_path_old else
    stop("strategy A needs BLASTr 0.1.7 at cfg$pkg_path_old (git worktree at 60ce902)")
  suppressMessages(pkgload::load_all(pkg_old, quiet = TRUE))
  withr::local_envvar(c(BLAST_USAGE_REPORT = "false"))
  res <- parallel_blast(
    query_seqs = unname(queries),
    db_path = cfg$db,
    total_cores = cfg$P,
    num_threads = cfg$T,
    perc_id = cfg$perc_id,
    perc_qcov_hsp = cfg$qcov,
    num_alignments = cfg$num_alignments,
    mt_mode = as.character(if (is.null(cfg$mt_mode) || is.na(cfg$mt_mode)) 2 else cfg$mt_mode),
    verbose = "silent"
  )
  long <- blastr_long(res, cfg$num_alignments)
  ex <- attr(res, "BLASTr_metadata")$exit_codes
  log <- data.frame(chunk = NA, n_queries = length(queries), start = as.numeric(t0),
    end = as.numeric(Sys.time()), status = max(ex$exit_code), stderr = "")
} else {
  chunks <- switch(cfg$strategy,
    D_raw = split_k(queries, 1L),
    B_single = list(queries),
    C_batch = split_k(queries, cfg$k)
  )
  P <- if (cfg$strategy == "B_single") 1L else cfg$P
  r <- run_chunks(chunks, P, cfg$T)
  log <- r$log
  h <- r$hits
  if (is.null(h)) {
    long <- data.frame(seq_id = character(), rank = integer(), subject = character(),
      staxid = character(), identity = numeric(), qcovhsp = numeric(), evalue = numeric(), bitscore = numeric())
  } else {
    names(h)[1:15] <- c("seq_id", "subject", "identity", "length", "mismatch", "gaps", "qstart", "qend",
      "sstart", "send", "evalue", "bitscore", "qcovhsp", "staxid", "stitle")
    h <- h[order(h$seq_id), ]
    h$rank <- stats::ave(seq_along(h$seq_id), h$seq_id, FUN = seq_along)
    long <- h[, c("seq_id", "rank", "subject", "staxid", "identity", "qcovhsp", "evalue", "bitscore")]
  }
}
wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

utils::write.table(long, file.path(cfg$out_dir, "hits.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
utils::write.table(log, file.path(cfg$out_dir, "jobs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(format(wall, digits = 10), file.path(cfg$out_dir, "inner_wall.txt"))

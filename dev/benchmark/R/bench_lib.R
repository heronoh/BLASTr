# Benchmark driver helpers. Source from the BLASTr repo root:
#   source("dev/benchmark/R/bench_lib.R")

bench_dir <- normalizePath("dev/benchmark", mustWork = TRUE)
results_dir <- file.path(bench_dir, "results")
dir.create(file.path(results_dir, "runs"), recursive = TRUE, showWarnings = FALSE)

DBS <- list(
  mini = file.path(bench_dir, "db", "mini_multimarker"),
  core_nt = "/data/databases/core_nt_nov26/core_nt",
  core_nt_lgc = paste(
    "/data/databases/core_nt_nov26/core_nt",
    "/data/databases/LGC12Sdb/jun26/LGC12Sdb_complete_noGaps-2026-06-02.fasta"
  )
)

default_cfg <- list(
  phase = "p0", db_name = "mini", query_set = "S100", n_queries = NULL,
  strategy = "C_batch", P = 1L, T = 1L, k = NA_integer_, mt_mode = NA_integer_,
  num_alignments = 3L, perc_id = 80L, qcov = 80L, cache = "warm", rep = 1L
)

parse_time_v <- function(path) {
  l <- readLines(path, warn = FALSE)
  get <- function(key) {
    x <- grep(key, l, fixed = TRUE, value = TRUE)
    if (!length(x)) return(NA_character_)
    trimws(sub(".*: ", "", x[[1]]))
  }
  el <- get("Elapsed (wall clock) time")
  parts <- as.numeric(strsplit(el, ":")[[1]])
  wall <- sum(parts * 60^(rev(seq_along(parts)) - 1))
  data.frame(
    wall_s = wall,
    user_s = as.numeric(get("User time (seconds)")),
    sys_s = as.numeric(get("System time (seconds)")),
    max_rss_single_mb = as.numeric(get("Maximum resident set size (kbytes)")) / 1024,
    fs_inputs_blocks = as.numeric(get("File system inputs")),
    major_faults = as.numeric(get("Major (requiring I/O) page faults")),
    exit_status = as.integer(get("Exit status"))
  )
}

# Sum over processes of the last utime+stime seen for each pid (seconds).
# Catches blastn launched by orphaned mirai daemons (strategy A); processes
# shorter than the sampling interval are missed, so it is a lower bound.
sampled_cpu <- function(pids_path) {
  if (!file.exists(pids_path)) return(NA_real_)
  p <- utils::read.delim(pids_path, header = FALSE, col.names = c("pid", "comm", "ticks"))
  sum(tapply(p$ticks, p$pid, max)) / 100
}

summarise_ts <- function(path) {
  ts <- utils::read.delim(path)
  if (nrow(ts) < 2) {
    return(data.frame(
      sys_cpu_busy_pct = NA, sda_read_gb = NA, sda_util_pct = NA, load1_mean = NA,
      run_pss_peak_gb = NA, run_read_gb = NA, heron_running_max = NA, mem_avail_min_gb = NA,
      sampled_cpu_s = NA
    ))
  }
  d <- function(x) utils::tail(x, 1) - x[[1]]
  dt <- d(ts$t)
  data.frame(
    sys_cpu_busy_pct = 100 * d(ts$cpu_busy_jiffies) / d(ts$cpu_total_jiffies),
    sda_read_gb = d(ts$sda_read_sectors) * 512 / 1e9,
    sda_util_pct = 100 * d(ts$sda_io_ms) / 1000 / dt,
    load1_mean = mean(ts$load1),
    run_pss_peak_gb = max(ts$run_pss_kb) / 1024^2,
    run_read_gb = max(ts$run_read_bytes) / 1e9,
    heron_running_max = max(ts$heron_blast_running),
    mem_avail_min_gb = min(ts$mem_avail_kb) / 1024^2,
    sampled_cpu_s = sampled_cpu(paste0(path, ".pids"))
  )
}

run_config <- function(..., sample_interval = 2) {
  cfg <- utils::modifyList(default_cfg, list(...))
  cfg$db <- DBS[[cfg$db_name]]
  cfg$query_fasta <- file.path(bench_dir, "data", paste0(cfg$query_set, ".fasta"))
  cfg$pkg_path <- normalizePath(file.path(bench_dir, "..", ".."))
  # BLASTr 0.1.7 (one blastn per sequence) for strategy A: a worktree at
  # 60ce902, the last commit before the batched parallel_blast()
  #   git worktree add --detach ../BLASTr-v017 60ce902
  cfg$pkg_path_old <- normalizePath(file.path(cfg$pkg_path, "..", "BLASTr-v017"), mustWork = FALSE)
  cfg$run_id <- paste0(
    format(Sys.time(), "%Y%m%d-%H%M%S"), "_", cfg$phase, "_", cfg$strategy,
    "_P", cfg$P, "T", cfg$T, "_k", cfg$k, "_m", cfg$mt_mode, "_n", cfg$num_alignments, "_r", cfg$rep
  )
  cfg$out_dir <- file.path(results_dir, "runs", cfg$run_id)
  dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)
  cfg_path <- file.path(cfg$out_dir, "config.json")
  jsonlite::write_json(cfg, cfg_path, auto_unbox = TRUE, null = "null", na = "null", pretty = TRUE)

  ts_path <- file.path(cfg$out_dir, "timeseries.tsv")
  sampler <- processx::process$new(
    file.path(bench_dir, "R", "sampler.sh"),
    c(cfg$run_id, ts_path, as.character(sample_interval)),
    cleanup = TRUE
  )
  Sys.sleep(sample_interval + 0.5)

  time_path <- file.path(cfg$out_dir, "time.txt")
  res <- processx::run(
    "/usr/bin/time",
    c("-v", "-o", time_path, file.path(R.home("bin"), "Rscript"),
      file.path(bench_dir, "R", "run_one.R"), cfg_path),
    env = c("current", BLASTR_BENCH_RUN = cfg$run_id),
    error_on_status = FALSE, stderr_to_stdout = TRUE
  )
  writeLines(res$stdout, file.path(cfg$out_dir, "run_one.log"))

  Sys.sleep(sample_interval + 0.5)
  file.create(paste0(ts_path, ".stop"))
  sampler$wait(10000)
  sampler$kill()

  hits <- utils::read.delim(file.path(cfg$out_dir, "hits.tsv"))
  jobs <- utils::read.delim(file.path(cfg$out_dir, "jobs.tsv"))
  n_q <- length(grep("^>", readLines(cfg$query_fasta)))
  if (!is.null(cfg$n_queries)) n_q <- min(n_q, cfg$n_queries)

  row <- cbind(
    data.frame(
      run_id = cfg$run_id, time = format(Sys.time()), host = Sys.info()[["nodename"]],
      phase = cfg$phase, db = cfg$db_name, query_set = cfg$query_set, n_queries = n_q,
      strategy = cfg$strategy, P = cfg$P, T = cfg$T, threads_total = cfg$P * cfg$T,
      k = cfg$k, mt_mode = cfg$mt_mode, num_alignments = cfg$num_alignments,
      cache = cfg$cache, rep = cfg$rep,
      n_jobs = nrow(jobs), failed_jobs = sum(jobs$status != 0),
      n_hits = nrow(hits), n_queries_with_hit = length(unique(hits$seq_id))
    ),
    parse_time_v(time_path),
    summarise_ts(ts_path)
  )
  # /usr/bin/time misses orphaned mirai daemons (strategy A): use the larger value
  row$cpu_s <- pmax(row$user_s + row$sys_s, row$sampled_cpu_s, na.rm = TRUE)
  row$cpu_eff_pct <- 100 * row$cpu_s / (row$wall_s * row$threads_total)
  row$queries_per_h <- 3600 * n_q / row$wall_s

  runs_path <- file.path(results_dir, "runs.tsv")
  utils::write.table(row, runs_path,
    sep = "\t", quote = FALSE, row.names = FALSE,
    append = file.exists(runs_path), col.names = !file.exists(runs_path)
  )
  if (res$status != 0 || any(jobs$status != 0)) {
    warning("Run ", cfg$run_id, " had errors; see ", cfg$out_dir)
  }
  row
}

# Wait until heron's BLAST is still stopped and load is acceptable.
server_ok <- function(max_load = 45) {
  st <- system("ps -u heron -o stat=,comm= | awk '$2==\"blastn\" && $1 !~ /^T/' | wc -l", intern = TRUE)
  load1 <- as.numeric(strsplit(readLines("/proc/loadavg"), " ")[[1]][1])
  ok <- as.integer(st) == 0 && load1 < max_load
  attr(ok, "detail") <- sprintf("heron running blastn=%s load1=%.1f", st, load1)
  ok
}

# Drop a database's files from the page cache (no root: posix_fadvise
# DONTNEED; pages mapped by a running process stay). Returns the % still cached.
evict_db <- function(db = DBS$core_nt) {
  files <- Sys.glob(paste0(db, "*"))
  py <- "import os,sys
for f in sys.argv[1:]:
    fd=os.open(f,os.O_RDONLY); os.posix_fadvise(fd,0,0,os.POSIX_FADV_DONTNEED); os.close(fd)"
  system2("python3", c("-c", shQuote(py), shQuote(files)))
  cached_pct(db)
}

cached_pct <- function(db = DBS$core_nt) {
  files <- Sys.glob(paste0(db, "*"))
  out <- system2("fincore", c("-b", "-n", "-o", "RES,SIZE", shQuote(files)), stdout = TRUE)
  v <- matrix(as.numeric(unlist(strsplit(trimws(out), "\\s+"))), ncol = 2, byrow = TRUE)
  100 * sum(v[, 1]) / sum(v[, 2])
}

# Other users' blastn on the same machine (cold runs need the DB to ourselves)
others_blast <- function() {
  as.integer(system("ps -eo user=,comm= | awk '$2==\"blastn\" && $1!=\"noreh\"' | wc -l", intern = TRUE))
}

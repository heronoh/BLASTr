# Phase cold: same strategies as phase 1, each run starting with core_nt
# evicted from the page cache (first read comes from the RAID).
# Requires the DB to ourselves: waits while any other user runs blastn.
# Run from repo root:
#   nohup nice -n 5 Rscript dev/benchmark/R/06_phase_cold.R 2 > dev/benchmark/results/phase_cold.log 2>&1 &
source("dev/benchmark/R/bench_lib.R")
set.seed(3)

args <- commandArgs(trailingOnly = TRUE)
reps <- if (length(args)) as.integer(args[[1]]) else 2L

g <- function(strategy, P, T, k = NA, mt_mode = NA) data.frame(strategy, P, T, k, mt_mode)
grid <- rbind(
  g("B_single", 1, 32, mt_mode = 2),
  g("B_single", 1, 8, mt_mode = 2),
  g("C_batch", 4, 8, k = 25),
  g("D_raw", 32, 1),
  g("A_blastr", 32, 1)
)
plan <- grid[rep(seq_len(nrow(grid)), each = reps), ]
plan$rep <- rep(seq_len(reps), times = nrow(grid))
# first run = B T32 (fast) to measure the cold read; the rest shuffled
plan <- rbind(plan[1, ], plan[-1, ][sample(nrow(plan) - 1), ])
cat("Cold runs:", nrow(plan), "\n")

cold_log <- file.path(results_dir, "cold_log.tsv")
for (i in seq_len(nrow(plan))) {
  x <- plan[i, ]
  repeat {
    load1 <- as.numeric(strsplit(readLines("/proc/loadavg"), " ")[[1]][1])
    if (others_blast() == 0 && load1 < 50) break
    message(format(Sys.time()), " waiting: others_blast=", others_blast(), " load1=", load1)
    Sys.sleep(120)
  }
  before <- evict_db()
  a <- list(
    phase = "pcold", db_name = "core_nt", query_set = "S100", strategy = x$strategy,
    P = as.integer(x$P), T = as.integer(x$T), k = as.integer(x$k), mt_mode = as.integer(x$mt_mode),
    num_alignments = 3L, cache = "cold", rep = as.integer(x$rep)
  )
  r <- tryCatch(do.call(run_config, c(a, sample_interval = 5)), error = function(e) {
    message("ERROR: ", conditionMessage(e)); NULL
  })
  if (!is.null(r)) {
    utils::write.table(data.frame(run_id = r$run_id, cached_pct_start = before, cached_pct_end = cached_pct()),
      cold_log, sep = "\t", quote = FALSE, row.names = FALSE,
      append = file.exists(cold_log), col.names = !file.exists(cold_log))
    message(sprintf("[%d/%d] %s %s cached_start=%.1f%% wall=%.0fs cpu=%.2fh sda=%.0fGB util=%.0f%% hits=%d",
      i, nrow(plan), format(Sys.time(), "%H:%M"), r$run_id, before, r$wall_s, r$cpu_s / 3600,
      r$sda_read_gb, r$sda_util_pct, r$n_hits))
  }
}
message("Phase cold done")

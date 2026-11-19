# Phase 1: strategy screening on core_nt with WARM page cache.
# Constraints (2026-11-19): shared server, heron's BLAST paused (must stay
# untouched), renato using ~30 cores -> at most 40 threads per run, no cache
# eviction. Pauses automatically if heron's blastn resumes or load is high.
# Run from repo root:
#   nohup Rscript dev/benchmark/R/04_phase1.R > dev/benchmark/results/phase1.log 2>&1 &
source("dev/benchmark/R/bench_lib.R")
set.seed(2)

args <- commandArgs(trailingOnly = TRUE)
reps <- if (length(args)) as.integer(args[[1]]) else 2L

g <- function(strategy, P, T, k = NA, mt_mode = NA, num_alignments = 3, query_set = "S100") {
  data.frame(strategy, P, T, k, mt_mode, num_alignments, query_set)
}
grid <- rbind(
  # current mode (BLASTr) and its raw equivalent
  g("A_blastr", 32, 1),
  g("D_raw", 32, 1),
  g("D_raw", 16, 2),
  # one blastn with everything
  g("B_single", 1, 32, mt_mode = 1),
  g("B_single", 1, 32, mt_mode = 2),
  g("B_single", 1, 16, mt_mode = 2),
  # batches: same 32 threads split in different ways
  g("C_batch", 32, 1, k = 4),
  g("C_batch", 8, 4, k = 13),
  g("C_batch", 4, 8, k = 25),
  g("C_batch", 8, 4, k = 5),
  g("C_batch", 16, 2, k = 7)
)

plan <- grid[rep(seq_len(nrow(grid)), each = reps), ]
plan$rep <- rep(seq_len(reps), times = nrow(grid))
plan <- plan[sample(nrow(plan)), ]
cat("Phase 1 runs:", nrow(plan), "\n")

for (i in seq_len(nrow(plan))) {
  x <- plan[i, ]
  while (!server_ok(max_load = 50)) {
    message(format(Sys.time()), " waiting: ", attr(server_ok(), "detail"))
    Sys.sleep(120)
  }
  a <- list(
    phase = "p1", db_name = "core_nt", query_set = x$query_set, strategy = x$strategy,
    P = as.integer(x$P), T = as.integer(x$T), k = as.integer(x$k), mt_mode = as.integer(x$mt_mode),
    num_alignments = as.integer(x$num_alignments), rep = as.integer(x$rep)
  )
  r <- tryCatch(do.call(run_config, c(a, sample_interval = 5)), error = function(e) {
    message("ERROR: ", conditionMessage(e)); NULL
  })
  if (!is.null(r)) {
    message(sprintf("[%d/%d] %s %s wall=%.0fs cpu=%.2fh hits=%d heron_running=%d load=%.0f",
      i, nrow(plan), format(Sys.time(), "%H:%M"), r$run_id, r$wall_s, r$cpu_s / 3600, r$n_hits,
      r$heron_running_max, r$load1_mean))
  }
}
message("Phase 1 done")

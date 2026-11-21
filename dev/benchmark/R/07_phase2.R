# Phase 2 (warm core_nt): what the auto strategy of parallel_blast() needs.
#  a) scaling in N (S100/S500/S2000) for one blastn (B) vs P parallel batches (C)
#  b) threads of a single blastn (T 8..48, mt_mode 2)
#  c) cost and hits of num_alignments 1..50 (new default 10)
#  d) one BLASTr run on S500 to confirm the per-query cost is linear in N
# At most 48 threads per run (shared server). Run from repo root:
#   nohup nice -n 5 Rscript dev/benchmark/R/07_phase2.R 2 > dev/benchmark/results/phase2.log 2>&1 &
source("dev/benchmark/R/bench_lib.R")
set.seed(4)

args <- commandArgs(trailingOnly = TRUE)
reps <- if (length(args)) as.integer(args[[1]]) else 2L

g <- function(strategy, P, T, k = NA, mt_mode = NA, num_alignments = 3, query_set = "S100", r = reps) {
  data.frame(strategy, P, T, k, mt_mode, num_alignments, query_set, r)
}
n_of <- c(S100 = 100, S500 = 500, S2000 = 2000)
grid <- list()
for (qs in c("S500", "S2000")) {
  n <- n_of[[qs]]
  grid[[length(grid) + 1]] <- g("B_single", 1, 32, mt_mode = 2, query_set = qs)
  for (pt in list(c(2, 16), c(4, 8), c(8, 4))) {
    grid[[length(grid) + 1]] <- g("C_batch", pt[1], pt[2], k = ceiling(n / pt[1]), query_set = qs)
  }
  grid[[length(grid) + 1]] <- g("C_batch", 4, 8, k = 25, query_set = qs)
}
for (t in c(8, 48)) grid[[length(grid) + 1]] <- g("B_single", 1, t, mt_mode = 2, query_set = "S500")
for (na in c(1, 5, 10, 25, 50)) grid[[length(grid) + 1]] <- g("B_single", 1, 32, mt_mode = 2, num_alignments = na)
grid[[length(grid) + 1]] <- g("C_batch", 32, 1, k = 4, num_alignments = 10) # exact per-query results, 10 hits
grid[[length(grid) + 1]] <- g("A_blastr", 32, 1, query_set = "S500", r = 1)
grid <- do.call(rbind, grid)

plan <- grid[rep(seq_len(nrow(grid)), times = grid$r), ]
plan$rep <- stats::ave(seq_len(nrow(plan)), paste(plan$strategy, plan$P, plan$T, plan$k, plan$mt_mode,
  plan$num_alignments, plan$query_set), FUN = seq_along)
plan <- plan[sample(nrow(plan)), ]
cat("Phase 2 runs:", nrow(plan), "\n")

for (i in seq_len(nrow(plan))) {
  x <- plan[i, ]
  repeat {
    load1 <- as.numeric(strsplit(readLines("/proc/loadavg"), " ")[[1]][1])
    if (others_blast() == 0 && load1 < 50) break
    message(format(Sys.time()), " waiting: others_blast=", others_blast(), " load1=", load1)
    Sys.sleep(120)
  }
  a <- list(
    phase = "p2", db_name = "core_nt", query_set = x$query_set, strategy = x$strategy,
    P = as.integer(x$P), T = as.integer(x$T), k = as.integer(x$k), mt_mode = as.integer(x$mt_mode),
    num_alignments = as.integer(x$num_alignments), rep = as.integer(x$rep)
  )
  r <- tryCatch(do.call(run_config, c(a, sample_interval = 5)), error = function(e) {
    message("ERROR: ", conditionMessage(e)); NULL
  })
  if (!is.null(r)) {
    message(sprintf("[%d/%d] %s %s wall=%.0fs cpu=%.2fh hits=%d with_hit=%d load=%.0f",
      i, nrow(plan), format(Sys.time(), "%H:%M"), r$run_id, r$wall_s, r$cpu_s / 3600, r$n_hits,
      r$n_queries_with_hit, r$load1_mean))
  }
}
message("Phase 2 done")

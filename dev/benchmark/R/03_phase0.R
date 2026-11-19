# Phase 0: mini multi-marker DB (warm, small). Validates harness, measures the
# fixed per-call overhead and checks result agreement between strategies.
# Run from repo root: nohup Rscript dev/benchmark/R/03_phase0.R > dev/benchmark/results/phase0.log 2>&1 &
source("dev/benchmark/R/bench_lib.R")
set.seed(1)

grid <- rbind(
  expand.grid(strategy = "D_raw", P = c(1, 8, 32), T = 1, k = NA, mt_mode = NA, num_alignments = 3, query_set = "S100"),
  expand.grid(strategy = "A_blastr", P = c(1, 8, 32), T = 1, k = NA, mt_mode = NA, num_alignments = 3, query_set = "S100"),
  expand.grid(strategy = "B_single", P = 1, T = c(1, 8, 32), k = NA, mt_mode = c(0, 1, 2), num_alignments = 3, query_set = "S100"),
  expand.grid(strategy = "C_batch", P = c(8, 32), T = NA, k = c(5, 10, 25, 50), mt_mode = NA, num_alignments = 3, query_set = "S100"),
  # result agreement / cost vs number of alignments
  expand.grid(strategy = c("D_raw", "B_single"), P = NA, T = NA, k = NA, mt_mode = NA,
    num_alignments = c(1, 5, 10, 25, 50), query_set = "S100"),
  # scaling in N
  expand.grid(strategy = c("D_raw", "A_blastr", "B_single", "C_batch"), P = NA, T = NA, k = NA, mt_mode = NA,
    num_alignments = 3, query_set = "S500")
)
grid$T[grid$strategy == "C_batch" & grid$P == 8] <- 4
grid$T[grid$strategy == "C_batch" & grid$P == 32] <- 1
# defaults for the agreement / scaling rows
na_p <- is.na(grid$P)
grid$P[na_p & grid$strategy %in% c("D_raw", "A_blastr")] <- 32
grid$T[is.na(grid$T) & grid$strategy %in% c("D_raw", "A_blastr")] <- 1
grid$P[na_p & grid$strategy == "B_single"] <- 1
grid$T[is.na(grid$T) & grid$strategy == "B_single"] <- 32
grid$mt_mode[is.na(grid$mt_mode) & grid$strategy == "B_single"] <- 0
grid$P[na_p & grid$strategy == "C_batch"] <- 8
grid$T[is.na(grid$T) & grid$strategy == "C_batch"] <- 4
grid$k[is.na(grid$k) & grid$strategy == "C_batch"] <- 25

reps <- 3
plan <- grid[rep(seq_len(nrow(grid)), each = reps), ]
plan$rep <- rep(seq_len(reps), times = nrow(grid))
plan <- plan[sample(nrow(plan)), ]
cat("Phase 0 runs:", nrow(plan), "\n")

for (i in seq_len(nrow(plan))) {
  g <- plan[i, ]
  while (!server_ok()) {
    message(format(Sys.time()), " waiting: ", attr(server_ok(), "detail"))
    Sys.sleep(60)
  }
  args <- list(
    phase = "p0", db_name = "mini", query_set = g$query_set, strategy = as.character(g$strategy),
    P = as.integer(g$P), T = as.integer(g$T), k = as.integer(g$k), mt_mode = as.integer(g$mt_mode),
    num_alignments = as.integer(g$num_alignments), rep = as.integer(g$rep)
  )
  r <- tryCatch(do.call(run_config, args), error = function(e) {
    message("ERROR: ", conditionMessage(e)); NULL
  })
  if (!is.null(r)) {
    message(sprintf("[%d/%d] %s wall=%.2fs cpu=%.1fs hits=%d", i, nrow(plan), r$run_id, r$wall_s, r$cpu_s, r$n_hits))
  }
}
message("Phase 0 done")

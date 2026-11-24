# Phase 3 (warm core_nt): does a second process help beyond ~32 threads?
#   B 1x32 m2 | B 1x64 m2 | C 2x32 m2 (seqs split in 2 batches, one per process)
#   on S500 and S2000, 2 reps each = 12 runs (~1 h).
# Uses 64 of the 80 threads, so it only starts at night (22:00-07:00) once the
# other analyses are over (BLASTR_BENCH_WAIT_FOR, see bench_lib.R), with no
# other user's blastn and a low load. Outside the window it waits for the
# next night; it gives up after max_days.
# Run from repo root:
#   nohup Rscript dev/benchmark/R/08_phase_threads.R > dev/benchmark/results/phase_threads.log 2>&1 &
source("dev/benchmark/R/bench_lib.R")
set.seed(5)

max_days <- 3
night <- function() {
  h <- as.integer(format(Sys.time(), "%H"))
  h >= 22L || h < 7L
}
load1 <- function() as.numeric(strsplit(readLines("/proc/loadavg"), " ")[[1]][1])
server_free <- function() {
  busy <- c(
    busy_processes(),
    others_blastn = others_blast() > 0
  )
  list(ok = !any(busy) && load1() < 12, busy = names(busy)[busy], load = load1())
}

g <- function(strategy, P, T, k, query_set) {
  data.frame(strategy, P, T, k, mt_mode = 2L, num_alignments = 3L, query_set)
}
n_of <- c(S500 = 500, S2000 = 2000)
grid <- do.call(rbind, lapply(names(n_of), function(qs) rbind(
  g("B_single", 1, 32, NA, qs),
  g("B_single", 1, 64, NA, qs),
  g("C_batch", 2, 32, n_of[[qs]] / 2, qs)
)))
plan <- grid[rep(seq_len(nrow(grid)), each = 2L), ]
plan$rep <- rep(1:2, nrow(grid))
plan <- plan[sample(nrow(plan)), ]
cat("Phase threads runs:", nrow(plan), "\n")

deadline <- Sys.time() + max_days * 86400
wait_free <- function() {
  last <- ""
  repeat {
    if (Sys.time() > deadline) stop("gave up: server not free at night within ", max_days, " days")
    s <- server_free()
    if (night() && s$ok) return(invisible())
    msg <- if (!night()) "outside the night window" else
      paste0("busy: ", paste(c(s$busy, if (s$load >= 12) "load1>=12"), collapse = ", "))
    if (msg != last) message(format(Sys.time()), " waiting (", msg, ")")
    last <- msg
    Sys.sleep(60)
  }
}

# Warm-up (not analysed): brings core_nt back to the page cache, since the
# first run of phase 2 was partly cold (350 s instead of ~240 s).
wait_free()
invisible(tryCatch(run_config(phase = "p3warm", db_name = "core_nt", query_set = "S100",
  strategy = "B_single", P = 1L, T = 32L, mt_mode = 2L, num_alignments = 3L, rep = 1L,
  sample_interval = 5), error = function(e) message("warm-up ERROR: ", conditionMessage(e))))
message(format(Sys.time()), " warm-up done")

for (i in seq_len(nrow(plan))) {
  x <- plan[i, ]
  wait_free()
  a <- list(
    phase = "p3", db_name = "core_nt", query_set = x$query_set, strategy = x$strategy,
    P = as.integer(x$P), T = as.integer(x$T), k = as.integer(x$k), mt_mode = 2L,
    num_alignments = 3L, rep = as.integer(x$rep)
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
message("Phase threads done")

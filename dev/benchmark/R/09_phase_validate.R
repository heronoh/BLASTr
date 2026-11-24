# Phase 4 (warm core_nt): validate the released parallel_blast() (strategy E,
# automatic strategy, search_alignments 50, check_n 10) against one raw
# blastn (B) under the same conditions, and against BLASTr 0.1.7 (A).
#   S100 / S500 / S2000: E with 32 cores, B 1x32 mt_mode 2; 2 reps each
#   S2000: E with 64 cores (-> 2 processes x 32 threads); 2 reps
#   S100: A (BLASTr 0.1.7, 32 processes x 1 thread); 1 rep
# num_alignments 10 (the new default) everywhere. ~1 h.
# Starts only when the server is idle: none of the analyses in
# BLASTR_BENCH_WAIT_FOR running (see bench_lib.R), no blastn of anyone,
# load1 < 10. Gives up after max_days.
# Run from repo root:
#   nohup nice -n 5 Rscript dev/benchmark/R/09_phase_validate.R > dev/benchmark/results/phase_validate.log 2>&1 &
source("dev/benchmark/R/bench_lib.R")
set.seed(6)

max_days <- 7
load1 <- function() as.numeric(strsplit(readLines("/proc/loadavg"), " ")[[1]][1])
# "[b]in/blastn": the bracket keeps pgrep from matching the shell that runs it
running <- function(pattern) {
  out <- suppressWarnings(system2("pgrep", c("-f", shQuote(pattern)), stdout = TRUE))
  length(out) > 0L
}
server_free <- function() {
  busy <- c(
    busy_processes(),
    blastn = running("[b]in/blastn ")
  )
  list(ok = !any(busy) && load1() < 10, busy = c(names(busy)[busy], if (load1() >= 10) "load1>=10"))
}

deadline <- Sys.time() + max_days * 86400
wait_free <- function() {
  last <- ""
  repeat {
    if (Sys.time() > deadline) stop("gave up: server not idle within ", max_days, " days")
    s <- server_free()
    if (s$ok) return(invisible())
    msg <- paste(s$busy, collapse = ", ")
    if (msg != last) message(format(Sys.time()), " waiting (busy: ", msg, ")")
    last <- msg
    Sys.sleep(60)
  }
}

g <- function(strategy, P, T, query_set, mt_mode = NA, r = 2L) {
  data.frame(strategy, P, T, mt_mode, query_set, r)
}
grid <- rbind(
  do.call(rbind, lapply(c("S100", "S500", "S2000"), function(qs) rbind(
    g("E_blastr", 1, 32, qs),
    g("B_single", 1, 32, qs, mt_mode = 2)
  ))),
  g("E_blastr", 2, 32, "S2000"),
  g("A_blastr", 32, 1, "S100", r = 1L)
)
plan <- grid[rep(seq_len(nrow(grid)), times = grid$r), ]
plan$rep <- stats::ave(seq_len(nrow(plan)), paste(plan$strategy, plan$P, plan$T, plan$query_set), FUN = seq_along)
plan <- plan[sample(nrow(plan)), ]
cat("Phase validate runs:", nrow(plan), "\n")

# Warm-up (not analysed): core_nt back in the page cache.
wait_free()
invisible(tryCatch(run_config(phase = "p4warm", db_name = "core_nt", query_set = "S100",
  strategy = "B_single", P = 1L, T = 32L, mt_mode = 2L, num_alignments = 10L, rep = 1L,
  sample_interval = 5), error = function(e) message("warm-up ERROR: ", conditionMessage(e))))
message(format(Sys.time()), " warm-up done")

for (i in seq_len(nrow(plan))) {
  x <- plan[i, ]
  wait_free()
  a <- list(
    phase = "p4", db_name = "core_nt", query_set = x$query_set, strategy = x$strategy,
    P = as.integer(x$P), T = as.integer(x$T), mt_mode = as.integer(x$mt_mode),
    num_alignments = 10L, rep = as.integer(x$rep)
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
message("Phase validate done")

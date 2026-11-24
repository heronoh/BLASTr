# Summaries and figures from results/runs.tsv (+ per-run hits.tsv).
# Run from repo root: Rscript dev/benchmark/R/05_analyze.R [phase ...]
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})
`%||%` <- function(a, b) if (is.null(a)) b else a
args <- commandArgs(trailingOnly = TRUE)
res_dir <- "dev/benchmark/results"
fig_dir <- file.path(res_dir, "figs")
dir.create(fig_dir, showWarnings = FALSE)

runs <- utils::read.delim(file.path(res_dir, "runs.tsv"))
# older rows name this column after the watched user
names(runs) <- sub("^[a-z]+_running_max$", "watched_running_max", names(runs))
if (length(args)) runs <- filter(runs, phase %in% args)

runs <- runs |>
  mutate(
    config = case_when(
      strategy == "B_single" ~ sprintf("B single T%d m%d", T, mt_mode),
      strategy == "C_batch" ~ sprintf("C batch k%d %dx%d", k, P, T),
      strategy == "E_blastr" ~ sprintf("E BLASTr auto %d cores", P * T),
      TRUE ~ sprintf("%s %dx%d", strategy, P, T)
    ),
    strategy_label = recode(strategy,
      A_blastr = "A: BLASTr 0.1.7 (1 query/processo)",
      E_blastr = "E: BLASTr novo (estratégia automática)",
      D_raw = "D: blastn puro (1 query/processo)",
      B_single = "B: 1 blastn com todas",
      C_batch = "C: lotes de k queries"
    )
  )

theme_set(theme_minimal(base_size = 11))
save <- function(p, name, w = 9, h = 5.5) {
  ggsave(file.path(fig_dir, paste0(name, ".png")), p, width = w, height = h, dpi = 150, bg = "white")
}

# 1. Summary table ------------------------------------------------------------
summ <- runs |>
  group_by(phase, db, query_set, n_queries, strategy, config, P, T, k, mt_mode, num_alignments) |>
  summarise(
    reps = n(),
    wall_med = median(wall_s), wall_min = min(wall_s), wall_max = max(wall_s),
    cpu_med = median(cpu_s), qph_med = median(queries_per_h),
    pss_gb = median(run_pss_peak_gb), sda_gb = median(sda_read_gb),
    hits_med = median(n_hits), with_hit_med = median(n_queries_with_hit),
    with_hit_min = min(n_queries_with_hit), with_hit_max = max(n_queries_with_hit),
    load = median(load1_mean), watched_running = max(watched_running_max),
    .groups = "drop"
  ) |>
  arrange(phase, db, query_set, num_alignments, wall_med)
utils::write.table(summ, file.path(res_dir, "summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

for (ph in unique(runs$phase)) {
  for (dbn in unique(runs$db[runs$phase == ph])) {
    d <- filter(runs, phase == ph, db == dbn, num_alignments == 3, query_set == "S100")
    if (!nrow(d)) next
    # 2. Strategy comparison (S100, 3 alignments)
    p <- ggplot(d, aes(x = reorder(config, wall_s, median), y = wall_s, colour = strategy_label)) +
      geom_point(size = 2, alpha = 0.8) +
      stat_summary(fun = median, geom = "crossbar", width = 0.5, linewidth = 0.3) +
      coord_flip() +
      labs(
        x = NULL, y = "Tempo total (s)", colour = NULL,
        title = sprintf("100 seqs, 3 alinhamentos - banco %s (%s)", dbn, ph),
        subtitle = "pontos = repeticoes; barra = mediana"
      ) +
      theme(legend.position = "bottom") +
      guides(colour = guide_legend(nrow = 2))
    save(p, sprintf("%s_%s_strategies_wall", ph, dbn))

    p <- ggplot(d, aes(x = wall_s, y = cpu_s / 3600, colour = strategy_label)) +
      geom_point(aes(size = run_pss_peak_gb), alpha = 0.7) +
      scale_x_log10() + scale_y_log10() +
      labs(x = "Tempo total (s, log)", y = "CPU (h, log)", size = "Pico PSS (GB)", colour = NULL,
        title = sprintf("Tempo x CPU por configuracao - %s (%s)", dbn, ph)) +
      theme(legend.position = "right")
    save(p, sprintf("%s_%s_pareto", ph, dbn))
  }

  # 3. Number of alignments: cost and seqs with hit
  d <- filter(runs, phase == ph, query_set == "S100", strategy %in% c("D_raw", "B_single", "C_batch", "A_blastr"))
  d_al <- d |>
    group_by(db, strategy_label, config, num_alignments) |>
    filter(n_distinct(d$num_alignments[d$config == first(config)]) > 2) |>
    summarise(wall = median(wall_s), cpu = median(cpu_s), with_hit = median(n_queries_with_hit), .groups = "drop")
  if (nrow(d_al)) {
    p <- ggplot(d_al, aes(num_alignments, with_hit, colour = config)) +
      geom_line() + geom_point() +
      scale_x_log10(breaks = c(1, 3, 5, 10, 25, 50)) +
      facet_wrap(~db) +
      labs(x = "num_alignments", y = "Seqs com pelo menos 1 hit (de 100)", colour = NULL,
        title = "Seqs com hit x numero de alinhamentos pedidos",
        subtitle = "com -qcov_hsp_perc 80, a lista e cortada antes do filtro de cobertura")
    save(p, sprintf("%s_alignments_with_hit", ph))
    p <- ggplot(d_al, aes(num_alignments, wall, colour = config)) +
      geom_line() + geom_point() +
      scale_x_log10(breaks = c(1, 3, 5, 10, 25, 50)) +
      facet_wrap(~db, scales = "free_y") +
      labs(x = "num_alignments", y = "Tempo total (s)", colour = NULL, title = "Custo x numero de alinhamentos")
    save(p, sprintf("%s_alignments_wall", ph))
  }

  # 4. Scaling in N
  d <- filter(runs, phase == ph, num_alignments == 3)
  d_n <- d |>
    group_by(db, strategy_label, config) |>
    filter(n_distinct(n_queries) > 1) |>
    ungroup()
  if (nrow(d_n)) {
    p <- ggplot(d_n, aes(n_queries, queries_per_h, colour = config)) +
      stat_summary(fun = median, geom = "line") + geom_point(alpha = 0.6) +
      scale_x_log10() + scale_y_log10() + facet_wrap(~db) +
      labs(x = "Numero de seqs (log)", y = "Seqs por hora (log)", colour = NULL, title = "Vazao x tamanho do conjunto")
    save(p, sprintf("%s_scaling_n", ph))
  }

  # 5. P x T heatmap for batches
  d <- filter(runs, phase == ph, strategy == "C_batch", num_alignments == 3, query_set == "S100")
  if (n_distinct(paste(d$P, d$T)) > 3) {
    dh <- d |> group_by(db, P, T, k) |> summarise(wall = median(wall_s), .groups = "drop")
    p <- ggplot(dh, aes(factor(P), factor(T), fill = wall)) +
      geom_tile() + geom_text(aes(label = round(wall)), size = 3) +
      facet_grid(db ~ k, labeller = label_both) +
      scale_fill_viridis_c(direction = -1) +
      labs(x = "Processos (P)", y = "Threads por processo (T)", fill = "Tempo (s)", title = "Lotes: tempo por P x T x k")
    save(p, sprintf("%s_heatmap_PxT", ph), w = 11)
  }
}

# 6. Result agreement vs reference (D_raw, P max, same db/set/alignments; ------
#    without D_raw in the group, the single blastn B with most threads) -------
read_hits <- function(run_id) {
  f <- file.path(res_dir, "runs", run_id, "hits.tsv")
  h <- utils::read.delim(f, colClasses = "character")
  h$subject <- sub("^.*\\|([^|]+)\\|?$", "\\1", h$subject) # gi|..|emb|ACC| -> ACC
  h
}
agree <- list()
for (key in unique(paste(runs$phase, runs$db, runs$query_set, runs$num_alignments))) {
  d <- filter(runs, paste(phase, db, query_set, num_alignments) == key)
  ref <- d |> filter(strategy == "D_raw") |> arrange(desc(P), rep) |> slice(1)
  if (!nrow(ref)) ref <- d |> filter(strategy == "B_single") |> arrange(desc(T), rep) |> slice(1)
  if (!nrow(ref)) next
  rh <- read_hits(ref$run_id)
  ids <- sub("^>", "", grep("^>", readLines(file.path("dev/benchmark/data", paste0(ref$query_set, ".fasta"))), value = TRUE))
  ref_top <- rh |> filter(rank == "1") |> select(seq_id, ref_top = subject)
  ref_set <- split(rh$subject, rh$seq_id)
  for (i in seq_len(nrow(d))) {
    h <- read_hits(d$run_id[i])
    top <- h |> filter(rank == "1") |> select(seq_id, top = subject)
    s <- split(h$subject, h$seq_id)
    cmp <- data.frame(seq_id = ids) |>
      left_join(ref_top, by = "seq_id") |>
      left_join(top, by = "seq_id")
    same_set <- vapply(ids, function(q) setequal(ref_set[[q]] %||% character(), s[[q]] %||% character()), logical(1))
    agree[[length(agree) + 1]] <- data.frame(
      run_id = d$run_id[i], phase = d$phase[i], db = d$db[i], query_set = d$query_set[i],
      num_alignments = d$num_alignments[i], config = d$config[i], rep = d$rep[i],
      same_hit_status_pct = 100 * mean(is.na(cmp$ref_top) == is.na(cmp$top)),
      same_top1_pct = 100 * mean(
        (is.na(cmp$ref_top) & is.na(cmp$top)) |
          (!is.na(cmp$ref_top) & !is.na(cmp$top) & cmp$ref_top == cmp$top)
      ),
      same_subject_set_pct = 100 * mean(same_set)
    )
  }
}
if (length(agree)) {
  agree <- bind_rows(agree)
  utils::write.table(agree, file.path(res_dir, "agreement.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  a <- agree |>
    group_by(phase, db, query_set, num_alignments, config) |>
    summarise(across(ends_with("_pct"), min), .groups = "drop")
  utils::write.table(a, file.path(res_dir, "agreement_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}

# 7. Searches of the projects that use BLASTr (BLASTr_logs) -------------------
# Each project writes <dir>/BLASTr_logs/blastr_runs.tsv (parallel_blast(log_dir =)).
# Paths: BLASTR_BENCH_LOGS (":"-separated), or the project folders inside
# BLASTR_BENCH_PROJECTS (default ~/prjcts).
log_files <- strsplit(Sys.getenv("BLASTR_BENCH_LOGS", ""), ":", fixed = TRUE)[[1]]
if (!length(log_files)) {
  prj <- path.expand(Sys.getenv("BLASTR_BENCH_PROJECTS", "~/prjcts"))
  log_files <- unique(c(
    Sys.glob(file.path(prj, "*", "BLASTr_logs", "blastr_runs.tsv")),
    Sys.glob(file.path(prj, "*", "mtbrcdR_logs", "BLASTr_logs", "blastr_runs.tsv"))
  ))
}
log_files <- log_files[file.exists(log_files)]
if (length(log_files)) {
  plog <- bind_rows(lapply(log_files, function(f) {
    x <- utils::read.delim(f, colClasses = "character", check.names = FALSE)
    x$project <- basename(sub("/(mtbrcdR_logs/)?BLASTr_logs/blastr_runs.tsv$", "", f))
    x$log_file <- f
    x
  }))
  num <- c("n_input", "n_seqs", "len_median", "db_gb", "processes", "threads", "n_batches",
           "wall_s", "cpu_s", "cpu_h", "cpu_efficiency", "pct_with_hits", "pct_best_id", "num_alignments")
  for (k in intersect(num, names(plog))) plog[[k]] <- suppressWarnings(as.numeric(plog[[k]]))
  if (!"function" %in% names(plog)) plog[["function"]] <- "parallel_blast"
  plog[["function"]][is.na(plog[["function"]]) | plog[["function"]] == ""] <- "parallel_blast"
  plog$s_per_seq <- plog$wall_s / plog$n_seqs
  utils::write.table(plog, file.path(res_dir, "project_logs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  plog_summ <- plog |>
    group_by(project, `function`) |>
    summarise(
      runs = n(), seqs = sum(n_seqs, na.rm = TRUE), db_gb = max(db_gb, na.rm = TRUE),
      wall_h = sum(wall_s, na.rm = TRUE) / 3600, cpu_h = sum(cpu_h, na.rm = TRUE),
      s_per_seq_med = median(s_per_seq, na.rm = TRUE),
      strategies = paste(unique(paste0(processes, "x", threads, " b", n_batches)), collapse = "; "),
      .groups = "drop"
    )
  utils::write.table(plog_summ, file.path(res_dir, "project_logs_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  print(as.data.frame(plog_summ))

  # seconds per sequence of the real searches, with the benchmark (core_nt in
  # RAM) for reference: new BLASTr (E), and BLASTr 0.1.7 (A)
  bench <- runs |>
    filter(db == "core_nt", strategy %in% c("A_blastr", "E_blastr"), cache != "cold") |>
    transmute(source = strategy_label, n_seqs = n_queries, s_per_seq = wall_s / n_queries)
  pts <- bind_rows(
    plog |> transmute(source = paste0(project, ": ", `function`), n_seqs, s_per_seq),
    bench
  ) |> filter(is.finite(s_per_seq), s_per_seq > 0)
  p <- ggplot(pts, aes(n_seqs, s_per_seq, colour = source)) +
    geom_point(size = 2.5, alpha = 0.8) +
    scale_x_log10() + scale_y_log10() +
    labs(x = "sequências na busca", y = "segundos por sequência (log)", colour = NULL,
         title = "Buscas reais dos projetos (BLASTr_logs) e benchmark no core_nt") +
    theme(legend.position = "bottom", legend.direction = "vertical")
  save(p, "project_logs_s_per_seq", h = 6.5)
}

cat("Summary rows:", nrow(summ), " figs:", length(list.files(fig_dir)), "\n")

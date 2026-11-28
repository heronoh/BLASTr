# Multi-marker benchmark, step 3: the tie rules (tie_rules_lib.R) on the
# queries of 17_ searched by 18_, as parallel_blast() keeps the hits by
# default (identity >= 80, 100 subjects). Decision procedure fixed before the
# results (hub 2026-11-24-2113, with the review's adjustments of the same time).
#
# Scenarios:
#   in          species in the database (only the template entry left out)
#   in_noident  the same, only queries with no full-length identical hit of
#               their species (a twin of the template, e.g. RefSeq copy or the
#               specimen's own barcode, would make "in" an identical match)
#   out         every hit of the species removed (any species called is wrong)
# Cuts (review of hub 2026-11-24-2228): the call as blast_lca() gives it
# ("none"), and as the mtbrcdR reports it: rank no deeper than the Similarity
# score of the best hit allows (species >= 98 or 97, genus >= 95) and an
# open-nomenclature species ("X sp.") kept at the genus.
# Outputs (results/, one set per cut, suffix "" / "_cut98" / "_cut97"):
#   primer_bench_summary.tsv    % species right / wrong, genus right / wrong by
#                               group, marker, marker x group, primer x group
#   primer_bench_paired.tsv     each rule minus the current one (and the gap
#                               rule minus obitag_floor), same queries, 95%
#                               intervals by bootstrap over templates
#   primer_bench_decision.tsv   per rule and set: significant gain / loss in
#                               species right (in), species wrong (in, out);
#                               passes = no significant loss and at least one
#                               significant gain; tie-break: fewest species
#                               wrong out of the database
#   primer_bench_worse_markers.tsv  markers where a rule is significantly worse
#   primer_bench_grid.tsv       winner of right - p * wrong, mixing "in" (or,
#                               added after the first results, "in_noident")
#                               with a fraction f of absent species (f 0 to
#                               0.8, p 1 to 3)
#   primer_bench_flags.tsv      species calls of obitag_floor by signal and how
#                               many are wrong
# Rules added after the first results (21_primer_bench_extra_rules.R) are read
# from results/primer_bench/rules_res_extra*.rds when they exist.
# Without cut suffix only:
#   primer_bench_truncated.tsv  % of queries whose tied set reaches the 100th hit
#   primer_bench_full_list.tsv  rbcL with every hit of the search (up to 500)
#                               against the default 100
# Per-query results are kept in results/primer_bench/rules_res.rds.
# Run from repo root: Rscript dev/benchmark/R/19_primer_bench_eval.R
suppressMessages(devtools::load_all(".", quiet = TRUE))
suppressPackageStartupMessages(library(dplyr))
source("dev/benchmark/R/tie_rules_lib.R")
set.seed(19)

out_dir <- "dev/benchmark/results/primer_bench"
res_dir <- "dev/benchmark/results"
tax_db <- path.expand("~/.cache/R/BLASTr/taxonomy/taxonomy_2026-11-12.sqlite") # taxonomy of core_nt_nov26
rules_used <- c(
  "current", "current+sp", "cov_k0", "common_prop", "obitag", "obitag_floor", "obitag_floor+sp", "obitag_floor+gap"
)
full_list_markers <- "rbcL"
full_list_rules <- c("current", "obitag_floor", "obitag_floor+sp", "obitag_floor+gap")
cuts <- list(none = NULL, cut98 = c(species = 98, genus = 95), cut97 = c(species = 97, genus = 95))
min_power <- 50L # fewer queries in a set: "no power to tell"
q <- readRDS(file.path(out_dir, "queries.rds"))

# one truth per distinct sequence (sequences shared by two species left out);
# the template of a sequence is its bootstrap cluster
b <- bench_truth(q)
truth <- b$truth
q <- b$q
qq <- b$qq
full_setting <- function(h) {
  h <- h[h$pident >= 80, ]
  h$blast_rank <- as.integer(stats::ave(seq_len(nrow(h)), h$qid, FUN = seq_along))
  h
}

res_rds <- file.path(out_dir, "rules_res.rds")
if (file.exists(res_rds)) {
  cached <- readRDS(res_rds)
} else {
  hs <- bench_hits(out_dir, truth, tax_db)
  s <- hs$s
  twin <- unique(s$qid[(s$species == s$true_species) %in% TRUE & s$nident == s$qlen & s$length == s$qlen])
  default <- run_scenarios(s, hs$tax, tax_db, default_setting, 100L, rules_used)
  fl_qid <- match(qq$seq[qq$marker %in% full_list_markers], truth$seq)
  full <- run_scenarios(s[s$qid %in% fl_qid, ], hs$tax, tax_db, full_setting, 500L, full_list_rules)
  cached <- list(twin = twin, default = default, full = full)
  saveRDS(cached, res_rds)
  rm(s, hs)
}
for (extra_rds in list.files(out_dir, "^rules_res_extra.*[.]rds$", full.names = TRUE)) {
  extra <- readRDS(extra_rds)
  for (sc in names(cached$default$runs)) {
    cached$default$runs[[sc]] <- rbind(cached$default$runs[[sc]], extra$runs[[sc]][names(cached$default$runs[[sc]])])
  }
  cached$default$flags <- rbind(cached$default$flags, extra$flags)
  message("extra rules: ", paste(unique(extra$runs[["in"]]$rule), collapse = ", "))
}
twin <- cached$twin
message("queries with a full-length identical hit of their species: ", length(twin), " of ", nrow(truth))

# one row per primer and query, with the scenario "in_noident"
as_rows <- function(runs) {
  runs$in_noident <- runs[["in"]][!runs[["in"]]$qid %in% twin, ]
  runs$in_noident$scenario <- "in_noident"
  res <- do.call(rbind, runs)
  res$setting <- "default"
  res$true_species <- truth$species[res$qid]
  res$true_genus <- truth$genus[res$qid]
  res$seq <- truth$seq[res$qid]
  res$cluster <- truth$acc[res$qid]
  x <- inner_join(res, qq, by = "seq", relationship = "many-to-many")
  x$marker_group <- paste(x$marker, x$group)
  x$primer_group <- paste(x$primer, x$group)
  x
}

# the call as the mtbrcdR reports it: no deeper than the score of the best hit
# allows; "X sp." at the genus
apply_cut <- function(x, cut) {
  if (is.null(cut)) return(x)
  k <- x$rank %in% "species" & !(x$score >= cut[["species"]]) %in% TRUE
  x$rank[k] <- "genus"
  x$species[k] <- NA_character_
  k <- x$rank %in% c("species", "genus") & !(x$score >= cut[["genus"]]) %in% TRUE
  x$rank[k] <- "above genus"
  x$species[k] <- NA_character_
  x$genus[k] <- NA_character_
  k <- x$rank %in% "species" & grepl(open_nomenclature_pattern, x$species, perl = TRUE)
  x$rank[k] <- "genus"
  x$species[k] <- NA_character_
  x
}

boot_cluster <- function(v, cl, B = 1000L) {
  sm <- tapply(v, cl, sum)
  nn <- tapply(v, cl, length)
  b <- vapply(seq_len(B), function(i) {
    k <- sample.int(length(sm), replace = TRUE)
    sum(sm[k]) / sum(nn[k])
  }, numeric(1L))
  c(mean = 100 * sum(sm) / sum(nn), lo = 100 * unname(stats::quantile(b, 0.025)), hi = 100 * unname(stats::quantile(b, 0.975)))
}

# paired differences of every rule against `against`, bootstrap over templates
paired_diff <- function(x, against) {
  base <- x[x$rule == against, c("scenario", "primer", "seq", "sp_ok", "sp_wrong", "gen_wrong")]
  names(base)[4:6] <- c("ok0", "wrong0", "gwrong0")
  pd <- inner_join(x[x$rule != against, ], base, by = c("scenario", "primer", "seq"))
  bind_rows(lapply(c("all", "group", "marker"), function(k) {
    pd |>
      mutate(set = if (k == "all") "all" else as.character(.data[[k]])) |>
      group_by(scenario, rule, set) |>
      group_modify(function(d, g) {
        dr <- d$sp_ok - d$ok0
        dw <- d$sp_wrong - d$wrong0
        r <- rbind(
          right = boot_cluster(dr, d$cluster), wrong = boot_cluster(dw, d$cluster),
          genus_wrong = boot_cluster(d$gen_wrong - d$gwrong0, d$cluster),
          `right-1wrong` = boot_cluster(dr - dw, d$cluster), `right-2wrong` = boot_cluster(dr - 2 * dw, d$cluster)
        )
        tibble::tibble(
          n = nrow(d), templates = length(unique(d$cluster)), measure = rownames(r),
          mean = round(r[, "mean"], 2), lo = round(r[, "lo"], 2), hi = round(r[, "hi"], 2)
        )
      }) |>
      ungroup()
  })) |>
    mutate(against = against, .before = 1L)
}

sets <- c("group", "marker", "marker_group", "primer_group")
x0 <- as_rows(cached$default$runs)
fl0 <- cached$default$flags

for (cn in names(cuts)) {
  sfx <- if (cn == "none") "" else paste0("_", cn)
  wr <- function(d, name) write.table(d, file.path(res_dir, paste0("primer_bench_", name, sfx, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  x <- apply_cut(x0, cuts[[cn]])
  x$sp_ok <- x$rank %in% "species" & (x$species == x$true_species) %in% TRUE
  x$sp_wrong <- x$rank %in% "species" & !x$sp_ok
  x$gen_wrong <- x$rank %in% c("species", "genus") & !(x$genus == x$true_genus) %in% TRUE

  out <- summarise_rules(x, sets)
  out$power <- ifelse(out$n < min_power, "no power to tell", "")
  wr(out, "summary")

  if (cn == "none") {
    tr <- bind_rows(lapply(c(sets, "all"), function(k) {
      x |>
        mutate(set = if (k == "all") "all" else as.character(.data[[k]])) |>
        group_by(scenario, rule, set) |>
        summarise(n = n(), truncated = round(100 * mean(truncated), 1), .groups = "drop")
    }))
    wr(tr, "truncated")
  }

  paired <- bind_rows(paired_diff(x, "current"), paired_diff(x[x$rule %in% c("obitag_floor", "obitag_floor+gap"), ], "obitag_floor"))
  wr(paired, "paired")

  # decision: no significant loss and at least one significant gain
  pc <- paired[paired$against == "current", ]
  sig <- function(sc, ms, better) {
    d <- pc[pc$scenario == sc & pc$measure == ms, c("rule", "set", "lo", "hi")]
    d$gain <- if (better == "up") d$lo > 0 else d$hi < 0
    d$loss <- if (better == "up") d$hi < 0 else d$lo > 0
    d$crit <- paste(sc, ms)
    d
  }
  crit <- bind_rows(
    sig("in", "right", "up"), sig("in", "wrong", "down"),
    sig("in_noident", "right", "up"), sig("in_noident", "wrong", "down"),
    sig("out", "wrong", "down")
  )
  decision <- crit |>
    filter(set %in% c("all", "vertebrates", "invertebrates", "plants")) |>
    group_by(rule, set) |>
    summarise(
      gains = paste(crit[gain], collapse = "; "), losses = paste(crit[loss], collapse = "; "),
      passes = !any(loss) && any(gain), .groups = "drop"
    )
  wout <- out[out$scenario == "out" & out$set == "all", c("rule", "species_wrong")]
  decision$out_wrong_all <- wout$species_wrong[match(decision$rule, wout$rule)]
  wr(decision, "decision")

  worse <- crit |>
    filter(!set %in% c("all", "vertebrates", "invertebrates", "plants"), loss) |>
    select(rule, set, crit)
  wr(worse, "worse_markers")

  # winner of right - p * wrong over a grid of absent fraction f and weight p
  grid <- expand.grid(
    in_scenario = c("in", "in_noident"), f = seq(0, 0.8, by = 0.1), p = c(1, 1.5, 2, 2.5, 3),
    set = c("all", "vertebrates", "invertebrates", "plants"), stringsAsFactors = FALSE
  )
  grid$winner <- NA_character_
  grid$margin <- NA_real_
  for (i in seq_len(nrow(grid))) {
    a <- out[out$set == grid$set[i], ]
    u <- vapply(unique(a$rule), function(r) {
      ri <- a[a$rule == r & a$scenario == grid$in_scenario[i], ]
      ro <- a[a$rule == r & a$scenario == "out", ]
      (1 - grid$f[i]) * (ri$species_right - grid$p[i] * ri$species_wrong) - grid$f[i] * grid$p[i] * ro$species_wrong
    }, numeric(1L))
    o <- sort(u, decreasing = TRUE)
    grid$winner[i] <- names(o)[1L]
    grid$margin[i] <- round(o[1L] - o[2L], 2)
  }
  wr(grid, "grid")

  # signals of obitag_floor on its species calls
  ofl <- x[x$rule == "obitag_floor" & x$rank %in% "species" & x$scenario %in% c("in", "out"), ]
  ofl$flag <- "none"
  if (!is.null(fl0) && nrow(fl0)) {
    fl <- fl0[fl0$rule == "obitag_floor", ]
    k <- match(paste(ofl$scenario, ofl$qid), paste(fl$scenario, fl$query))
    ofl$flag[!is.na(k)] <- fl$flag[k[!is.na(k)]]
  }
  fsum <- ofl |>
    group_by(scenario, group, flag) |>
    summarise(species_calls = n(), wrong = round(100 * mean(sp_wrong), 1), .groups = "drop")
  wr(fsum, "flags")

  message("== cut: ", cn)
  print(out[out$set %in% c("all", "vertebrates", "invertebrates", "plants"), ], row.names = FALSE)
  print(as.data.frame(decision), row.names = FALSE)
}

# rbcL: the default 100 hits against every hit of the search (up to 500)
xf <- as_rows(cached$full$runs)
xf <- xf[xf$marker %in% full_list_markers, ]
xd <- x0[x0$marker %in% full_list_markers & x0$rule %in% full_list_rules, ]
fl_out <- bind_rows(lapply(names(cuts), function(cn) {
  bind_rows(
    mutate(summarise_rules(apply_cut(xd, cuts[[cn]]), "marker"), hits = "100"),
    mutate(summarise_rules(apply_cut(xf, cuts[[cn]]), "marker"), hits = "500")
  ) |>
    filter(set != "all") |>
    mutate(cut = cn)
}))
trf <- bind_rows(
  summarise(group_by(xd, scenario, rule), truncated = round(100 * mean(truncated), 1), .groups = "drop") |> mutate(hits = "100"),
  summarise(group_by(xf, scenario, rule), truncated = round(100 * mean(truncated), 1), .groups = "drop") |> mutate(hits = "500")
)
fl_out <- left_join(fl_out, trf, by = c("scenario", "rule", "hits"))
write.table(fl_out, file.path(res_dir, "primer_bench_full_list.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(fl_out[, c("cut", "hits", "scenario", "rule", "n", "species_right", "species_wrong", "genus_right", "genus_wrong", "median_tied", "truncated")], row.names = FALSE)

message("queries: ", length(unique(x0$seq)), " sequences, ", nrow(qq), " primer-query pairs, ", length(unique(x0$cluster)), " templates")
message("finished")

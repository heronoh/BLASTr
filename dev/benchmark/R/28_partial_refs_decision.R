# Decision among the variants of the rule by references (weight of the query
# bases a hit does not cover, in the choice of the best hits and in the
# window), by the criterion fixed before the results (hub 2026-11-25-0246),
# on the runs of 27_primer_bench_partial_refs.R.
#
# Against the delivered rule (refs_b1_w1), same queries (one row per primer
# and query), calls cut at 98 by the Similarity score of the best hit tied by
# each rule; bootstrap over templates (1000).
#   1. no significant loss (95%) at any coverage level (1, 0.97, 0.9, 0.8), in
#      all queries and in each group (vertebrates, invertebrates, plants):
#      species right "in", species wrong "in", species wrong "out", genus wrong
#      "out";
#   2. a significant gain (97.5%, two measures) in the main measure, the
#      paired difference averaged over the four levels, all queries: genus
#      right "out" or species right "in".
# More than one passing: least species wrong "out" (mean over the levels),
# then the largest window weight.
# Output (results/): primer_bench_partial_decision.tsv (per variant: losses,
# gains, passes, chosen, cells too wide) and primer_bench_partial_paired.tsv
# (every paired difference, with `wide` when its interval allows a worsening
# larger than the variant's main gain and `no_power` under 50 rows).
# Run from repo root: Rscript dev/benchmark/R/28_partial_refs_decision.R
suppressMessages(devtools::load_all(".", quiet = TRUE))
suppressPackageStartupMessages(library(dplyr))
source("dev/benchmark/R/tie_rules_lib.R")
set.seed(28)

out_dir <- "dev/benchmark/results/primer_bench"
res_dir <- "dev/benchmark/results"
base_rule <- "refs_b1_w1"
cut98 <- c(species = 98, genus = 95)
groups3 <- c("vertebrates", "invertebrates", "plants")

b <- bench_truth(readRDS(file.path(out_dir, "queries.rds")))
truth <- b$truth
# TEST_DIR: the quick check of 27_ (TEST_N)
test_dir <- Sys.getenv("TEST_DIR")
if (nzchar(test_dir)) out_dir <- res_dir <- test_dir
files <- list.files(out_dir, "^rules_res_partial_.*[.]rds$", full.names = TRUE)
stopifnot(nzchar(test_dir) || length(files) == 8L)
x <- do.call(rbind, lapply(files, readRDS))

# calls as mtbrcdR reports them, cut by the score of each rule's best tied hit
k <- x$rank %in% "species" & !(x$score_tied >= cut98[["species"]]) %in% TRUE
x$rank[k] <- "genus"
x$species[k] <- NA_character_
k <- x$rank %in% c("species", "genus") & !(x$score_tied >= cut98[["genus"]]) %in% TRUE
x$rank[k] <- "above genus"
x$species[k] <- NA_character_
x$genus[k] <- NA_character_
k <- x$rank %in% "species" & grepl(open_nomenclature_pattern, x$species, perl = TRUE)
x$rank[k] <- "genus"
x$species[k] <- NA_character_
x$sp_ok <- x$rank %in% "species" & (x$species == truth$species[x$qid]) %in% TRUE
x$sp_wrong <- x$rank %in% "species" & !x$sp_ok
x$gen_ok <- x$rank %in% c("species", "genus") & (x$genus == truth$genus[x$qid]) %in% TRUE
x$gen_wrong <- x$rank %in% c("species", "genus") & !x$gen_ok
x$seq <- truth$seq[x$qid]
x$cluster <- truth$acc[x$qid]
x <- inner_join(x, b$qq, by = "seq", relationship = "many-to-many")

measures <- c("sp_ok", "sp_wrong", "gen_ok", "gen_wrong")
base <- x[x$rule == base_rule, c("coverage", "scenario", "primer", "seq", measures)]
names(base)[match(measures, names(base))] <- paste0(measures, "0")
pd <- inner_join(x[x$rule != base_rule, ], base, by = c("coverage", "scenario", "primer", "seq"))
for (m in measures) pd[[paste0("d_", m)]] <- as.numeric(pd[[m]]) - as.numeric(pd[[paste0(m, "0")]])

# mean of v and its bootstrap interval over clusters (templates)
boot_ci <- function(v, cl, level, B = 1000L) {
  if (!length(v)) return(c(mean = NA_real_, lo = NA_real_, hi = NA_real_))
  sm <- tapply(v, cl, sum)
  nn <- tapply(v, cl, length)
  bs <- vapply(seq_len(B), function(i) {
    k <- sample.int(length(sm), replace = TRUE)
    sum(sm[k]) / sum(nn[k])
  }, numeric(1L))
  a <- (1 - level) / 2
  c(mean = 100 * sum(sm) / sum(nn), lo = 100 * unname(stats::quantile(bs, a)),
    hi = 100 * unname(stats::quantile(bs, 1 - a)))
}

# 1. constraints: per coverage level, all and each group
cons <- list(c("in", "sp_ok", "down"), c("in", "sp_wrong", "up"), c("out", "sp_wrong", "up"), c("out", "gen_wrong", "up"))
c_rows <- list()
for (cc in cons) {
  d0 <- pd[pd$scenario == cc[1], ]
  for (set in c("all", groups3)) {
    d1 <- if (set == "all") d0 else d0[d0$group == set, ]
    r <- d1 |>
      group_by(rule, coverage) |>
      group_modify(function(d, g) {
        ci <- boot_ci(d[[paste0("d_", cc[2])]], d$cluster, 0.95)
        tibble::tibble(n = nrow(d), mean = ci[["mean"]], lo = ci[["lo"]], hi = ci[["hi"]])
      }) |>
      ungroup()
    r$scenario <- cc[1]
    r$measure <- cc[2]
    r$set <- set
    r$kind <- "constraint"
    r$loss <- (if (cc[3] == "down") r$hi < 0 else r$lo > 0) %in% TRUE
    c_rows[[length(c_rows) + 1L]] <- r
  }
}
cons_tab <- bind_rows(c_rows)

# 2. main measure: averaged over the coverage levels, all queries, 97.5%
mains <- list(c("out", "gen_ok"), c("in", "sp_ok"))
m_rows <- list()
for (mm in mains) {
  d0 <- pd[pd$scenario == mm[1], ]
  r <- d0 |>
    group_by(rule) |>
    group_modify(function(d, g) {
      ci <- boot_ci(d[[paste0("d_", mm[2])]], d$cluster, 0.975)
      tibble::tibble(n = nrow(d), mean = ci[["mean"]], lo = ci[["lo"]], hi = ci[["hi"]])
    }) |>
    ungroup()
  r$coverage <- NA_real_
  r$scenario <- mm[1]
  r$measure <- mm[2]
  r$set <- "all"
  r$kind <- "main"
  r$gain <- (r$lo > 0) %in% TRUE
  m_rows[[length(m_rows) + 1L]] <- r
}
main_tab <- bind_rows(m_rows)

# description (review, not a criterion): a constraint cell whose interval
# allows a worsening larger than the variant's main gain is marked `wide`;
# fewer than 50 rows, `no_power`
ref_gain <- main_tab |> group_by(rule) |> summarise(ref_gain = max(0, mean), .groups = "drop")
cons_tab <- cons_tab |>
  left_join(ref_gain, by = "rule") |>
  mutate(worst = ifelse(measure == "sp_ok", -lo, hi), wide = (worst > ref_gain) %in% TRUE, no_power = n < 50L) |>
  select(-worst)

# tie-break: species wrong "out", mean over the levels (not paired)
tb <- x |>
  filter(scenario == "out") |>
  group_by(rule) |>
  summarise(sp_wrong_out = round(100 * mean(sp_wrong), 2), .groups = "drop")
window_w <- function(r) c(`1` = 1, `11` = 1 / 11, `0` = 0)[sub("^.*_w", "", r)]

dec <- cons_tab |>
  group_by(rule) |>
  summarise(losses = sum(loss), loss_where = paste(unique(paste(scenario, measure, set, coverage)[loss]), collapse = "; "),
            wide = sum(wide), no_power = sum(no_power), .groups = "drop") |>
  left_join(main_tab |> group_by(rule) |>
              summarise(gains = sum(gain), gain_where = paste(paste(scenario, measure)[gain], collapse = "; "),
                        .groups = "drop"), by = "rule") |>
  left_join(tb, by = "rule") |>
  filter(rule != "current") |>
  mutate(passes = losses == 0 & gains > 0)
ok <- dec[dec$passes, ]
dec$chosen <- FALSE
if (nrow(ok)) {
  ok <- ok[order(ok$sp_wrong_out, -window_w(ok$rule)), ]
  dec$chosen[dec$rule == ok$rule[1]] <- TRUE
}
dec <- dec[order(match(dec$rule, rule_order)), ]

paired <- bind_rows(cons_tab, main_tab) |> mutate(across(c(mean, lo, hi), ~ round(.x, 2)))
write.table(dec, file.path(res_dir, "primer_bench_partial_decision.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(paired, file.path(res_dir, "primer_bench_partial_paired.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(as.data.frame(dec[, c("rule", "losses", "wide", "no_power", "gains", "gain_where", "sp_wrong_out", "passes", "chosen")]), row.names = FALSE)
print(as.data.frame(main_tab[, c("rule", "scenario", "measure", "mean", "lo", "hi", "gain")]), row.names = FALSE, digits = 3)
message(if (nrow(ok)) paste("chosen:", ok$rule[1]) else "no variant passes: at_best_refs stays as delivered")

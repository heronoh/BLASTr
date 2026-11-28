# Multi-marker benchmark: tie rules added after the first results, on the same
# hits (no new search). 19_primer_bench_eval.R reads them with the others.
#   obitag_floor_ident      the floor of obitag_floor also when the best hit is
#                           identical to the query (HO, 2026-11-24, after the
#                           long rbcL amplicon: references within 0.5% of an
#                           identical best hit join the tied set)
#   obitag_floor_ident+gap  the same with the gap of obitag_floor+gap
# Output: results/primer_bench/rules_res_extra.rds
# Run from repo root: Rscript dev/benchmark/R/21_primer_bench_extra_rules.R
suppressMessages(devtools::load_all(".", quiet = TRUE))
suppressPackageStartupMessages(library(dplyr))
source("dev/benchmark/R/tie_rules_lib.R")

out_dir <- "dev/benchmark/results/primer_bench"
tax_db <- path.expand("~/.cache/R/BLASTr/taxonomy/taxonomy_2026-11-12.sqlite") # taxonomy of core_nt_nov26
extra_rules <- c("obitag_floor_ident", "obitag_floor_ident+gap")

b <- bench_truth(readRDS(file.path(out_dir, "queries.rds")))
hs <- bench_hits(out_dir, b$truth, tax_db)
extra <- run_scenarios(hs$s, hs$tax, tax_db, default_setting, 100L, extra_rules)
saveRDS(extra, file.path(out_dir, "rules_res_extra.rds"))
message("finished")

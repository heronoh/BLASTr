# Rank calibration on EM291 plant barcodes confirmed by barcodes (client
# data: read from the barcodes project, only aggregated numbers written).
# Truth: the expected genus (specialist, confirmed by both markers).
source("dev/benchmark/R/11_rank_calibration_eval.R")
f <- "/home/noreh/prjcts/barcodes/analyses/eDNA_vale_dmx/EM291_os1e2/experimentos/blast_qcov80/calibracao--barcodes_confirmados.tsv"
x <- read.delim(f, colClasses = "character")
out_dir <- "dev/benchmark/results/rank_calibration"
hits_rds <- file.path(out_dir, "em291_hits.rds")
if (!file.exists(hits_rds)) {
  s <- parallel_blast_hsps(
    unique(x$sequence), "/data/databases/core_nt_nov26/core_nt",
    total_cores = 32L, max_target_seqs = 100L, max_hsps = 1L, perc_id = 70, perc_qcov_hsp = 50,
    fields = c("saccver", "pident", "length", "qcovhsp", "bitscore", "staxid"),
    verbose = "silent"
  )
  saveRDS(s, hits_rds)
}
hits <- readRDS(hits_rds)
et <- suppressWarnings(get_tax(x$expected_taxid))
truth <- data.frame(
  Sequence = x$sequence, set = x$marcador, species = NA_character_,
  genus = et$genus[match(x$expected_taxid, et$taxid)], family = et$family[match(x$expected_taxid, et$taxid)]
)
truth <- truth[!duplicated(truth$Sequence), ]
res <- eval_ranking(hits, truth)
sm <- summarise_ranking(res)
write.table(sm, "dev/benchmark/results/rank_calibration_em291.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
print(sm[, c("set", "criterion", "n", "top1_genus", "top1_family", "lca_wrong", "lca_genus_ok", "median_tie")], row.names = FALSE)
dis <- disagreements(res)
cat("\nsequences where the criteria disagree on the genus:", length(unique(dis$Sequence)), "\n")
d <- summarise_ranking(dis)
print(d[d$set == "all", c("criterion", "n", "top1_genus", "lca_wrong", "lca_genus_ok")], row.names = FALSE)

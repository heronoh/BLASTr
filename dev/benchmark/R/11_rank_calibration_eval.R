# Evaluate the ranking criteria on queries of known taxonomy.
# eval_ranking(hits, truth):
#   hits:  Sequence, saccver, pident, length, qcovhsp, bitscore, staxid
#          (parallel_blast_hsps(), max_hsps 1), BLAST order kept in `hit`
#   truth: Sequence, set, species, genus, family (NA when unknown)
# For each criterion, (a) the first hit and (b) the lowest common rank of the
# hits tied with the best one are compared with the truth.
suppressMessages(devtools::load_all(".", quiet = TRUE))

criteria <- function() {
  w <- c(1, 2, 3, 5, 10, 20, 50)
  c(
    list(identity = list(key = function(h) h$pident, tol = function(k, h) k >= max(k) - 0.5)),
    stats::setNames(lapply(w, function(wi) list(
      key = function(h) (wi * h$pident + h$qcovhsp) / (wi + 1),
      tol = function(k, h) k >= max(k) - 0.5
    )), paste0("w", w)),
    list(bitscore = list(key = function(h) h$bitscore, tol = function(k, h) k >= 0.98 * max(k))),
    # barcodes: best by bitscore; tied = bitscore within 2% or identity >= the best one's
    list(combined = list(key = function(h) h$bitscore, tol = function(k, h) k >= 0.98 * max(k) | h$pident >= h$pident[1]))
  )
}

eval_ranking <- function(hits, truth, taxonomy_db = NULL) {
  tax <- suppressWarnings(get_tax(hits$staxid, taxonomy_db = taxonomy_db))
  hits <- dplyr::left_join(hits, tax[, c("taxid", "species", "genus", "family")], by = c(staxid = "taxid"))
  hits <- hits[order(hits$Sequence, hits$hit), ]
  by_q <- split(hits, hits$Sequence)
  by_q <- by_q[names(by_q) %in% truth$Sequence]
  rows <- list()
  for (cn in names(criteria())) {
    cr <- criteria()[[cn]]
    r <- lapply(names(by_q), function(sq) {
      h <- by_q[[sq]]
      k <- cr$key(h)
      o <- order(-k, -h$pident, -h$qcovhsp, -h$bitscore, h$hit)
      h <- h[o, ]
      k <- k[o]
      tie <- h[cr$tol(k, h), ]
      lca <- if (length(unique(tie$species)) == 1L && !anyNA(tie$species)) "species" else
        if (length(unique(tie$genus)) == 1L && !anyNA(tie$genus)) "genus" else
          if (length(unique(tie$family)) == 1L && !anyNA(tie$family)) "family" else "above"
      data.frame(
        Sequence = sq, criterion = cn,
        top_species = h$species[1], top_genus = h$genus[1], top_family = h$family[1],
        top_id = h$pident[1], top_qcov = h$qcovhsp[1],
        lca_rank = lca, lca_species = if (lca == "species") tie$species[1] else NA,
        lca_genus = if (lca %in% c("species", "genus")) tie$genus[1] else NA,
        lca_family = if (lca %in% c("species", "genus", "family")) tie$family[1] else NA,
        n_tie = nrow(tie)
      )
    })
    rows[[cn]] <- do.call(rbind, r)
  }
  res <- dplyr::left_join(do.call(rbind, rows), truth, by = "Sequence", suffix = c("", "_true"))
  res$criterion <- factor(res$criterion, names(criteria()))
  res
}

ok <- function(pred, true) ifelse(is.na(true), NA, !is.na(pred) & pred == true)

summarise_ranking <- function(res) {
  res$sp1 <- ok(res$top_species, res$species)
  res$gen1 <- ok(res$top_genus, res$genus)
  res$fam1 <- ok(res$top_family, res$family)
  # LCA: right if every rank it gives matches the truth; resolved at species / genus
  res$lca_wrong <- (!is.na(res$lca_species) & !is.na(res$species) & res$lca_species != res$species) |
    (!is.na(res$lca_genus) & !is.na(res$genus) & res$lca_genus != res$genus) |
    (!is.na(res$lca_family) & !is.na(res$family) & res$lca_family != res$family)
  res$lca_sp_ok <- ok(res$lca_species, res$species)
  res$lca_gen_ok <- ok(res$lca_genus, res$genus)
  agg <- function(d) data.frame(
    n = nrow(d),
    top1_species = round(100 * mean(d$sp1, na.rm = TRUE), 1),
    top1_genus = round(100 * mean(d$gen1, na.rm = TRUE), 1),
    top1_family = round(100 * mean(d$fam1, na.rm = TRUE), 1),
    lca_wrong = round(100 * mean(d$lca_wrong, na.rm = TRUE), 1),
    lca_species_ok = round(100 * mean(d$lca_sp_ok, na.rm = TRUE), 1),
    lca_genus_ok = round(100 * mean(d$lca_gen_ok, na.rm = TRUE), 1),
    median_tie = stats::median(d$n_tie)
  )
  sets <- split(res, list(res$set, res$criterion), drop = TRUE)
  out <- do.call(rbind, lapply(names(sets), function(k) cbind(
    set = sets[[k]]$set[1], criterion = as.character(sets[[k]]$criterion[1]), agg(sets[[k]])
  )))
  all <- do.call(rbind, lapply(split(res, res$criterion), function(d) cbind(set = "all", criterion = as.character(d$criterion[1]), agg(d))))
  rbind(out, all)
}

# Queries where the criteria disagree on the genus of the first hit.
disagreements <- function(res) {
  g <- tapply(res$top_genus, res$Sequence, function(x) length(unique(x)) > 1L)
  res[res$Sequence %in% names(g)[g], ]
}

if (sys.nframe() == 0L) {
  d <- "dev/benchmark/results/rank_calibration"
  hits <- readRDS(file.path(d, "hits.rds"))
  ref <- readRDS(file.path(d, "references.rds"))
  truth <- data.frame(Sequence = ref$query, set = ref$marker, species = ref$species, genus = ref$genus, family = ref$family)
  res <- eval_ranking(hits, truth)
  sm <- summarise_ranking(res)
  write.table(sm, file.path("dev/benchmark/results", "rank_calibration_public.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  dis <- disagreements(res)
  dsm <- summarise_ranking(dis)
  write.table(dsm, file.path("dev/benchmark/results", "rank_calibration_public_disagree.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  print(sm[sm$set == "all", ], row.names = FALSE)
  cat("\nqueries where the criteria disagree on the genus:", length(unique(dis$Sequence)), "\n")
  print(dsm[dsm$set == "all", ], row.names = FALSE)
}

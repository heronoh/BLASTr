# Shared code of the tie-rule calibrations (15_tie_rules.R on public
# references, 16_tie_rules_real_seqs.R on real sequences). See 15_tie_rules.R for the rules.
# Needs the package loaded (devtools::load_all()) and dplyr.

# voucher codes in the GenBank titles ("... voucher LBPV-31886 cytochrome ...")
parse_voucher <- function(title) {
  m <- regmatches(title, regexec("voucher ((?:[A-Z0-9][A-Za-z0-9:_./-]*(?: |$)){1,3})", title, perl = TRUE))
  v <- vapply(m, function(x) if (length(x) > 1L) x[2L] else NA_character_, character(1L))
  v <- sub("\\s*\\b(12S|16S|18S|28S|5\\.8S|COI|CO1|COX1|COXI|ITS[12]?|ND[0-9]|CYTB)\\b.*$", "", v, perl = TRUE)
  v <- toupper(gsub("[^A-Za-z0-9]", "", v))
  ifelse(is.na(v) | !nzchar(v), NA_character_, v)
}
# query positions (in query coordinates) of the differences of each hit, from
# btop: a substitution or a gap in the subject at its query position, an
# insertion in the subject at the next query position
btop_diffs <- function(btop, qstart) {
  tok <- regmatches(btop, gregexpr("[0-9]+|[A-Za-z*-]{2}", btop))
  mapply(function(t, q0) {
    isnum <- grepl("^[0-9]", t)
    n <- suppressWarnings(as.integer(ifelse(isnum, t, "0")))
    adv <- ifelse(isnum, n, ifelse(substr(t, 1L, 1L) == "-", 0L, 1L))
    pos <- q0 + c(0L, cumsum(adv))[seq_along(t)]
    as.integer(pos[!isnum])
  }, tok, as.integer(qstart), SIMPLIFY = FALSE, USE.NAMES = FALSE)
}
n_in <- function(p, a, b) sum(p >= a & p <= b)
`%||%` <- function(a, b) if (is.null(a)) b else a

# What each hit has where it differs from the query, from btop: "p<pos>" =
# the subject base (or "-") at that query position; "i<pos>" = the bases the
# subject inserts before that query position. Two hits that both match the
# query at a position agree there.
btop_states <- function(btop, qstart) {
  tok <- regmatches(btop, gregexpr("[0-9]+|[A-Za-z*-]{2}", btop))
  mapply(function(t, q0) {
    isnum <- grepl("^[0-9]", t)
    n <- suppressWarnings(as.integer(ifelse(isnum, t, "0")))
    qgap <- !isnum & substr(t, 1L, 1L) == "-"
    adv <- ifelse(isnum, n, ifelse(qgap, 0L, 1L))
    pos <- q0 + c(0L, cumsum(adv))[seq_along(t)]
    ch <- toupper(substr(t, 2L, 2L))
    k <- !isnum & !qgap
    v <- if (any(k)) stats::setNames(ch[k], paste0("p", pos[k])) else stats::setNames(character(0L), character(0L))
    if (any(qgap)) {
      ins <- tapply(ch[qgap], pos[qgap], paste, collapse = "")
      v <- c(v, stats::setNames(as.character(ins), paste0("i", names(ins))))
    }
    v
  }, tok, as.integer(qstart), SIMPLIFY = FALSE, USE.NAMES = FALSE)
}

# Differences between two hits (references) on the query region lo..hi,
# from their states.
pair_diff <- function(a, b, lo, hi) {
  inr <- function(v) {
    if (!length(v)) return(v)
    nm <- names(v)
    p <- as.integer(substring(nm, 2L))
    v[ifelse(startsWith(nm, "p"), p >= lo & p <= hi, p > lo & p <= hi)]
  }
  a <- inr(a)
  b <- inr(b)
  keys <- union(names(a), names(b))
  if (!length(keys)) return(0L)
  sum(is.na(a[keys]) | is.na(b[keys]) | a[keys] != b[keys])
}

# Hits of 14_tie_rules_search.R-like searches (fields saccver, stitle, pident,
# length, nident, qstart, qend, qlen, qcovhsp, bitscore, staxid, btop) with
# the query id `qid` of each row, the taxonomy table `tax` and the true
# species of each query (`true_species[qid]`): ordered in BLAST order, with
# taxonomy, informative flags, voucher codes and the positions of the
# differences.
prep_hits <- function(s, qid, tax, true_species) {
  s$qid <- qid
  s <- s[order(s$qid, s$hit), ]
  s$blast_rank <- as.integer(stats::ave(seq_len(nrow(s)), s$qid, FUN = seq_along))
  num <- c("pident", "length", "nident", "qstart", "qend", "qlen", "qcovhsp", "bitscore")
  s[num] <- lapply(s[num], as.numeric)
  s$staxid <- trimws(sub(";.*", "", as.character(s$staxid)))
  k <- match(s$staxid, tax$taxid)
  s$species <- tax$species[k]
  s$genus <- tax$genus[k]
  s$named <- !is.na(s$species) & !grepl(open_nomenclature_pattern, s$species, perl = TRUE)
  s$open <- !is.na(s$species) & !s$named & !is.na(s$genus)
  s$d_aln <- s$length - s$nident
  s$d_cov <- s$d_aln + (s$qlen - (s$qend - s$qstart + 1))
  s$uninformative_category <- uninformative_category(s$saccver, s$stitle, uninformative_references())
  s$informative <- is.na(s$uninformative_category)
  s$true_species <- true_species[s$qid]
  s$voucher <- parse_voucher(s$stitle)
  s$diffs <- btop_diffs(s$btop, s$qstart)
  s$states <- btop_states(s$btop, s$qstart)
  chk <- lengths(s$diffs) == s$d_aln
  message("btop differences = length - nident in ", round(100 * mean(chk), 2), "% of the hits")
  s
}

# rules ---------------------------------------------------------------------
# tied set of one query (rows of `h` already ranked by rank_blast_hits)
tie_sets <- function(h) {
  pool <- !is.na(h$rank)
  qmin <- function(d) stats::ave(ifelse(pool, d, Inf), h$query, FUN = min)
  out <- list(
    current = h$at_best,
    aln_k0 = pool & h$d_aln <= qmin(h$d_aln),
    aln_k1 = pool & h$d_aln <= qmin(h$d_aln) + 1,
    cov_k0 = pool & h$d_cov <= qmin(h$d_cov),
    cov_k1 = pool & h$d_cov <= qmin(h$d_cov) + 1
  )
  # common region with the first hit
  common <- function(p) {
    res <- logical(nrow(h))
    for (rows in split(seq_len(nrow(h)), h$query)) {
      rows <- rows[pool[rows]]
      if (!length(rows)) next
      a <- rows[1L] # first by the Similarity score
      for (r in rows) {
        lo <- max(h$qstart[a], h$qstart[r])
        hi <- min(h$qend[a], h$qend[r])
        L <- hi - lo + 1
        if (L < 0.5 * h$qlen[a]) {
          res[r] <- TRUE
          next
        }
        res[r] <- n_in(h$diffs[[r]], lo, hi) <= n_in(h$diffs[[a]], lo, hi) + floor(p * L)
      }
    }
    res
  }
  out$common_k0 <- common(0)
  out$common_prop <- common(0.005)
  out
}

# open nomenclature votes only when strictly closer than every named species
# tied with it; `closer(i, j)`: hit i closer to the query than hit j
sp_rule <- function(h, tied, closer) {
  for (rows in split(seq_len(nrow(h)), h$query)) {
    t <- rows[tied[rows]]
    nm <- t[h$named[t]]
    op <- t[h$open[t]]
    if (!length(nm) || !length(op)) next
    for (o in op) if (!all(vapply(nm, function(n) closer(o, n), logical(1L)))) tied[o] <- FALSE
  }
  tied
}
# open nomenclature with the voucher code of a named hit, at the same distance
voucher_rule <- function(h, tied, same) {
  hits <- 0L
  for (rows in split(seq_len(nrow(h)), h$query)) {
    t <- rows[tied[rows]]
    nm <- t[h$named[t] & !is.na(h$voucher[t])]
    op <- t[h$open[t] & !is.na(h$voucher[t])]
    for (o in op) {
      m <- nm[h$voucher[nm] == h$voucher[o]]
      if (length(m) && any(vapply(m, function(n) same(o, n), logical(1L)))) {
        tied[o] <- FALSE
        hits <- hits + 1L
      }
    }
  }
  attr(tied, "dropped") <- hits
  tied
}

key_of <- function(h) round(similarity_score(h$pident, h$qcovhsp), 6)
common_cmp <- function(h, i, j) {
  lo <- max(h$qstart[i], h$qstart[j])
  hi <- min(h$qend[i], h$qend[j])
  c(n_in(h$diffs[[i]], lo, hi), n_in(h$diffs[[j]], lo, hi))
}

# OBITools obitag rule (obitools4.metabarcoding.org/obitools/obitag): the
# best hits BH (fewest differences to the query over the whole query, d) and
# every hit at least as similar to a BH as the query is: differences to the
# BH on the region both cover / its length <= d / query length. A hit whose
# common region with the BH is under half the query cannot be compared and
# is kept. With `gap`, a species is called only if the query is not farther
# from the BH than the farthest reference of that species is (its spread in
# the hits); a species with one reference (or identical ones) is called only
# on an identical match; otherwise the call goes to the genus. Returns the
# tied set, with the queries to cap at the genus in attribute "cap".
# Weight of the query bases a hit does not cover (added after the real runs,
# hub 2026-11-25-0218, -0222, -0244): the distance of a hit to the query is
# (differences + w * uncovered) / (aligned query bases + w * uncovered); w = 1
# is the whole distance (as delivered), w = 0 the rate on the aligned region.
# `best_w`: weight to choose the best hits ("S": the first by the Similarity
# score, the rank of BLASTr); `window_w`: weight of the window.
weighted_distance <- function(h, rows, w) {
  span <- h$qend[rows] - h$qstart[rows] + 1
  unc <- h$qlen[rows] - span
  (h$d_aln[rows] + w * unc) / (span + w * unc)
}
obitag_sets <- function(h, gap = FALSE, floor = 0, floor_identical = FALSE, best_w = 1, window_w = 1) {
  pool <- !is.na(h$rank)
  res <- logical(nrow(h))
  cap <- integer(0L)
  flags <- data.frame(query = integer(0L), flag = character(0L))
  for (rows in split(seq_len(nrow(h)), h$query)) {
    rows <- rows[pool[rows]]
    if (!length(rows)) next
    if (identical(best_w, "S")) {
      bh <- rows[h$rank[rows] == min(h$rank[rows])]
    } else {
      dq <- weighted_distance(h, rows, best_w)
      bh <- rows[dq <= min(dq) + 1e-12]
    }
    t0 <- min(weighted_distance(h, bh, window_w))
    # window: the query's distance to the best hit; with `floor`, never under
    # it unless the match is identical (with `floor_identical`, not even then)
    t <- if (t0 > 0 || floor_identical) max(t0, floor) else 0
    sig <- vapply(h$states[bh], function(v) paste(names(v), v, collapse = ";"), character(1L))
    bh_u <- bh[!duplicated(sig)]
    inc <- rows %in% bh
    for (b in bh_u) {
      for (j in which(!inc)) {
        r <- rows[j]
        lo <- max(h$qstart[b], h$qstart[r])
        hi <- min(h$qend[b], h$qend[r])
        L <- hi - lo + 1
        inc[j] <- L < 0.5 * h$qlen[b] || pair_diff(h$states[[b]], h$states[[r]], lo, hi) / L <= t + 1e-9
      }
    }
    res[rows[inc]] <- TRUE
    if (t0 > 0) {
      sp <- unique(h$species[bh])
      if (length(sp) == 1L && !is.na(sp) && isTRUE(h$named[bh[1L]])) {
        same <- rows[h$species[rows] %in% sp]
        spread <- 0
        for (b in bh_u) for (r in same) {
          lo <- max(h$qstart[b], h$qstart[r])
          hi <- min(h$qend[b], h$qend[r])
          L <- hi - lo + 1
          if (L >= 0.5 * h$qlen[b]) spread <- max(spread, pair_diff(h$states[[b]], h$states[[r]], lo, hi) / L)
        }
        # two signals: one reference (or identical ones) only; or the query
        # farther from the best hit than the references of the species are
        if (spread == 0) {
          flags <- rbind(flags, data.frame(query = h$query[rows[1L]], flag = "single_reference"))
        } else if (t0 > spread + 1e-9) {
          flags <- rbind(flags, data.frame(query = h$query[rows[1L]], flag = "beyond_spread"))
        }
        if (gap && t0 > spread + 1e-9) cap <- c(cap, h$query[rows[1L]])
      }
    }
  }
  attr(res, "cap") <- cap
  attr(res, "flags") <- flags
  res
}

all_rules <- function(h, which = NULL) {
  if (!is.null(which)) return(some_rules(h, which))
  t <- tie_sets(h)
  key <- key_of(h)
  t$`current+sp` <- sp_rule(h, t$current, function(i, j) key[i] > key[j])
  t$`aln_k0+sp` <- sp_rule(h, t$aln_k0, function(i, j) h$d_aln[i] < h$d_aln[j])
  t$`cov_k0+sp` <- sp_rule(h, t$cov_k0, function(i, j) h$d_cov[i] < h$d_cov[j])
  t$`common_k0+sp` <- sp_rule(h, t$common_k0, function(i, j) {
    d <- common_cmp(h, i, j)
    d[1L] < d[2L]
  })
  t$`common_prop+sp` <- sp_rule(h, t$common_prop, function(i, j) {
    d <- common_cmp(h, i, j)
    d[1L] < d[2L]
  })
  t$`current+voucher` <- voucher_rule(h, t$current, function(i, j) key[i] == key[j])
  t$`common_prop+voucher` <- voucher_rule(h, t$common_prop, function(i, j) {
    d <- common_cmp(h, i, j)
    d[1L] == d[2L]
  })
  t$obitag <- obitag_sets(h)
  t$`obitag+gap` <- obitag_sets(h, gap = TRUE)
  t$`current+gap` <- structure(t$current, cap = attr(t$`obitag+gap`, "cap"))
  t
}

# Only the rules in `which` (faster for large benchmarks).
some_rules <- function(h, which) {
  key <- key_of(h)
  pool <- !is.na(h$rank)
  qmin <- function(d) stats::ave(ifelse(pool, d, Inf), h$query, FUN = min)
  t <- list()
  if ("current" %in% which || "current+sp" %in% which) t$current <- h$at_best
  if ("current+sp" %in% which) t$`current+sp` <- sp_rule(h, t$current, function(i, j) key[i] > key[j])
  if ("cov_k0" %in% which) t$cov_k0 <- pool & h$d_cov <= qmin(h$d_cov)
  if ("common_prop" %in% which) t$common_prop <- tie_sets(h)$common_prop
  if ("obitag" %in% which) t$obitag <- obitag_sets(h)
  if (any(c("obitag_floor", "obitag_floor+gap", "obitag_floor+sp") %in% which)) {
    of <- obitag_sets(h, gap = TRUE, floor = 0.005)
    t$obitag_floor <- structure(as.logical(of), flags = attr(of, "flags"))
    if ("obitag_floor+gap" %in% which) t$`obitag_floor+gap` <- structure(as.logical(of), cap = attr(of, "cap"))
    if ("obitag_floor+sp" %in% which) {
      t$`obitag_floor+sp` <- sp_rule(h, as.logical(of), function(i, j) h$d_cov[i] < h$d_cov[j])
    }
  }
  if (any(c("obitag_floor_ident", "obitag_floor_ident+gap") %in% which)) {
    oi <- obitag_sets(h, gap = TRUE, floor = 0.005, floor_identical = TRUE)
    t$obitag_floor_ident <- structure(as.logical(oi), flags = attr(oi, "flags"))
    if ("obitag_floor_ident+gap" %in% which) t$`obitag_floor_ident+gap` <- structure(as.logical(oi), cap = attr(oi, "cap"))
  }
  # refs_b<best weight>_w<window weight>: obitag_floor_ident with the weights
  # of the uncovered bases (1, 11 = 1/11 or S = Similarity score, 0)
  wt <- c(`1` = 1, `11` = 1 / 11, `0` = 0)
  for (r in grep("^refs_b(1|11|S|0)_w(1|11|0)$", which, value = TRUE)) {
    bw <- sub("^refs_b([^_]+)_.*$", "\\1", r)
    ww <- sub("^.*_w", "", r)
    o <- obitag_sets(h, floor = 0.005, floor_identical = TRUE,
                     best_w = if (bw == "S") "S" else wt[[bw]], window_w = wt[[ww]])
    t[[r]] <- structure(as.logical(o), flags = attr(o, "flags"))
  }
  t[intersect(which, names(t))]
}

# evaluation ----------------------------------------------------------------
# Ancestor of blast_lca() (consensus 0.9, named_first) for one tied set.
evaluate <- function(h, tied, rule, tax, tax_db, cap = integer(0L)) {
  x <- tibble::tibble(
    Sequence = as.character(h$query), staxid = h$staxid, at_best = tied,
    informative = h$informative, rank = h$rank, blast_rank = h$blast_rank
  )
  l <- suppressWarnings(blast_lca(x, taxonomy = tax))
  lt <- suppressWarnings(get_tax(unique(stats::na.omit(l$`LCA ranked taxid`)), taxonomy_db = tax_db))
  m <- match(l$`LCA ranked taxid`, lt$taxid)
  out <- data.frame(
    qid = as.integer(l$Sequence), rule = rule, rank = l$`LCA ranked rank`,
    species = ifelse(l$`LCA ranked rank` %in% "species", lt$species[m], NA_character_),
    genus = ifelse(l$`LCA ranked rank` %in% c("species", "genus"), lt$genus[m], NA_character_),
    n_tied = l$`N tied`
  )
  # Similarity score of the best tied informative hit (the score mtbrcdR cuts
  # each rule by since 657887b, hub 2026-11-25-0230)
  sc <- similarity_score(h$pident, h$qcovhsp)
  k <- tied %in% TRUE & h$informative %in% TRUE
  st <- tapply(sc[k], h$query[k], max)
  out$score_tied <- as.numeric(st[as.character(out$qid)])
  # species calls capped at the genus
  k <- out$qid %in% cap & out$rank %in% "species"
  out$rank[k] <- "genus"
  out$species[k] <- NA_character_
  out
}

# Hits as parallel_blast() keeps them by default: identity >= 80, the first
# 100 subjects in BLAST order.
default_setting <- function(h) {
  h <- h[h$pident >= 80, ]
  h$blast_rank <- as.integer(stats::ave(seq_len(nrow(h)), h$qid, FUN = seq_along))
  h[h$blast_rank <= 100L, ]
}

# Ranks the hits as BLASTr does and evaluates every rule. `limit`: number of
# hits kept, to flag tied sets reaching the last one.
run_rules <- function(h, tax, tax_db, which = NULL, limit = 100L) {
  h$query <- h$qid
  h$identity <- h$pident
  h <- rank_blast_hits(h, rank_by = "similarity", tolerance = 0.5)
  rules <- all_rules(h, which)
  dropped <- vapply(rules, function(t) if (is.null(attr(t, "dropped"))) NA_integer_ else attr(t, "dropped"), integer(1L))
  res <- do.call(rbind, lapply(names(rules), function(r) {
    evaluate(h, as.logical(rules[[r]]), r, tax, tax_db, cap = attr(rules[[r]], "cap") %||% integer(0L))
  }))
  attr(res, "voucher_dropped") <- dropped[!is.na(dropped)]
  # tied set reaching the last hit kept (the list may cut it, RE-02)
  last <- stats::ave(h$blast_rank, h$query, FUN = max)
  full <- last >= limit
  res$truncated <- FALSE
  for (r in names(rules)) {
    tq <- unique(h$query[as.logical(rules[[r]]) & full & h$blast_rank == last])
    res$truncated[res$rule == r & res$qid %in% tq] <- TRUE
  }
  fl <- lapply(names(rules), function(r) {
    f <- attr(rules[[r]], "flags")
    if (is.null(f) || !nrow(f)) return(NULL)
    cbind(rule = r, f)
  })
  attr(res, "flags") <- do.call(rbind, fl)
  res
}

rule_order <- c(
  "current", "current+sp", "current+voucher", "aln_k0", "aln_k0+sp", "aln_k1", "cov_k0", "cov_k0+sp",
  "cov_k1", "common_k0", "common_k0+sp", "common_prop", "common_prop+sp", "common_prop+voucher",
  "obitag", "obitag+gap", "current+gap", "obitag_floor", "obitag_floor+sp", "obitag_floor+gap",
  "obitag_floor_ident", "obitag_floor_ident+gap",
  as.vector(outer(c("refs_b1", "refs_b11", "refs_bS", "refs_b0"), c("_w1", "_w11", "_w0"), paste0))
)

# % of the queries with the right / wrong species and genus, by `sets`
# (columns of `res` giving the groups; "all" is added).
summarise_rules <- function(res, sets, by = c("scenario", "setting")) {
  long <- dplyr::bind_rows(c(
    lapply(sets, function(k) dplyr::mutate(res, set = as.character(.data[[k]]))),
    list(dplyr::mutate(res, set = "all"))
  ))
  long |>
    dplyr::mutate(
      at_sp = rank %in% "species",
      at_gen = rank %in% c("species", "genus"),
      sp_ok = at_sp & (species == true_species) %in% TRUE,
      gen_ok = at_gen & (genus == true_genus) %in% TRUE
    ) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(by, "rule", "set")))) |>
    dplyr::summarise(
      n = dplyr::n(),
      species_right = round(100 * mean(sp_ok), 1),
      species_wrong = round(100 * mean(at_sp & !sp_ok), 1),
      genus_right = round(100 * mean(gen_ok), 1),
      genus_wrong = round(100 * mean(at_gen & !gen_ok), 1),
      median_tied = stats::median(n_tied, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::across(dplyr::all_of(c(by, "set"))), match(rule, rule_order)) |>
    as.data.frame()
}

# primer benchmark (17_ to 21_) ---------------------------------------------
# One truth per distinct sequence (sequences shared by two species left out);
# `qq`: one row per primer and query.
bench_truth <- function(q) {
  truth <- q |>
    dplyr::group_by(seq) |>
    dplyr::summarise(
      n_sp = dplyr::n_distinct(species), species = dplyr::first(species), genus = dplyr::first(genus),
      acc = dplyr::first(acc), .groups = "drop"
    ) |>
    dplyr::filter(n_sp == 1L)
  q <- q[q$seq %in% truth$seq, ]
  qq <- q[, c("primer", "marker", "group", "seq")]
  qq <- qq[!duplicated(paste(qq$primer, qq$seq)), ]
  list(truth = truth, q = q, qq = qq)
}

# Hits of the primer benchmark ready for the rules.
bench_hits <- function(out_dir, truth, tax_db) {
  s <- readRDS(file.path(out_dir, "hits.rds"))
  s <- s[s$Sequence %in% truth$seq, ]
  s$staxid <- trimws(sub(";.*", "", as.character(s$staxid)))
  tax <- suppressWarnings(get_tax(unique(s$staxid), taxonomy_db = tax_db))
  list(s = prep_hits(s, match(s$Sequence, truth$seq), tax, truth$species), tax = tax)
}

# Similarity score of the best informative hit of each query (the score the
# mtbrcdR cuts by rank).
best_score <- function(h) {
  h <- h[h$informative, ]
  sc <- tapply(similarity_score(h$pident, h$qcovhsp), h$qid, max)
  data.frame(qid = as.integer(names(sc)), score = as.numeric(sc))
}

# The rules in `which` with the species in the database ("in") and with every
# hit of the species removed ("out"); `setting` picks the hits kept.
run_scenarios <- function(s, tax, tax_db, setting, limit, which) {
  runs <- list()
  flags <- list()
  for (sc in c("in", "out")) {
    h <- s
    if (sc == "out") h <- h[!(h$species %in% h$true_species), ]
    h <- setting(h)
    r <- run_rules(h, tax, tax_db, which = which, limit = limit)
    r$scenario <- sc
    if (!is.null(attr(r, "flags"))) flags[[sc]] <- cbind(scenario = sc, attr(r, "flags"))
    r <- merge(r, best_score(h), by = "qid", all.x = TRUE)
    runs[[sc]] <- r
    message(format(Sys.time()), " ", sc, " done")
  }
  list(runs = runs, flags = do.call(rbind, flags))
}

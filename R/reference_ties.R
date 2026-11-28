# Tie rule by references, after the one of OBITools' obitag: the hits tied
# with the best one are the references as close to the best reference as the
# query is, compared with each other on the region they share. It needs the
# `btop` of each hit (what each reference has where it differs from the
# query). See `at_best_refs` in parallel_blast() and `ties` in blast_lca().

# Smallest window of the rule, as a fraction of the query length: with the
# query identical or closer than this to the best reference, references up to
# this far from the best one are tied with it (0.5%, about the 0.5 point of
# Similarity score of the default tie rule).
refs_floor_default <- 0.005

# Two references are compared only on the query region both cover; when it is
# shorter than this fraction of the query, the reference is tied (it cannot
# be told apart from the best one).
refs_min_common <- 0.5

# What each hit has where it differs from the query, from its btop and query
# start: a character vector named "p<pos>" (the subject base, or "-", at that
# query position) and "i<pos>" (bases the subject inserts before that query
# position), with the attributes `pos` and `ins`. Two hits that both match the
# query at a position agree there. NULL when the btop is missing.
btop_states <- function(btop, qstart) {
  out <- vector("list", length(btop))
  ok <- !is.na(btop) & nzchar(btop) & !is.na(qstart)
  if (!any(ok)) return(out)
  tok <- regmatches(btop[ok], gregexpr("[0-9]+|[A-Za-z*-]{2}", btop[ok]))
  out[ok] <- mapply(function(t, q0) {
    isnum <- grepl("^[0-9]", t)
    n <- suppressWarnings(as.integer(ifelse(isnum, t, "0")))
    qgap <- !isnum & substr(t, 1L, 1L) == "-"
    adv <- ifelse(isnum, n, ifelse(qgap, 0L, 1L))
    pos <- q0 + c(0L, cumsum(adv))[seq_along(t)]
    ch <- toupper(substr(t, 2L, 2L))
    k <- !isnum & !qgap
    # paste0() of an empty vector would give one name
    v <- if (any(k)) stats::setNames(ch[k], paste0("p", pos[k])) else stats::setNames(character(0L), character(0L))
    p <- as.integer(pos[k])
    ins <- rep(FALSE, sum(k))
    if (any(qgap)) {
      gi <- tapply(ch[qgap], pos[qgap], paste, collapse = "")
      v <- c(v, stats::setNames(as.character(gi), paste0("i", names(gi))))
      p <- c(p, as.integer(names(gi)))
      ins <- c(ins, rep(TRUE, length(gi)))
    }
    structure(v, pos = p, ins = ins)
  }, tok, as.integer(qstart[ok]), SIMPLIFY = FALSE, USE.NAMES = FALSE)
  out
}

# Differences between two hits on the query region lo..hi, from their states.
state_diff <- function(a, b, lo, hi) {
  inr <- function(v) {
    if (!length(v)) return(v)
    p <- attr(v, "pos")
    v[ifelse(attr(v, "ins"), p > lo & p <= hi, p >= lo & p <= hi)]
  }
  a <- inr(a)
  b <- inr(b)
  keys <- union(names(a), names(b))
  if (!length(keys)) return(0L)
  sum(is.na(a[keys]) | is.na(b[keys]) | a[keys] != b[keys])
}

# Distance of each hit to the query: differences in the alignment plus query
# bases outside it, over the query length.
query_distance <- function(length, nident, qstart, qend, qlen) {
  ((length - nident) + (qlen - (qend - qstart + 1))) / qlen
}

# The tie rule by references for the hits of several queries. `pool`: the
# hits that compete (the ranked ones). Returns a list: `tied` (logical; NA for
# the queries with a hit of the pool without btop), `best` (the hits closest
# to the query) and `distance` (to the query).
reference_ties <- function(query, pool, btop, qstart, qend, length, nident, qlen,
                           floor = refs_floor_default, min_common = refs_min_common) {
  n <- base::length(query)
  tied <- rep(FALSE, n)
  best <- rep(FALSE, n)
  d <- query_distance(as.numeric(length), as.numeric(nident), as.numeric(qstart), as.numeric(qend), as.numeric(qlen))
  pool <- pool %in% TRUE & !is.na(d)
  for (rows in split(seq_len(n), factor(query, levels = unique(query)))) {
    rows <- rows[pool[rows]]
    if (!base::length(rows)) next
    if (anyNA(btop[rows]) || !all(nzchar(btop[rows]))) {
      tied[rows] <- NA
      next
    }
    st <- btop_states(btop[rows], qstart[rows])
    t0 <- min(d[rows])
    t <- max(t0, floor)
    b <- d[rows] <= t0 + 1e-12
    best[rows[b]] <- TRUE
    sig <- vapply(st[b], function(v) paste(names(v), v, collapse = ";"), character(1L))
    bu <- which(b)[!duplicated(sig)]
    inc <- b
    for (i in bu) {
      for (j in which(!inc)) {
        lo <- max(qstart[rows[i]], qstart[rows[j]])
        hi <- min(qend[rows[i]], qend[rows[j]])
        L <- hi - lo + 1
        inc[j] <- L < min_common * qlen[rows[i]] || state_diff(st[[i]], st[[j]], lo, hi) / L <= t + 1e-9
      }
    }
    tied[rows] <- inc
  }
  list(tied = tied, best = best, distance = d)
}

# Signal on the species of the best reference, for the tie rule by
# references: when the query differs from the best reference, and the best
# references are of one named species, `"single reference"` when that species
# has only references identical to the best one, `"beyond spread"` when the
# query is farther from the best reference than the references of the species
# are from each other; NA otherwise. Returns a data frame with `query`,
# `Species flag` and `Flagged species`.
reference_species_flag <- function(query, pool, btop, qstart, qend, length, nident, qlen, species, named,
                                   min_common = refs_min_common) {
  qs <- unique(query)
  flag <- rep(NA_character_, base::length(qs))
  fsp <- rep(NA_character_, base::length(qs))
  d <- query_distance(as.numeric(length), as.numeric(nident), as.numeric(qstart), as.numeric(qend), as.numeric(qlen))
  pool <- pool %in% TRUE & !is.na(d)
  by_q <- split(seq_along(query), factor(query, levels = qs))
  for (g in seq_along(qs)) {
    rows <- by_q[[g]][pool[by_q[[g]]]]
    if (!base::length(rows) || anyNA(btop[rows])) next
    t0 <- min(d[rows])
    if (t0 <= 0) next
    b <- rows[d[rows] <= t0 + 1e-12]
    sp <- unique(species[b])
    if (base::length(sp) != 1L || is.na(sp) || !isTRUE(named[b[1L]])) next
    same <- rows[species[rows] %in% sp]
    st <- btop_states(btop[c(b, same)], qstart[c(b, same)])
    sb <- st[seq_along(b)]
    ss <- st[base::length(b) + seq_along(same)]
    sig <- vapply(sb, function(v) paste(names(v), v, collapse = ";"), character(1L))
    bu <- which(!duplicated(sig))
    spread <- 0
    for (i in bu) {
      for (j in seq_along(same)) {
        lo <- max(qstart[b[i]], qstart[same[j]])
        hi <- min(qend[b[i]], qend[same[j]])
        L <- hi - lo + 1
        if (L >= min_common * qlen[b[i]]) spread <- max(spread, state_diff(sb[[i]], ss[[j]], lo, hi) / L)
      }
    }
    if (spread == 0) {
      flag[g] <- "single reference"
    } else if (t0 > spread + 1e-9) {
      flag[g] <- "beyond spread"
    }
    if (!is.na(flag[g])) fsp[g] <- sp
  }
  data.frame(query = qs, `Species flag` = flag, `Flagged species` = fsp, check.names = FALSE)
}

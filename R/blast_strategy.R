#' @title Plan the BLAST search strategy
#'
#' @description
#' Decide how [parallel_blast()] splits a search: how many `blastn`
#' processes run at once, how many threads each one uses, the BLAST+
#' multithreading mode and how many query sequences go in each batch.
#'
#' Each `blastn` process reads the whole database once per batch, so the
#' cost of a search grows with the number of batches and of concurrent
#' processes, not with the number of sequences. The rules below come from
#' the BLASTr benchmark (`dev/benchmark/` in the source repository):
#'
#' * **Small databases** (under `small_db_gb`): reading the database is
#'   cheap. One single-threaded process per core, each with an equal share
#'   of the sequences; results are identical to one search per sequence.
#' * **Large databases that fit in RAM**: few processes with many threads
#'   each, splitting the database between the threads (`mt_mode = 2`).
#'   Threads per process are limited to the number of database volumes and
#'   to `max_threads` (32); the remaining cores run more processes, each
#'   with its own share of the sequences. On `core_nt` (92 volumes), with
#'   32 cores one `blastn` searched 500 sequences in 170 s, against 164 s
#'   for 2 processes x 16 threads and 243 s for 8 x 4 (2000 sequences: 239 s,
#'   311 s, 367 s). With 64 cores, one `blastn` with 64 threads was slower
#'   than with 32 (152 s vs 120 s for 500 sequences, 265 s vs 206 s for
#'   2000), while 2 processes x 32 threads took 87 s and 182 s.
#' * **Large databases that do not fit in RAM**: one process and one batch,
#'   so the database is read from disk only once (on `core_nt`, each extra
#'   read cost about 6 minutes). The ETA is then unknown until the end.
#'
#' Otherwise, batches are kept at or under `max_batch_size` sequences and
#' there is at least one batch per process. Each batch reads the whole
#' database again (about 70 s on `core_nt` in RAM, with 32 threads; 2000
#' sequences in 80 batches took 1190 s instead of 345 s in 4), so batches
#' are only as many as needed for the ETA to show up during long searches.
#'
#' @param n_seqs Number of unique query sequences.
#' @param db_path Path(s) to the BLAST database(s). Several databases are
#'   searched together; their sizes are added.
#' @param total_cores Number of CPU cores (threads) to use in total.
#' @param ram_gb RAM available for the search, in GB. `NULL` uses the total
#'   RAM of the machine.
#' @param num_threads,batch_size,mt_mode Optional fixed values; `NULL`
#'   chooses automatically.
#' @param small_db_gb Size under which a database is treated as small.
#' @param max_threads Largest number of threads for one `blastn` process on
#'   a large database; more threads than this did not make BLAST faster.
#' @param max_batch_size Largest number of sequences in a single batch.
#'
#' @returns A list of class `blastr_strategy` with `processes`, `threads`,
#'   `mt_mode`, `batch_size`, `n_batches`, the database facts (`db_path`,
#'   `db_gb`, `db_volumes`, `ram_gb`, `fits_ram`) and `reason`.
#'
#' @examples
#' # needs BLAST+ (see install_dependencies())
#' if (has_blast()) {
#'   fasta <- system.file("extdata", "minimal_db_blast.fasta", package = "BLASTr")
#'   db <- file.path(tempdir(), "minimal_db_blast")
#'   make_blast_db(fasta, db_path = db)
#'   seqs <- readLines(system.file("extdata", "seqs_test.txt", package = "BLASTr"))
#'   # a small database: one single-threaded process per core
#'   blast_strategy(n_seqs = 2000, db_path = db, total_cores = 8)
#' }
#'
#' @family search functions
#' @export
blast_strategy <- function(
  n_seqs,
  db_path,
  total_cores = 1L,
  ram_gb = NULL,
  num_threads = NULL,
  batch_size = NULL,
  mt_mode = NULL,
  small_db_gb = 1,
  max_batch_size = 1000L,
  max_threads = 32L
) {
  total_cores <- max(1L, as.integer(total_cores))
  n_seqs <- max(1L, as.integer(n_seqs))
  db <- blast_db_info(db_path)
  if (is.null(ram_gb)) {
    ram_gb <- system_ram_gb()
  }
  fits_ram <- is.na(ram_gb) || db$gb <= 0.9 * ram_gb

  if (db$gb < small_db_gb) {
    threads <- 1L
    processes <- total_cores
    mode <- 0L
    reason <- "small database: one single-threaded process per core"
  } else if (isTRUE(fits_ram)) {
    threads <- max(1L, min(total_cores, db$volumes, max_threads))
    processes <- max(1L, total_cores %/% threads)
    mode <- 2L
    reason <- "large database in RAM: threads split the database (mt_mode 2)"
  } else {
    threads <- max(1L, min(total_cores, max_threads))
    processes <- 1L
    mode <- 2L
    reason <- paste(
      "database larger than the RAM: one process and one batch,",
      "so the disk is read only once (no ETA until the end)"
    )
  }

  if (!is.null(num_threads)) {
    threads <- max(1L, min(as.integer(num_threads), total_cores))
    processes <- max(1L, total_cores %/% threads)
  }
  if (!is.null(mt_mode)) {
    mode <- as.integer(mt_mode)
  }
  processes <- min(processes, n_seqs)

  if (is.null(batch_size) && !isTRUE(fits_ram) && db$gb >= small_db_gb) {
    batch_size <- n_seqs
  }
  plan <- batch_plan(n_seqs, processes, batch_size, max_batch_size)
  batch_size <- plan$batch_size
  n_batches <- plan$n_batches
  processes <- min(processes, n_batches)

  structure(
    list(
      processes = as.integer(processes),
      threads = as.integer(threads),
      mt_mode = as.integer(mode),
      batch_size = batch_size,
      n_batches = n_batches,
      n_seqs = n_seqs,
      db_path = db$path,
      db_gb = db$gb,
      db_volumes = db$volumes,
      ram_gb = ram_gb,
      fits_ram = fits_ram,
      reason = reason
    ),
    class = "blastr_strategy"
  )
}

#' @export
print.blastr_strategy <- function(x, ...) {
  if (is.na(x$threads)) {
    cli::cli_text(
      "{x$n_seqs} sequences in {x$n_batches} batch{?es} of up to ",
      "{x$batch_size}; {x$processes} single-threaded process{?es}."
    )
  } else {
    cli::cli_text(
      "{x$n_seqs} sequences in {x$n_batches} batch{?es} of up to ",
      "{x$batch_size}; {x$processes} process{?es} x {x$threads} thread{?s} ",
      "(mt_mode {x$mt_mode})."
    )
  }
  cli::cli_text(
    "Database: {format(round(x$db_gb, 2), nsmall = 2)} GB in ",
    "{x$db_volumes} volume{?s}; RAM: ",
    "{if (is.na(x$ram_gb)) 'unknown' else paste(round(x$ram_gb), 'GB')}. ",
    "{x$reason}."
  )
  invisible(x)
}

# Paths, total size (GB) and number of volumes of one or more BLAST
# databases, smallest first. A single string with spaces that is not itself
# a database is split, as in `-db "db1 db2"`.
blast_db_info <- function(db_path) {
  db_path <- as.character(db_path)
  if (
    length(db_path) == 1L &&
      grepl(" ", db_path, fixed = TRUE) &&
      !db_files_exist(db_path)
  ) {
    db_path <- strsplit(trimws(db_path), "\\s+")[[1L]]
  }
  info <- lapply(db_path, function(p) {
    files <- db_files(p)
    data.frame(
      path = p,
      gb = sum(file.size(files), na.rm = TRUE) / 1e9,
      volumes = max(1L, length(grep("\\.[np]sq$", files)))
    )
  })
  info <- do.call(rbind, info)
  info <- info[order(info$gb), , drop = FALSE]
  list(
    path = info$path,
    gb = sum(info$gb),
    volumes = sum(info$volumes)
  )
}

db_files <- function(db_path) {
  dir <- dirname(db_path)
  if (!dir.exists(dir)) {
    return(character(0L))
  }
  base <- basename(db_path)
  hits <- list.files(dir, full.names = TRUE)
  hits[startsWith(basename(hits), paste0(base, "."))]
}

db_files_exist <- function(db_path) {
  any(grepl(
    "\\.(nal|pal|nin|pin|nsq|psq|ndb|pdb)$",
    db_files(db_path)
  ))
}

# The value given to `-db`: databases smallest first, quoted when needed.
db_arg <- function(db_paths) {
  q <- ifelse(grepl(" ", db_paths, fixed = TRUE), paste0("\"", db_paths, "\""), db_paths)
  paste(q, collapse = " ")
}

system_ram_gb <- function() {
  ram <- tryCatch(ps::ps_system_memory()$total, error = function(e) NA_real_)
  ram / 1024^3
}

# Batches: at least one per process, a multiple of the processes, none
# larger than max_batch_size; or the batch size given.
batch_plan <- function(n_seqs, processes, batch_size = NULL, max_batch_size = 1000L) {
  if (is.null(batch_size)) {
    n_batches <- max(processes, ceiling(n_seqs / max_batch_size))
    n_batches <- processes * ceiling(n_batches / processes)
    batch_size <- ceiling(n_seqs / n_batches)
  }
  batch_size <- max(1L, as.integer(batch_size))
  list(batch_size = batch_size, n_batches = as.integer(ceiling(n_seqs / batch_size)))
}

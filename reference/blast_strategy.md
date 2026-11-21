# Plan the BLAST search strategy

Decide how
[`parallel_blast()`](https://heronoh.github.io/BLASTr/reference/parallel_blast.md)
splits a search: how many `blastn` processes run at once, how many
threads each one uses, the BLAST+ multithreading mode and how many query
sequences go in each batch.

Each `blastn` process reads the whole database once per batch, so the
cost of a search grows with the number of batches and of concurrent
processes, not with the number of sequences. The rules below come from
the BLASTr benchmark (`dev/benchmark/` in the source repository):

- **Small databases** (under `small_db_gb`): reading the database is
  cheap. One single-threaded process per core, each with an equal share
  of the sequences; results are identical to one search per sequence.

- **Large databases that fit in RAM**: few processes with many threads
  each, splitting the database between the threads (`mt_mode = 2`).
  Threads per process are limited to the number of database volumes and
  to `max_threads` (32); the remaining cores run more processes, each
  with its own share of the sequences. On `core_nt` (92 volumes), with
  32 cores one `blastn` searched 500 sequences in 170 s, against 164 s
  for 2 processes x 16 threads and 243 s for 8 x 4 (2000 sequences: 239
  s, 311 s, 367 s). With 64 cores, one `blastn` with 64 threads was
  slower than with 32 (152 s vs 120 s for 500 sequences, 265 s vs 206 s
  for 2000), while 2 processes x 32 threads took 87 s and 182 s.

- **Large databases that do not fit in RAM**: one process and one batch,
  so the database is read from disk only once (on `core_nt`, each extra
  read cost about 6 minutes). The ETA is then unknown until the end.

Otherwise, batches are kept at or under `max_batch_size` sequences and
there is at least one batch per process. Each batch reads the whole
database again (about 70 s on `core_nt` in RAM, with 32 threads; 2000
sequences in 80 batches took 1190 s instead of 345 s in 4), so batches
are only as many as needed for the ETA to show up during long searches.

## Usage

``` r
blast_strategy(
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
)
```

## Arguments

- n_seqs:

  Number of unique query sequences.

- db_path:

  Path(s) to the BLAST database(s). Several databases are searched
  together; their sizes are added.

- total_cores:

  Number of CPU cores (threads) to use in total.

- ram_gb:

  RAM available for the search, in GB. `NULL` uses the total RAM of the
  machine.

- num_threads, batch_size, mt_mode:

  Optional fixed values; `NULL` chooses automatically.

- small_db_gb:

  Size under which a database is treated as small.

- max_batch_size:

  Largest number of sequences in a single batch.

- max_threads:

  Largest number of threads for one `blastn` process on a large
  database; more threads than this did not make BLAST faster.

## Value

A list of class `blastr_strategy` with `processes`, `threads`,
`mt_mode`, `batch_size`, `n_batches`, the database facts (`db_path`,
`db_gb`, `db_volumes`, `ram_gb`, `fits_ram`) and `reason`.

## Examples

``` r
if (FALSE) { # \dontrun{
blast_strategy(
  n_seqs = 2000,
  db_path = "/data/databases/core_nt/core_nt",
  total_cores = 32,
  ram_gb = 256
)
} # }
```

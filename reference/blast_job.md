# Background BLAST jobs

[`parallel_blast()`](https://heronoh.github.io/BLASTr/reference/parallel_blast.md)
with `background = TRUE` returns a job and the search goes on in a
separate process, so R stays free (and the search goes on even if the R
session ends). Use these functions to follow the job and get its
results:

- `blast_job()` opens a job from its directory, e.g. in a new R session.

- `blast_job_status()` shows the progress and the estimated time left.

- `blast_job_result()` returns the results like
  [`parallel_blast()`](https://heronoh.github.io/BLASTr/reference/parallel_blast.md),
  waiting for the job to finish if `wait = TRUE`.

- `blast_job_cancel()` stops the job.

- `blast_job_resume()` restarts a job whose scheduler was stopped (e.g.
  the machine rebooted); finished batches are kept.

## Usage

``` r
blast_job(job_dir)

blast_job_status(job)

blast_job_cancel(job)

blast_job_resume(job)

blast_job_result(job, wait = TRUE, progress = rlang::is_interactive())
```

## Arguments

- job_dir:

  Directory of an existing job.

- job:

  A job returned by
  [`parallel_blast()`](https://heronoh.github.io/BLASTr/reference/parallel_blast.md)
  or `blast_job()`, or the path to its directory.

- wait:

  Wait for the job to finish? If `FALSE` and the job is not finished,
  `blast_job_result()` returns `NULL`.

- progress:

  Show a progress bar while waiting.

## Value

`blast_job()` and `blast_job_resume()` return the job;
`blast_job_status()` a list of class `blastr_job_status`;
`blast_job_result()` a tibble (see
[`parallel_blast()`](https://heronoh.github.io/BLASTr/reference/parallel_blast.md));
`blast_job_cancel()` the job, invisibly.

## Examples

``` r
if (FALSE) { # \dontrun{
job <- parallel_blast(seqs, db_path, total_cores = 32, background = TRUE)
blast_job_status(job)
# ... later, or in another R session:
job <- blast_job(job$dir)
res <- blast_job_result(job)
} # }
```

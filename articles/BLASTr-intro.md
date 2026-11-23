# BLASTr-intro

``` r

library(BLASTr)
```

## Searching in the background

Large searches (thousands of sequences against `core_nt`) can take
hours. With `background = TRUE`,
[`parallel_blast()`](https://heronoh.github.io/BLASTr/reference/parallel_blast.md)
starts the search in a separate process and returns a job at once, so
the R session stays free. The job keeps its batches and results in a
directory (`job_dir`; by default a new directory in the BLASTr cache),
so it goes on even if the R session ends, and can be opened again from
another session.

``` r

job <- parallel_blast(
  query_seqs = seqs,
  db_path = "/data/databases/core_nt_2026-11-17/core_nt",
  total_cores = 32,
  background = TRUE,
  job_dir = "~/blast_jobs/run_01" # a fixed place, to find it again
)

blast_job_status(job) # batches done, progress and estimated time left
```

To get the results,
[`blast_job_result()`](https://heronoh.github.io/BLASTr/reference/blast_job.md)
waits for the end (with a progress bar in interactive sessions) and
returns the same table as
[`parallel_blast()`](https://heronoh.github.io/BLASTr/reference/parallel_blast.md);
with `wait = FALSE` it returns `NULL` while the job is not finished.

``` r

blast_results <- blast_job_result(job)

# later, or in another R session
job <- blast_job("~/blast_jobs/run_01")
blast_job_status(job)
blast_results <- blast_job_result(job, wait = FALSE)
```

[`blast_job_cancel()`](https://heronoh.github.io/BLASTr/reference/blast_job.md)
stops a job. If the process running it stopped (for example, the machine
was restarted),
[`blast_job_status()`](https://heronoh.github.io/BLASTr/reference/blast_job.md)
reports the job as stopped and
[`blast_job_resume()`](https://heronoh.github.io/BLASTr/reference/blast_job.md)
restarts it: the batches already finished are kept and only the others
are searched again.

``` r

blast_job_cancel(job)
blast_job_resume("~/blast_jobs/run_01")
```

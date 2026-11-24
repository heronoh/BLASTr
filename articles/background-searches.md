# Long searches in the background

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
  db_path = "~/blast_databases/core_nt_2026-11-17/core_nt",
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

## Session information

``` r

sessionInfo()
#> R version 4.6.2 (2026-10-29)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 26.04.1 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.32.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] BLASTr_0.1.7
#> 
#> loaded via a namespace (and not attached):
#>  [1] compiler_4.6.2  fastmap_1.2.0   cli_3.6.6       tools_4.6.2    
#>  [5] htmltools_0.5.9 otel_0.3.0      yaml_2.3.12     rmarkdown_2.32 
#>  [9] knitr_1.52      jsonlite_2.0.0  xfun_0.62       digest_0.6.39  
#> [13] lifecycle_1.0.5 rlang_1.3.0     evaluate_1.0.5
```

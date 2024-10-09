get_tax_by_taxID(
  organisms_taxIDs = 13662, parse_result = F
  # ,
  # parse_result = ,
  # total_cores =
  # total_cores =
)




# parallel ----
source("/home/heron/prjcts/BLASTr/R/parallel_get_tax.R")

parallel_get_tax(
  organisms_taxIDs = c(
    "2978354", "10113", 0000, 2
    # ,
    # "02","3", "123124155","1131"
  ),
  retry_times = 2, parse_result = TRUE

  # ,
  # parse_result = ,
  # total_cores =
)


# organisms_taxIDs <- c("2978354", "10113", 0000,2)
organisms_taxIDs <- c("2978354", "13662", 0000, 2)
organisms_taxIDs <- "13662"
env_name <- "entrez-env"
verbose <- TRUE
parse_result <- FALSE
total_cores <- 2
retry_times <- 2

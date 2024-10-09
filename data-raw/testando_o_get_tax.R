source("/home/heron/prjcts/BLASTr/R/get_tax_by_taxID.R")


library(BLASTr)

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
  retry_times = 2, parse_result = T

  # ,
  # parse_result = ,
  # total_cores =
)


# organisms_taxIDs <- c("2978354", "10113", 0000,2)
organisms_taxIDs <- c("2978354", "13662", 0000, 2)
organisms_taxIDs <- "13662"
env_name <- "entrez-env"
verbose <- T
parse_result <- F
total_cores <- 2
retry_times <- 2

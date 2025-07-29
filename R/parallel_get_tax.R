#' Retrieve Taxonomic Ranks for a List of NCBI Taxonomy Tax IDs in Parallel
#'
#' Retrieves taxonomy ranks for a list of NCBI Taxonomy Tax IDs using parallel processing. The function queries the NCBI database and can retry fetching taxonomy information multiple times if initial attempts fail.
#'
#' @param organisms_taxIDs A character vector of NCBI Taxonomy Tax IDs to retrieve taxonomy information for.
#' @param parse_result Logical indicating whether to parse the taxonomy information into a tibble (`TRUE`, default) or return the raw output as returned by `efetch` (`FALSE`).
#' @param total_cores Integer specifying the number of cores to use for parallel processing. Defaults to `1`.
#' @param retry_times Integer specifying the number of times to retry fetching taxonomy information if it fails. Defaults to `10`.
#' @param verbose Logical indicating whether to print verbose messages during the process. Default is `FALSE`.
#' @param env_name Character string specifying the name of the conda environment where `efetch` is installed. Default is `"blastr-entrez-env"`.
#'
#' @return A tibble containing the taxonomic ranks for the given Tax IDs.
#'
#' @examples
#' \dontrun{
#' # Retrieve taxonomy for multiple Tax IDs using 4 cores
#' tax_ids <- c("9606", "10090", "10116") # Human, Mouse, Rat
#' tax_info <- parallel_get_tax(tax_ids, total_cores = 4)
#'
#' # Retrieve unparsed taxonomy result
#' tax_info_unparsed <- parallel_get_tax(tax_ids, parse_result = FALSE)
#'
#' # Increase the number of retry attempts
#' tax_info <- parallel_get_tax(tax_ids, retry_times = 20)
#'
#' # Enable verbose output
#' tax_info <- parallel_get_tax(tax_ids, verbose = TRUE)
#' }
#' @export
parallel_get_tax <- function(
    organisms_taxIDs,
    parse_result = TRUE,
    total_cores = 1,
    retry_times = 10,
    verbose = FALSE,
    env_name = "blastr-entrez-env") {
  organisms_taxIDs <- unique(organisms_taxIDs)
  check_cmd("efetch", env_name = env_name)

  # Set up parallel plan
  # future::plan(future::multisession, workers = total_cores)
  if (total_cores > 1L) {
    future::plan(future::multisession(),
      workers = total_cores
    )
  }
  if (rlang::is_true(parse_result)) {
    results <- tibble::tibble(
      "Sci_name" = character(0L),
      "query_taxID" = character(0L),
      # "Division (NCBI)" = character(0L),
      "Superkingdom (NCBI)" = character(0L),
      "Kingdom (NCBI)" = character(0L),
      "Phylum (NCBI)" = character(0L),
      "Subphylum (NCBI)" = character(0L),
      "Class (NCBI)" = character(0L),
      "Subclass (NCBI)" = character(0L),
      "Order (NCBI)" = character(0L),
      "Suborder (NCBI)" = character(0L),
      "Family (NCBI)" = character(0L),
      "Subfamily (NCBI)" = character(0L),
      "Genus (NCBI)" = character(0L)
    )
  }


  if (rlang::is_false(parse_result)) {
    # create empty tibble for binding
    results <- tibble::tibble(
      "Rank" = character(0L),
      "ScientificName" = character(0L),
      "query_taxID" = character(0L),
      "Sci_name" = character(0L)
    )
  }


  res_taxid <- character(0L)
  retry_count <- 0L
  while (isTRUE(retry_count < retry_times) && isFALSE(all(organisms_taxIDs %in% res_taxid))) {
    message(paste0("retrying ", retry_count, " of ", retry_times))
    # message(organisms_taxIDs)
    if (total_cores > 1L) {
      results_temp <- furrr::future_map_dfr(
        .x = organisms_taxIDs,
        .f = get_tax_by_taxID,
        parse_result = parse_result,
        env_name = env_name,
        verbose = verbose,
        .options = furrr::furrr_options(seed = NULL)
      )
    } else {
      results_temp <- purrr::map_dfr(
        .x = organisms_taxIDs,
        .f = get_tax_by_taxID,
        parse_result = parse_result,
        env_name = env_name,
        verbose = verbose
      )
    }
    retry_count <- retry_count + 1L

    results <- dplyr::bind_rows(results, results_temp) |>
      dplyr::distinct()

    # taxids that were found
    res_taxid <- unique(results$query_taxID)
    # missing taxids
    organisms_taxIDs <- organisms_taxIDs[!(organisms_taxIDs %in% results$query_taxID)]
  }
  # message for problematic taxIDs
  if (length(organisms_taxIDs[!(organisms_taxIDs %in% results$query_taxID)]) != 0) {
    message(paste0(
      "The following taxIDs could not be retrieved even after ",
      retry_times,
      " attempts:\n",
      organisms_taxIDs[!(organisms_taxIDs %in% results$query_taxID
      )]
    ))
  }
  # Reset future plan to default
  if (total_cores > 1L) {
    future::plan(future::sequential)
  }

  return(results)
}

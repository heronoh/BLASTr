#' @title Get taxonomy ranks for a list of _NCBI Taxonomy Tax ID_
#'
#' @description Get taxonomy ranks for a list of _NCBI Taxonomy Tax ID_ in parallel.
#' @param organisms_taxIDs Vector of _NCBI Taxonomy Tax ID_ to retrieve taxonomy for.
#' @param parse_result Should the taxonomy be returned as the _efetch_ returns it or should it be parsed into a tibble.
#' @param total_cores Number of threads to use. Defaults to 1.
#' @param parse_result Should the taxonomy be returned as the `efetch` returns it or should it be parsed into a tibble.
#'
#' @return Vector of Tax Ids.
#'
#' @export
parallel_get_tax <- function(organisms_taxIDs,
                             parse_result = TRUE, # default value for parse_result
                             total_cores = 1,
                             retry_times = 10,
                             verbose = FALSE,
                             env_name = "entrez-env") {
  organisms_taxIDs <- unique(organisms_taxIDs)
  check_cmd("efetch", env_name = env_name)

  # Set up parallel plan
  # future::plan(future::multisession, workers = total_cores)
  if (total_cores > 1L) {
    future::plan(future::multisession(),
      workers = total_cores
    )
  }

  results <- tibble::tibble(
    "Sci_name" = character(0L),
    "query_taxID" = character(0L),
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
  res_taxid <- character(0L)
  retry_count <- 0L
  while (isTRUE(retry_count < retry_times) && isFALSE(all(organisms_taxIDs %in% res_taxid))) {
    # message(retry_count)
    # message(organisms_taxIDs)
    if (total_cores > 1L) {
      results_temp <- furrr::future_map_dfr(
        .x = organisms_taxIDs,
        .f = get_tax_by_taxID,
        parse_result = parse_result,
        env_name = env_name,
        verbose = verbose,
        .options = furrr::furrr_options(seed = NULL),
        .progress = TRUE
      )
    } else {
      results_temp <- purrr::map_dfr(
        .x = organisms_taxIDs,
        .f = get_tax_by_taxID,
        parse_result = parse_result,
        env_name = env_name,
        verbose = verbose,
        .progress = TRUE
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

  return(results)
}

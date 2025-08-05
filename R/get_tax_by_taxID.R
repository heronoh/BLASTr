#' Retrieve Taxonomic Ranks Using NCBI Taxonomy Tax IDs
#'
#' Retrieves complete taxonomy information for given NCBI Taxonomy Tax IDs by querying the NCBI database using the `efetch` command.
#'
#' @param organisms_taxIDs A character vector of NCBI Taxonomy Tax IDs for which to retrieve taxonomy information.
#' @param parse_result Logical indicating whether to parse the taxonomy information into a tibble (`TRUE`, default) or return the raw output as returned by `efetch` (`FALSE`).
#' @param verbose Logical indicating whether to print verbose messages during the process. Default is `FALSE`.
#' @param env_name Character string specifying the name of the conda environment where `efetch` is installed. Default is `"blastr-entrez-env"`.
#'
#' @return A tibble containing the taxonomic ranks for the given Tax IDs.
#'
#' @examples
#' \dontrun{
#' # Retrieve taxonomy for a single Tax ID
#' tax_info <- get_tax_by_taxID("9606") # Human
#'
#' # Retrieve taxonomy for multiple Tax IDs
#' tax_ids <- c("9606", "10090", "10116") # Human, Mouse, Rat
#' tax_info <- get_tax_by_taxID(tax_ids)
#'
#' # Get unparsed taxonomy result
#' raw_tax_info <- get_tax_by_taxID("9606", parse_result = FALSE)
#'
#' # Enable verbose output
#' tax_info <- get_tax_by_taxID("9606", verbose = TRUE)
#' }
#' @export
get_tax_by_taxID <- function(
    organisms_taxIDs,
    parse_result = TRUE,
    verbose = FALSE,
    env_name = "blastr-entrez-env") {
  `%>%` <- dplyr::`%>%`
  .data <- rlang::.data

  organisms_taxIDs <- as.character(organisms_taxIDs)

  check_cmd("efetch", env_name = env_name)

  # run entrez command
  organism_xml <- condathis::run(
    "efetch",
    "-db", "taxonomy",
    "-id", organisms_taxIDs,
    "-format", "xml",
    env_name = env_name,
    verbose = verbose,
    error = "continue"
  )

  # Create empty tibble to return in case of error
  organism_tbl_parsed_empty <- tibble::tibble(
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

  # testing if efetch returned a valid result
  if (isFALSE(stringr::str_detect(organism_xml$stdout, "TaxId"))) {
    message(paste0("------------------------> unable to retrieve taxonomy for: ", organisms_taxIDs, "\t"))
    return(organism_tbl_parsed_empty)
  }

  # testing integrity of xml output
  xml_teste <- tryCatch(
    {
      xml_file <- xml2::read_xml(organism_xml$stdout)
      message("XML file is valid!")
    },
    error = function(e) {
      message(paste("Error in XML file:", e$message))
    }
  )

  if (isTRUE(stringr::str_detect(xml_teste, "^Error"))) {
    message(paste0("------------------------> unable to retrieve taxonomy for: ", organisms_taxIDs, "\t"))
    return(organism_tbl_parsed_empty)
  }

  # processing valid xml
  organism_list <- organism_xml$stdout |>
    xml2::read_xml() |>
    xml2::as_list()

  message(paste0("taxonomy succesfully retrieved for: ", organisms_taxIDs))

  organism_tbl <- organism_list$TaxaSet$Taxon$LineageEx |>
    unname() |>
    purrr::map(function(x) {
      tibble::tibble(
        Rank = x$Rank[[1]],
        ScientificName = x$ScientificName[[1]]
      )
    }) |>
    purrr::list_rbind() |>
    dplyr::mutate("query_taxID" = organisms_taxIDs) |>
    dplyr::mutate(
      "Sci_name" = unlist(organism_list$TaxaSet$Taxon$ScientificName)
    )

  if (rlang::is_true(parse_result)) {
    temp_names_tbl <- tibble::tibble(
      "Sci_name" = character(0L),
      "query_taxID" = character(0L),
      # "division" = character(0L),
      "superkingdom" = character(0L),
      "kingdom" = character(0L),
      "phylum" = character(0L),
      "subphylum" = character(0L),
      "class" = character(0L),
      "subclass" = character(0L),
      "order" = character(0L),
      "suborder" = character(0L),
      "family" = character(0L),
      "subfamily" = character(0L),
      "genus" = character(0L)
    )

    organism_tbl_parsed <- organism_tbl |>
      dplyr::filter(!(.data$Rank %in% c("no rank", "clade"))) |>
      dplyr::distinct() |>
      dplyr::filter(.data$Rank %in% c(
        "superkingdom", "kingdom",
        "phylum", "subphylum", "class",
        "subclass", "order", "suborder",
        "family", "subfamily", "genus"
      )) |>
      tidyr::pivot_wider(
        id_cols = c("query_taxID", "Sci_name"),
        names_from = "Rank",
        values_from = c("ScientificName")
      ) |>
      dplyr::bind_rows(
        temp_names_tbl
      )
    organism_tbl_parsed <- organism_tbl_parsed |>
      dplyr::relocate(
        "Sci_name", "query_taxID",
        "superkingdom", "kingdom",
        "phylum", "subphylum", "class", "subclass", "order",
        "suborder", "family", "subfamily", "genus"
      ) |>
      dplyr::rename(
        "Superkingdom (NCBI)" = "superkingdom",
        "Kingdom (NCBI)" = "kingdom",
        "Phylum (NCBI)" = "phylum",
        "Subphylum (NCBI)" = "subphylum",
        "Class (NCBI)" = "class",
        "Subclass (NCBI)" = "subclass",
        "Order (NCBI)" = "order",
        "Suborder (NCBI)" = "suborder",
        "Family (NCBI)" = "family",
        "Subfamily (NCBI)" = "subfamily",
        "Genus (NCBI)" = "genus"
      ) %>%
      dplyr::filter(dplyr::if_any(dplyr::everything(), ~ !base::is.na(.)))
    organism_tbl_final <- organism_tbl_parsed
  }
  if (rlang::is_false(parse_result)) {
    organism_tbl_unparsed_empty <- tibble::tibble(
      "Rank" = character(0L),
      "ScientificName" = character(0L),
      "query_taxID" = character(0L),
      "Sci_name" = character(0L)
    )
    organism_tbl_unparsed <- organism_tbl |>
      dplyr::bind_rows(
        organism_tbl_unparsed_empty
      )
    organism_tbl_final <- organism_tbl_unparsed
  }
  return(organism_tbl_final)
}

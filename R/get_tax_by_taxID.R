#' @title Get taxonomy ranks using _NCBI Taxonomy Tax ID_
#'
#' @description Recover complete taxonomy for a given _NCBI Taxonomy Tax ID_
#'
#' @param organisms_taxIDs Vector of _NCBI Taxonomy Tax ID_ to retrieve taxonomy for.
#' @param parse_result Should the taxonomy be returned as the `efetch` returns it or should it be parsed into a tibble.
#' @param total_cores Number of threads to use. Defaults to 1.
#' @param show_command Should the `efetch` command on command line be printed for debugging. Defaults to TRUE.
#'
#' @inheritParams parallel_get_tax
#'
#' @return A tibble with all the taxonomic ranks for the corresponding taxID.
#'
#' @export
get_tax_by_taxID <- function(organisms_taxIDs,
                             parse_result = TRUE,
                             verbose = FALSE,
                             env_name = "entrez-env") {
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
    verbose = verbose
  )


  if (isFALSE(stringr::str_detect(organism_xml$stdout, "TaxId"))) {
    message(paste0("------------------------> unable to retrieve taxonomy for: ", organisms_taxIDs, "\t"))
    return(organism_tbl_parsed_empty)
  }

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
    dplyr::mutate("Sci_name" = unlist(organism_list$TaxaSet$Taxon$ScientificName))

  # create empty tibble to return in case of error (required for parallel_get_tax() to work)
  organism_tbl_parsed_empty <- tibble::tibble(
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



  if (rlang::is_true(parse_result)) {
    temp_names_tbl <- tibble::tibble(
      "Sci_name" = character(0L),
      "query_taxID" = character(0L),
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

    organism_tbl_parsed <- organism_tbl %>%
      dplyr::filter(!(.data$Rank %in% c("no rank", "clade"))) %>%
      # dplyr::select(-c("taxID")) %>%
      dplyr::distinct() %>%
      # dplyr::filter(Rank %in% c("kingdom","phylum","class","order","family")) %>%
      dplyr::filter(.data$Rank %in% c(
        "superkingdom", "kingdom",
        "phylum", "subphylum", "class",
        "subclass", "order", "suborder",
        "family", "subfamily", "genus"
      )) %>%
      tidyr::pivot_wider(
        id_cols = c("query_taxID", "Sci_name"),
        names_from = "Rank",
        values_from = c("ScientificName")
      ) |>
      dplyr::bind_rows(
        temp_names_tbl
      )

    organism_tbl_parsed <- organism_tbl_parsed %>%
      dplyr::relocate(
        "Sci_name", "query_taxID", "superkingdom", "kingdom",
        "phylum", "subphylum", "class", "subclass", "order",
        "suborder", "family", "subfamily", "genus"
      ) %>%
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
    # create empty tibble for binding
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

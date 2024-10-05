#' @title Get taxonomy ranks using _NCBI Taxonomy Tax ID_
#'
#' @description Recover complete taxonomy for a given _NCBI Taxonomy Tax ID_
#'
#' @param organism_taxID Vector of _NCBI Taxonomy Tax ID_ to retrieve taxonomy for.
#' @param parse_result Should the taxonomy be returned as the `efetch` returns it or should it be parsed into a tibble.
#' @param total_cores Number of threads to use. Defaults to 1.
#' @param show_command Should the `efetch` command on command line be printed for debugging. Defaults to TRUE.

#'
#' @return A tibble with all the taxonomic ranks for the corresponding taxID.
#'
#' @export
get_tax_by_taxID <- function(organism_taxID,
                             parse_result = TRUE,
                             total_cores = 1,
                             verbose = FALSE,
                             env_name = "entrez-env") {
  `%>%` <- dplyr::`%>%`
  .data <- rlang::.data

  # This space ensures the taxID won't have it first digit chopped of
  # organism_taxID_parsed <- paste0(" ", organism_taxID)

  # run entrez command
  organism_xml <- condathis::run(
    "efetch",
    "-db", "taxonomy",
    "-id", organism_taxID,
    "-format", "xml",
    env_name = env_name,
    verbose = verbose
  )

  # entrez_json <- organism_xml$stdout
  entrez_xml <- organism_xml$stdout |>
    xml2::read_xml()


  if (length(organism_xml$stdout) == 0) {
    message(paste0("------------------------> unable to retrieve taxonomy for: ", organism_taxID))

    return(tibble::tibble())
  } else {
    if (isTRUE(xml2::xml_length(entrez_xml) == 1)) {
      organism_list <- entrez_xml |>
        xml2::as_list()

      message(paste0("taxonomy succesfully retrieved for: ", organism_taxID))

      organism_tbl <- organism_list$TaxaSet$Taxon$LineageEx |>
        unname() |>
        purrr::map(function(x) {
          tibble::tibble(
            Rank = x$Rank[[1]],
            ScientificName = x$ScientificName[[1]]
          )
        }) |>
        purrr::list_rbind() |>
        dplyr::mutate("query_taxID" = organism_taxID) |>
        dplyr::mutate("Sci_name" = unlist(organism_df$TaxaSet$Taxon$ScientificName))

      if (isFALSE(parse_result)) {
        return(organism_tbl)
      } else {
        organism_tbl_parsed_empty <- tibble::tibble(
          "Sci_name" = NA,
          "query_taxID" = NA,
          "superkingdom" = NA,
          "kingdom" = NA,
          "phylum" = NA,
          "subphylum" = NA,
          "class" = NA,
          "subclass" = NA,
          "order" = NA,
          "suborder" = NA,
          "family" = NA,
          "subfamily" = NA,
          "genus" = NA
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
          )

        organism_tbl_parsed <- dplyr::bind_rows(
          organism_tbl_parsed_empty,
          organism_tbl_parsed
        ) %>%
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



        return(organism_tbl_parsed)
      }
    } else {
      return(tibble::tibble())
    }
  }
}

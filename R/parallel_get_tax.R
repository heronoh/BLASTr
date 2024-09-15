#' @title Get taxonomy ranks for a list of _NCBI Taxonomy Tax ID_
#'
#' @description Get taxonomy ranks for a list of _NCBI Taxonomy Tax ID_ in parallel.
#' @param taxIDs Vector of _NCBI Taxonomy Tax ID_ to retrieve taxonomy for.
#' @param parse_result Should the taxonomy be returned as the _efetch_ returns it or should it be parsed into a tibble.
#' @param total_cores Number of threads to use. Defaults to 1.
#'
#' @param perc_id Lowest identity percentage cutoff.
#'   Passed on to BLAST+ _-perc_identity_.
#' @param perc_qcov_hsp Lowest query coverage per HSP percentage cutoff.
#' Passed on to BLAST+ _-qcov_hsp_perc_.
#' @param num_alignments Number of alignments to retrieve from BLAST. Max = 6.
#' @param blast_type One of the available BLAST+ search engines,
#'   one of: "blastn", "blastp", "blastx", "tblastn", "tblastx".
#'
#' @return Unformatted BLAST results.
#'   For results formatted as tibble, please use `BLASTr::get_blast_results()`
#'
#' @examples
#' blast_res <- BLASTr::run_blast(
#'   asv = "CTAGCCATAAACTTAAATGAAGCTATACTAAACTCGTTCGCCAG
#'     AGTACTACAAGCGAAAGCTTAAAACTCATAGGACTTGGCGGTGTTTCAGACCCAC",
#'   db_path = "/data/databases/nt/nt",
#'   perc_id = 80,
#'   num_thread = 1,
#'   perc_qcov_hsp = 80,
#'   num_alignments = 2
#' )
#'
#' @export
#'
#'

parallel_get_tax <- function(organisms_taxIDs,
                             parse_result = TRUE, # default value for parse_result
                             total_cores = 1) {
  `%>%` <- dplyr::`%>%`

  # Set up parallel plan
  future::plan(multisession, workers = total_cores)




# Helper function to process one organism_taxID
get_single_tax_by_taxID <- function(organism_taxID,
                                    parse_result = TRUE,
                                    total_cores=1) {
  `%>%` <- dplyr::`%>%`




  #This space ensures the taxID won't have it first digit chopped of
  organism_taxID_parsed <- paste0(" ",organism_taxID)

  #print(class(organism_taxID_parsed))

  #Generate entrez _efetch_ command
  entrez_cmd <- paste0("efetch -db taxonomy -id ${organism_taxID_parsed} -format xml -json")

  print(entrez_cmd)

  organism_xml <- shell_exec(cmd = entrez_cmd)
  # organism_xml <- system2(command = "efetch",args = c("-db", "taxonomy", "-id", organism_taxID, "-format", "xml", "-json"))

  if (length(organism_xml$stdout) == 0) {

    message(paste0("------------------------> unable to retrieve taxonomy for: ",organism_taxID))

    return(tibble::tibble("Sci_name" = NA,
                          "query_taxID" = NA,
                          "Superkingdom (NCBI)" = NA,
                          "Kingdom (NCBI)" = NA,
                          "Phylum (NCBI)" = NA,
                          "Subphylum (NCBI)" = NA,
                          "Class (NCBI)" = NA,
                          "Subclass (NCBI)" = NA,
                          "Order (NCBI)" = NA,
                          "Suborder (NCBI)" = NA,
                          "Family (NCBI)" = NA,
                          "Subfamily (NCBI)" = NA,
                          "Genus (NCBI)" = NA))
  }else{

    if (jsonlite::validate(organism_xml$stdout)) {
      organism_df <- organism_xml$stdout %>%
        jsonlite::fromJSON()

      message(paste0("taxonomy succesfully retrieved for: ",organism_taxID))

      organism_tbl <- organism_df$TaxaSet$Taxon$LineageEx$Taxon %>%
        dplyr::as_tibble() %>%
        dplyr::mutate("query_taxID" = organism_taxID) %>%
        dplyr::mutate("Sci_name" = organism_df$TaxaSet$Taxon$ScientificName)

      if (isFALSE(parse_result)) {
        return(organism_tbl)

      }else{

        organism_tbl_parsed_empty <- tibble::tibble("Sci_name" = NA,
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
                                                    "genus" = NA)

        organism_tbl_parsed <- organism_tbl %>%
          dplyr::filter(!Rank %in% c("no rank","clade")) %>%
          dplyr::select(-c("TaxId")) %>%
          unique() %>%
          # dplyr::filter(Rank %in% c("kingdom","phylum","class","order","family")) %>%
          dplyr::filter(Rank %in% c("superkingdom","kingdom",
                                    "phylum","subphylum","class",
                                    "subclass","order","suborder",
                                    "family","subfamily","genus")) %>%
          tidyr::pivot_wider(
            id_cols = c("query_taxID","Sci_name"),
            names_from = "Rank",
            values_from = c("ScientificName"))

        organism_tbl_parsed <- dplyr::bind_rows(organism_tbl_parsed_empty,
                                                organism_tbl_parsed) %>%
          dplyr::relocate("Sci_name","query_taxID","superkingdom","kingdom",
                          "phylum","subphylum","class","subclass","order",
                          "suborder","family","subfamily","genus") %>%
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
            "Genus (NCBI)" = "genus") %>%
          dplyr::filter(dplyr::if_any(dplyr::everything(), ~!base::is.na(.)))



        return(organism_tbl_parsed)
      }


    }else{

      return(tibble::tibble())

    }
  }



}
# Apply the function in parallel using future_map_dfr
results <- furrr::future_map_dfr(organisms_taxIDs,
                                 get_single_tax_by_taxID,

                                 .progress = TRUE,
                                 .options = furrr::furrr_options(seed = NULL))

return(results)
}



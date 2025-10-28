# Get `taxid` from header of fasta file

# TODO: @luciorq WIP - This do not actually work yet.
get_taxid_by_accession <- function(accession) {
  cmd <- paste0(
    "esearch -db nucleotide -query ", accession,
    " | elink -target taxonomy",
    " | efetch -format xml"
  )
  result <- condathis::run(
    cmd = "bash",
    "-c", cmd,
    env_name = "blastr-entrez-env",
    error = "cancel",
    verbose = "full"
  )
  taxid_str <- xml2::read_xml(result$stdout) |>
    xml2::xml_find_first("//TaxId") |>
    xml2::xml_text()

  return(taxid_str)
}


# AY882416.1 Homo sapiens isolate 38_U6a(Tor41) mitochondrion, complete genome
get_taxid_by_accession("AY882416")


# AP011979.1 Gymnotus carapo mitochondrial DNA, almost complete genome
get_taxid_by_accession("AP011979.1")

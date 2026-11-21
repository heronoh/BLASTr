# A tiny NCBI taxdump: root > cellular organisms > Eukaryota (superkingdom)
# > Metazoa > Chordata > Primates (order) > Hominidae > Homo > Homo sapiens,
# plus Pan troglodytes; 9999 was merged into 9606.
tiny_taxdump <- function(dir) {
  dmp <- function(...) paste0(paste(..., sep = "\t|\t"), "\t|")
  writeLines(c(
    dmp(1, 1, "no rank"), dmp(131567, 1, "cellular root"), dmp(2759, 131567, "superkingdom"),
    dmp(33208, 2759, "kingdom"), dmp(7711, 33208, "phylum"), dmp(9443, 7711, "order"),
    dmp(9604, 9443, "family"), dmp(9605, 9604, "genus"), dmp(9606, 9605, "species"),
    dmp(9596, 9604, "genus"), dmp(9598, 9596, "species")
  ), file.path(dir, "nodes.dmp"))
  writeLines(c(
    dmp(1, "root", "", "scientific name"), dmp(131567, "cellular organisms", "", "scientific name"),
    dmp(2759, "Eukaryota", "", "scientific name"), dmp(33208, "Metazoa", "", "scientific name"),
    dmp(7711, "Chordata", "", "scientific name"), dmp(9443, "Primates", "", "scientific name"),
    dmp(9604, "Hominidae", "", "scientific name"), dmp(9605, "Homo", "", "scientific name"),
    dmp(9606, "Homo sapiens", "", "scientific name"), dmp(9606, "human", "", "genbank common name"),
    dmp(9596, "Pan", "", "scientific name"), dmp(9598, "Pan troglodytes", "", "scientific name")
  ), file.path(dir, "names.dmp"))
  writeLines(dmp(9999, 9606), file.path(dir, "merged.dmp"))
  dir
}


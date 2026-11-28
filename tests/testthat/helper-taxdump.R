# A tiny NCBI taxdump: root > cellular organisms > Eukaryota (superkingdom)
# > Metazoa > Chordata > Primates (order) > Hominoidea (superfamily, not one of
# the ranks of get_tax()) > Hominidae > Homo > Homo sapiens, plus Pan
# troglodytes and Hylobatidae > Hylobates > Hylobates lar; 9999 was merged
# into 9606.
tiny_taxdump <- function(dir) {
  dmp <- function(...) paste0(paste(..., sep = "\t|\t"), "\t|")
  writeLines(c(
    dmp(1, 1, "no rank"), dmp(131567, 1, "cellular root"), dmp(2759, 131567, "superkingdom"),
    dmp(33208, 2759, "kingdom"), dmp(7711, 33208, "phylum"), dmp(9443, 7711, "order"),
    dmp(314295, 9443, "superfamily"), dmp(9604, 314295, "family"), dmp(9605, 9604, "genus"), dmp(9606, 9605, "species"),
    dmp(9596, 9604, "genus"), dmp(9598, 9596, "species"),
    dmp(9577, 314295, "family"), dmp(9578, 9577, "genus"), dmp(9580, 9578, "species")
  ), file.path(dir, "nodes.dmp"))
  writeLines(c(
    dmp(1, "root", "", "scientific name"), dmp(131567, "cellular organisms", "", "scientific name"),
    dmp(2759, "Eukaryota", "", "scientific name"), dmp(33208, "Metazoa", "", "scientific name"),
    dmp(7711, "Chordata", "", "scientific name"), dmp(9443, "Primates", "", "scientific name"),
    dmp(314295, "Hominoidea", "", "scientific name"), dmp(9604, "Hominidae", "", "scientific name"), dmp(9605, "Homo", "", "scientific name"),
    dmp(9606, "Homo sapiens", "", "scientific name"), dmp(9606, "human", "", "genbank common name"),
    dmp(9596, "Pan", "", "scientific name"), dmp(9598, "Pan troglodytes", "", "scientific name"),
    dmp(9577, "Hylobatidae", "", "scientific name"), dmp(9578, "Hylobates", "", "scientific name"),
    dmp(9580, "Hylobates lar", "", "scientific name")
  ), file.path(dir, "names.dmp"))
  writeLines(dmp(9999, 9606), file.path(dir, "merged.dmp"))
  dir
}


# A taxonomy database of tiny_taxdump(), removed with the calling test.
lca_db <- function() {
  dir <- tiny_taxdump(withr::local_tempdir(.local_envir = parent.frame()))
  make_taxonomy_db(file.path(withr::local_tempdir(.local_envir = parent.frame()), "tax.sqlite"), taxdump = dir) |>
    suppressMessages()
}

testthat::test_that("UNITE and SILVA classifications become names to look up", {
  u <- unite_levels(c(
    "k__Fungi;p__Basidiomycota;f__Entolomataceae;g__Entoloma;s__Entoloma_sp",
    "k__Fungi;f__Thelephoraceae;g__Thelephoraceae_gen_Incertae_sedis;s__Thelephoraceae_sp",
    "k__Viridiplantae;g__Bupleurum;s__Bupleurum_fruticosum"
  ))
  testthat::expect_equal(u[[1]], c("Fungi", "Basidiomycota", "Entolomataceae", "Entoloma", NA))
  testthat::expect_equal(u[[2]], c("Fungi", "Thelephoraceae", "Thelephoraceae", NA))
  testthat::expect_equal(u[[3]], c("Viridiplantae", "Bupleurum", "Bupleurum fruticosum"))
  s <- silva_levels("Bacteria;Pseudomonadota;Klebsiella;uncultured bacterium")
  testthat::expect_equal(s[[1]], c("Bacteria", "Pseudomonadota", "Klebsiella", NA))
})

testthat::test_that("names_to_taxids takes the deepest name found and resolves homonyms", {
  dir <- tiny_taxdump(withr::local_tempdir())
  dmp <- function(...) paste0(paste(..., sep = "\t|\t"), "\t|")
  # a second (fake) genus "Homo", in Hylobatidae
  cat(dmp(50000, 9577, "genus"), "\n", sep = "", file = file.path(dir, "nodes.dmp"), append = TRUE)
  cat(dmp(50000, "Homo", "", "scientific name"), "\n", sep = "", file = file.path(dir, "names.dmp"), append = TRUE)
  db <- make_taxonomy_db(file.path(withr::local_tempdir(), "tax.sqlite"), taxdump = dir) |>
    suppressMessages()
  levels <- list(
    c("Primates", "Hominidae", "Homo", "Homo sapiens"),
    c("Primates", "Hominidae", "Homo", NA),
    c("Primates", "Hylobatidae", "Homo", NA),
    c("Primates", "Hominidae", "Homo", "Homo sapiens"),
    c("Nothing", NA)
  )
  res <- names_to_taxids(levels, taxonomy_db = db)
  testthat::expect_equal(res$taxid, c("9606", "9605", "50000", "9606", NA))
  testthat::expect_equal(res$name, c("Homo sapiens", "Homo", "Homo", "Homo sapiens", NA))
})

testthat::test_that("read_fasta_records joins lines; make_unique_ids dedups", {
  f <- withr::local_tempfile(fileext = ".fasta")
  writeLines(c(">a x y", "acg", "UU-", ">b", "TTT"), f)
  x <- read_fasta_records(f)
  testthat::expect_equal(x$header, c("a x y", "b"))
  testthat::expect_equal(x$seq, c("ACGUU", "TTT"))
  testthat::expect_equal(make_unique_ids(c("A1", "B 2", "A1")), c("A1", "B_2", "A1_2"))
})

testthat::test_that("the online taxonomy table has the columns of get_tax()", {
  xml <- xml2::read_xml(paste0(
    "<TaxaSet><Taxon><TaxId>9606</TaxId><ScientificName>Homo sapiens</ScientificName>",
    "<Rank>species</Rank><AkaTaxIds><TaxId>9999</TaxId></AkaTaxIds><LineageEx>",
    "<Taxon><TaxId>131567</TaxId><ScientificName>cellular organisms</ScientificName><Rank>no rank</Rank></Taxon>",
    "<Taxon><TaxId>2759</TaxId><ScientificName>Eukaryota</ScientificName><Rank>superkingdom</Rank></Taxon>",
    "<Taxon><TaxId>9605</TaxId><ScientificName>Homo</ScientificName><Rank>genus</Rank></Taxon>",
    "</LineageEx></Taxon></TaxaSet>"
  ))
  res <- online_tax_table(c("9606", "9999", "1234"), parse_taxa_xml(xml), tax_ranks_default)
  dir <- tiny_taxdump(withr::local_tempdir())
  db <- make_taxonomy_db(file.path(withr::local_tempdir(), "tax.sqlite"), taxdump = dir) |>
    suppressMessages()
  testthat::expect_equal(names(res), names(get_tax("9606", taxonomy_db = db)))
  testthat::expect_equal(res$found, c(TRUE, TRUE, FALSE))
  testthat::expect_equal(res$merged_into, c(NA, "9606", NA))
  testthat::expect_equal(res$domain[1:2], c("Eukaryota", "Eukaryota"))
  testthat::expect_equal(res$genus[1], "Homo")
  testthat::expect_equal(res$lineage[1], "Eukaryota; Homo; Homo sapiens")
  testthat::expect_equal(res$lineage_taxids[1], "2759; 9605; 9606")
})

testthat::test_that("get_tax_by_taxID and get_tax_by_name ask the NCBI", {
  testthat::skip_on_cran()
  testthat::skip_if_offline("eutils.ncbi.nlm.nih.gov")
  res <- get_tax_by_taxID(c("9606", "10090"))
  testthat::expect_equal(res$sci_name, c("Homo sapiens", "Mus musculus"))
  testthat::expect_equal(res$genus, c("Homo", "Mus"))
  nm <- get_tax_by_name(c("Homo sapiens", "human"))
  testthat::expect_equal(nm$taxid, c("9606", "9606"))
  testthat::expect_equal(nm$matched_by, c("scientific name", "other name"))
})

testthat::test_that("list_ncbi_blast_dbs lists the NCBI databases", {
  testthat::skip_on_cran()
  testthat::skip_if_offline("ftp.ncbi.nlm.nih.gov")
  dbs <- list_ncbi_blast_dbs()
  testthat::expect_true(all(c("core_nt", "16S_ribosomal_RNA") %in% dbs$db))
  testthat::expect_true(all(c("db", "type", "size_gb", "last_updated") %in% names(dbs)))
})

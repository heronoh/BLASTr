testthat::test_that("uninformative_references reads the shipped table", {
  u <- uninformative_references()
  testthat::expect_named(u, c("field", "pattern", "category", "note"))
  testthat::expect_true(all(u$field %in% c("header", "accession")))
  testthat::expect_match(attr(u, "md5"), "^[0-9a-f]{32}$")
  testthat::expect_identical(attr(uninformative_references(u), "md5"), attr(u, "md5"))
  testthat::expect_error(uninformative_references(data.frame(field = "x", pattern = "a", category = "c")))
  testthat::expect_error(uninformative_references(data.frame(field = "header", pattern = "(a", category = "c")))
})

testthat::test_that("uninformative_category matches headers and accessions", {
  u <- uninformative_references()
  hd <- c(
    "Picea abies chloroplast rbcL", "uncultured bacterium clone X", "MAG: Bacteroidales bacterium",
    "Chordata environmental sample clone 1", "Cloning vector pUC19", "Synthetic construct GFP",
    "Spermatophyta environmental sample", "Homo sapiens chimeric gene", "Unidentified fungus",
    "PREDICTED: Nomascus leucogenys uncharacterized LOC115832432", "Environmental protection agency strain"
  )
  sub <- c(paste0("A", 1:9), "XR_004027927.1", "A11")
  cat <- uninformative_category(sub, hd, u)
  testthat::expect_equal(cat, c(
    NA, "uncultured", "metagenome", "environmental", "synthetic", "synthetic",
    "environmental", "chimeric", "unidentified", "bad reference", NA
  ))
  # accession with another version, and gi|...|acc| ids
  testthat::expect_equal(uninformative_category(c("XR_004027927.2", "gi|1|ref|XR_004027927.1|"), c("x", "y"), u), rep("bad reference", 2))
  testthat::expect_true(all(is.na(uninformative_category("A1", "uncultured", NULL))))
})

testthat::test_that("hits are ranked among the informative ones", {
  u <- uninformative_references()
  hits <- tibble::tibble(
    query = c(rep("q1", 4), rep("q2", 2)),
    subject = paste0("s", 1:6),
    `subject header` = c("uncultured bacterium", "Picea abies", "Picea glauca", "Pinus taeda", "MAG: x", "uncultured y"),
    staxid = c("77133", "3329", "3330", "3352", "1", "2"),
    identity = c(100, 99, 99, 95, 98, 97), qcovhsp = 100, bitscore = c(500, 490, 490, 450, 480, 470),
    blast_rank = c(1:4, 1:2)
  )
  long <- blast_long(hits, seq_ids = c("q1", "q2"), seqs = c("AAA", "CCC"), num_alignments = 2, uninformative = u)
  q1 <- long[long$Sequence == "AAA", ]
  testthat::expect_equal(q1$subject, c("s2", "s3", "s1"))
  testthat::expect_equal(q1$rank, c(1L, 1L, NA))
  testthat::expect_equal(q1$at_best, c(TRUE, TRUE, FALSE))
  testthat::expect_equal(q1$informative, c(TRUE, TRUE, FALSE))
  testthat::expect_equal(q1$uninformative_category, c(NA, NA, "uncultured"))
  # no informative hit: the uninformative ones are ranked among themselves
  q2 <- long[long$Sequence == "CCC", ]
  testthat::expect_equal(q2$rank, c(1L, 2L))
  testthat::expect_equal(q2$at_best, c(TRUE, FALSE))
  s <- blast_summary(long, search_alignments = 50L)
  testthat::expect_equal(s$`Check ID`, c("check ID", "check ID"))
  testthat::expect_equal(s$`N hits at best`, c(2L, 1L))
  testthat::expect_equal(s$`N uninformative`, c(1L, 2L))
  testthat::expect_equal(s$`Best uninformative score`, c(100, round(similarity_score(98, 100), 6)))
  # a long result without flags, flagged again
  plain <- blast_long(hits, seq_ids = c("q1", "q2"), seqs = c("AAA", "CCC"))
  testthat::expect_equal(plain$subject[plain$Sequence == "AAA"][1], "s1")
  again <- as_blast_long(plain, uninformative = u)
  testthat::expect_equal(again$subject[again$Sequence == "AAA"][1], "s2")
})

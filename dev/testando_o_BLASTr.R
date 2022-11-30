library(devtools)
library(BLASTr)
load_all(path = "/home/heron/prjcts/omics/metaseqs/BLASTr")

ASVs_test <- c("CTAGCCATAAACTTAAATGAAGCTATACTAAACTCGTTCGCCAGAGTACTACAAGCGAAAGCTTAAAACTCATAGGACTTGGCGGTGTTTCAGACCCAC",
               "CTAGCCATAAACTTAAATGAAGCTATACTAAACTCGTTCGCCAGAGTACTACAAGTGAAAGCTTAAAACTCATAGGACTTGGCGGTGTTTCAGACCCAC",
               "GCCAAATTTGTGTTTTGTCCTTCGTTTTTAGTTAATTGTTACTGGCAAATGACTAACGACAAATGATAAATTACTAATAC")




# blaste_res_1ASV <- get_blast_results(asv = "CTAGCCATAAACTTAAATGAAGCTATACTAAACTCGTTCGCCAGAGTACTACAAGCGAAAGCTTAAAACTCATAGGACTTGGCGGTGTTTCAGACCCAC",num_thread = 10)


future::availableCores()


# options(BLASTr.dbapth = "/data/databases/nt/nt")

blast_res <- BLASTr::parallel_blast(asvs = ASVs_test[1],
               db_path = "/data/databases/nt/nt",
               out_file = "/home/heron/prjcts/omics/BLASTr_run/blast_out.csv",
               out_RDS = "/home/heron/prjcts/omics/BLASTr_run/blast_out.RDS",
               # blast_cmd = "blastn -db {db_path} -outfmt '6 std qcovhsp' -max_hsps 1 -perc_identity 80 -qcov_hsp_perc 80 -num_threads {as.character(num_thread)} -num_alignments {as.character(num_alignments)}",
               total_cores = 77,
               perc_ID = 80,
               num_thread = 1,
               perc_qcov_hsp = 80,
               num_alignments = 4
               # ,
               # blast_type = "blastn"
               )


blast_res1 <- BLASTr::run_blast(asv = ASVs_test[3],
               db_path = "/data/databases/nt/nt",
               perc_ID = 80,
               num_thread = 1,
               perc_qcov_hsp = 80,
               num_alignments = 2,
               blast_type = "blastn")





teste <- readRDS(file = "/home/heron/prjcts/omics/metaseqs/blast_out.RDS")
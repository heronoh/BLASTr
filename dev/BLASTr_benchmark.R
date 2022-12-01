
# Testing
#
library(tictoc)
library(BLASTr)

#load test ASVs ----
ASVs_test <- c("CTAGCCATAAACTTAAATGAAGCTATACTAAACTCGTTCGCCAGAGTACTACAAGCGAAAGCTTAAAACTCATAGGACTTGGCGGTGTTTCAGACCCAC",
               "CTAGCCATAAACTTAAATGAAGCTATACTAAACTCGTTCGCCAGAGTACTACAAGTGAAAGCTTAAAACTCATAGGACTTGGCGGTGTTTCAGACCCAC",
               "GCCAAATTTGTGTTTTGTCCTTCGTTTTTAGTTAATTGTTACTGGCAAATGACTAACGACAAATGATAAATTACTAATAC",
               "AACATTGTATTTTGTCTTTGGGGCCTGGGCAGGTGCAGTAGGAACTTCACTTAGAATAATTATTCGTACTGAGCTTGGGCATCCAGGAAGACTTATCGGGGATGATCAAATCTATAATGTAATTGTTACAGCACATGCATTTGTGATAATTTTTTTTATAGTAATACCTATTATGATT",
               "ACTATACCTATTATTCGGCGCATGAGCTGGAGTCCTAGGCACAGCTCTAAGCCTCCTTATTCGAGCCGAGCTGGGCCAGCCAGGCAACCTTCTAGGTAACGACCACATCTACAACGTTATCGTCACAGCCCATGCATTTGTAATAATCTTCTTCATAGTAATACCCATCATAATCGGAGGCTTTGGCAACTGACTAGTTCCCCTAATAATCGGTGCCCCCGATATG",
               "TTAGCCATAAACATAAAAGTTCACATAACAAGAACTTTTGCCCGAGAACTACTAGCAACAGCTTAAAACTCAAAGGACTTGGCGGTGCTTTATATCCAC",
               "AACTCTTTACTTTATCTTTGGTACTTGAGCGGGTCTGATAGGATCTTCACTGAGAGCTCTAATCCGTCTCGAATTAGGACAACCAGGTTCTCTTCTAGGAAACGATCAAGTATATAACACTATTGTAACTGCC",
               "TTAGCCCTAAACACAGATAATTACATAAACAAAATTATTCGCCAGAGTACTACCAGCAACAGCTTAAAACTCAAAGGACTTGGCGGTGCTTTATATCCTT",
               "TTAGCCCTAAACACAGATAATTACATAAACAAAATTATTCGCCAGAGTACTACTAGCAACAGCTTAAAACTCAAAGGACTTGGCGGTGCTTTATATCCTT",
               "TTAGCCCTAAACATAGATAGTTTTACAACAAAATAATTCGCCAGAGGACTACTAGCAATAGCTTAAAACTCAAAGGACTTGGCGGTGCTTTATATCCCT",
               "CTAGCCCTAAATCTAGATACCTCCCATCACACATGTATCCGCCTGAGAACTACGAGCACAAACGCTTAAAACTCTAAGGACTTGGCGGTGCCCCAAACCCAC",
               "TTAGCCTTAAACCTAGATGTTAAATTTACAAGCAACATCCGCCAGGGTACTACAAGCGCTAGCTTAAAACCCAAAGGACTTGACGGTGTCTCAGACCCAC",
               "CACCCTCTACTTAGTATTCGGTGCCTGAGCCGGAATAGTTGGCACAGCCCTTAGCCTTCTAATTCGGGCAGAGCTGGCCCAACCTGGCGCCCTCCTAGGTGATGATCAAATTTACAACGTTATCGTTACTGCT",
               "CACCCTTTACCTAGTTTTCGGTGCATGAGCCGGAATAGTGGGCACAGCTTTAAGCCTCCTAATCCGAGCAGAATTAAGCCAGCCCGGTTCACTCCTGGGGGACGACCAGATTTACAACGTAATCGTAACAGCA",
               "TTAGCCATAAACATAAAAGTTCACATAACAAGAACTTTCGCCCGAGAACTACTAGCAACAGCTTAAAACTCAAAGGACTTGGCGGTGCTTTATATCCAC",
               "AACTTTATATTTAATATTTGCTTTATTTTCAGGATTATTAGGTACTGCTATGTCAGTATTAATAAGATTAGAATTAAGTGGACCAGGAGTGCAATTCATATCAAATAACCAATTATATAATAGTATTGTTACAGCT",
               "AACTTTATACCTAAATTTTTGGAACCTGAGCTGGTCTTCTAGGAGCCTCTCTGAGGGCTCTAATCCGATTGGAATTAGGGCAGCCCGGATCCCTTATAGAAAATGATCAAATTTACAACACAATCGTAACAGCT",
               "CACCCTCTATTTAGTGTTTGGTGCCTGGGCCGGCATAGTCGGAACAGCCCTAAGTCTCCTAATCCGAGCAGAGCTCGGCCAACCTGGAACCCTCCTTGGTGACGACCAAGTATATAATGTCATCGTTACCGCC",
               "TTAGCCCTAAACCCAGATATTTAACAAACAAAAATATCCGCCAGAGAACTACGAGCAAACGCTTAAAACTCTAAGGACTTGGCGGTGCCCTAAACCCCC",
               "TTAGCCATAAACCTAAATAATTAAATTTAACAAAACTATTTGCCAGAGAACTACTAGCCACAGCTTAAAACTCAAAGGACTTGGCGGTACTTTATATCCAT",
               "TTAGCCTTAAACCCCAATAACTCCACCAACAAAACTATTCGCCAGAACACTACAAGCAATAGCTTAAAACTCAAAGGACCTGGCGGTGCTTTATATCCGT",
               "CACCCTCTACCTAGTATTTGGTGCATGAGCCGGAATAGTTGGAACCGCCCTCAGCCTCTTAATTCGAGCAGAGCTCAGTCAACCCGGATCCCTACTGGGCGACGATCAGATCTATAATGTTATCGTTACTGCA",
               "CTAGCCCTAAACCCAAATAGTTACATAACAAAACTATTCGCCAGAGTACTACTCGCAACTGCCTAAAACTCAAAGGACTTGGCGGTGCTTCACATCCAC",
               "TTAGCCCTAAATCTTGATGCTATTTCATACTAAAGCATCCGCCTGAGAACTACGAGCACAAACGCTTAAAACTCTAAGGACTTGGCGGTGCCCTAAACCCAC",
               "CTAGCCGTAAACCAAAACAATTAATAAACAAAATTGTTCGCCAGAGTACTACTAGCAACAGCTTAAAACTCAAAGGACTTGGCGGTGCTTTATATCCTT",
               "CTGGCCCTAAATCTAGATGCTTATGCTACCAAAGCATCCGCCCGAGGACTACGAGCACAAACGCTTAAAACTCTAAGGACTTGGCGGTGCCCCAAACCCAC",
               "CACCCTCTACTTAAGTATTTGGTGCTTGAGCCGGAATAGTTGGTACAGCCCTCAGCCTCCTAATTCGAGCAGAGCTCAGCCAACCCGGAGCCCTTCTAGGCGATGATCAAATTTATAATGTTATCGTTACCGCA",
               "CACCCTCTACTTAGTATTTGGTGCTTGAGCCGGAATAGTTGGCACAGCCCTAAGTCTGCTAATTCGAGCAGAATTAGCCCAACCCGGCGCCCTTCTGGGGGATGACCAAATTTATAATGTTATTGTTACTGCC",
               "CACCCTCTACTTAGTATTTGGTGCCTGAGCCGGTATAGTTGGCACAGCCCTCAGCCTCCTAATTCGAGCGGAGCTTAGCCAACCCGGAGCCCTACTAGGTGATGATCAGATTTATAATGTTATCGTTACCGCA",
               "AACCCTATACTTAATCTTCGGCGCATGGGCCGGAATAGTGGGAACCGCCCTAAGCCTATTAATTCGGGCTGAACTCAGCCAACCGGGGTCTCTCCTTGGCGATGATCAAATTTACAACGTCATTGTAACCGCT",
               "CACCCTTTACCTAGTTTTCGGTGCATGAGCCGGAATAGTGGGCACAGCTTTAAGCCTCCTAATCCGAGCAGAATTAAGCCAGCCCGGTTCACTCCTGGGGGACGACCAGATTTACAACGTAATCGTAACAGCA",
               "CACCCTCTATATAGTATTTGGTGCCTGAGCCGGAATGATTGGCACGGCCCTAAGCCTGCTAATTCGAGCAGAGCTCAGCCAACCAGGGGCCCTTCTGGGGGACGACCAAATCTACAACGTAATCGTCACCGCC",
               "CACCCTCTACCTAGTATTTGGTGCTTGAGCCGGCATAGTTGGTACAGCCCTTAGCCTCTTAATCCGAGCGGAGCTTAGCCAACCCGGGTCCCTACTGGGCGATGACCAGATTTATAATGTTATCGTTACTGCA",
               "TTGGCCCTAAATCTAGACACTTATCTCCTACCAAAGTGTCCGCCCGAGAACTACGAGCACGAACGCTTAAAACTCTAAGGACTTGGCGGTGCTCCAAACCCAC",
               "CTAGCCCTAAATCTTGATGTTCCCCAAACACAAACATCCGCCCGAGAACTACGAGCACAAACGCTTAAAACTCTAGGGACTTGGCGGTGCCCTAAACCCAC",
               "CTGGCCCTAAATCTTGATACTTATTATACCAAAGTATCCGCCAGAGAACTACGAGCACAAACGCTTAAAACTCTAAGGACTTGGCGGTGCCCTAAACCCAC",
               "TTAGCCCTAAACTTAGATAGTTATCCTAAACAAAACTATCCGCCAGAGAACTACTAGCAATAGCTTAAAACCCAAAGGACTTGGCGGTGCTTTACATCCCT",
               "TTAGCCCTAAACTTAGATAGTTACCCTAAACAAAGCTATCCGCCAGAGAACTACTAGCAATAGCTTAAAACTCAAAGGACTTGGCGGTGCTTTACATCCCT",
               "TTAGCCCTAAACATAGATAATTTTAACATAACAAAATTGCATGCCAGAGAACTACTAGCAACAGCTTAAAACTCAAAGGACTTGGCGGTGCTTTATATCCCT",
               "TTAGCCCTAAACCTAGATAGTTAGCTCAAACAAAACTATCCGCCAGAGAACTACTAGCAACAGCTTAAAACTCAAAGGGCTTGGCGGTGCTTTATATCCCT",
               "CTGGCCCTAAATCTCGATGCTTACATAACCAAAGCATCCGCCTGAGGACTACGAGCACAAACGCTTAAAACTCTAAGGACTTGGCGGTGCCCCAAACCCAC",
               "TCAGTCATAAACCTAGACGTCCTAATACAGCAAGACGTCCGCCCGGGTACTACGAGCGCTTGGCTTGAAACCCAAAGGACCTGACGGTGCCTCAGACCCCC",
               "CTAGCCCTAAACCTAAATAATCGACCAACAAGATTATTCGCCAGAGTACTACTAGCAACAGCCTAAAACTCAAAGGACTTGACGGTGCTTCATATCCAT",
               "TTAGTTCTAAACCCAAATAGTTCAACCAACAAAACTATTCACCAGGGTACTACAAGCAACAGCTTAAGACTCAAAAGACTTGGCAGTGCTTTATATCCCT",
               "TTAGCCTTAAACTCAGATGCTACAAATACAAATAACATCCGCCAGGGTACTACAAGCGCTAGCTTAAAACCCAAAGGACTTGACGGTGTCTCAGACCCAC",
               "TTAGCCCTTAACAAAGGTGTACCCCACACATCCTACACCCGCCAGAAGATTACGAGCCAAGCTTAAGACTCAAAGGACTTGACAGCACTTCAAATCCAC",
               "CTAGCCCTAAACTAAAACAGTTCTCAAACAAAACTGTTCGCCAGAGTACTACTAGCAATAGCTTAAAACTCAAAGGACTTGGCGGTGCTTTACATCCTT",
               "CACCCTATATTTAATCTTTGGTGCCTGAGCCGGAATAGTCGGTACAGCCCTTAGCCTTCTCATTCGAGCAGAGCTTGCCCAACCCGGCGCCCTATTGGGTGATGACCAAATTTATAATGTTATTGTCACTGCC",
               "CACCCTTTACCTAGTATTCGGTGCCTGAGCTGGGATAGTTGGCACAGCTCTTAGCCTCCTCATCCGAGCAGAATTAAGTCAACCTGGCTCCCTGCTAGGCGATGATCAAATTTACAATGTTATCGTAACCGCA",
               "CACCCTCTACCTGGTATTTGGTGCATGAGCCGGAATAGTTGGAACCGCCCTCAGCCTATTAATTCGAGCAGAGCTTAGCCAACCCGGATCCCTACTGGGCGACGATCAGATCTATAATGTTATCGTTACTGCG")




# Sequential versions ----

# sequencial com 80 threads ----
tictoc::tic("Sequential - Purrr - 80 threads")
teste80 <- purrr::map_dfr(ASVs_test, get_blast_results, num_threads = 80,
                          db_path = "/data/databases/nt/nt",
                          # out_file = "/home/heron/prjcts/omics/BLASTr_run/blast_out.csv",
                          # out_RDS = "/home/heron/prjcts/omics/BLASTr_run/blast_out.RDS",
                          perc_ID = 80,
                          perc_qcov_hsp = 80,
                          num_alignments = 4,
                          blast_type = "blastn")
tictoc::toc()

# sequencial com 60 threads ----
tictoc::tic("Sequential - Purrr - 60 threads")
teste60 <- purrr::map_dfr(ASVs_test, get_blast_results, num_threads = 60,
                          db_path = "/data/databases/nt/nt",
                          # out_file = "/home/heron/prjcts/omics/BLASTr_run/blast_out.csv",
                          # out_RDS = "/home/heron/prjcts/omics/BLASTr_run/blast_out.RDS",
                          perc_ID = 80,
                          perc_qcov_hsp = 80,
                          num_alignments = 4,
                          blast_type = "blastn")
tictoc::toc()

# sequencial com 40 threads ----
tictoc::tic("Sequential - Purrr - 40 threads")
teste40 <- purrr::map_dfr(ASVs_test, get_blast_results, num_threads = 40,
                          db_path = "/data/databases/nt/nt",
                          # out_file = "/home/heron/prjcts/omics/BLASTr_run/blast_out.csv",
                          # out_RDS = "/home/heron/prjcts/omics/BLASTr_run/blast_out.RDS",
                          perc_ID = 80,
                          perc_qcov_hsp = 80,
                          num_alignments = 4,
                          blast_type = "blastn")
tictoc::toc()

# sequencial com 20 threads ----
tictoc::tic("Sequential - Purrr - 20 threads")
teste20 <- purrr::map_dfr(ASVs_test, get_blast_results, num_threads = 20,
                          db_path = "/data/databases/nt/nt",
                          # out_file = "/home/heron/prjcts/omics/BLASTr_run/blast_out.csv",
                          # out_RDS = "/home/heron/prjcts/omics/BLASTr_run/blast_out.RDS",
                          perc_ID = 80,
                          perc_qcov_hsp = 80,
                          num_alignments = 4,
                          blast_type = "blastn")
tictoc::toc()

# sequencial com 10 threads ----
tictoc::tic("Sequential - Purrr - 10 threads")
teste10 <- purrr::map_dfr(ASVs_test, get_blast_results, num_threads = 10,
                          db_path = "/data/databases/nt/nt",
                          # out_file = "/home/heron/prjcts/omics/BLASTr_run/blast_out.csv",
                          # out_RDS = "/home/heron/prjcts/omics/BLASTr_run/blast_out.RDS",
                          perc_ID = 80,
                          perc_qcov_hsp = 80,
                          num_alignments = 4,
                          blast_type = "blastn")
tictoc::toc()

# sequencial com 5 threads ----
tictoc::tic("Sequential - Purrr 5 threads")
teste5 <- purrr::map_dfr(ASVs_test, get_blast_results, num_threads = 5,
                         db_path = "/data/databases/nt/nt",
                         # out_file = "/home/heron/prjcts/omics/BLASTr_run/blast_out.csv",
                         # out_RDS = "/home/heron/prjcts/omics/BLASTr_run/blast_out.RDS",
                         perc_ID = 80,
                         perc_qcov_hsp = 80,
                         num_alignments = 4,
                         blast_type = "blastn")
tictoc::toc()
# sequencial com 2 threads ----
tictoc::tic("Sequential - Purrr 2 threads")
teste2 <- purrr::map_dfr(ASVs_test, get_blast_results, num_threads = 2, db_path = "/data/databases/nt/nt",
                         # out_file = "/home/heron/prjcts/omics/BLASTr_run/blast_out.csv",
                         # out_RDS = "/home/heron/prjcts/omics/BLASTr_run/blast_out.RDS",
                         perc_ID = 80,
                         perc_qcov_hsp = 80,
                         num_alignments = 4,
                         blast_type = "blastn")
tictoc::toc()
# sequencial com 1 threads ----
tictoc::tic("Sequential - Purrr 1 threads")
teste1 <- purrr::map_dfr(ASVs_test, get_blast_results, num_threads = 1, db_path = "/data/databases/nt/nt",
                         perc_ID = 80,
                         perc_qcov_hsp = 80,
                         num_alignments = 4,
                         blast_type = "blastn")
tictoc::toc() # Sequential - Purrr 1 threads: 133.009 sec elapsed


# testing with furrr on parallel----

cores_to_be_used <- future::availableCores() - 2 # Usar todos os cores -2 = 78
future::plan(future::multisession(workers = cores_to_be_used))

# paralela com 80 threads ----
tictoc::tic("Parallel - Furrr 80 threads")
teste_F80 <- furrr::future_map_dfr(ASVs_test,
                                   get_blast_results,
                                   num_thread = 80,
                                   .options = furrr::furrr_options(seed = NULL),
                                   db_path = "/data/databases/nt/nt",
                                   perc_ID = 80,
                                   perc_qcov_hsp = 80,
                                   num_alignments = 4,
                                   blast_type = "blastn")
tictoc::toc()


# paralela com 60 threads ----
tictoc::tic("Parallel - Furrr 60 threads")
teste_F60 <- furrr::future_map_dfr(ASVs_test,
                                   get_blast_results,
                                   num_thread = 60,
                                   .options = furrr::furrr_options(seed = NULL),
                                   db_path = "/data/databases/nt/nt",
                                   perc_ID = 80,
                                   perc_qcov_hsp = 80,
                                   num_alignments = 4,
                                   blast_type = "blastn")
tictoc::toc()


# paralela com 40 threads ----
tictoc::tic("Parallel - Furrr 40 threads")
teste_F40 <- furrr::future_map_dfr(ASVs_test,
                                   get_blast_results,
                                   num_thread = 40,
                                   .options = furrr::furrr_options(seed = NULL),
                                   db_path = "/data/databases/nt/nt",
                                   perc_ID = 80,
                                   perc_qcov_hsp = 80,
                                   num_alignments = 4,
                                   blast_type = "blastn")
tictoc::toc()


# paralela com 20 threads ----
tictoc::tic("Parallel - Furrr 20 threads")
teste_F20 <- furrr::future_map_dfr(ASVs_test,
                                   get_blast_results,
                                   num_thread = 20,
                                   .options = furrr::furrr_options(seed = NULL),
                                   db_path = "/data/databases/nt/nt",
                                   perc_ID = 80,
                                   perc_qcov_hsp = 80,
                                   num_alignments = 4,
                                   blast_type = "blastn")
tictoc::toc()


# paralela com 10 threads ----
tictoc::tic("Parallel - Furrr 10 threads")
teste_F10 <- furrr::future_map_dfr(ASVs_test,
                                   get_blast_results,
                                   num_thread = 10,
                                   .options = furrr::furrr_options(seed = NULL),
                                   db_path = "/data/databases/nt/nt",
                                   perc_ID = 80,
                                   perc_qcov_hsp = 80,
                                   num_alignments = 4,
                                   blast_type = "blastn")
tictoc::toc()


# paralela com 5 threads ----
tictoc::tic("Parallel - Furrr 5 threads")
teste_F5 <- furrr::future_map_dfr(ASVs_test,
                                  get_blast_results,
                                  num_thread = 5,
                                  .options = furrr::furrr_options(seed = NULL),
                                  db_path = "/data/databases/nt/nt",
                                  perc_ID = 80,
                                  perc_qcov_hsp = 80,
                                  num_alignments = 4,
                                  blast_type = "blastn")
tictoc::toc()


# paralela com 2 threads ----
tictoc::tic("Parallel - Furrr 2 threads")
teste_F2 <- furrr::future_map_dfr(ASVs_test,
                                  get_blast_results,
                                  num_thread = 2,
                                  .options = furrr::furrr_options(seed = NULL),
                                  db_path = "/data/databases/nt/nt",
                                  perc_ID = 80,
                                  perc_qcov_hsp = 80,
                                  num_alignments = 4,
                                  blast_type = "blastn")
tictoc::toc()



# paralela com 1 threads ----
tictoc::tic("Parallel - Furrr 1 threads")
teste_F1 <- furrr::future_map_dfr(ASVs_test,
                                  get_blast_results,
                                  num_thread = 1,
                                  .options = furrr::furrr_options(seed = NULL),
                                  db_path = "/data/databases/nt/nt",
                                  perc_ID = 80,
                                  perc_qcov_hsp = 80,
                                  num_alignments = 4,
                                  blast_type = "blastn")
tictoc::toc()







base::save.image("~/prjcts/omics/BLASTr_benchmark/BLASTr_benchmark.RData")













#plots ----


library(dplyr)
library(tidyr)
library(ggplot2)
library(artists)
options(scipen = 33)
options(pillar.sigfig = 7)


bench_res <- read.csv(file = "/home/heron/prjcts/omics/BLASTr_benchmark/BLASTr_bench_res.csv",sep = ";") %>% as_tibble()



bench_plot <- bench_res %>%
  unite(Num_ASVs, Type, col = "group",remove = F,sep = " ASVs ") %>%
  ggplot(aes(x=num_threads,
             y=seconds/Num_ASVs,
             # y=seconds,
             col=group,
             Group=group,
             shape=Package,
             linetype=Package))+
  geom_point(size=3)+
  geom_path()+
  scale_y_log10()+
  scale_color_manual(values = c("#2A6BBF", "#3FACAA", "#F2C84B", "#D99036")) +
  theme_light(base_size = 12)+
  ylab(label = "Seconds per ASV")+
  xlab(label = "Number of threads per BLAST search")+
  geom_rect(inherit.aes=FALSE, aes(xmin=7.5, xmax=12.5,
                                   ymin=0, ymax=Inf),
            color="transparent",
            fill="orange",
            alpha=0.005) +
  guides(colour = guide_legend("Group")) +
  xlim(0,85) +
  scale_x_continuous(breaks=seq(0, 80, 5))

ggsave(bench_plot,
       filename = "~/prjcts/omics/BLASTr/dev/bench_plot.pdf",
       device = "pdf",
       units = "in",
       width = 10,
       height = 6,dpi = 300)




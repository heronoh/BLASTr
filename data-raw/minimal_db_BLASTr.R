# select minimal DB for installation testing ----

# load libs ----
library(DECIPHER)
library(Biostrings)
library(tidyr)
library(tibble)
library(stringr)


# select IDs from the test results ----

# TODO: @heronoh How teste20 was obtained?

teste20$`1_subject` |>
  unique() |>
  paste0(collapse = ",")



# get original seqs from original nt ----
# blastdbcmd -db /data/databases/nt/nt -entry "KX902240.1,KJ710708.1" -outfmt %f

# blastdbcmd -db /data/databases/nt/nt -entry "AP011979.1,CP030121.1,MG967958.1,MN849868.1,KX381515.1,MH611404.1,JN319671.1,MT048504.1,KU290564.1,MN013407.1,KF209672.1,HM405192.1,GU701891.1,KF479445.1,KU374630.1,MF123486.1,LT599484.1,KR020500.1,KX381421.1,JN989232.1,MT048590.1,AY894142.1,KX180157.1,HM640211.1,MK506873.1,JN988778.1,MH537267.1,KU494528.1,MZ051800.1,KF562408.1,U77157.1,JQ235619.1,KJ710708.1,KX786344.1,KP202276.1,KX381634.1,KX381509.1,KX902240.1,MW652802.1,KY987555.1,XR_004027927.1,MG755611.1,MN885912.1,KT818547.1,MZ051986.1,FJ749051.1,KX868689.1" -outfmt %f > dev/minimal_db_BLASTr.fasta


# create sql conection to store seqs ----
db_BLASTr <- dbConnect(
  SQLite(), "/home/heron/prjcts/BLASTr/dev/minimal_db/minimal_db_BLASTr.sql"
)

# get sequences for the minimal test database ----
Seqs2DB(
  seqs = "/home/heron/prjcts/BLASTr/dev/minimal_db/minimal_db_BLASTr.fasta",
  type = "FASTA", dbFile = db_BLASTr, identifier = "BLASTr"
)

dna_BLASTr <- SearchDB(
  dbFile = db_BLASTr, identifier = "BLASTr", nameBy = "description"
)


# trim too long sequences ----
dna_BLASTr[2] <- subseq(
  x = dna_BLASTr[2],
  start = 885700,
  end = 885900
)

dna_BLASTr[2] <- subseq(
  x = dna_BLASTr[2],
  start = 885700,
  end = 885900
)

dna_BLASTr |> names()

# update tbl with STARTs and ENDs larger than match ----
BLASTr_tbl <- teste20 |>
  select(
    c("1_subject header", "1_subject", "1_subject start", "1_subject end")
  ) |>
  mutate(
    START = if_else(
      `1_subject start` > `1_subject end`, `1_subject end`, `1_subject start`
    ),
    END = if_else(
      `1_subject start` > `1_subject end`, `1_subject start`, `1_subject end`
    )
  ) |>
  filter(!is.na(`1_subject header`)) |>
  unique()




# get subsequence of the original DB seq to make shortest DB ----
dna_BLASTr[str_detect(
  string = names(dna_BLASTr),
  pattern = BLASTr_tbl$`1_subject`[8]
)]

dna_BLASTr[2]


for (seq in 1:nrow(BLASTr_tbl)) {
  dna_BLASTr[str_detect(
    string = names(dna_BLASTr),
    pattern = BLASTr_tbl$`1_subject`[seq]
  )] <- subseq(
    x = dna_BLASTr[str_detect(
      string = names(dna_BLASTr),
      pattern = BLASTr_tbl$`1_subject`[seq]
    )],
    start = if_else(
      (BLASTr_tbl$START[seq] - 30 >= 1), (BLASTr_tbl$START[seq] - 30), 1
    ),
    end = if_else(
      (BLASTr_tbl$END[seq] + 30) >= BLASTr_tbl$END[seq],
      BLASTr_tbl$END[seq],
      (BLASTr_tbl$END[seq] + 30)
    )
  )
}

dna_BLASTr |> seqlengths()
# save shortest DB into fasta

Biostrings::writeXStringSet(
  x = dna_BLASTr,
  filepath = "/home/heron/prjcts/BLASTr/dev/minimal_db/shortest_minimal_db_BLASTr.fasta",
  format = "fasta"
)


# format blast db ----

# /prjcts/BLASTr/dev/minimal_db$ makeblastdb -in shortest_minimal_db_BLASTr.fasta -dbtype "nucl" -parse_seqids -hash_index

# test! ----





ASVs_test <- readr::read_lines("dev/asv_test.txt")



# paralela com 2 threads ----
tictoc::tic("Parallel - Furrr 2 threads")

ASVs_test

teste_F2local_short <- furrr::future_map_dfr(
  get_blast_results,
  num_thread = 2,
  .options = furrr::furrr_options(seed = NULL),
  # db_path = "/data/databases/nt/nt",
  # db_path = "/home/heron/prjcts/BLASTr/dev/minimal_db/shortest_minimal_db_BLASTr.fasta",
  db_path = fs::path_package("BLASTr", "inst", "ext", "minimal_db"),
  perc_id = 80,
  perc_qcov_hsp = 80,
  num_alignments = 4,
  blast_type = "blastn"
)
tictoc::toc()

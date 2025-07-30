# BLASTr

An R package for parallel taxonomic classification of metabarcoding sequences.

<!-- badges: start -->
<!--[![R-CMD-check](https://github.com/heronoh/BLASTr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/heronoh/BLASTr/actions/workflows/R-CMD-check.yaml)-->
<!-- badges: end -->

The BLASTr package is a powerful tool for performing BLAST operations from within the R environment.
Despite being initially developed for metagenomic applications, the BLASTr package is flexible for other applications.
It can be used with any installed NCBI BLAST+ search strategy, configuration, and any database on UNIX and Windows platforms through Windows Subsystem for Linux.
This makes it a versatile tool that can be used in various applications and contexts requiring sequence identification on tabular outputs.
Additionally, the package includes documentation and functions for parsing and analyzing the results of BLAST searches, making it easier for users to extract useful information from their BLAST results.
Overall, the BLASTr package is a valuable tool for bioinformaticians and researchers who need to perform BLAST operations from within R.

## Requirements

The BLASTr package requires the NCBI BLAST+ to be installed.
The easiest way to perform its installation is on the UNIX command line.

``` bash
#install NCBI-BLAST+

sudo apt update
sudo apt install ncbi-blast+

#chech installation using one of its applicaitons
blastn -help

# identify the executable complete path
which blastn

```

## Installation

<!--
### R-Universe

The main way to install `BLASTr` is through:

```r
install.packages('BLASTr', repos = "https://heronoh.r-universe.dev")
```
-->

### Development version

You can install the development version of BLASTr from [GitHub](https://github.com/heronoh/BLASTr) with:

```r
# install.packages("remotes")
remotes::install_github("heronoh/BLASTr")
```

## Documentation

The main place to look for help and documentation is <heronoh.github.io/BLASTr/>

## Database configuration

### Obtaining NCBI complete databases

Identifications can be performed using NCBI complete databases, such as NT, which are readily available to download and update.
This can be performed using the *BLAST+* _script_ [update_blastdb.pl](https://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl).

``` bash
#set a folder to download the desired database (for example, the nt database)
BLAST_DB_PATH="/data/database/blast/nt"

#create dir
mkdir -p  "${BLAST_DB_PATH}"

#enter dir
cd "${BLAST_DB_PATH}"

#suggestion: use screen or tmux to emulate a terminal. The downloads usually takes long.
#          tmux: https://tmuxcheatsheet.com/
#          screen: https://kapeli.com/cheat_sheets/screen.docset/Contents/Resources/Documents/index

#user BLAST+ executable to download/update db files
update_blastdb --passive --decompress nt

#set permissions to enable usage by all users
chown root "${BLAST_DB_PATH}"/*
chmod 755 "${BLAST_DB_PATH}"/*
```

Another option is to download it directly fom the [NCBI ftp site](https://ftp.ncbi.nlm.nih.gov/blast/db/). It can be parallelized and is the best choice when you want to download only the new files.

```bash
#suggestion: use screen or tmux to emulate a terminal. The downloads usually takes long.
#          tmux: https://tmuxcheatsheet.com/
#          screen: https://kapeli.com/cheat_sheets/screen.docset/Contents/Resources/Documents/index

# download volumes and md5 check files
seq -w 000 150 | parallel wget https://ftp.ncbi.nlm.nih.gov/blast/db/nt.{}.tar.gz -t 0 --show-progress
seq -w 000 150 | parallel wget https://ftp.ncbi.nlm.nih.gov/blast/db/nt.{}.tar.gz.md5 -t 0 --show-progress
     # where 000 is the first volume and 150, the last (up to now).

ls *5 | parallel md5sum -c {} >> check.txt
sort check.txt > check_sort.txt

ls *tar.gz | parallel tar -xvzf {}
```

### Formating a custom database

In case you prefer, any FASTA file can be formated as a BLAST+ database, using the *BLAST+* _script_ *_makeblastdb_*.

``` bash
#set the path to your fasta file (replace the example below)
DB_FILE="/data/database/my_db.fasta"

#check parameters and usage
makeblastdb

#format your db
makeblastdb -in "${DB_FILE}" -dbtype "nucl" -parse_seqids -hash_index
```


## Testing

The package installation can be with a mock BLAST formatted database provided in this [link](https://drive.google.com/file/d/1Qy4w4KIGSTiGjx-J4BrcyN6Y5wtRYBHl/view?usp=sharing).
Alternatively, you can download the unformated mock database (a fasta file) and format it using the ncbi-blast+ functionalities, as you would do for any other custom database.
You can obtain the fasta [here](https://drive.google.com/file/d/1WKIwq7RySleuSotWjZ4OJUXqbXELaE7q/view?usp=sharing).

```bash
#set the path to your fasta file (replace the example below)
DB_FILE="/data/database/shortest_minimal_db_BLASTr.fasta"

#check parameters and usage
makeblastdb

#format your db
makeblastdb -in "${DB_FILE}" -dbtype "nucl" -parse_seqids -hash_index

```

For the testing, execute the following commands on your R console.

```r
library(BLASTr)


#here are 8 ASVs to be tested with the mock blast DB

ASVs_test <- c(
  "CTAGCCATAAACTTAAATGAAGCTATACTAAACTCGTTCGCCAGAGTACTACAAGCGAAAGCTTAAAACTCATAGGACTTGGCGGTGTTTCAGACCCAC",
  "CTAGCCATAAACTTAAATGAAGCTATACTAAACTCGTTCGCCAGAGTACTACAAGTGAAAGCTTAAAACTCATAGGACTTGGCGGTGTTTCAGACCCAC",
  "GCCAAATTTGTGTTTTGTCCTTCGTTTTTAGTTAATTGTTACTGGCAAATGACTAACGACAAATGATAAATTACTAATAC",
  "AACATTGTATTTTGTCTTTGGGGCCTGGGCAGGTGCAGTAGGAACTTCACTTAGAATAATTATTCGTACTGAGCTTGGGCATCCAGGAAGACTTATCGGGGATGATCAAATCTATAATGTAATTGTTACAGCACATGCATTTGTGATAATTTTTTTTATAGTAATACCTATTATGATT",
  "ACTATACCTATTATTCGGCGCATGAGCTGGAGTCCTAGGCACAGCTCTAAGCCTCCTTATTCGAGCCGAGCTGGGCCAGCCAGGCAACCTTCTAGGTAACGACCACATCTACAACGTTATCGTCACAGCCCATGCATTTGTAATAATCTTCTTCATAGTAATACCCATCATAATCGGAGGCTTTGGCAACTGACTAGTTCCCCTAATAATCGGTGCCCCCGATATG",
  "TTAGCCATAAACATAAAAGTTCACATAACAAGAACTTTTGCCCGAGAACTACTAGCAACAGCTTAAAACTCAAAGGACTTGGCGGTGCTTTATATCCAC"
)

blast_res <- BLASTr::parallel_blast(
  asvs = ASVs_test,                                                 # vector of sequences to be searched
  db_path = "/data/database/shortest_minimal_db_BLASTr.fasta",      # path to a formated blast database
  out_file = NULL,                                                  # path to a .csv file to be created with results (on an existing folder)
  out_RDS = NULL,                                                   # path to a .RDS file to be created with results (on an existing folder)
  perc_id = 80,                                                     # minimum identity percentual cutoff
  perc_qcov_hsp = 80,                                               # minimum percentual coverage of query sequence by subject sequence cutoff
  num_threads = 1,                                                  # number of threads/cores to run each blast on
  total_cores = 8,                                                  # number of tota threads/cores to alocate all blast searches
  num_alignments = 3,                                               # maximum number of alignments/matches to retrieve results for each query sequence
  blast_type = "blastn"                                             # blast search engine to use
)

# check identificaitons results

blast_res

#or

View(blast_res)
```


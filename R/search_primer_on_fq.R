library(tibble)

search_primers_on_fq <- function(primer_seqs, fastq_paths) {
  # Define IUPAC codes as a named list
  iupac_codes <- list(
    A = "A", C = "C", G = "G", T = "T",
    W = "[AT]", S = "[CG]",
    M = "[AC]", K = "[GT]",
    R = "[AG]", Y = "[CT]",
    B = "[CGT]", D = "[AGT]",
    H = "[ACT]", V = "[ACG]",
    N = "[ACGT]"
  )


  # Convert each DNA sequence into a parsed regular expression
  parsed_primers <- sapply(primer_seqs, function(primer) {
    paste0(sapply(strsplit(primer, NULL)[[1]],
                  function(x) iupac_codes[[x]]),
           collapse = "")
  }, USE.NAMES = TRUE)

  # Initialize a list to store results
  results <- list()

  # Process each file and primer combination
  for (fastq_path in fastq_paths) {
    # Determine if the file is compressed
    is_gz <- grepl("\\.gz$", fastq_path)

    # Command to count total reads
    if (is_gz) {
      total_reads_cmd <- sprintf("zgrep -c '^@' %s", fastq_path)
    } else {
      total_reads_cmd <- sprintf("grep -c '^@' %s", fastq_path)
    }

    # Execute total reads command
    total_reads <- tryCatch(
      as.numeric(system(total_reads_cmd, intern = TRUE)),
      warning = function(w) 0,
      error = function(e) 0
    )

    # Process each primer
    for (primer_name in names(primer_seqs)) {
      primer_seq <- primer_seqs[[primer_name]]
      parsed_primer <- parsed_primers[[primer_name]]

      # Command to count primer matches
      if (is_gz) {
        primer_match_cmd <- sprintf("zgrep -cE '%s' %s", parsed_primer, fastq_path)
      } else {
        primer_match_cmd <- sprintf("grep -cE '%s' %s", parsed_primer, fastq_path)
      }

      # Execute primer match command
      primer_matches <- tryCatch(
        as.numeric(system(primer_match_cmd, intern = TRUE)),
        warning = function(w) 0,
        error = function(e) 0
      )

      # Calculate percentage
      percentage <- if (total_reads > 0) (primer_matches / total_reads) * 100 else 0

      # Append results
      results <- append(results, list(tibble::tibble(
        file_name = fastq_path,
        primer_name = primer_name,
        primer_sequence = primer_seq,
        parsed_primer = parsed_primer,
        total_reads = total_reads,
        primer_matches = primer_matches,
        percentage = percentage
      )))
    }
  }

  # Combine all results into a single tibble
  result_tibble <- bind_rows(results)

  return(result_tibble)
}

# Example usage
# primer_seqs <- c(FORWARD = "GGDACWGGWTGAACWGTWTAYCCHCC", REVERSE = "CCCTACGGGNGGCWGCAG")
# fastq_paths <- c("file1.fastq.gz", "file2.fastq.gz")
# results <- search_primers_on_fq(primer_seqs, fastq_paths)
# print(results)

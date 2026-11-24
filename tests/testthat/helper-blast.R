# Tests that run BLAST+ need BLAST+ 2.16.0 or newer, from the conda environment
# of BLASTr (install_dependencies()) or the system; they never install it.
skip_if_no_blast <- function() {
  testthat::skip_if_not(
    has_blast(),
    "BLAST+ 2.16.0 or newer not found (conda environment of BLASTr or system)"
  )
}

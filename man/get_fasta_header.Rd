% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_fasta_header.R
\name{get_fasta_header}
\alias{get_fasta_header}
\title{Get header from SubjectID}
\usage{
get_fasta_header(id, db_path, env_name = "blast-env", verbose = FALSE)
}
\arguments{
\item{id}{SubjectID from BLAST results or any NCBI Nucleotide identifier.}

\item{db_path}{Complete path do formatted BLAST database.}

\item{env_name}{The name of the conda environment with the parameter (i.e. "blast-env")}

\item{verbose}{Should condathis::run() internal command be shown?}
}
\value{
Complete identifier for the SubjectID as is on the database.
}
\description{
Retrieve complete sequence header from BLAST DB based on SubjectIDs
}
\examples{
\dontrun{
BLASTr::get_fasta_header(id = "AP011979.1", db_path = "/data/databases/nt/nt")
}
}

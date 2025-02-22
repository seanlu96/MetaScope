% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metascope_blast.R
\name{add_in_taxa_ncbi}
\alias{add_in_taxa_ncbi}
\title{Adds in taxa if input used NCBI database}
\usage{
add_in_taxa_ncbi(metascope_id_in, accession, BPPARAM)
}
\arguments{
\item{metascope_id_in}{MetaScope ID file with NCBI taxa qnames}

\item{BPPARAM}{An optional BiocParallelParam instance determining the
parallel back-end to be used during evaluation.}
}
\value{
data.frame or tibble of taxonomy information
}
\description{
Returns MetaScope Table with NCBI taxa in separate columns
}
\keyword{internal}

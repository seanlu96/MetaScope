% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convert_animalcules_patho.R
\name{convert_animalcules_patho}
\alias{convert_animalcules_patho}
\title{Create a multi-assay experiment from PathoScope 2.0 output for usage with
animalcules}
\usage{
convert_animalcules_patho(
  patho_counts,
  annot_path,
  which_annot_col,
  end_string = "-sam-report.tsv"
)
}
\arguments{
\item{patho_counts}{Character string, a directory filepath to the counts ID CSVs output by
\code{metascope_id()}.}

\item{annot_path}{The filepath to the CSV annotation file for the samples.}

\item{which_annot_col}{The name of the column of the annotation file
containing the sample IDs. These should be the same as the
\code{meta_counts} root filenames.}

\item{end_string}{The end string used at the end of the metascope_id files.
Default is "-sam-report.tsv".}
}
\value{
Returns a MultiAssay Experiment file of combined sample counts data.
The MultiAssay Experiment will have a counts assay ("MGX").
}
\description{
This function serves as a legacy integration method for usage with
PathoScope 2.0 outputs. Upon completion of the PathoScope 2.0 pipeline, users can analyze and visualize
abundances in their samples using the animalcules package. After running this function, the user should save the returned MAE
to an RDS file using a function like \code{saveRDS} to upload the output into
the \code{animalcules} package.
}

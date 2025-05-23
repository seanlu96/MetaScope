#' Make a Bowtie2 index
#'
#' This function is a wrapper for the \code{Rbowtie2::bowtie2_build} function.
#' It will create either small (.bt2) or large Bowtie2 indexes (.bt2l)
#' depending on the combined size of the reference fasta files.
#'
#' @param ref_dir The path to the directory that contains the reference files
#'   either uncompressed or compressed (.gz). NOTE: This directory should
#'   contain only the reference fasta files to be indexed.
#' @param lib_dir The path to the directory where Bowtie2 index files should
#'   be created.
#' @param lib_name The basename of the index file to be created (without the
#'   .bt2 or .bt2l extension)
#' @param bowtie2_build_options Optional: Options that can be passed to the
#'   mk_bowtie_index() function. All options should be passed as one string. To
#'   see all the available options that can be passed to the function use
#'   Rbowtie2::bowtie2_build_usage(). NOTE: Do not specify threads here.
#' @param threads The number of threads available to the function. Default is 1
#'   thread.
#' @param overwrite Whether existing files should be overwritten. Default is
#'   FALSE.
#' @return Creates the Bowtie2 indexes of the supplied reference .fasta files.
#'   Returns the path to the directory containing these files.
#'
#' @export
#'
#' @examples
#' #### Create a bowtie index from the example reference library
#'
#' ## Create a temporary directory to store the reference library
#' ref_temp <- tempfile()
#' dir.create(ref_temp)
#'
#' tmp_accession <- system.file("extdata", "example_accessions.sql", package = "MetaScope")
#'
#' ## Download reference genome
#' download_refseq('Bovismacovirus', reference = FALSE, representative = FALSE,
#'                 out_dir = ref_temp, compress = TRUE, patho_out = FALSE,
#'                 caching = TRUE, accession_path = tmp_accession)
#'
#' ## Create the reference library index files in the current directory
#' mk_bowtie_index(ref_dir = ref_temp, lib_dir = ref_temp,
#'                 lib_name = "target", threads = 1, overwrite = FALSE)
#'
#' ## Remove temporary directory
#' unlink(ref_temp, recursive = TRUE)
#'

mk_bowtie_index <- function(ref_dir, lib_dir, lib_name, bowtie2_build_options,
                            threads = 1, overwrite = FALSE) {
    # If user doesn't specify their own build options, then use default options
    if (missing(bowtie2_build_options)) {
        bowtie2_build_options <- paste("--threads", threads)
    } else {
        bowtie2_build_options <- paste(bowtie2_build_options,
                                       "--threads", threads)
    }
    ref_dir <- dir(ref_dir, full.names = TRUE)
    # Convert user specified path to absolute path for debugging purposes
    lib_dir <- tools::file_path_as_absolute(lib_dir)
    # Call the bowtie2_build function from Rbowtie2
    Rbowtie2::bowtie2_build(references = ref_dir,
                            bt2Index = file.path(lib_dir, lib_name),
                            ... = bowtie2_build_options,
                            overwrite = overwrite)
    message("Index building complete")
    return(tools::file_path_as_absolute(lib_dir))
}

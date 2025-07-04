getExtension <- function(file){
  ex <- strsplit(basename(file), split="\\.")[[1]]
  return(ex[-1])
}

# download_accessions.R
#' Download indexes required for MetaScope ID and MetaBlast modules
#'
#' This is a necessary step for all samples utilizing NCBI and SILVA databases
#' in the MetaScope pipeline. As specified by the user,
#' \code{download_accessions} will automatically download the NCBI accessions
#' database, the SILVA taxonomy database, and or the NCBI Blast 16S database and
#' prepare consolidated databases for downstream use with the MetaID and
#' MetaBLAST modules. This package relies on the \code{taxonomizr} package.
#'
#' @param ind_dir Character string. Directory filepath where indices should be
#'   saved. Required.
#' @param tmp_dir Character path to directory for storing temp files. (Useful
#'   to avoid redownloading) Defaults to \code{file.path(ind_dir, "tmp")}
#' @param remove_tmp_dir Delete tmp_dir after downloads are complete? Defaults
#'   to \code{TRUE}
#' @param NCBI_accessions_database Logical. Download taxonomizr NCBI accessions
#'   database? Defaults to \code{TRUE}.
#' @param NCBI_accessions_name Character string. Filename (with or without
#'   extension) to save taxonomizr NCBI accessions database. Defaults to
#'   \code{"accessionTaxa.sql"}.
#' @param silva_taxonomy_database Logical. Download SILVA taxonomy database?
#'   Defaults to \code{TRUE}.
#' @param silva_taxonomy_name Character string. Filename (with or without
#'   extension) to save SILVA taxonomy database.
#'   Defaults to the file supplied with the package,
#'   \code{"all_silva_headers.rds"}.
#' @param blast_16S_database Logical. Download BLAST 16S database?
#'   Defaults to \code{TRUE}.
#' @param blast_16S_name Character string. Filename (without extension) to save \
#'   BLAST 16S database.
#'   Defaults to the file supplied with the package,
#'   \code{"16S_ribosomal_RNA"}.
#'
#' @return Exports database(s) with names and to location specified by the user.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   download_accessions(
#'     ind_dir = "C:/Users/JohnSmith/Research",
#'     tmp_dir = file.path(ind_dir, "tmp"),
#'     remove_tmp_dir = TRUE,
#'     NCBI_accessions_database = TRUE,
#'     NCBI_accessions_name = "accessionTaxa.sql",
#'     silva_taxonomy_database = TRUE,
#'     silva_taxonomy_name = "all_silva_headers.rds")
#' }
#'

download_accessions <- function(ind_dir,
                                tmp_dir = file.path(ind_dir, "tmp"),
                                remove_tmp_dir = TRUE,
                                NCBI_accessions_database = TRUE,
                                NCBI_accessions_name = "accessionTaxa",
                                silva_taxonomy_database = TRUE,
                                silva_taxonomy_name = "all_silva_headers",
                                blast_16S_database = TRUE,
                                blast_16S_name = "16S_ribosomal_RNA") {
  if (!file.exists(ind_dir)) {
    stop("ind_dir does not exist")
  }

  if (!file.exists(tmp_dir)) {
    dir.create(tmp_dir)
  }
  if (NCBI_accessions_database) {
    message("Downloading NCBI accessions database")
    # Check if NCBI_accessions_name has extension
    if (!identical(getExtension(NCBI_accessions_name), "sql")) {
      NCBI_accessions_name <- paste0(NCBI_accessions_name, ".sql")
    }
    taxonomizr::prepareDatabase(file.path(ind_dir, NCBI_accessions_name),
                                tmpDir = tmp_dir)
  }
  if (silva_taxonomy_database) {
    message("Downloading SILVA taxonomy database")
    location <- "https://github.com/wejlab/metascope-docs/raw/main/all_silva_headers.rds"
    # CHeck if SILVA_taxonomy_name has extension
    if (!identical(getExtension(silva_taxonomy_name), "rds")) {
      NCBI_accessions_name <- paste0(silva_taxonomy_name, ".rds")
    }
    destination <- paste(ind_dir, silva_taxonomy_name, sep = "/")
    utils::download.file(location, destination)
  }

  if (remove_tmp_dir) {
    unlink(tmp_dir)
  }

  message("Download completed!")
}


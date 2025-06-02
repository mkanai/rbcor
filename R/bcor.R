#' Read BCOR File
#'
#' @description Opens a BCOR (binary correlation) file and returns a bcor object.
#' Supports both standard BCOR format (magic: "bcor1.1") and extended format
#' (magic: "bcor1.1x") which stores diagonal values explicitly.
#'
#' @param filename Path to the BCOR file
#' @param read_header Whether to read the header and metadata immediately (default: TRUE)
#' @param packed_threshold Size threshold for using packed symmetric matrices (default: 1000)
#'
#' @return A bcor object containing:
#' \itemize{
#'   \item{ptr}{External pointer to the C++ bcor object}
#'   \item{filename}{Path to the BCOR file}
#'   \item{nSNPs}{Number of SNPs in the file}
#'   \item{nSamples}{Number of samples used to compute correlations}
#'   \item{is_extended}{Whether the file uses extended format with diagonal values}
#' }
#'
#' @export
#' @examples
#' \dontrun{
#' bcor <- read_bcor("path/to/file.bcor")
#'
#' # Access metadata
#' meta <- bcor$get_meta()
#'
#' # Read full correlation matrix
#' corr_full <- bcor$read_corr()
#'
#' # Read subset of correlations
#' corr_subset <- bcor$read_corr(snps = c(1, 10, 20))
#'
#' # Read rectangular submatrix
#' corr_rect <- bcor$read_corr(snps = 1:10, snps2 = 11:20)
#'
#' # Read sparse matrix with threshold
#' corr_sparse <- bcor$read_corr(snps = 1:100, sparse = TRUE, threshold = 0.1)
#'
#' # Read packed symmetric matrix (memory efficient for large matrices)
#' corr_packed <- bcor$read_corr(packed = TRUE)
#'
#' # Get diagonal values (useful for extended format files)
#' diag_vals <- bcor$get_diagonal()
#' }
read_bcor <- function(filename, read_header = TRUE, packed_threshold = 1000) {
  if (!file.exists(filename)) {
    stop("File does not exist: ", filename)
  }

  bcor_obj <- bcor_open(filename, read_header)
  bcor_obj$packed_threshold <- as.integer(packed_threshold)
  class(bcor_obj) <- c("bcor", class(bcor_obj))

  # Add methods to the object
  bcor_obj$get_meta <- function() {
    bcor_get_meta(bcor_obj$ptr)
  }

  bcor_obj$get_diagonal <- function() {
    bcor_get_diagonal(bcor_obj$ptr)
  }

  bcor_obj$read_corr <- function(snps = NULL, snps2 = NULL, sparse = FALSE, packed = FALSE, threshold = 0.0) {
    # Validate inputs
    if (!is.null(snps)) {
      if (!is.numeric(snps) || any(snps < 1) || any(snps > bcor_obj$nSNPs)) {
        stop("SNP indices must be between 1 and ", bcor_obj$nSNPs)
      }
    }

    if (!is.null(snps2)) {
      if (!is.numeric(snps2) || any(snps2 < 1) || any(snps2 > bcor_obj$nSNPs)) {
        stop("SNP indices must be between 1 and ", bcor_obj$nSNPs)
      }
    }

    # Determine matrix dimensions for packed logic
    if (!is.null(snps2)) {
      # Rectangular matrix - always return regular matrix
      if (is.null(snps)) {
        stop("Must specify snps when using snps2")
      }

      if (packed) {
        warning("Rectangular matrices cannot use packed storage. Returning regular matrix.")
      }

      if (sparse) {
        return(bcor_read_corr_sparse2(bcor_obj$ptr, as.integer(snps), as.integer(snps2), threshold))
      } else {
        return(bcor_read_corr2(bcor_obj$ptr, as.integer(snps), as.integer(snps2)))
      }
    }

    # Square matrix or subset
    actual_size <- if (is.null(snps)) bcor_obj$nSNPs else length(snps)

    # Handle sparse vs packed conflicts
    if (sparse && packed) {
      warning("Both sparse=TRUE and packed=TRUE specified. Using sparse=TRUE.")
      packed <- FALSE
    }

    # Check if packed storage should be used
    use_packed <- packed && actual_size >= bcor_obj$packed_threshold

    if (packed && !use_packed) {
      warning(
        "Matrix size (", actual_size, "x", actual_size,
        ") below packed threshold (", bcor_obj$packed_threshold,
        "). Returning regular matrix."
      )
    }

    # Read correlations
    if (sparse) {
      bcor_read_corr_sparse(bcor_obj$ptr, snps, threshold)
    } else if (use_packed) {
      bcor_read_corr_packed(bcor_obj$ptr, snps)
    } else {
      bcor_read_corr(bcor_obj$ptr, snps)
    }
  }

  return(bcor_obj)
}


#' Print bcor Object
#'
#' @param x A bcor object
#' @param ... Additional arguments (ignored)
#'
#' @method print bcor
#' @export
print.bcor <- function(x, ...) {
  cat("BCOR file:", x$filename, "\n")
  cat("Number of SNPs:", x$nSNPs, "\n")
  cat("Number of samples:", x$nSamples, "\n")
  if (x$is_extended) {
    cat("Format: Extended BCOR (with diagonal values)\n")
  }
}

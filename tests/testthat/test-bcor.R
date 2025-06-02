test_that("read_bcor can open a BCOR file", {
  bcor_file <- system.file("extdata", "data.bcor", package = "rbcor")

  expect_true(file.exists(bcor_file))

  bcor <- read_bcor(bcor_file)

  expect_s3_class(bcor, "bcor")
  expect_type(bcor, "list")
  expect_true("ptr" %in% names(bcor))
  expect_true("filename" %in% names(bcor))
  expect_true("nSNPs" %in% names(bcor))
  expect_true("nSamples" %in% names(bcor))

  expect_equal(bcor$filename, bcor_file)
  expect_type(bcor$nSNPs, "integer")
  expect_type(bcor$nSamples, "integer")

  # Check specific values from the test file
  expect_equal(bcor$nSNPs, 55L)
  expect_equal(bcor$nSamples, 5363L)
})

test_that("read_bcor fails with non-existent file", {
  expect_error(
    read_bcor("non_existent_file.bcor"),
    "File does not exist"
  )
})

test_that("get_meta returns valid metadata", {
  bcor_file <- system.file("extdata", "data.bcor", package = "rbcor")
  bcor <- read_bcor(bcor_file)

  meta <- bcor$get_meta()

  expect_s3_class(meta, "data.frame")
  expect_equal(nrow(meta), bcor$nSNPs)

  # Check required columns
  required_cols <- c("rsid", "position", "chromosome", "allele1", "allele2")
  expect_true(all(required_cols %in% names(meta)))

  # Check data types
  expect_type(meta$rsid, "character")
  expect_type(meta$position, "integer")
  expect_type(meta$chromosome, "character")
  expect_type(meta$allele1, "character")
  expect_type(meta$allele2, "character")

  # Check specific values
  # First SNP
  expect_equal(meta$rsid[1], "rs1")
  expect_equal(meta$position[1], 1L)
  expect_equal(meta$chromosome[1], "01")
  expect_equal(meta$allele1[1], "A")
  expect_equal(meta$allele2[1], "G")

  # Last SNP
  expect_equal(meta$rsid[bcor$nSNPs], "rs55")
  expect_equal(meta$position[bcor$nSNPs], 55L)
  expect_equal(meta$chromosome[bcor$nSNPs], "01")
  expect_equal(meta$allele1[bcor$nSNPs], "A")
  expect_equal(meta$allele2[bcor$nSNPs], "G")

  # Check some middle values
  expect_equal(meta$rsid[10], "rs10")
  expect_equal(meta$position[10], 10L)

  # Check that all chromosomes are "01"
  expect_true(all(meta$chromosome == "01"))

  # Check that all alleles are A/G
  expect_true(all(meta$allele1 == "A"))
  expect_true(all(meta$allele2 == "G"))

  # Check that positions are sequential
  expect_equal(meta$position, 1:bcor$nSNPs)
})

test_that("get_meta fails with invalid input", {
  # Test with non-bcor object
  expect_error(
    {
      fake_bcor <- list()
      fake_bcor$get_meta()
    },
    "attempt to apply non-function"
  )
})

test_that("read_corr reads full correlation matrix", {
  bcor_file <- system.file("extdata", "data.bcor", package = "rbcor")
  bcor <- read_bcor(bcor_file)

  # Read full matrix
  corr_full <- bcor$read_corr()

  expect_true(is.matrix(corr_full))
  expect_equal(nrow(corr_full), bcor$nSNPs)
  expect_equal(ncol(corr_full), bcor$nSNPs)

  # Check matrix properties
  # Should be symmetric
  expect_equal(corr_full, t(corr_full), tolerance = 1e-6)

  # Diagonal should be 1
  expect_equal(diag(corr_full), rep(1, bcor$nSNPs), tolerance = 1e-6)

  # Values should be between -1 and 1
  expect_true(all(corr_full >= -1 & corr_full <= 1))
})

test_that("read_corr matches expected LD file", {
  bcor_file <- system.file("extdata", "data.bcor", package = "rbcor")
  ld_file <- system.file("extdata", "data.ld", package = "rbcor")

  bcor <- read_bcor(bcor_file)
  corr_bcor <- bcor$read_corr()

  # Read LD file
  ld_data <- as.matrix(read.table(ld_file))

  # Compare dimensions
  expect_equal(dim(corr_bcor), dim(ld_data))

  # Compare values (allowing for small numerical differences)
  # Remove dimnames to avoid mismatch
  expect_equal(unname(corr_bcor), unname(ld_data), tolerance = 1e-5)
})

test_that("read_corr can read subsets of SNPs", {
  bcor_file <- system.file("extdata", "data.bcor", package = "rbcor")
  bcor <- read_bcor(bcor_file)

  # Read subset
  snps_subset <- c(1, 5, 10)
  corr_subset <- bcor$read_corr(snps = snps_subset)

  expect_true(is.matrix(corr_subset))
  expect_equal(nrow(corr_subset), bcor$nSNPs)
  expect_equal(ncol(corr_subset), length(snps_subset))

  # Read full matrix for comparison
  corr_full <- bcor$read_corr()

  # Check that subset matches corresponding columns from full matrix
  expect_equal(corr_subset, corr_full[, snps_subset], tolerance = 1e-6)
})

test_that("read_corr can read rectangular submatrices", {
  bcor_file <- system.file("extdata", "data.bcor", package = "rbcor")
  bcor <- read_bcor(bcor_file)

  # Read rectangular submatrix
  snps1 <- c(1, 3, 5)
  snps2 <- c(2, 4, 6, 8)
  corr_rect <- bcor$read_corr(snps = snps1, snps2 = snps2)

  expect_true(is.matrix(corr_rect))
  expect_equal(nrow(corr_rect), length(snps1))
  expect_equal(ncol(corr_rect), length(snps2))

  # Read full matrix for comparison
  corr_full <- bcor$read_corr()

  # Check that rectangular subset matches
  expect_equal(corr_rect, corr_full[snps1, snps2], tolerance = 1e-6)
})

test_that("read_corr validates SNP indices", {
  bcor_file <- system.file("extdata", "data.bcor", package = "rbcor")
  bcor <- read_bcor(bcor_file)

  # Test invalid indices
  expect_error(
    bcor$read_corr(snps = c(0, 1)),
    "SNP indices must be between 1 and"
  )
  expect_error(
    bcor$read_corr(snps = c(1, bcor$nSNPs + 1)),
    "SNP indices must be between 1 and"
  )
  expect_error(
    bcor$read_corr(snps2 = c(0, 1)),
    "SNP indices must be between 1 and"
  )

  # Test missing snps with snps2
  expect_error(
    bcor$read_corr(snps2 = c(1, 2)),
    "Must specify snps when using snps2"
  )
})

test_that("read_corr returns sparse matrices", {
  bcor_file <- system.file("extdata", "data.bcor", package = "rbcor")
  bcor <- read_bcor(bcor_file)

  # Read sparse matrix with threshold
  threshold <- 0.5
  corr_sparse <- bcor$read_corr(sparse = TRUE, threshold = threshold)

  expect_s4_class(corr_sparse, "sparseMatrix")
  expect_equal(nrow(corr_sparse), bcor$nSNPs)
  expect_equal(ncol(corr_sparse), bcor$nSNPs)

  # Convert to dense for checking
  corr_dense <- as.matrix(corr_sparse)

  # Check that values below threshold are zero
  small_values <- abs(corr_dense) < threshold & corr_dense != 0
  expect_equal(sum(small_values), 0)

  # Check diagonal is preserved
  expect_equal(diag(corr_dense), rep(1, bcor$nSNPs))
})

test_that("read_corr sparse subset works correctly", {
  bcor_file <- system.file("extdata", "data.bcor", package = "rbcor")
  bcor <- read_bcor(bcor_file)

  snps_subset <- c(1, 5, 10, 15)
  threshold <- 0.3
  corr_sparse_subset <- bcor$read_corr(
    snps = snps_subset,
    sparse = TRUE, threshold = threshold
  )

  expect_s4_class(corr_sparse_subset, "sparseMatrix")
  expect_equal(nrow(corr_sparse_subset), bcor$nSNPs)
  expect_equal(ncol(corr_sparse_subset), length(snps_subset))

  # Read dense version for comparison
  corr_dense_subset <- bcor$read_corr(snps = snps_subset)

  # Convert sparse to dense
  corr_sparse_dense <- as.matrix(corr_sparse_subset)

  # Check that thresholding was applied correctly
  for (i in 1:nrow(corr_dense_subset)) {
    for (j in 1:ncol(corr_dense_subset)) {
      if (abs(corr_dense_subset[i, j]) <= threshold) {
        expect_equal(corr_sparse_dense[i, j], 0)
      } else {
        expect_equal(corr_sparse_dense[i, j], corr_dense_subset[i, j],
          tolerance = 1e-6
        )
      }
    }
  }
})

test_that("print.bcor works correctly", {
  bcor_file <- system.file("extdata", "data.bcor", package = "rbcor")
  bcor <- read_bcor(bcor_file)

  output <- capture.output(print(bcor))

  expect_length(output, 3)
  expect_match(output[1], "BCOR file:")
  expect_match(output[2], "Number of SNPs:")
  expect_match(output[3], "Number of samples:")
})

test_that("read_corr handles edge cases", {
  bcor_file <- system.file("extdata", "data.bcor", package = "rbcor")
  bcor <- read_bcor(bcor_file)

  # Single SNP
  corr_single <- bcor$read_corr(snps = 1)
  expect_equal(dim(corr_single), c(bcor$nSNPs, 1))
  expect_equal(corr_single[1, 1], 1)

  # All SNPs explicitly
  all_snps <- 1:bcor$nSNPs
  corr_all <- bcor$read_corr(snps = all_snps)
  corr_full <- bcor$read_corr()
  expect_equal(corr_all, corr_full)

  # Rectangular with single SNP
  corr_rect_single <- bcor$read_corr(snps = 1, snps2 = c(2, 3, 4))
  expect_equal(dim(corr_rect_single), c(1, 3))
})

test_that("packed threshold configuration works correctly", {
  bcor_file <- system.file("extdata", "data.bcor", package = "rbcor")

  # Test default threshold
  bcor_default <- read_bcor(bcor_file)
  expect_equal(bcor_default$packed_threshold, 1000L)

  # Test custom threshold
  bcor_custom <- read_bcor(bcor_file, packed_threshold = 50)
  expect_equal(bcor_custom$packed_threshold, 50L)
})

test_that("read_corr with packed=TRUE creates dspMatrix when above threshold", {
  bcor_file <- system.file("extdata", "data.bcor", package = "rbcor")
  bcor <- read_bcor(bcor_file, packed_threshold = 10) # Low threshold for testing

  # Should create dspMatrix since 55x55 >= 10
  corr_packed <- bcor$read_corr(packed = TRUE)

  expect_s4_class(corr_packed, "dspMatrix")
  expect_equal(dim(corr_packed), c(bcor$nSNPs, bcor$nSNPs))

  # Check it's symmetric packed storage
  expect_equal(corr_packed@uplo, "U")

  # Verify values match regular matrix
  corr_regular <- bcor$read_corr(packed = FALSE)
  expect_equal(as.matrix(corr_packed), corr_regular, tolerance = 1e-10)
})

test_that("read_corr with packed=TRUE falls back to regular matrix below threshold", {
  bcor_file <- system.file("extdata", "data.bcor", package = "rbcor")
  bcor <- read_bcor(bcor_file) # Default threshold 1000

  # Should return regular matrix since 55x55 < 1000
  expect_warning(
    corr_result <- bcor$read_corr(packed = TRUE),
    "Matrix size \\(55x55\\) below packed threshold \\(1000\\)\\. Returning regular matrix\\."
  )

  expect_true(is.matrix(corr_result))
  expect_false(inherits(corr_result, "dspMatrix"))
})

test_that("packed storage works with SNP subsets", {
  bcor_file <- system.file("extdata", "data.bcor", package = "rbcor")
  bcor <- read_bcor(bcor_file, packed_threshold = 10)

  # Test subset that meets threshold
  snps_subset <- 1:15 # 15x15 >= 10
  corr_packed_subset <- bcor$read_corr(snps = snps_subset, packed = TRUE)

  expect_s4_class(corr_packed_subset, "dspMatrix")
  expect_equal(dim(corr_packed_subset), c(length(snps_subset), length(snps_subset)))

  # Verify values against the rectangular subset using snps2 for square matrix
  corr_regular_subset <- bcor$read_corr(snps = snps_subset, snps2 = snps_subset, packed = FALSE)
  expect_equal(as.matrix(corr_packed_subset), corr_regular_subset, tolerance = 1e-10)
})

test_that("packed matrix supports matrix operations", {
  bcor_file <- system.file("extdata", "data.bcor", package = "rbcor")
  bcor <- read_bcor(bcor_file, packed_threshold = 10)

  corr_packed <- bcor$read_corr(packed = TRUE)
  corr_regular <- bcor$read_corr(packed = FALSE)

  # Test indexing
  expect_equal(corr_packed[1, 5], corr_regular[1, 5])
  expect_equal(corr_packed[10, 20], corr_regular[10, 20])

  # Test matrix multiplication
  x <- rep(1, ncol(corr_packed))
  result_packed <- corr_packed %*% x
  result_regular <- corr_regular %*% x
  expect_equal(as.vector(result_packed), as.vector(result_regular))
})

test_that("sparse and packed conflict handling", {
  bcor_file <- system.file("extdata", "data.bcor", package = "rbcor")
  bcor <- read_bcor(bcor_file, packed_threshold = 10)

  # Should warn and use sparse
  expect_warning(
    corr_result <- bcor$read_corr(sparse = TRUE, packed = TRUE, threshold = 0.5),
    "Both sparse=TRUE and packed=TRUE specified\\. Using sparse=TRUE\\."
  )

  expect_s4_class(corr_result, "sparseMatrix")
})

test_that("rectangular matrices cannot use packed storage", {
  bcor_file <- system.file("extdata", "data.bcor", package = "rbcor")
  bcor <- read_bcor(bcor_file, packed_threshold = 10)

  # Should warn and return regular matrix
  expect_warning(
    corr_result <- bcor$read_corr(snps = 1:5, snps2 = 6:10, packed = TRUE),
    "Rectangular matrices cannot use packed storage\\. Returning regular matrix\\."
  )

  expect_true(is.matrix(corr_result))
  expect_equal(dim(corr_result), c(5, 5))
})

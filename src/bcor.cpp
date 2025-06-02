// [[Rcpp::depends(RcppArmadillo)]]
#include "bcor.h"
#include <Rcpp.h>
#include <cmath>
#include <stdexcept>

// Define static constexpr for linking
constexpr size_t Bcor::BUFFER_SIZE;

Bcor::Bcor(const std::string &filename, bool read_header) : fname(filename)
{
    fh.open(fname, std::ios::binary);
    if (!fh.is_open())
    {
        throw std::runtime_error("Cannot open file: " + fname +
                                 " (check file permissions and path)");
    }

    try
    {
        validateFileIntegrity();
        readHeader();
        if (read_header)
        {
            readMeta();
        }
    }
    catch (...)
    {
        if (fh.is_open())
            fh.close();
        throw;
    }
}

Bcor::~Bcor()
{
    if (fh.is_open())
    {
        fh.close();
    }
}

void Bcor::validateFileIntegrity() const
{
    fh.seekg(0, std::ios::end);
    uint64_t actual_size = fh.tellg();
    fh.seekg(0, std::ios::beg);

    if (actual_size < 32)
    { // Minimum header size
        throw std::runtime_error("File '" + fname + "' is too small to be a valid bcor file");
    }
}

void Bcor::checkFilePosition(const std::string &operation) const
{
    if (fh.fail())
    {
        throw std::runtime_error("File I/O error during " + operation +
                                 " in file '" + fname + "'");
    }
}

void Bcor::readHeader()
{
    char magic[8];
    fh.read(magic, 7);
    magic[7] = '\0';
    checkFilePosition("header read");

    if (std::string(magic) != "bcor1.1")
    {
        throw std::runtime_error("File '" + fname + "' is not from LDStore!");
    }

    fh.read(reinterpret_cast<char *>(&fsize), sizeof(uint64_t));
    checkFilePosition("file size read");
    fh.read(reinterpret_cast<char *>(&nSamples), sizeof(uint32_t));
    checkFilePosition("nSamples read");
    fh.read(reinterpret_cast<char *>(&nSNPs), sizeof(uint32_t));
    checkFilePosition("nSNPs read");
    fh.read(reinterpret_cast<char *>(&compression), sizeof(uint8_t));
    checkFilePosition("compression read");
    fh.read(reinterpret_cast<char *>(&corr_block_offset), sizeof(uint64_t));
    checkFilePosition("correlation block offset read");

    // Validate compression type
    if (compression > 3)
    {
        throw std::runtime_error("Invalid compression type: " + std::to_string(compression));
    }
}

void Bcor::readMeta()
{
    rsid.resize(nSNPs);
    position.resize(nSNPs);
    chromosome.resize(nSNPs);
    allele1.resize(nSNPs);
    allele2.resize(nSNPs);

    for (uint32_t snp = 0; snp < nSNPs; ++snp)
    {
        uint32_t L_buffer, index;
        uint16_t L_rsid, L_chromosome;
        uint32_t L_allele1, L_allele2;

        fh.read(reinterpret_cast<char *>(&L_buffer), sizeof(uint32_t));
        fh.read(reinterpret_cast<char *>(&index), sizeof(uint32_t));

        fh.read(reinterpret_cast<char *>(&L_rsid), sizeof(uint16_t));
        rsid[snp].resize(L_rsid);
        fh.read(&rsid[snp][0], L_rsid);

        fh.read(reinterpret_cast<char *>(&position[snp]), sizeof(uint32_t));

        fh.read(reinterpret_cast<char *>(&L_chromosome), sizeof(uint16_t));
        chromosome[snp].resize(L_chromosome);
        fh.read(&chromosome[snp][0], L_chromosome);

        fh.read(reinterpret_cast<char *>(&L_allele1), sizeof(uint32_t));
        allele1[snp].resize(L_allele1);
        fh.read(&allele1[snp][0], L_allele1);

        fh.read(reinterpret_cast<char *>(&L_allele2), sizeof(uint32_t));
        allele2[snp].resize(L_allele2);
        fh.read(&allele2[snp][0], L_allele2);

        if (L_buffer != (20 + L_rsid + L_chromosome + L_allele1 + L_allele2))
        {
            throw std::runtime_error("Meta information block corrupted");
        }
    }
}

// Optimized buffered I/O methods
template <typename T>
T Bcor::readBuffered(uint64_t offset) const
{
    if (offset < buffer_start_offset || offset + sizeof(T) > buffer_end_offset)
    {
        fillBuffer(offset);
    }

    size_t buffer_pos = offset - buffer_start_offset;
    T value;
    std::memcpy(&value, &read_buffer[buffer_pos], sizeof(T));
    return value;
}

void Bcor::fillBuffer(uint64_t target_offset) const
{
    buffer_start_offset = target_offset;
    fh.seekg(target_offset);

    size_t to_read = std::min(BUFFER_SIZE,
                              static_cast<size_t>(fsize - target_offset));
    read_buffer.resize(to_read);
    fh.read(reinterpret_cast<char *>(read_buffer.data()), to_read);
    buffer_end_offset = buffer_start_offset + fh.gcount();

    if (fh.gcount() != static_cast<std::streamsize>(to_read) && !fh.eof())
    {
        throw std::runtime_error("Failed to read expected amount of data from file");
    }
}

// Input validation methods
void Bcor::validateSNPIndices(const std::vector<int> &snps) const
{
    for (int snp : snps)
    {
        if (snp < 0 || snp >= static_cast<int>(nSNPs))
        {
            throw std::runtime_error("SNP index " + std::to_string(snp) +
                                     " out of range [0, " + std::to_string(nSNPs - 1) + "]");
        }
    }
}

template <typename Container>
void Bcor::validateAndConvertIndices(const Container &r_indices,
                                     std::vector<int> &cpp_indices) const
{
    cpp_indices.reserve(r_indices.size());
    for (auto idx : r_indices)
    {
        int cpp_idx = idx - 1; // Convert from 1-based to 0-based
        if (cpp_idx < 0 || cpp_idx >= static_cast<int>(nSNPs))
        {
            throw std::runtime_error("SNP index " + std::to_string(idx) +
                                     " out of range [1, " + std::to_string(nSNPs) + "]");
        }
        cpp_indices.push_back(cpp_idx);
    }
}

uint64_t Bcor::getIndex(uint32_t snp_x, uint32_t snp_y) const
{
    if (snp_x > snp_y)
    {
        return static_cast<uint64_t>(nSNPs) * (nSNPs - 1) / 2 -
               static_cast<uint64_t>(nSNPs - snp_y) * ((nSNPs - snp_y) - 1) / 2 +
               snp_x - snp_y - 1;
    }
    else
    {
        return static_cast<uint64_t>(nSNPs) * (nSNPs - 1) / 2 -
               static_cast<uint64_t>(nSNPs - snp_x) * ((nSNPs - snp_x) - 1) / 2 +
               snp_y - snp_x - 1;
    }
}

// Template-based optimized correlation reading
template <uint8_t comp_type>
double Bcor::readCorrPairTyped(uint32_t snp_x, uint32_t snp_y) const
{
    using Traits = CompressionTraits<comp_type>;

    uint64_t index = getIndex(snp_x, snp_y);
    uint64_t offset = corr_block_offset + index * Traits::bytes;

    typename Traits::type val = readBuffered<typename Traits::type>(offset);

    if (val == Traits::na_value)
    {
        return NA_REAL;
    }
    return std::ldexp(static_cast<double>(val), -1 * (8 * Traits::bytes - 2)) - 1.0;
}

// Runtime dispatcher for correlation reading
double Bcor::readCorrPair(uint32_t snp_x, uint32_t snp_y, bool seek) const
{
    switch (compression)
    {
    case 0:
        return readCorrPairTyped<0>(snp_x, snp_y);
    case 1:
        return readCorrPairTyped<1>(snp_x, snp_y);
    case 2:
        return readCorrPairTyped<2>(snp_x, snp_y);
    case 3:
        return readCorrPairTyped<3>(snp_x, snp_y);
    default:
        throw std::runtime_error("Unknown compression type");
    }
}

// Legacy methods (kept for compatibility)
uint64_t Bcor::getIntNA(uint8_t n_bytes) const
{
    switch (n_bytes)
    {
    case 1:
        return 208;
    case 2:
        return 53248;
    case 4:
        return 3489660928;
    case 8:
        return 14987979559889010688ULL;
    default:
        throw std::runtime_error("Only 1, 2, 4 and 8 bytes are supported!");
    }
}

double Bcor::convertIntToFloat(uint64_t x, uint8_t n_bytes) const
{
    if (x == getIntNA(n_bytes))
    {
        return NA_REAL;
    }
    else
    {
        return std::ldexp(static_cast<double>(x), -1 * (8 * n_bytes - 2));
    }
}

Rcpp::DataFrame Bcor::getMeta()
{
    if (rsid.empty())
    {
        fh.seekg(0);
        readHeader();
        readMeta();
    }

    // Convert position to IntegerVector
    Rcpp::IntegerVector pos_vec(position.size());
    for (size_t i = 0; i < position.size(); ++i)
    {
        pos_vec[i] = static_cast<int>(position[i]);
    }

    return Rcpp::DataFrame::create(
        Rcpp::Named("rsid") = rsid,
        Rcpp::Named("position") = pos_vec,
        Rcpp::Named("chromosome") = chromosome,
        Rcpp::Named("allele1") = allele1,
        Rcpp::Named("allele2") = allele2,
        Rcpp::_["stringsAsFactors"] = false);
}

arma::mat Bcor::readCorr(const std::vector<int> &snps)
{
    if (snps.empty())
    {
        // Read full matrix
        arma::mat corr(nSNPs, nSNPs, arma::fill::eye);

        // Clear buffer to ensure sequential reading works optimally
        buffer_start_offset = UINT64_MAX;
        buffer_end_offset = 0;

        for (uint32_t snp_x = 0; snp_x < nSNPs - 1; ++snp_x)
        {
            for (uint32_t snp_y = snp_x + 1; snp_y < nSNPs; ++snp_y)
            {
                double val = readCorrPair(snp_x, snp_y, false);
                corr(snp_y, snp_x) = val;
                corr(snp_x, snp_y) = val;
            }
        }
        return corr;
    }
    else
    {
        // Validate SNP indices using consolidated method
        validateSNPIndices(snps);

        // Read subset
        arma::mat corr(nSNPs, snps.size());

        for (size_t i = 0; i < snps.size(); ++i)
        {
            int snp_x = snps[i];
            for (uint32_t snp_y = 0; snp_y < nSNPs; ++snp_y)
            {
                if (snp_x == static_cast<int>(snp_y))
                {
                    corr(snp_y, i) = 1.0;
                }
                else
                {
                    corr(snp_y, i) = readCorrPair(snp_x, snp_y, true);
                }
            }
        }
        return corr;
    }
}

arma::mat Bcor::readCorr(const std::vector<int> &snps1, const std::vector<int> &snps2)
{
    if (snps1.empty() || snps2.empty())
    {
        throw std::runtime_error("Both SNP lists must be non-empty");
    }

    // Validate indices using consolidated methods
    validateSNPIndices(snps1);
    validateSNPIndices(snps2);

    // Pre-allocate result matrix
    arma::mat corr(snps1.size(), snps2.size());

    for (size_t i = 0; i < snps1.size(); ++i)
    {
        for (size_t j = 0; j < snps2.size(); ++j)
        {
            if (snps1[i] == snps2[j])
            {
                corr(i, j) = 1.0;
            }
            else
            {
                corr(i, j) = readCorrPair(snps1[i], snps2[j], true);
            }
        }
    }

    return corr;
}

arma::sp_mat Bcor::readCorrSparse(const std::vector<int> &snps, double threshold)
{
    arma::mat dense = readCorr(snps);

    // Convert to sparse matrix, keeping only values above threshold
    arma::sp_mat sparse(dense.n_rows, dense.n_cols);

    for (size_t i = 0; i < dense.n_rows; ++i)
    {
        for (size_t j = 0; j < dense.n_cols; ++j)
        {
            if (std::abs(dense(i, j)) > threshold)
            {
                sparse(i, j) = dense(i, j);
            }
        }
    }

    return sparse;
}

// Rcpp export functions

// [[Rcpp::export]]
Rcpp::List bcor_open(std::string filename, bool read_header = true)
{
    Bcor *bcor = new Bcor(filename, read_header);

    Rcpp::XPtr<Bcor> ptr(bcor, true);

    return Rcpp::List::create(
        Rcpp::Named("ptr") = ptr,
        Rcpp::Named("filename") = bcor->getFname(),
        Rcpp::Named("nSNPs") = static_cast<int>(bcor->getNumOfSNPs()),
        Rcpp::Named("nSamples") = static_cast<int>(bcor->getNumOfSamples()));
}

// [[Rcpp::export]]
Rcpp::DataFrame bcor_get_meta(SEXP bcor_ptr)
{
    Rcpp::XPtr<Bcor> bcor(bcor_ptr);
    return bcor->getMeta();
}

// [[Rcpp::export]]
arma::mat bcor_read_corr(SEXP bcor_ptr, Rcpp::Nullable<Rcpp::IntegerVector> snps = R_NilValue)
{
    Rcpp::XPtr<Bcor> bcor(bcor_ptr);

    if (snps.isNull())
    {
        std::vector<int> empty_vec;
        return bcor->readCorr(empty_vec);
    }
    else
    {
        std::vector<int> snps_cpp;
        bcor->validateAndConvertIndices(Rcpp::IntegerVector(snps), snps_cpp);
        return bcor->readCorr(snps_cpp);
    }
}

// [[Rcpp::export]]
arma::mat bcor_read_corr2(SEXP bcor_ptr, Rcpp::IntegerVector snps1, Rcpp::IntegerVector snps2)
{
    Rcpp::XPtr<Bcor> bcor(bcor_ptr);

    std::vector<int> snps1_cpp, snps2_cpp;
    bcor->validateAndConvertIndices(snps1, snps1_cpp);
    bcor->validateAndConvertIndices(snps2, snps2_cpp);

    return bcor->readCorr(snps1_cpp, snps2_cpp);
}

// [[Rcpp::export]]
arma::sp_mat bcor_read_corr_sparse(SEXP bcor_ptr, Rcpp::Nullable<Rcpp::IntegerVector> snps = R_NilValue,
                                   double threshold = 0.0)
{
    Rcpp::XPtr<Bcor> bcor(bcor_ptr);

    if (snps.isNull())
    {
        std::vector<int> empty_vec;
        return bcor->readCorrSparse(empty_vec, threshold);
    }
    else
    {
        std::vector<int> snps_cpp;
        bcor->validateAndConvertIndices(Rcpp::IntegerVector(snps), snps_cpp);
        return bcor->readCorrSparse(snps_cpp, threshold);
    }
}

// [[Rcpp::export]]
arma::sp_mat bcor_read_corr_sparse2(SEXP bcor_ptr, Rcpp::IntegerVector snps1,
                                    Rcpp::IntegerVector snps2, double threshold = 0.0)
{
    Rcpp::XPtr<Bcor> bcor(bcor_ptr);

    std::vector<int> snps1_cpp, snps2_cpp;
    bcor->validateAndConvertIndices(snps1, snps1_cpp);
    bcor->validateAndConvertIndices(snps2, snps2_cpp);

    return bcor->readCorrSparse(snps1_cpp, snps2_cpp, threshold);
}

// [[Rcpp::export]]
Rcpp::S4 bcor_read_corr_packed(SEXP bcor_ptr, Rcpp::Nullable<Rcpp::IntegerVector> snps = R_NilValue)
{
    Rcpp::XPtr<Bcor> bcor(bcor_ptr);

    arma::mat dense;
    if (snps.isNull())
    {
        std::vector<int> empty_vec;
        dense = bcor->readCorr(empty_vec);
    }
    else
    {
        std::vector<int> snps_cpp;
        bcor->validateAndConvertIndices(Rcpp::IntegerVector(snps), snps_cpp);
        // For packed storage, we need a square symmetric matrix
        // So we read the correlation matrix for the subset of SNPs vs themselves
        dense = bcor->readCorr(snps_cpp, snps_cpp);
    }

    // Create dspMatrix (symmetric packed matrix)
    int n = dense.n_rows;
    int packed_size = n * (n + 1) / 2;

    // Extract upper triangular part (column-major order for R)
    Rcpp::NumericVector x(packed_size);
    int idx = 0;
    for (int j = 0; j < n; ++j)
    {
        for (int i = 0; i <= j; ++i)
        {
            x[idx++] = dense(i, j);
        }
    }

    // Create dspMatrix S4 object
    Rcpp::S4 dsp("dspMatrix");
    dsp.slot("x") = x;
    dsp.slot("Dim") = Rcpp::IntegerVector::create(n, n);
    dsp.slot("uplo") = Rcpp::CharacterVector::create("U");

    return dsp;
}

arma::sp_mat Bcor::readCorrSparse(const std::vector<int> &snps1, const std::vector<int> &snps2, double threshold)
{
    arma::mat dense = readCorr(snps1, snps2);

    // Convert to sparse matrix, keeping only values above threshold
    arma::sp_mat sparse(dense.n_rows, dense.n_cols);

    for (size_t i = 0; i < dense.n_rows; ++i)
    {
        for (size_t j = 0; j < dense.n_cols; ++j)
        {
            if (std::abs(dense(i, j)) > threshold)
            {
                sparse(i, j) = dense(i, j);
            }
        }
    }

    return sparse;
}

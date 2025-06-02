#ifndef BCOR_H
#define BCOR_H

#include <RcppArmadillo.h>
#include <string>
#include <fstream>
#include <vector>
#include <cstring>
#include <algorithm>
#include <numeric>

// Template specializations for compression types
template <uint8_t compression_type>
struct CompressionTraits;

template <>
struct CompressionTraits<0>
{
    using type = uint16_t;
    static constexpr uint8_t bytes = 2;
    static constexpr uint64_t na_value = 53248;
};

template <>
struct CompressionTraits<1>
{
    using type = uint32_t;
    static constexpr uint8_t bytes = 4;
    static constexpr uint64_t na_value = 3489660928;
};

template <>
struct CompressionTraits<2>
{
    using type = uint64_t;
    static constexpr uint8_t bytes = 8;
    static constexpr uint64_t na_value = 14987979559889010688ULL;
};

template <>
struct CompressionTraits<3>
{
    using type = uint8_t;
    static constexpr uint8_t bytes = 1;
    static constexpr uint64_t na_value = 208;
};

class Bcor
{
private:
    std::string fname;
    mutable std::ifstream fh;

    // Header information
    uint64_t fsize;
    uint32_t nSamples;
    uint32_t nSNPs;
    uint8_t compression;
    uint64_t corr_block_offset;
    bool is_extended_format; // True if magic is "bcor1.x"

    // Meta information
    std::vector<std::string> rsid;
    std::vector<uint32_t> position;
    std::vector<std::string> chromosome;
    std::vector<std::string> allele1;
    std::vector<std::string> allele2;

    // Diagonal values for extended format
    mutable std::vector<double> diagonal_values;
    mutable bool diagonal_loaded = false;

    // Buffered I/O for performance
    static constexpr size_t BUFFER_SIZE = 64 * 1024; // 64KB buffer
    mutable std::vector<uint8_t> read_buffer;
    mutable uint64_t buffer_start_offset = UINT64_MAX;
    mutable uint64_t buffer_end_offset = 0;

    // Private methods
    void readHeader();
    void readMeta();
    void validateFileIntegrity() const;
    void checkFilePosition(const std::string &operation) const;

    // Optimized I/O methods
    template <typename T>
    T readBuffered(uint64_t offset) const;
    void fillBuffer(uint64_t target_offset) const;

    // Template-based correlation reading
    template <uint8_t comp_type>
    double readCorrPairTyped(uint32_t snp_x, uint32_t snp_y) const;

    template <uint8_t comp_type>
    double readDiagonalTyped(uint32_t snp_idx) const;

    void loadDiagonalValues() const;
    double getDiagonalValue(uint32_t snp_idx) const;

    uint64_t getIndex(uint32_t snp_x, uint32_t snp_y) const;
    double readCorrPair(uint32_t snp_x, uint32_t snp_y, bool seek = true) const;

    // Legacy methods (kept for compatibility)
    double convertIntToFloat(uint64_t x, uint8_t n_bytes) const;
    uint64_t getIntNA(uint8_t n_bytes) const;

public:
    Bcor(const std::string &filename, bool read_header = true);
    ~Bcor();

    // Getters
    std::string getFname() const { return fname; }
    uint64_t getFsize() const { return fsize; }
    uint32_t getNumOfSNPs() const { return nSNPs; }
    uint32_t getNumOfSamples() const { return nSamples; }
    bool isExtendedFormat() const { return is_extended_format; }

    // Get diagonal values (for extended format)
    std::vector<double> getDiagonalValues() const;

    // Get meta information as R data frame
    Rcpp::DataFrame getMeta();

    // Input validation (public for Rcpp exports)
    void validateSNPIndices(const std::vector<int> &snps) const;
    template <typename Container>
    void validateAndConvertIndices(const Container &r_indices,
                                   std::vector<int> &cpp_indices) const;

    // Read correlations
    arma::mat readCorr(const std::vector<int> &snps);
    arma::mat readCorr(const std::vector<int> &snps1, const std::vector<int> &snps2);
    arma::sp_mat readCorrSparse(const std::vector<int> &snps, double threshold = 0.0);
    arma::sp_mat readCorrSparse(const std::vector<int> &snps1, const std::vector<int> &snps2, double threshold = 0.0);
};

#endif // BCOR_H
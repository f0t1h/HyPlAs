/**
 * @file subprograms.hpp
 * @brief Declarations for HyPlAs C++ subprograms
 */

#ifndef HYPLAS_SUBPROGRAMS_HPP
#define HYPLAS_SUBPROGRAMS_HPP

#include <string>
#include <vector>

namespace hyplas {

// ---------------------------------------------------------------------------
// innotin
// ---------------------------------------------------------------------------

struct InnotinParams {
    std::string main_fastq;
    std::vector<std::string> subset_fastqs;
};

int run_innotin(const InnotinParams& params);

// ---------------------------------------------------------------------------
// select-missing-reads
// ---------------------------------------------------------------------------

struct SelectMissingReadsParams {
    std::string paf_path;
    std::string gaf_path;
    std::string fastq_path;
    std::string output_path;
    std::string prediction_path;
};

int run_select_missing_reads(const SelectMissingReadsParams& params);

// ---------------------------------------------------------------------------
// split-plasmid-reads
// ---------------------------------------------------------------------------

struct SplitPlasmidReadsParams {
    std::string gaf_path;
    std::string fastq_path;
    std::string prediction_path;
    std::string plasmid_out_path;
    std::string chr_out_path;
    std::string unknown_neither_path;
    std::string unknown_both_path;
    std::string unmapped_path;
    bool output_chromosomal = false;
};

int run_split_plasmid_reads(const SplitPlasmidReadsParams& params);

} // namespace hyplas

#endif // HYPLAS_SUBPROGRAMS_HPP

/**
 * @file hyplas_common.hpp
 * @brief Common utilities for HyPlAs C++ tools
 * 
 * This header contains shared code for TSV parsing, contig classification,
 * and safe I/O operations used across HyPlAs tools.
 */

#ifndef HYPLAS_COMMON_HPP
#define HYPLAS_COMMON_HPP

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <system_error>
#include <unordered_map>
#include <zlib.h>

#include "kseq.h"
#include "mio.hpp"
#include "mview.hpp"

// Initialize kseq for gzFile - must be outside namespace
KSEQ_INIT(gzFile, gzread)

namespace hyplas {

// Build-time configurable debug logging
#ifndef HYPLAS_DEBUG
#define HYPLAS_DEBUG 0
#endif

#if HYPLAS_DEBUG
#define HYPLAS_LOG(...) std::fprintf(stderr, __VA_ARGS__)
#else
#define HYPLAS_LOG(...) ((void)0)
#endif

// Type aliases using mview
using MmapView = mview::memory_view<mio::mmap_source>;
using LiteView = MmapView::liteview;

/**
 * @brief Classification type for contigs
 */
enum class ContigType {
    Plasmid,
    Chromosome,
    Unknown
};

/**
 * @brief Parse a TSV file containing contig classifications
 * 
 * @param path Path to the TSV file
 * @return Map of contig names to their classification types
 */
[[nodiscard]] inline std::unordered_map<std::string, ContigType> 
parse_plasmid_tsv(const std::string& path) {
    std::unordered_map<std::string, ContigType> plasmid_contigs;
    std::error_code error;
    
    mio::mmap_source tsv_mmap = mio::make_mmap_source(path, error);
    if (error) {
        std::fprintf(stderr, "Error mapping file: %s, exiting...\n", error.message().c_str());
        std::exit(EXIT_FAILURE);
    }

    MmapView tsv_view{&tsv_mmap};
    tsv_view.skip_next('\n');  // Skip header

    while (tsv_view.e < tsv_view.cap) {
        tsv_view.extend_until("\t\n");
        char next_char = tsv_view.at_end();
        ContigType contig_type = ContigType::Plasmid;
        std::string contig_name = static_cast<std::string>(tsv_view);

        if (next_char == '\t') {
            tsv_view.skip_next('\t');
            tsv_view.extend_until("\t\n");
            
            if (tsv_view == "plasmid") {
                contig_type = ContigType::Plasmid;
            } else if (tsv_view == "chromosome") {
                contig_type = ContigType::Chromosome;
            } else {
                contig_type = ContigType::Unknown;
            }
        }

        HYPLAS_LOG("%s\t%s\n", contig_name.c_str(), 
                   static_cast<std::string>(tsv_view).c_str());
        plasmid_contigs[contig_name] = contig_type;
        tsv_view.skip_next('\n');
    }

    return plasmid_contigs;
}

/**
 * @brief Safely open a gzip pipe for writing
 * 
 * @param output_path Path to the output file
 * @return FILE pointer to the pipe, or nullptr on error
 * 
 * @note The caller is responsible for calling pclose() on the returned FILE*.
 *       This function validates the path to prevent shell injection.
 */
[[nodiscard]] inline FILE* open_gzip_write_pipe(const std::string& output_path) {
    // Basic validation to prevent shell injection
    for (char c : output_path) {
        if (c == ';' || c == '|' || c == '&' || c == '$' || 
            c == '`' || c == '\n' || c == '\r') {
            std::fprintf(stderr, "Error: Invalid character in output path\n");
            return nullptr;
        }
    }
    
    // Use quotes around the path for safety
    std::string cmd = "gzip - > \"" + output_path + "\"";
    return popen(cmd.c_str(), "w");
}

} // namespace hyplas

#endif // HYPLAS_COMMON_HPP

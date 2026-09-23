/**
 * @file fastx.hpp
 * @brief FASTQ/alignment read operations and gzip output
 *
 * Houses the memory-mapped alignment parser shared by GAF and PAF inputs, an
 * RAII gzip writer, and the read-selection subprograms (innotin,
 * select-missing-reads, split-plasmid-reads) formerly shipped as hyplas-utils.
 */

#ifndef HYPLAS_FASTX_HPP
#define HYPLAS_FASTX_HPP

#include <filesystem>
#include <string>
#include <gtl/phmap.hpp>
#include <vector>

#include <zlib.h>

#include "mio.hpp"
#include "mview.hpp"

namespace hyplas {

/// Zero-copy cursor over an mmap'd (or in-memory) alignment buffer.
using MmapView = mview::byte_view;
using LiteView = mview::byte_view::span_ref;

/**
 * @brief One logical alignment record spanning consecutive lines that share a
 *        read id, as found in minigraph GAF and minimap2 PAF output.
 *
 * Construction advances @p view past the consumed lines. Column 6 (the path /
 * target field) is split on '<' and '>' into target tokens; a leading
 * orientation marker is skipped when present. Reaching end-of-input yields a
 * record with @ref is_terminal set.
 */
struct AlignmentEntry {
    bool is_terminal = false;
    LiteView id{};
    std::vector<LiteView> contigs;

    explicit AlignmentEntry(MmapView& view);
};

/**
 * @brief RAII writer for a gzip-compressed FASTQ stream.
 *
 * Replaces the previous popen("gzip - > file") pattern: no shell, no injection
 * surface. Default-constructed instances are closed; open() reports failure by
 * returning false.
 */
class GzWriter {
public:
    GzWriter() = default;

    ~GzWriter() { close(); }

    GzWriter(const GzWriter&) = delete;
    GzWriter& operator=(const GzWriter&) = delete;

    GzWriter(GzWriter&& other) noexcept : fp_(other.fp_) { other.fp_ = nullptr; }
    GzWriter& operator=(GzWriter&& other) noexcept {
        if (this != &other) {
            close();
            fp_ = other.fp_;
            other.fp_ = nullptr;
        }
        return *this;
    }

    /// @brief Open @p path for gzip writing. Returns false on failure.
    bool open(const std::filesystem::path& path) {
        close();
        fp_ = gzopen(path.string().c_str(), "wb");
        if (fp_) {
            // gzprintf uses an internal stack buffer (BUFSIZ, typically 8 KB).
            // Long reads (e.g. >4 KB) overflow it and are silently truncated.
            // Set a 4 MB buffer so reads up to ~1 MB are written correctly.
            gzbuffer(fp_, 1 << 22);
        }
        return fp_ != nullptr;
    }

    void close() {
        if (fp_) {
            gzclose(fp_);
            fp_ = nullptr;
        }
    }

    [[nodiscard]] bool is_open() const { return fp_ != nullptr; }

    /// @brief Write a FASTQ record. A null comment is emitted as empty.
    void write_fastq(const char* name, const char* comment,
                     const char* seq, const char* qual) {
        gzprintf(fp_, "@%s %s\n%s\n+\n%s\n",
                 name, comment ? comment : "", seq, qual);
    }

private:
    gzFile fp_ = nullptr;
};

// ---------------------------------------------------------------------------
// Read-selection subprogram parameters
// ---------------------------------------------------------------------------

struct InnotinParams {
    std::vector<std::string> main_fastqs;
    std::vector<std::string> subset_fastqs;
    std::string output_path;
};

struct SelectMissingReadsParams {
    std::string paf_path;
    std::vector<std::string> fastq_paths;
    std::string output_path;
};

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

/**
 * @brief Emit (as FASTA) the reads from @p InnotinParams::main_fastqs whose
 *        ids are NOT present in any subset FASTQ file. Inputs may be plain or
 *        gzip-compressed. Returns EXIT_SUCCESS/EXIT_FAILURE.
 */
[[nodiscard]] int innotin(const InnotinParams& params);

/**
 * @brief Emit reads from @p SelectMissingReadsParams::fastq_paths that appear
 *        as targets in the PAF. Missing PAF yields an empty gzip output and
 *        success.
 */
[[nodiscard]] int select_missing_reads(const SelectMissingReadsParams& params);

/**
 * @brief Split reads by their GAF alignment into plasmid / chromosome-and-
 *        plasmid / neither / unmapped gzip FASTQ streams.
 */
[[nodiscard]] int split_plasmid_reads(const SplitPlasmidReadsParams& params);

/**
 * @brief Copy the reads of @p source_fastq whose ids are in @p ids into a
 *        gzip FASTQ file. Returns 0 on success, 1 on I/O failure.
 */
[[nodiscard]] int write_reads_by_id(const std::filesystem::path& source_fastq,
                                    const gtl::flat_hash_set<std::string>& ids,
                                    const std::filesystem::path& output_fastq_gz);

} // namespace hyplas

#endif // HYPLAS_FASTX_HPP

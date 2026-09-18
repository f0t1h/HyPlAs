/**
 * @file fastx.cpp
 * @brief FASTQ/alignment read operations
 */

#include "fastx.hpp"

#include "error.hpp"
#include "contig_classification.hpp"

#include <cstdio>
#include <cstdlib>
#include <string>
#include <string_view>
#include <gtl/phmap.hpp>

#include <zlib.h>

#include "kseq.h"

// Initialize kseq for gzFile — must be at file scope, in exactly one TU.
KSEQ_INIT(gzFile, gzread)

namespace hyplas {

// ---------------------------------------------------------------------------
// Unified GAF/PAF record parser
// ---------------------------------------------------------------------------

AlignmentEntry::AlignmentEntry(MmapView& view) {
    if (view.e >= view.cap) {
        is_terminal = true;
        return;
    }

    bool first = true;
    do {
        MmapView row = view.line();
        row.extend_until('\t');

        LiteView current_id = static_cast<LiteView>(row);
        if (row != id) {
            if (!first) {
                return;
            }
            id = current_id;
        }
        first = false;

        row.skip_next_n('\t', 5);
        row.extend_until('\t');

        MmapView align_view = row.focus();

        // GAF paths begin with an orientation marker ('<' or '>'); PAF target
        // names do not. Skip the marker only when present so both formats work.
        if (align_view.s < align_view.cap &&
            (align_view[0] == '<' || align_view[0] == '>')) {
            align_view.skip_next("<>");
        }

        while (align_view.e < align_view.cap) {
            align_view.extend_until("<>");
            contigs.push_back(static_cast<LiteView>(align_view));
            align_view.skip_next("<>");
        }

        view.skip_next('\n');
    } while (view.e < view.cap);
}

// ---------------------------------------------------------------------------
// innotin: emit reads not present in any subset file (as FASTA)
// ---------------------------------------------------------------------------

int innotin(const InnotinParams& params) {
    gtl::flat_hash_set<std::string> subset_ids;

    // Collect read IDs from plain or gzip-compressed subset files.
    for (const auto& subset_path : params.subset_fastqs) {
        gzFile subset_fp = gzopen(subset_path.c_str(), "r");
        if (!subset_fp) {
            std::fprintf(stderr, "Error opening subset FASTQ file: %s\n",
                         subset_path.c_str());
            return EXIT_FAILURE;
        }

        kseq_t* seq = kseq_init(subset_fp);
        while (kseq_read(seq) >= 0) {
            subset_ids.emplace(seq->name.s);
        }
        kseq_destroy(seq);
        gzclose(subset_fp);
    }

    // Open output: file or stdout.
    FILE* out_fp = stdout;
    if (!params.output_path.empty()) {
        out_fp = std::fopen(params.output_path.c_str(), "w");
        if (!out_fp) {
            std::fprintf(stderr, "Error opening output file: %s\n",
                         params.output_path.c_str());
            return EXIT_FAILURE;
        }
    }

    // Stream each main input without materializing a concatenated FASTQ.
    for (const auto& main_path : params.main_fastqs) {
        gzFile main_fp = gzopen(main_path.c_str(), "r");
        if (!main_fp) {
            std::fprintf(stderr, "Error opening main FASTQ file: %s\n",
                         main_path.c_str());
            if (out_fp != stdout) std::fclose(out_fp);
            return EXIT_FAILURE;
        }

        kseq_t* seq = kseq_init(main_fp);
        while (kseq_read(seq) >= 0) {
            if (!subset_ids.contains(seq->name.s)) {
                std::fprintf(out_fp, ">%s", seq->name.s);
                if (seq->comment.s && seq->comment.l > 0) {
                    std::fprintf(out_fp, " %s", seq->comment.s);
                }
                std::fprintf(out_fp, "\n%s\n", seq->seq.s);
            }
        }
        kseq_destroy(seq);
        gzclose(main_fp);
    }

    if (out_fp != stdout) {
        std::fclose(out_fp);
    }
    return EXIT_SUCCESS;
}

// ---------------------------------------------------------------------------
// select-missing-reads: select reads that overlap with plasmid reads
// ---------------------------------------------------------------------------

int select_missing_reads(const SelectMissingReadsParams& params) {
    std::error_code error;

    mio::mmap_source paf_mmap = mio::make_mmap_source(params.paf_path, error);
    if (error) {
        std::fprintf(stderr, "Error mapping PAF file: %s\n", error.message().c_str());
        std::fprintf(stderr, "Creating empty output since there are no mappings\n");

        GzWriter empty;
        empty.open(params.output_path);  // best-effort empty gzip file
        return EXIT_SUCCESS;
    }

    GzWriter plasmid_out;
    if (!plasmid_out.open(params.output_path)) {
        std::fprintf(stderr, "Error opening output file\n");
        return EXIT_FAILURE;
    }

    MmapView paf_view{paf_mmap};
    gtl::flat_hash_map<std::string, int> reads_to_use;

    for (AlignmentEntry entry{paf_view}; !entry.is_terminal; entry = AlignmentEntry{paf_view}) {
        for (const LiteView& v : entry.contigs) {
            MmapView contig_view = paf_view.sub(v);
            reads_to_use[static_cast<std::string>(contig_view)] = 1;
        }
    }

    for (const auto& fastq_path : params.fastq_paths) {
        gzFile fastq_fp = gzopen(fastq_path.c_str(), "r");
        if (!fastq_fp) {
            std::fprintf(stderr, "Error opening FASTQ file: %s\n", fastq_path.c_str());
            return EXIT_FAILURE;
        }

        kseq_t* seq = kseq_init(fastq_fp);
        while (kseq_read(seq) >= 0) {
            auto it = reads_to_use.find(seq->name.s);
            if (it != reads_to_use.end() && it->second > 0) {
                plasmid_out.write_fastq(seq->name.s, seq->comment.s,
                                        seq->seq.s, seq->qual.s);
                --it->second;
            }
        }
        kseq_destroy(seq);
        gzclose(fastq_fp);
    }

    return EXIT_SUCCESS;
}

// ---------------------------------------------------------------------------
// split-plasmid-reads: split reads based on alignment to plasmid/chromosome
// ---------------------------------------------------------------------------

int split_plasmid_reads(const SplitPlasmidReadsParams& params) {
    std::error_code error;
    mio::mmap_source gaf_mmap = mio::make_mmap_source(params.gaf_path, error);
    if (error) {
        std::fprintf(stderr, "Error mapping GAF file: %s\n", error.message().c_str());
        return EXIT_FAILURE;
    }

    gzFile fastq_fp = gzopen(params.fastq_path.c_str(), "r");
    if (!fastq_fp) {
        std::fprintf(stderr, "Error opening FASTQ file: %s\n", params.fastq_path.c_str());
        return EXIT_FAILURE;
    }

    GzWriter plasmid_out;
    GzWriter unknown_neither_out;
    GzWriter unknown_both_out;
    GzWriter unmapped_out;
    if (!plasmid_out.open(params.plasmid_out_path) ||
        !unknown_neither_out.open(params.unknown_neither_path) ||
        !unknown_both_out.open(params.unknown_both_path) ||
        !unmapped_out.open(params.unmapped_path)) {
        std::fprintf(stderr, "Error opening output files\n");
        gzclose(fastq_fp);
        return EXIT_FAILURE;
    }

    GzWriter chr_out;
    if (params.output_chromosomal) {
        if (!chr_out.open(params.chr_out_path)) {
            std::fprintf(stderr, "Error opening chromosome output file\n");
            gzclose(fastq_fp);
            return EXIT_FAILURE;
        }
    }

    gtl::flat_hash_map<std::string, ContigType> plasmid_contigs;
    try {
        plasmid_contigs = parse_prediction_tsv(params.prediction_path);
    } catch (const HyplasError& e) {
        std::fprintf(stderr, "%s\n", e.what());
        gzclose(fastq_fp);
        return EXIT_FAILURE;
    }

    kseq_t* seq = kseq_init(fastq_fp);
    MmapView gaf_view{gaf_mmap};

    for (AlignmentEntry entry{gaf_view}; !entry.is_terminal; entry = AlignmentEntry{gaf_view}) {
        int64_t l = kseq_read(seq);

        MmapView gid = gaf_view.sub(entry.id);

        while (l >= 0) {
            if (gid != std::string_view(seq->name.s)) {
                unmapped_out.write_fastq(seq->name.s, seq->comment.s,
                                         seq->seq.s, seq->qual.s);
                l = kseq_read(seq);
            } else {
                break;
            }
        }

        int from_chromosome = 0;
        int from_plasmid = 0;

        for (const LiteView& v : entry.contigs) {
            MmapView contig_view = gaf_view.sub(v);
            std::string contig_name = static_cast<std::string>(contig_view);

            auto it = plasmid_contigs.find(contig_name);
            if (it != plasmid_contigs.end()) {
                if (it->second == ContigType::Chromosome) {
                    ++from_chromosome;
                } else if (it->second == ContigType::Plasmid) {
                    ++from_plasmid;
                }
            }
        }

        const char* comment = seq->comment.s ? seq->comment.s : "";

        if (from_chromosome == 0 && from_plasmid == 0) {
            unknown_neither_out.write_fastq(seq->name.s, comment, seq->seq.s, seq->qual.s);
        } else if (from_plasmid > 0 && from_chromosome == 0) {
            plasmid_out.write_fastq(seq->name.s, comment, seq->seq.s, seq->qual.s);
        } else if (from_chromosome > 0 && from_plasmid > 0) {
            unknown_both_out.write_fastq(seq->name.s, comment, seq->seq.s, seq->qual.s);
        } else if (params.output_chromosomal && from_chromosome > 0) {
            chr_out.write_fastq(seq->name.s, comment, seq->seq.s, seq->qual.s);
        }
    }

    kseq_destroy(seq);
    gzclose(fastq_fp);

    return EXIT_SUCCESS;
}

// ---------------------------------------------------------------------------
// write_reads_by_id: copy selected reads into a gzip FASTQ file
// ---------------------------------------------------------------------------

void write_reads_by_id(const std::filesystem::path& source_fastq,
                       const gtl::flat_hash_set<std::string>& ids,
                       const std::filesystem::path& output_fastq_gz) {
    gzFile fp = gzopen(source_fastq.string().c_str(), "r");
    if (!fp) {
        throw HyplasError("cannot open FASTQ: " + source_fastq.string());
    }

    GzWriter out;
    if (!out.open(output_fastq_gz)) {
        gzclose(fp);
        throw HyplasError("cannot open gzip output: " + output_fastq_gz.string());
    }

    kseq_t* seq = kseq_init(fp);
    while (kseq_read(seq) >= 0) {
        if (ids.count(seq->name.s)) {
            out.write_fastq(seq->name.s, seq->comment.s, seq->seq.s, seq->qual.s);
        }
    }
    kseq_destroy(seq);
    gzclose(fp);
}

} // namespace hyplas

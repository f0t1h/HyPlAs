/**
 * @file gfa.hpp
 * @brief GFA/FASTA graph transforms used by the HyPlAs pipeline
 *
 * These are pure functions over files. On unrecoverable errors they throw
 * HyplasError; the caller decides whether that triggers soft-fail recovery.
 */

#ifndef HYPLAS_GFA_HPP
#define HYPLAS_GFA_HPP

#include <cstddef>
#include <filesystem>
#include <ostream>
#include <string>
#include <gtl/phmap.hpp>
#include <vector>

namespace hyplas {

/**
 * @brief Rewrite a GFA, dropping empty segments and bypassing them with links.
 *
 * For every empty (sequence-less) segment, new 0M links are created that
 * connect each of its predecessors directly to each of its successors.
 */
void fix_gfa_empty_segments(const std::filesystem::path& input,
                            const std::filesystem::path& output);

/**
 * @brief Remove k-mer overlaps from a GFA (Unicycler-style asymmetric trim).
 *
 * SPAdes emits uniform overlaps (e.g. 53M for k=53) that minigraph cannot
 * consume. This trims segment sequences and rewrites all links as 0M.
 */
void remove_gfa_overlaps(const std::filesystem::path& input,
                         const std::filesystem::path& output);

/**
 * @brief Write segment sequences from a GFA as FASTA records.
 * @param min_length Skip segments shorter than this many bases.
 */
void extract_fasta_from_gfa(const std::filesystem::path& gfa,
                            const std::filesystem::path& fasta,
                            std::size_t min_length = 200);

/**
 * @brief Write a sub-GFA containing only the named segments and links between
 *        them. Path (P) lines are dropped; header/other lines pass through.
 */
void write_component_gfa(const std::vector<std::string>& segments,
                         const std::filesystem::path& source_gfa,
                         const std::filesystem::path& output_gfa);

/**
 * @brief Append FASTA records whose header contains "circular" to @p out.
 *
 * Records are keyed by contig name (first whitespace-delimited header token);
 * a name already present in @p written is skipped, and newly written names are
 * inserted into it.
 */
void append_circular_by_header(std::ostream& out,
                               const std::filesystem::path& fasta,
                               gtl::flat_hash_set<std::string>& written);

/**
 * @brief Append circular, plasmid-classified SR contigs to @p out as FASTA.
 *
 * A contig qualifies when it has a self-loop link in @p gfa_path and is
 * classified "plasmid" in @p prediction_tsv. The header gets a " circular"
 * suffix. Names already in @p written are skipped.
 */
void append_circular_sr_plasmids(std::ostream& out,
                                 const std::filesystem::path& gfa_path,
                                 const std::filesystem::path& fasta_path,
                                 const std::filesystem::path& prediction_tsv,
                                 gtl::flat_hash_set<std::string>& written);

} // namespace hyplas

#endif // HYPLAS_GFA_HPP

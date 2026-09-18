/**
 * @file gfa.cpp
 * @brief GFA/FASTA graph transform implementations
 */

#include "gfa.hpp"

#include "error.hpp"
#include "log.hpp"
#include "contig_classification.hpp"
#include "mview.hpp"
#include "stage.hpp"

#include <array>
#include <cstdint>
#include <fstream>
#include <string_view>
#include <utility>
#include <gtl/btree.hpp>
#include <gtl/phmap.hpp>

namespace hyplas {

void fix_gfa_empty_segments(const std::filesystem::path& input,
                            const std::filesystem::path& output) {
    std::string buf = read_file(input);
    std::ofstream out(output);
    if (!out) {
        throw HyplasError("cannot open GFA output for empty segment fixing: " + output.string());
    }

    gtl::flat_hash_set<std::string> empty_segments;
    gtl::flat_hash_map<std::string, std::vector<std::pair<std::string, std::string>>> incoming;
    gtl::flat_hash_map<std::string, std::vector<std::pair<std::string, std::string>>> outgoing;

    std::vector<std::string_view> deferred_lines;

    // First pass: identify empty segments and collect links
    mview::for_each_line(buf, [&](std::string_view line) {
        if (line.empty()) return;

        if (line[0] == 'S') {
            // Segment line: S<tab>name<tab>sequence
            std::array<std::string_view, 4> f{};
            std::size_t nf = mview::split_fields(line, '\t', f);
            std::string_view name = nf > 1 ? f[1] : std::string_view{};
            std::string_view seq  = nf > 2 ? f[2] : std::string_view{};
            if (seq.empty()) {
                empty_segments.insert(std::string(name));
            } else {
                out << line << '\n';
            }
        } else if (line[0] == 'L') {
            deferred_lines.push_back(line);
        } else {
            out << line << '\n';
        }
    });

    // Second pass: process links
    for (std::string_view lline : deferred_lines) {
        std::array<std::string_view, 6> f{};
        std::size_t nf = mview::split_fields(lline, '\t', f);
        std::string from_name  (nf > 1 ? f[1] : std::string_view{});
        std::string from_orient(nf > 2 ? f[2] : std::string_view{});
        std::string to_name    (nf > 3 ? f[3] : std::string_view{});
        std::string to_orient  (nf > 4 ? f[4] : std::string_view{});

        if (empty_segments.contains(from_name)) {
            outgoing[from_name].emplace_back(to_name, to_orient);
        } else if (empty_segments.contains(to_name)) {
            incoming[to_name].emplace_back(from_name, from_orient);
        } else {
            out << lline << '\n';
        }
    }

    // Create new links bypassing empty segments
    for (const auto& seg : empty_segments) {
        for (const auto& [in_name, in_orient] : incoming[seg]) {
            for (const auto& [out_name, out_orient] : outgoing[seg]) {
                out << "L\t" << in_name << '\t' << in_orient
                    << '\t' << out_name << '\t' << out_orient << "\t0M\n";
            }
        }
    }
}

void remove_gfa_overlaps(const std::filesystem::path& input,
                         const std::filesystem::path& output) {
    // Graph-aware overlap removal equivalent to Unicycler's remove_all_overlaps().
    // SPAdes uses uniform k-mer overlap (e.g., 53M for k=53).
    //
    // Algorithm:
    // 1. Parse GFA to get segments and links
    // 2. Build forward/reverse link maps
    // 3. Group edges into two sets based on constraints:
    //    - Complement edges must be in opposite groups
    //    - Edges connecting to same side of segment must be in same group
    // 4. Apply asymmetric trimming based on grouping
    // 5. Output GFA with 0M overlaps

    using SignedSeg = int64_t;  // Signed segment ID (negative = reverse complement)
    using Edge = std::pair<SignedSeg, SignedSeg>;

    struct Segment {
        std::string name;
        std::string sequence;
        std::string tags;
    };

    struct Link {
        std::string from_name;
        bool from_is_forward;
        std::string to_name;
        bool to_is_forward;
        std::string tags;  // Optional GFA tags (MQ, NM, RC, FC, KC, ID, etc.)
    };

    std::string buf = read_file(input);

    // Parse GFA
    gtl::flat_hash_map<std::string, int64_t> name_to_id;
    gtl::flat_hash_map<int64_t, std::string> id_to_name;
    std::vector<Segment> segments;
    std::vector<Link> links;
    std::vector<std::string_view> header_lines;
    int overlap = 0;
    int64_t next_id = 1;

    mview::for_each_line(buf, [&](std::string_view line) {
        if (line.empty()) return;

        if (line[0] == 'S') {
            std::array<std::string_view, 4> f{};
            std::size_t nf = mview::split_fields(line, '\t', f);
            std::string name(nf > 1 ? f[1] : std::string_view{});
            std::string seq (nf > 2 ? f[2] : std::string_view{});
            std::string tags(nf > 3 ? f[3] : std::string_view{});

            name_to_id[name] = next_id;
            id_to_name[next_id] = name;
            segments.emplace_back(name, seq, tags);
            next_id++;
        } else if (line[0] == 'L') {
            std::array<std::string_view, 7> f{};
            std::size_t nf = mview::split_fields(line, '\t', f);
            std::string from_name  (nf > 1 ? f[1] : std::string_view{});
            std::string from_orient(nf > 2 ? f[2] : std::string_view{});
            std::string to_name    (nf > 3 ? f[3] : std::string_view{});
            std::string to_orient  (nf > 4 ? f[4] : std::string_view{});
            std::string overlap_str(nf > 5 ? f[5] : std::string_view{});
            std::string link_tags  (nf > 6 ? f[6] : std::string_view{});

            links.emplace_back(Link{from_name, from_orient == "+", to_name, to_orient == "+", link_tags});

            // Parse overlap value
            if (overlap == 0 && !overlap_str.empty() && overlap_str.back() == 'M') {
                overlap_str.pop_back();
                try {
                    overlap = std::stoi(overlap_str);
                } catch (...) {}
            }
        } else if (line[0] == 'H') {
            header_lines.push_back(line);
        }
        // Drop P (path) lines and anything else — P lines reference full-assembly
        // paths that are invalid for a component sub-GFA and break downstream tools.
    });

    if (overlap == 0) {
        log("INFO", "GFA has no overlaps - writing filtered S/L lines");
        std::ofstream out(output);
        if (!out) {
            throw HyplasError("cannot create output GFA file: " + output.string());
        }
        for (const auto& h : header_lines) out << h << '\n';
        for (const auto& seg : segments) {
            out << "S\t" << seg.name << '\t' << seg.sequence;
            if (!seg.tags.empty()) out << '\t' << seg.tags;
            out << '\n';
        }
        for (const auto& link : links) {
            out << "L\t" << link.from_name << '\t' << (link.from_is_forward ? '+' : '-')
                << '\t' << link.to_name << '\t' << (link.to_is_forward ? '+' : '-')
                << "\t0M";
            if (!link.tags.empty()) out << '\t' << link.tags;
            out << '\n';
        }
        return;
    }

    log("INFO", "Removing " + std::to_string(overlap) + "bp overlaps from GFA");

    // Build edge set and link maps
    // Edge: (signed_from, signed_to) where negative means reverse complement
    gtl::btree_set<Edge> all_edges;
    gtl::flat_hash_map<SignedSeg, std::vector<SignedSeg>> forward_links;  // downstream
    gtl::flat_hash_map<SignedSeg, std::vector<SignedSeg>> reverse_links;  // upstream

    for (const auto& link : links) {
        auto from_it = name_to_id.find(link.from_name);
        auto to_it = name_to_id.find(link.to_name);
        if (from_it == name_to_id.end() || to_it == name_to_id.end()) continue;

        SignedSeg from_seg = link.from_is_forward ? from_it->second : -from_it->second;
        SignedSeg to_seg = link.to_is_forward ? to_it->second : -to_it->second;

        all_edges.insert({from_seg, to_seg});
        all_edges.insert({-to_seg, -from_seg});  // Complement edge

        forward_links[from_seg].emplace_back(to_seg);
        reverse_links[to_seg].emplace_back(from_seg);
        forward_links[-to_seg].emplace_back(-from_seg);
        reverse_links[-from_seg].emplace_back(-to_seg);
    }

    // Calculate trim amounts
    int large_half = (overlap + 1) / 2;
    int small_half = overlap / 2;

    // Build constraint maps using gtl::btree_map (Edge is comparable via std::pair)
    gtl::btree_map<Edge, gtl::btree_set<Edge>> must_match, must_differ;

    // Constraint 1: Complement edges must be in opposite groups
    for (const auto& edge : all_edges) {
        Edge rev_edge = {-edge.second, -edge.first};
        must_differ[edge].insert(rev_edge);
        must_differ[rev_edge].insert(edge);
    }

    // Constraint 2: Edges connecting to same side of segment must be in same group
    for (const auto& [seg, downstream] : forward_links) {
        if (downstream.size() > 1) {
            Edge first_edge = {seg, downstream[0]};
            for (size_t i = 1; i < downstream.size(); ++i) {
                Edge other_edge = {seg, downstream[i]};
                must_match[first_edge].insert(other_edge);
                must_match[other_edge].insert(first_edge);
            }
        }
    }
    for (const auto& [seg, upstream] : reverse_links) {
        if (upstream.size() > 1) {
            Edge first_edge = {upstream[0], seg};
            for (size_t i = 1; i < upstream.size(); ++i) {
                Edge other_edge = {upstream[i], seg};
                must_match[first_edge].insert(other_edge);
                must_match[other_edge].insert(first_edge);
            }
        }
    }

    // Constraint 3: Small segments (length == overlap) can't have large trim on both sides
    for (const auto& seg : segments) {
        if (static_cast<int>(seg.sequence.length()) == overlap) {
            auto it = name_to_id.find(seg.name);
            if (it == name_to_id.end()) continue;
            SignedSeg seg_id = it->second;

            for (SignedSeg signed_seg : {seg_id, -seg_id}) {
                auto fwd_it = forward_links.find(signed_seg);
                auto rev_it = reverse_links.find(signed_seg);
                if (fwd_it != forward_links.end() && rev_it != reverse_links.end()) {
                    for (SignedSeg downstream : fwd_it->second) {
                        Edge fwd_edge = {signed_seg, downstream};
                        for (SignedSeg upstream : rev_it->second) {
                            Edge rev_edge = {upstream, signed_seg};
                            must_match[fwd_edge].insert(rev_edge);
                            must_match[rev_edge].insert(fwd_edge);
                        }
                    }
                }
            }
        }
    }

    // Group edges using union-find style propagation
    gtl::btree_set<Edge> group_1, group_2;

    for (const auto& edge : all_edges) {
        if (group_1.contains(edge) || group_2.contains(edge)) {
            continue;
        }

        gtl::btree_set<Edge> new_group_1;
        gtl::btree_set<Edge> new_group_2;
        new_group_1.insert(edge);

        bool changed = true;
        while (changed) {
            changed = false;
            gtl::btree_set<Edge> temp_1;
            gtl::btree_set<Edge> temp_2;

            for (const auto& g1_edge : new_group_1) {
                for (const auto& match_edge : must_match[g1_edge]) {
                    if (!new_group_1.contains(match_edge)) {
                        temp_1.insert(match_edge);
                    }
                }
                for (const auto& differ_edge : must_differ[g1_edge]) {
                    if (!new_group_2.contains(differ_edge)) {
                        temp_2.insert(differ_edge);
                    }
                }
            }
            for (const auto& g2_edge : new_group_2) {
                for (const auto& match_edge : must_match[g2_edge]) {
                    if (!new_group_2.contains(match_edge)) {
                        temp_2.insert(match_edge);
                    }
                }
                for (const auto& differ_edge : must_differ[g2_edge]) {
                    if (!new_group_1.contains(differ_edge)) {
                        temp_1.insert(differ_edge);
                    }
                }
            }

            if (!temp_1.empty() || !temp_2.empty()) {
                changed = true;
                new_group_1.insert(temp_1.begin(), temp_1.end());
                new_group_2.insert(temp_2.begin(), temp_2.end());
            }
        }

        group_1.insert(new_group_1.begin(), new_group_1.end());
        group_2.insert(new_group_2.begin(), new_group_2.end());
    }

    // Determine trim amounts per segment
    // Group 1: trim more from end of start segment
    // Group 2: trim more from start of end segment
    gtl::btree_set<int64_t> large_trim_end;
    gtl::btree_set<int64_t> large_trim_start;

    for (const auto& edge : group_1) {
        SignedSeg start_seg = edge.first;
        if (start_seg > 0) {
            large_trim_end.insert(start_seg);
        } else {
            large_trim_start.insert(-start_seg);
        }
    }
    for (const auto& edge : group_2) {
        SignedSeg end_seg = edge.second;
        if (end_seg > 0) {
            large_trim_start.insert(end_seg);
        } else {
            large_trim_end.insert(-end_seg);
        }
    }

    // Write output GFA
    std::ofstream out(output);
    if (!out) {
        throw HyplasError("cannot create output GFA file: " + output.string());
    }

    // Write header lines
    for (const auto& h : header_lines) {
        out << h << '\n';
    }

    // Write segments with trimmed sequences
    for (const auto& seg : segments) {
        auto it = name_to_id.find(seg.name);
        int64_t seg_id = (it != name_to_id.end()) ? it->second : 0;

        int start_trim = large_trim_start.contains(seg_id) ? large_half : small_half;
        int end_trim = large_trim_end.contains(seg_id) ? large_half : small_half;

        std::string trimmed_seq = seg.sequence;
        int seq_len = static_cast<int>(trimmed_seq.length());
        if (seq_len > start_trim + end_trim) {
            trimmed_seq = trimmed_seq.substr(
                static_cast<size_t>(start_trim),
                static_cast<size_t>(seq_len - start_trim - end_trim));
        }

        out << "S\t" << seg.name << '\t' << trimmed_seq;
        if (!seg.tags.empty()) {
            out << '\t' << seg.tags;
        }
        out << '\n';
    }

    // Write links with 0M overlap (preserving optional tags)
    for (const auto& link : links) {
        out << "L\t" << link.from_name << '\t' << (link.from_is_forward ? '+' : '-')
            << '\t' << link.to_name << '\t' << (link.to_is_forward ? '+' : '-')
            << "\t0M";
        if (!link.tags.empty()) {
            out << '\t' << link.tags;
        }
        out << '\n';
    }
}

void extract_fasta_from_gfa(const std::filesystem::path& gfa,
                            const std::filesystem::path& fasta,
                            std::size_t min_length) {
    std::string buf = read_file(gfa);
    std::ofstream out(fasta);
    if (!out) {
        throw HyplasError("cannot open FASTA output: " + fasta.string());
    }

    mview::for_each_line(buf, [&](std::string_view line) {
        if (line.empty() || line[0] != 'S') return;

        std::array<std::string_view, 4> f{};
        std::size_t nf = mview::split_fields(line, '\t', f);
        std::string_view name = nf > 1 ? f[1] : std::string_view{};
        std::string_view seq  = nf > 2 ? f[2] : std::string_view{};

        if (seq.length() >= min_length) {
            out << '>' << name << '\n' << seq << '\n';
        }
    });
}

void write_component_gfa(const std::vector<std::string>& segments,
                         const std::filesystem::path& source_gfa,
                         const std::filesystem::path& output_gfa) {
    gtl::flat_hash_set<std::string> segs(segments.begin(), segments.end());

    std::string buf = read_file(source_gfa);
    std::ofstream out(output_gfa);
    if (!out) {
        throw HyplasError("cannot open GFA output for component extraction: " + output_gfa.string());
    }

    mview::for_each_line(buf, [&](std::string_view line) {
        if (line.empty()) { out << '\n'; return; }
        if (line[0] == 'S') {
            std::array<std::string_view, 3> f{};
            std::size_t nf = mview::split_fields(line, '\t', f);
            std::string name(nf > 1 ? f[1] : std::string_view{});
            if (segs.count(name)) out << line << '\n';
        } else if (line[0] == 'L') {
            std::array<std::string_view, 5> f{};
            std::size_t nf = mview::split_fields(line, '\t', f);
            std::string from(nf > 1 ? f[1] : std::string_view{});
            std::string to  (nf > 3 ? f[3] : std::string_view{});
            if (segs.count(from) && segs.count(to)) out << line << '\n';
        } else if (line[0] == 'P') {
            // Drop P (path) lines — they reference the full assembly and are
            // invalid for a component sub-GFA
        } else {
            // H lines and others — pass through
            out << line << '\n';
        }
    });
}

void append_circular_by_header(std::ostream& out,
                               const std::filesystem::path& fasta,
                               gtl::flat_hash_set<std::string>& written) {
    std::string buf = read_file(fasta);

    std::string header;
    std::string seq;
    bool is_circular = false;

    auto flush = [&]() {
        if (!is_circular || seq.empty()) return;
        std::string name = header.substr(1);
        name = name.substr(0, name.find(' '));
        if (written.insert(name).second) {
            out << header << '\n' << seq << '\n';
        }
    };

    mview::for_each_line(buf, [&](std::string_view line) {
        if (line.empty()) return;
        if (line[0] == '>') {
            flush();
            header = std::string(line);
            seq.clear();
            is_circular = (line.find("circular") != std::string_view::npos);
        } else {
            seq += line;
        }
    });
    flush();
}

void append_circular_sr_plasmids(std::ostream& out,
                                 const std::filesystem::path& gfa_path,
                                 const std::filesystem::path& fasta_path,
                                 const std::filesystem::path& prediction_tsv,
                                 gtl::flat_hash_set<std::string>& written) {
    // Circular contigs: those with a self-loop link (from == to) in the GFA.
    gtl::flat_hash_set<std::string> circular_contigs;
    {
        std::string gfa_buf = read_file(gfa_path);
        mview::for_each_line(gfa_buf, [&](std::string_view line) {
            if (line.empty() || line[0] != 'L') return;
            std::array<std::string_view, 5> f{};
            std::size_t nf = mview::split_fields(line, '\t', f);
            std::string_view from = nf > 1 ? f[1] : std::string_view{};
            std::string_view to   = nf > 3 ? f[3] : std::string_view{};
            if (from == to) circular_contigs.insert(std::string(from));
        });
    }
    if (circular_contigs.empty()) return;

    gtl::flat_hash_set<std::string> plasmid_contigs =
        plasmid_names(parse_prediction_tsv(prediction_tsv));

    std::string fasta_buf = read_file(fasta_path);

    std::string current_name;
    std::string current_seq;
    auto write_if_match = [&]() {
        if (current_name.empty() || current_seq.empty()) return;
        std::string contig_name = current_name.substr(0, current_name.find(' '));
        if (!written.count(contig_name)
            && circular_contigs.contains(contig_name)
            && plasmid_contigs.contains(contig_name)) {
            out << '>' << current_name << " circular\n" << current_seq << '\n';
            written.insert(contig_name);
        }
    };
    mview::for_each_line(fasta_buf, [&](std::string_view line) {
        if (line.empty()) return;
        if (line[0] == '>') {
            write_if_match();
            current_name = std::string(line.substr(1));
            current_seq.clear();
        } else {
            current_seq += line;
        }
    });
    write_if_match();
}

} // namespace hyplas

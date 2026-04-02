/**
 * @file pipeline.cpp
 * @brief HyPlAs pipeline stage implementations
 */

#include "pipeline.hpp"
#include "file_utils.hpp"
#include "process.hpp"
#include "hyplas_common.hpp"

#include <algorithm>
#include <cstdio>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <unordered_map>
#include <unordered_set>

namespace hyplas {

// ============================================================================
// Internal parameter structs for helper functions
// ============================================================================

struct InnotinParams {
    std::string main_fastq;
    std::vector<std::string> subset_fastqs;
    std::string output_path;
};

struct SelectMissingReadsParams {
    std::string paf_path;
    std::string fastq_path;
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

// Internal helper function declarations (implementations later in file)
static int run_innotin(const InnotinParams& params);
static int run_select_missing_reads(const SelectMissingReadsParams& params);
static int run_split_plasmid_reads(const SplitPlasmidReadsParams& params);

// ============================================================================
// Utility functions
// ============================================================================

int check_dependencies(bool use_spades) {
    std::fprintf(stderr, "Checking required tools%s:\n", 
                 use_spades ? " (SPAdes mode)" : "");
    bool all_ok = validate_tools(required_tools(use_spades), true);
    
    if (all_ok) {
        std::fprintf(stderr, "\nAll required tools found.\n");
        return 0;
    } else {
        std::fprintf(stderr, "\nSome tools are missing. Please install them and ensure they are in PATH.\n");
        return 1;
    }
}

// ============================================================================
// Pipeline class implementation
// ============================================================================

Pipeline::Pipeline(const PipelineConfig& config)
    : config_(config)
    , output_dir_(config.output_directory) {
}

void Pipeline::log(const std::string& level, const std::string& message) const {
    // Simple verbosity filtering
    static const std::unordered_map<std::string, int> levels = {
        {"DEBUG", 0}, {"INFO", 1}, {"WARNING", 2}, {"ERROR", 3}, {"CRITICAL", 4}
    };
    
    auto config_level = levels.find(config_.verbosity);
    auto msg_level = levels.find(level);
    
    if (config_level != levels.end() && msg_level != levels.end()) {
        if (msg_level->second >= config_level->second) {
            std::fprintf(stderr, "[%s] %s\n", level.c_str(), message.c_str());
        }
    }
}

// ============================================================================
// GFA/FASTA utilities
// ============================================================================

void Pipeline::fix_gfa_empty_segments(const std::filesystem::path& input,
                                       const std::filesystem::path& output) {
    std::ifstream in(input);
    std::ofstream out(output);
    
    if (!in || !out) {
        if (config_.soft_fail) soft_fail_exit();
        log("ERROR", "Cannot open GFA files for empty segment fixing");
        std::exit(EXIT_FAILURE);
    }

    std::unordered_set<std::string> empty_segments;
    std::unordered_map<std::string, std::vector<std::pair<std::string, std::string>>> incoming;
    std::unordered_map<std::string, std::vector<std::pair<std::string, std::string>>> outgoing;
    
    std::vector<std::string> deferred_lines;
    std::string line;
    
    // First pass: identify empty segments and collect links
    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }
        
        if (line[0] == 'S') {
            // Segment line: S<tab>name<tab>sequence
            std::istringstream iss(line);
            std::string type;
            std::string name;
            std::string seq;
            std::getline(iss, type, '\t');
            std::getline(iss, name, '\t');
            std::getline(iss, seq, '\t');
            
            if (seq.empty()) {
                empty_segments.insert(name);
            } else {
                out << line << '\n';
            }
        } else if (line[0] == 'L') {
            deferred_lines.push_back(line);
        } else {
            out << line << '\n';
        }
    }
    
    // Second pass: process links
    for (const auto& lline : deferred_lines) {
        std::istringstream iss(lline);
        std::string type, from_name, from_orient, to_name, to_orient, overlap;
        std::getline(iss, type, '\t');
        std::getline(iss, from_name, '\t');
        std::getline(iss, from_orient, '\t');
        std::getline(iss, to_name, '\t');
        std::getline(iss, to_orient, '\t');
        std::getline(iss, overlap, '\t');
        
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

void Pipeline::remove_gfa_overlaps(const std::filesystem::path& input,
                                   const std::filesystem::path& output) const {
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
    
    std::ifstream in(input);
    if (!in) {
        log("ERROR", "Cannot open GFA file for overlap removal: " + input.string());
        std::exit(EXIT_FAILURE);
    }
    
    // Parse GFA
    std::unordered_map<std::string, int64_t> name_to_id;
    std::unordered_map<int64_t, std::string> id_to_name;
    std::vector<Segment> segments;
    std::vector<Link> links;
    std::vector<std::string> header_lines;
    int overlap = 0;
    int64_t next_id = 1;
    
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        
        if (line[0] == 'S') {
            std::istringstream iss(line);
            std::string type, name, seq;
            std::getline(iss, type, '\t');
            std::getline(iss, name, '\t');
            std::getline(iss, seq, '\t');
            std::string tags;
            std::getline(iss, tags);
            
            name_to_id[name] = next_id;
            id_to_name[next_id] = name;
            segments.emplace_back(name, seq, tags);
            next_id++;
        } else if (line[0] == 'L') {
            std::istringstream iss(line);
            std::string type, from_name, from_orient, to_name, to_orient, overlap_str;
            std::getline(iss, type, '\t');
            std::getline(iss, from_name, '\t');
            std::getline(iss, from_orient, '\t');
            std::getline(iss, to_name, '\t');
            std::getline(iss, to_orient, '\t');
            std::getline(iss, overlap_str, '\t');
            std::string link_tags;
            std::getline(iss, link_tags);  // Capture optional tags (rest of line)
            
            links.emplace_back(Link{from_name, from_orient == "+", to_name, to_orient == "+", link_tags});
            
            // Parse overlap value
            if (overlap == 0 && !overlap_str.empty() && overlap_str.back() == 'M') {
                overlap_str.pop_back();
                try {
                    overlap = std::stoi(overlap_str);
                } catch (...) {}
            }
        } else {
            header_lines.push_back(line);
        }
    }
    
    if (overlap == 0) {
        log("INFO", "GFA has no overlaps - copying as-is");
        std::filesystem::copy_file(input, output,
                                   std::filesystem::copy_options::overwrite_existing);
        return;
    }
    
    log("INFO", "Removing " + std::to_string(overlap) + "bp overlaps from GFA");
    
    // Build edge set and link maps
    // Edge: (signed_from, signed_to) where negative means reverse complement
    std::set<Edge> all_edges;
    std::unordered_map<SignedSeg, std::vector<SignedSeg>> forward_links;  // downstream
    std::unordered_map<SignedSeg, std::vector<SignedSeg>> reverse_links;  // upstream
    
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
    
    // Build constraint maps using std::map (Edge is comparable via std::pair)
    std::map<Edge, std::set<Edge>> must_match, must_differ;
    
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
    std::set<Edge> group_1, group_2;
    
    for (const auto& edge : all_edges) {
        if (group_1.contains(edge) || group_2.contains(edge)) {
            continue;
        }
        
        std::set<Edge> new_group_1;
        std::set<Edge> new_group_2;
        new_group_1.insert(edge);
        
        bool changed = true;
        while (changed) {
            changed = false;
            std::set<Edge> temp_1;
            std::set<Edge> temp_2;
            
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
    std::set<int64_t> large_trim_end;
    std::set<int64_t> large_trim_start;
    
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
        log("ERROR", "Cannot create output GFA file: " + output.string());
        std::exit(EXIT_FAILURE);
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

void Pipeline::extract_fasta_from_gfa(const std::filesystem::path& gfa,
                                       const std::filesystem::path& fasta,
                                       size_t min_length) const {
    std::ifstream in(gfa);
    std::ofstream out(fasta);
    
    if (!in || !out) {
        log("ERROR", "Cannot open files for GFA->FASTA conversion");
        std::exit(EXIT_FAILURE);
    }
    
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] != 'S') continue;
        
        std::istringstream iss(line);
        std::string type, name, seq;
        std::getline(iss, type, '\t');
        std::getline(iss, name, '\t');
        std::getline(iss, seq, '\t');
        
        if (seq.length() >= min_length) {
            out << '>' << name << '\n' << seq << '\n';
        }
    }
}

void Pipeline::write_circular_contigs(const std::filesystem::path& assembly_fasta,
                                       int iteration) {
    auto final_path = output_dir_ / ("plasmids.final.it" + std::to_string(iteration) + ".fasta");
    
    std::ifstream in(assembly_fasta);
    std::ofstream out(final_path);
    
    if (!in || !out) {
        log("ERROR", "Cannot open files for circular contig extraction");
        std::exit(EXIT_FAILURE);
    }
    
    std::string line;
    std::string current_header;
    std::string current_seq;
    bool is_circular = false;
    
    auto write_if_circular = [&]() {
        if (is_circular && !current_seq.empty()) {
            out << current_header << '\n' << current_seq << '\n';
        }
    };
    
    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }
        
        if (line[0] == '>') {
            write_if_circular();
            current_header = line;
            current_seq.clear();
            is_circular = (line.find("circular") != std::string::npos);
        } else {
            current_seq += line;
        }
    }
    write_if_circular();
}

[[noreturn]] void Pipeline::soft_fail_exit() {
    log("WARNING", "Soft-fail: falling back to circular contigs from SR assembly");
    auto sr_fasta = output_dir_ / "unicycler_sr" / "assembly.fasta";
    write_circular_contigs(sr_fasta, 0);
    symlink_remaining_iterations(0);
    std::exit(0);
}

void Pipeline::symlink_remaining_iterations(int from_iteration) {
    auto source = output_dir_ / ("plasmids.final.it" + std::to_string(from_iteration) + ".fasta");
    
    for (int i = from_iteration + 1; i <= config_.propagate_rounds; ++i) {
        auto link = output_dir_ / ("plasmids.final.it" + std::to_string(i) + ".fasta");
        force_symlink(source.filename(), link);
    }
}

// ============================================================================
// Pipeline stages
// ============================================================================

std::filesystem::path Pipeline::run_unicycler_sr_assembly() {
    auto unicycler_sr_path = output_dir_ / "unicycler_sr";
    auto assembly_gfa = unicycler_sr_path / "assembly.gfa";
    auto assembly_fasta = unicycler_sr_path / "assembly.fasta";
    
    // Skip if output already exists (unless --force)
    if (!config_.force && file_readable(assembly_gfa) && file_readable(assembly_fasta)) {
        log("WARNING", "Unicycler SR output exists, skipping assembly. Use --force to rerun.");
        return unicycler_sr_path;
    }
    
    ensure_directory(unicycler_sr_path);
    
    std::vector<std::string> cmd = {
        "unicycler_hyplas_modified",
        "-o", unicycler_sr_path.string(),
        "-t", std::to_string(config_.threads),
        "-1", config_.short_reads[0].string(),
        "--min_component_size", "10"
    };
    
    if (config_.short_reads.size() > 1) {
        cmd.emplace_back("-2");
        cmd.emplace_back(config_.short_reads[1].string());
    }
    
    run_cmd(cmd, "unicycler SR assembly")
        .expect_file(assembly_gfa, Expect::NON_EMPTY)
        .expect_file(assembly_fasta, Expect::NON_EMPTY | Expect::FASTA)
        .or_die("unicycler SR assembly");
    
    return unicycler_sr_path;
}

std::filesystem::path Pipeline::run_spades_sr_assembly() {
    auto spades_path = output_dir_ / "spades_sr";
    auto unicycler_sr_path = output_dir_ / "unicycler_sr";
    auto spades_gfa = spades_path / "assembly_graph_with_scaffolds.gfa";
    
    ensure_directory(spades_path);
    ensure_directory(unicycler_sr_path);
    
    // Skip SPAdes if output already exists (unless --force)
    if (!config_.force && file_readable(spades_gfa)) {
        log("WARNING", "SPAdes output exists, skipping assembly. Use --force to rerun.");
        setup_from_spades_output();
        return unicycler_sr_path;
    }
    
    // Run SPAdes with default parameters: -k 53 --gfa11 --isolate -m 1024
    std::vector<std::string> cmd = {
        "spades.py",
        "-o", spades_path.string(),
        "-t", std::to_string(config_.threads),
        "-k", "53",
        "--gfa11",
        "--isolate",
        "-m", "1024",
        "-1", config_.short_reads[0].string()
    };
    
    if (config_.short_reads.size() > 1) {
        cmd.emplace_back("-2");
        cmd.emplace_back(config_.short_reads[1].string());
    }
    
    run_cmd(cmd, "SPAdes SR assembly")
        .expect_file(spades_gfa, Expect::NON_EMPTY)
        .or_die("SPAdes SR assembly");
    
    // SPAdes outputs assembly_graph_with_scaffolds.gfa in GFA 1.1 format
    // Set up unicycler_sr directory structure for compatibility with rest of pipeline
    setup_from_spades_output();
    
    return unicycler_sr_path;
}

void Pipeline::setup_from_spades_output() {
    auto spades_path = output_dir_ / "spades_sr";
    auto unicycler_sr_path = output_dir_ / "unicycler_sr";
    
    ensure_directory(unicycler_sr_path);
    
    // SPAdes GFA output location
    auto spades_gfa = spades_path / "assembly_graph_with_scaffolds.gfa";
    if (!file_readable(spades_gfa)) {
        // Try alternative location
        spades_gfa = spades_path / "assembly_graph.gfa";
    }
    
    if (!file_readable(spades_gfa)) {
        log("ERROR", "SPAdes GFA output not found at expected locations");
        std::exit(EXIT_FAILURE);
    }
    
    auto gfa_path = unicycler_sr_path / "assembly.gfa";
    auto fasta_path = unicycler_sr_path / "assembly.fasta";
    
    // Remove overlaps from SPAdes GFA (equivalent to Unicycler's overlap removal)
    // This is required for minigraph which doesn't support overlapping segments
    remove_gfa_overlaps(spades_gfa, gfa_path);
    
    // Also create 002_depth_filter.gfa for compatibility with rest of pipeline
    auto depth_filter_gfa = unicycler_sr_path / "002_depth_filter.gfa";
    std::filesystem::copy_file(gfa_path, depth_filter_gfa,
                               std::filesystem::copy_options::overwrite_existing);
    
    // Extract FASTA from overlap-removed GFA
    extract_fasta_from_gfa(gfa_path, fasta_path, 200);
    
    log("INFO", "SPAdes assembly set up in unicycler_sr directory");
}

void Pipeline::setup_from_existing_assembly() {
    auto unicycler_sr_path = output_dir_ / "unicycler_sr";
    ensure_directory(unicycler_sr_path);
    
    auto gfa_path = unicycler_sr_path / "002_depth_filter.gfa";
    auto fasta_path = unicycler_sr_path / "assembly.fasta";
    
    // Copy the provided assembly
    std::filesystem::copy_file(*config_.sr_assembly, gfa_path,
                               std::filesystem::copy_options::overwrite_existing);
    
    // Run unicycler with mock input to set up directory structure
    // (This matches the Python behavior)
    auto mock_fq = make_temp_file(".fq");
    run_cmd({
        "unicycler_hyplas_modified",
        "-s", mock_fq.string(),
        "-o", unicycler_sr_path.string()
    }, "unicycler setup from existing assembly")
        .or_die("unicycler setup from existing assembly");
    remove_if_exists(mock_fq);
    
    // Extract FASTA from GFA
    extract_fasta_from_gfa(gfa_path, fasta_path, 200);
}

std::filesystem::path Pipeline::run_platon_classifier() {
    auto unicycler_fasta = output_dir_ / "unicycler_sr" / "assembly.fasta";
    auto platon_path = output_dir_ / "classify";
    auto result_tsv = platon_path / "result.tsv";
    
    if (!config_.force && file_readable(result_tsv)) {
        log("WARNING", "Platon output exists, skipping. Use --force to rerun.");
        return platon_path;
    }
    
    run_cmd({
        "platon",
        "-c",
        "--db", config_.platon_db.string(),
        "--threads", std::to_string(config_.threads),
        "--prefix", "result",
        "--output", platon_path.string(),
        unicycler_fasta.string()
    }, "platon classification")
        .expect_file(result_tsv, Expect::NON_EMPTY)
        .or_die_if(!config_.soft_fail, "platon classification")
        .or_execute([this]{ soft_fail_exit(); });
    
    return platon_path;
}

std::filesystem::path Pipeline::process_platon_output(const std::filesystem::path& platon_dir) {
    auto result_tsv = platon_dir / "result.tsv";
    auto output_tsv = platon_dir / "result_p.tsv";
    
    constexpr double max_chr_rds = -7.9;
    constexpr double min_plasmid_rds = 0.7;
    
    std::ifstream in(result_tsv);
    std::ofstream out(output_tsv);
    
    if (!in || !out) {
        if (config_.soft_fail) soft_fail_exit();
        log("ERROR", "Cannot open Platon result files");
        std::exit(EXIT_FAILURE);
    }
    
    out << "ID\tPREDICTION\n";
    
    // Parse header to find column indices
    std::string header;
    std::getline(in, header);
    
    std::istringstream hss(header);
    std::string col;
    std::unordered_map<std::string, int> col_idx;
    int idx = 0;
    while (std::getline(hss, col, '\t')) {
        col_idx[col] = idx++;
    }
    
    // Required columns
    int id_col = col_idx.contains("ID") ? col_idx["ID"] : 0;
    int rds_col = col_idx.contains("RDS") ? col_idx["RDS"] : -1;
    int circular_col = col_idx.contains("Circular") ? col_idx["Circular"] : -1;
    int inc_col = col_idx.contains("Inc Type(s)") ? col_idx["Inc Type(s)"] : -1;
    int rep_col = col_idx.contains("# Replication") ? col_idx["# Replication"] : -1;
    int mob_col = col_idx.contains("# Mobilization") ? col_idx["# Mobilization"] : -1;
    int orit_col = col_idx.contains("# OriT") ? col_idx["# OriT"] : -1;
    int hits_col = col_idx.contains("# Plasmid Hits") ? col_idx["# Plasmid Hits"] : -1;
    int rrna_col = col_idx.contains("# rRNAs") ? col_idx["# rRNAs"] : -1;
    
    if (rds_col < 0) {
        if (config_.soft_fail) soft_fail_exit();
        log("ERROR", "RDS column not found in Platon output");
        std::exit(EXIT_FAILURE);
    }
    
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }
        
        std::vector<std::string> fields;
        std::istringstream lss(line);
        std::string field;
        while (std::getline(lss, field, '\t')) {
            fields.push_back(field);
        }
        
        if (fields.size() <= static_cast<size_t>(rds_col)) continue;
        
        int nfields = static_cast<int>(fields.size());
        auto at = [&](int i) -> const std::string& { return fields[static_cast<size_t>(i)]; };
        
        const std::string& id = at(id_col);
        double rds = std::stod(at(rds_col));
        
        bool is_chr = rds <= max_chr_rds;
        bool is_plasmid = rds >= min_plasmid_rds;
        
        // Additional plasmid indicators for ambiguous cases
        bool circular = (circular_col >= 0 && circular_col < nfields && 
                        at(circular_col) == "yes");
        bool has_inc = (inc_col >= 0 && inc_col < nfields && 
                       !at(inc_col).empty() && at(inc_col) != "0");
        
        int rep = (rep_col >= 0 && rep_col < nfields) ? 
                  std::stoi(at(rep_col)) : 0;
        int mob = (mob_col >= 0 && mob_col < nfields) ? 
                  std::stoi(at(mob_col)) : 0;
        int orit = (orit_col >= 0 && orit_col < nfields) ? 
                   std::stoi(at(orit_col)) : 0;
        int hits = (hits_col >= 0 && hits_col < nfields) ? 
                   std::stoi(at(hits_col)) : 0;
        int rrnas = (rrna_col >= 0 && rrna_col < nfields) ? 
                    std::stoi(at(rrna_col)) : 0;
        
        bool repmob = (rep + mob) > 0;
        bool hit = (rds > 0.5 && hits > 0 && rrnas == 0);
        
        if (is_plasmid) {
            out << id << "\tplasmid\n";
        } else if (!is_chr && !is_plasmid && (circular || has_inc || repmob || orit > 0 || hit)) {
            out << id << "\tplasmid\n";
        } else if (is_chr) {
            out << id << "\tchromosome\n";
        }
        // Ambiguous cases not in any category are omitted
    }
    
    return output_tsv;
}

std::filesystem::path Pipeline::run_minigraph_lr_to_sr() {
    auto sr_graph = output_dir_ / "unicycler_sr" / "assembly.gfa";
    auto sr_graph_fix = output_dir_ / "unicycler_sr" / "assembly_segfix.gfa";
    auto gaf_output = output_dir_ / "lr2assembly.gaf";
    
    if (!config_.force && file_readable(gaf_output)) {
        log("WARNING", "Minigraph output exists, skipping. Use --force to rerun.");
        return gaf_output;
    }
    
    // Fix empty segments
    fix_gfa_empty_segments(sr_graph, sr_graph_fix);
    
    RunOptions opts;
    opts.stdout_file = gaf_output;
    
    run_cmd({
        "minigraph",
        sr_graph_fix.string(),
        config_.long_reads->string(),
        "-t", std::to_string(config_.threads),
        "-x", "lr",
        "-c"
    }, "minigraph LR to SR assembly", opts)
        .expect_file(gaf_output, Expect::NON_EMPTY)
        .or_die_if(!config_.soft_fail, "minigraph LR to SR assembly")
        .or_execute([this]{ soft_fail_exit(); });
    
    return gaf_output;
}

ReadSelectionResult Pipeline::run_long_read_selection(
    const std::filesystem::path& prediction_tsv,
    const std::filesystem::path& graph_alignment) {
    
    auto plasmid_lr_path = output_dir_ / "plasmid_long_reads";
    ensure_directory(plasmid_lr_path);
    
    ReadSelectionResult result;
    result.plasmid_reads = plasmid_lr_path / "plasmid.fastq.gz";
    result.unknown_both = plasmid_lr_path / "unknown_both.fastq.gz";
    result.unknown_neither = plasmid_lr_path / "unknown_neither.fastq.gz";
    result.unmapped = plasmid_lr_path / "unmapped.fastq.gz";
    
    if (!config_.force && file_readable(result.plasmid_reads) && file_readable(result.unknown_both)) {
        log("WARNING", "Read selection outputs exist, skipping. Use --force to rerun.");
        return result;
    }
    
    SplitPlasmidReadsParams split_params;
    split_params.gaf_path = graph_alignment.string();
    split_params.fastq_path = config_.long_reads->string();
    split_params.prediction_path = prediction_tsv.string();
    split_params.plasmid_out_path = result.plasmid_reads.string();
    split_params.unknown_neither_path = result.unknown_neither.string();
    split_params.unknown_both_path = result.unknown_both.string();
    split_params.unmapped_path = result.unmapped.string();
    
    Result(run_split_plasmid_reads(split_params))
        .expect_file(result.plasmid_reads, Expect::GZIPPED | Expect::FASTQ)
        .expect_file(result.unknown_both, Expect::GZIPPED | Expect::FASTQ)
        .expect_file(result.unknown_neither, Expect::GZIPPED | Expect::FASTQ)
        .expect_file(result.unmapped, Expect::GZIPPED | Expect::FASTQ)
        .or_die_if(!config_.soft_fail, "split-plasmid-reads")
        .or_execute([this]{ soft_fail_exit(); });
    
    return result;
}

std::filesystem::path Pipeline::find_missing_long_reads(
    const std::vector<std::filesystem::path>& plasmid_files,
    const std::vector<std::filesystem::path>& unknown_files,
    int round) {
    
    auto prop_dir = output_dir_ / "prop_lr";
    ensure_directory(prop_dir);
    
    auto paf_output = prop_dir / ("lr.round." + std::to_string(round) + ".paf");
    
    if (!config_.force && file_readable(paf_output)) {
        log("WARNING", "Propagation round " + std::to_string(round) + " output exists, skipping.");
        return paf_output;
    }
    
    // Concatenate unknown reads to temp file
    auto temp_unknown = concat_gzipped(unknown_files);
    if (temp_unknown.empty()) {
        if (config_.soft_fail) soft_fail_exit();
        log("ERROR", "Failed to concatenate unknown reads");
        std::exit(EXIT_FAILURE);
    }

    // Run innotin to filter reads
    TempFile temp_filtered(".fasta");
    {
        InnotinParams innotin_params;
        innotin_params.main_fastq = temp_unknown.string();
        for (const auto& pf : plasmid_files) {
            innotin_params.subset_fastqs.emplace_back(pf.string());
        }
        innotin_params.output_path = temp_filtered.path();
        
        Result(run_innotin(innotin_params))
            .expect_file(temp_filtered.path(), Expect::FASTA)
            .or_die_if(!config_.soft_fail, "innotin")
            .or_execute([this]{ soft_fail_exit(); });
    }
    
    // Clean up temp unknown file
    remove_if_exists(temp_unknown);
    
    // Run minimap2
    std::vector<std::string> minimap_cmd = {
        "minimap2",
        temp_filtered.string()
    };
    for (const auto& pf : plasmid_files) {
        minimap_cmd.emplace_back(pf.string());
    }
    minimap_cmd.emplace_back("-o");
    minimap_cmd.emplace_back(paf_output.string());
    minimap_cmd.emplace_back("-t");
    minimap_cmd.emplace_back(std::to_string(config_.threads));
    
    run_cmd(minimap_cmd, "minimap2 propagation round " + std::to_string(round))
        .expect_file(paf_output)
        .or_die_if(!config_.soft_fail, "minimap2 propagation round " + std::to_string(round))
        .or_execute([this]{ soft_fail_exit(); });
    
    return paf_output;
}

std::filesystem::path Pipeline::extract_missing_long_reads(
    const std::filesystem::path& plasmid_alignment,
    const std::vector<std::filesystem::path>& unknown_reads) {
    
    // Output path: replace .paf with .fastq.gz
    auto output_path = plasmid_alignment;
    output_path.replace_extension(".fastq.gz");
    
    if (!config_.force && file_readable(output_path)) {
        log("WARNING", "Extracted reads exist, skipping.");
        return output_path;
    }
    
    // Concatenate unknown reads
    auto temp_unknown = concat_gzipped(unknown_reads);
    if (temp_unknown.empty()) {
        if (config_.soft_fail) soft_fail_exit();
        log("ERROR", "Failed to concatenate unknown reads");
        std::exit(EXIT_FAILURE);
    }
    
    SelectMissingReadsParams select_params;
    select_params.paf_path = plasmid_alignment.string();
    select_params.fastq_path = temp_unknown.string();
    select_params.output_path = output_path.string();
    
    Result(run_select_missing_reads(select_params))
        .expect_file(output_path, Expect::GZIPPED | Expect::FASTQ)
        .or_die_if(!config_.soft_fail, "select-missing-reads")
        .or_execute([this]{ soft_fail_exit(); });
    
    remove_if_exists(temp_unknown);
    
    return output_path;
}

std::filesystem::path Pipeline::run_unicycler_lr_assembly(
    const std::vector<std::filesystem::path>& plasmid_files,
    int iteration) {
    
    auto unicycler_sr_path = output_dir_ / "unicycler_sr";
    auto unicycler_lr_path = output_dir_ / ("unicycler_lr_" + std::to_string(iteration));
    auto assembly_fasta = unicycler_lr_path / "assembly.fasta";
    
    // Skip if output already exists (unless --force)
    if (!config_.force && file_readable(assembly_fasta)) {
        log("WARNING", "LR assembly iteration " + std::to_string(iteration) + 
            " exists, skipping. Use --force to rerun.");
        return assembly_fasta;
    }
    
    // Copy SR assembly directory as starting point
    copy_directory(unicycler_sr_path, unicycler_lr_path, true);
    
    // Remove files that will be regenerated
    remove_if_exists(unicycler_lr_path / "assembly.fasta");
    remove_if_exists(unicycler_lr_path / "assembly.gfa");
    
    // Concatenate plasmid reads
    auto temp_lr = concat_gzipped(plasmid_files);
    if (temp_lr.empty()) {
        if (config_.soft_fail) soft_fail_exit();
        log("ERROR", "Failed to concatenate plasmid reads");
        std::exit(EXIT_FAILURE);
    }
    
    std::vector<std::string> cmd = {
        "unicycler_hyplas_modified",
        "--verbosity", "1",
        "--keep", "3",
        "-o", unicycler_lr_path.string(),
        "-t", std::to_string(config_.threads),
        "-l", temp_lr.string()
    };
    
    // Handle short reads
    TempFile empty_fq1;
    TempFile empty_fq2;
    if (config_.short_reads.empty()) {
        // Create empty fastq files to satisfy unicycler
        empty_fq1 = TempFile(".fq");
        empty_fq2 = TempFile(".fq");
        cmd.emplace_back("-1");
        cmd.emplace_back(empty_fq1.string());
        cmd.emplace_back("-2");
        cmd.emplace_back(empty_fq2.string());
    } else {
        cmd.emplace_back("-1");
        cmd.emplace_back(config_.short_reads[0].string());
        if (config_.short_reads.size() > 1) {
            cmd.emplace_back("-2");
            cmd.emplace_back(config_.short_reads[1].string());
        }
    }
    
    run_cmd(cmd, "unicycler LR assembly iteration " + std::to_string(iteration))
        .expect_file(assembly_fasta, Expect::FASTA)
        .or_die_if(!config_.soft_fail, "unicycler LR assembly iteration " + std::to_string(iteration))
        .or_execute([this]{ soft_fail_exit(); });
    
    remove_if_exists(temp_lr);
    
    return assembly_fasta;
}

// ============================================================================
// Main pipeline execution
// ============================================================================

int Pipeline::run() {
    log("INFO", "Starting HyPlAs pipeline");
    log("INFO", "Output directory: " + output_dir_.string());
    
    ensure_directory(output_dir_);
    
    // 1. SR assembly (or use provided --sr-assembly)
    if (config_.sr_assembly) {
        log("INFO", "Using provided SR assembly: " + config_.sr_assembly->string());
        setup_from_existing_assembly();
    } else if (config_.use_spades) {
        log("INFO", "Running SPAdes SR assembly (--use-spades mode)");
        run_spades_sr_assembly();
    } else {
        log("INFO", "Running Unicycler SR assembly");
        run_unicycler_sr_assembly();
    }
    
    // 2. Classification
    log("INFO", "Running Platon classifier");
    auto platon_dir = run_platon_classifier();
    auto prediction_tsv = process_platon_output(platon_dir);
    
    // 3. Graph alignment (if long reads provided)
    if (!config_.long_reads) {
        log("INFO", "No long reads provided, skipping propagation");
        return 0;
    }
    
    log("INFO", "Running minigraph alignment");
    auto graph_alignment = run_minigraph_lr_to_sr();
    
    // 4. Initial read selection
    log("INFO", "Running initial read selection");
    auto reads = run_long_read_selection(prediction_tsv, graph_alignment);
    
    // 5. Check if any plasmid reads found
    size_t plasmid_lines = line_count(reads.plasmid_reads);
    if (plasmid_lines < 4) {
        log("INFO", "No plasmid long reads found");
        
        // Write only circular contigs from SR assembly
        auto sr_fasta = output_dir_ / "unicycler_sr" / "assembly.fasta";
        write_circular_contigs(sr_fasta, 0);
        
        // Symlink remaining iterations
        symlink_remaining_iterations(0);
        return 0;
    }
    
    // 6. Initial LR assembly
    log("INFO", "Running initial LR assembly");
    std::vector<std::filesystem::path> plasmid_files = {reads.plasmid_reads};
    auto lr_assembly = run_unicycler_lr_assembly(plasmid_files, 0);
    write_circular_contigs(lr_assembly, 0);
    
    // 7. Propagation loop
    std::vector<std::filesystem::path> unknown_files = {
        reads.unmapped, reads.unknown_both, reads.unknown_neither
    };
    
    for (int i = 0; i < config_.propagate_rounds; ++i) {
        log("INFO", "Propagation round " + std::to_string(i + 1));
        
        auto plasmid_alignment = find_missing_long_reads(plasmid_files, unknown_files, i);
        
        if (line_count(plasmid_alignment) == 0) {
            log("INFO", "No alignments in round " + std::to_string(i + 1) + ", stopping propagation");
            symlink_remaining_iterations(i);
            break;
        }
        
        auto new_reads = extract_missing_long_reads(
            plasmid_alignment, unknown_files);
        
        if (line_count(new_reads) == 0) {
            log("INFO", "No new reads in round " + std::to_string(i + 1) + ", stopping propagation");
            symlink_remaining_iterations(i);
            break;
        }
        
        plasmid_files.emplace_back(new_reads);
        lr_assembly = run_unicycler_lr_assembly(plasmid_files, i + 1);
        write_circular_contigs(lr_assembly, i + 1);
    }
    
    log("INFO", "Pipeline completed successfully");
    return 0;
}

// ============================================================================
// Internal helper functions (innotin, select-missing-reads, split-plasmid-reads)
// ============================================================================

namespace {

std::string extract_read_id(MmapView& view) {
    view.extend_until(" \t\n");
    return static_cast<std::string>(view);
}

} // anonymous namespace

int run_innotin(const InnotinParams& params) {
    std::unordered_set<std::string> subset_ids;
    
    // Collect all read IDs from subset files
    for (const auto& subset_path : params.subset_fastqs) {
        std::error_code error;
        mio::mmap_source subset_mmap = mio::make_mmap_source(subset_path, error);
        
        if (error) {
            std::fprintf(stderr, "Error mapping subset file '%s': %s\n", 
                        subset_path.c_str(), error.message().c_str());
            return EXIT_FAILURE;
        }
        
        MmapView view{&subset_mmap};
        
        while (view.s < view.cap) {
            view.skip_next('@');
            std::string read_id = extract_read_id(view);
            subset_ids.insert(read_id);
            
            // Skip to next record
            view.skip_next('\n');  // rest of header
            view.skip_next('\n');  // sequence
            view.skip_next('\n');  // + line
            view.skip_next('\n');  // quality
        }
    }
    
    // Process main file and output reads not in subsets
    gzFile main_fp = gzopen(params.main_fastq.c_str(), "r");
    if (!main_fp) {
        std::fprintf(stderr, "Error opening main FASTQ file: %s\n", 
                    params.main_fastq.c_str());
        return EXIT_FAILURE;
    }
    
    // Open output: file or stdout
    FILE* out_fp = stdout;
    if (!params.output_path.empty()) {
        out_fp = std::fopen(params.output_path.c_str(), "w");
        if (!out_fp) {
            std::fprintf(stderr, "Error opening output file: %s\n",
                        params.output_path.c_str());
            gzclose(main_fp);
            return EXIT_FAILURE;
        }
    }
    
    kseq_t* seq = kseq_init(main_fp);
    
    while (kseq_read(seq) >= 0) {
        std::string name{seq->name.s};
        
        if (subset_ids.find(name) == subset_ids.end()) {
            std::fprintf(out_fp, ">%s", seq->name.s);
            if (seq->comment.s && seq->comment.l > 0) {
                std::fprintf(out_fp, " %s", seq->comment.s);
            }
            std::fprintf(out_fp, "\n%s\n", seq->seq.s);
        }
    }
    
    kseq_destroy(seq);
    gzclose(main_fp);
    if (out_fp != stdout) {
        std::fclose(out_fp);
    }
    
    return EXIT_SUCCESS;
}

// ---------------------------------------------------------------------------
// select-missing-reads: Select reads that overlap with plasmid reads
// ---------------------------------------------------------------------------

namespace {

struct PafEntry {
    bool is_terminal = false;
    LiteView id{};
    std::vector<LiteView> contigs;

    explicit PafEntry(MmapView& paf_view) {
        if (paf_view.e >= paf_view.cap) {
            is_terminal = true;
            return;
        }

        bool first = true;
        do {
            paf_view.extend_until('\t');
            
            LiteView current_id = static_cast<LiteView>(paf_view);
            if (paf_view != id) {
                if (!first) {
                    return;
                }
                id = current_id;
            }
            first = false;

            paf_view.skip_next_n('\t', 5);
            paf_view.extend_until('\t');
            
            MmapView align_view = paf_view.focus();

            if (align_view[0] == '<' || align_view[0] == '>') {
                align_view.skip_next("<>");
            }

            while (align_view.e < align_view.cap) {
                align_view.extend_until("<>");
                contigs.push_back(static_cast<LiteView>(align_view));
                align_view.skip_next("<>");
            }

            paf_view.skip_next('\n');
        } while (paf_view.e < paf_view.cap);
    }
};

} // anonymous namespace

int run_select_missing_reads(const SelectMissingReadsParams& params) {
    std::error_code error;

    mio::mmap_source paf_mmap = mio::make_mmap_source(params.paf_path, error);
    if (error) {
        std::fprintf(stderr, "Error mapping PAF file: %s\n", error.message().c_str());
        std::fprintf(stderr, "Creating empty output since there are no mappings\n");
        
        FILE* empty_out = open_gzip_write_pipe(params.output_path);
        if (empty_out) pclose(empty_out);
        return EXIT_SUCCESS;
    }

    gzFile fastq_fp = gzopen(params.fastq_path.c_str(), "r");
    if (!fastq_fp) {
        std::fprintf(stderr, "Error opening FASTQ file: %s\n", params.fastq_path.c_str());
        return EXIT_FAILURE;
    }

    FILE* plasmid_out = open_gzip_write_pipe(params.output_path);
    if (!plasmid_out) {
        std::fprintf(stderr, "Error opening output file\n");
        gzclose(fastq_fp);
        return EXIT_FAILURE;
    }

    kseq_t* seq = kseq_init(fastq_fp);
    MmapView paf_view{&paf_mmap};

    std::unordered_map<std::string, int> reads_to_use;

    for (PafEntry entry{paf_view}; !entry.is_terminal; entry = PafEntry{paf_view}) {
        for (const LiteView& v : entry.contigs) {
            MmapView contig_view{&paf_mmap, v};
            reads_to_use[static_cast<std::string>(contig_view)] = 1;
        }
    }

    while (kseq_read(seq) >= 0) {
        std::string name{seq->name.s};
        auto it = reads_to_use.find(name);
        
        if (it != reads_to_use.end() && it->second > 0) {
            std::fprintf(plasmid_out, "@%s %s\n%s\n+\n%s\n",
                        seq->name.s,
                        seq->comment.s ? seq->comment.s : "",
                        seq->seq.s,
                        seq->qual.s);
            --it->second;
        }
    }

    kseq_destroy(seq);
    gzclose(fastq_fp);
    pclose(plasmid_out);

    return EXIT_SUCCESS;
}

// ---------------------------------------------------------------------------
// split-plasmid-reads: Split reads based on alignment to plasmid/chromosome
// ---------------------------------------------------------------------------

namespace {

struct GafEntry {
    bool is_terminal = false;
    LiteView id{};
    std::vector<LiteView> contigs;

    explicit GafEntry(MmapView& gaf_view) {
        if (gaf_view.e >= gaf_view.cap) {
            is_terminal = true;
            return;
        }

        bool first = true;
        do {
            gaf_view.extend_until('\t');
            
            LiteView current_id = static_cast<LiteView>(gaf_view);
            if (gaf_view != id) {
                if (!first) {
                    return;
                }
                id = current_id;
            }
            first = false;

            gaf_view.skip_next_n('\t', 5);
            gaf_view.extend_until('\t');
            
            MmapView align_view = gaf_view.focus();
            align_view.skip_next("<>");

            while (align_view.e < align_view.cap) {
                align_view.extend_until("<>");
                contigs.push_back(static_cast<LiteView>(align_view));
                align_view.skip_next("<>");
            }

            gaf_view.skip_next('\n');
        } while (gaf_view.e < gaf_view.cap);
    }
};

} // anonymous namespace

int run_split_plasmid_reads(const SplitPlasmidReadsParams& params) {
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

    FILE* plasmid_out = open_gzip_write_pipe(params.plasmid_out_path);
    FILE* unknown_neither_out = open_gzip_write_pipe(params.unknown_neither_path);
    FILE* unknown_both_out = open_gzip_write_pipe(params.unknown_both_path);
    FILE* unmapped_out = open_gzip_write_pipe(params.unmapped_path);

    if (!plasmid_out || !unknown_neither_out || !unknown_both_out || !unmapped_out) {
        std::fprintf(stderr, "Error opening output files\n");
        return EXIT_FAILURE;
    }

    FILE* chr_out = nullptr;
    if (params.output_chromosomal) {
        chr_out = open_gzip_write_pipe(params.chr_out_path);
        if (!chr_out) {
            std::fprintf(stderr, "Error opening chromosome output file\n");
            return EXIT_FAILURE;
        }
    }

    auto plasmid_contigs = parse_plasmid_tsv(params.prediction_path);

    kseq_t* seq = kseq_init(fastq_fp);
    MmapView gaf_view{&gaf_mmap};

    for (GafEntry entry{gaf_view}; !entry.is_terminal; entry = GafEntry{gaf_view}) {
        int64_t l = kseq_read(seq);

        MmapView gid{&gaf_mmap, entry.id};
        
        while (l >= 0) {
            if (gid != std::string_view(seq->name.s)) {
                std::fprintf(unmapped_out, "@%s %s\n%s\n+\n%s\n",
                            seq->name.s,
                            seq->comment.s ? seq->comment.s : "",
                            seq->seq.s,
                            seq->qual.s);
                l = kseq_read(seq);
            } else {
                break;
            }
        }

        int from_chromosome = 0;
        int from_plasmid = 0;

        for (const LiteView& v : entry.contigs) {
            MmapView contig_view{&gaf_mmap, v};
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
            std::fprintf(unknown_neither_out, "@%s %s\n%s\n+\n%s\n",
                        seq->name.s, comment, seq->seq.s, seq->qual.s);
        } else if (from_plasmid > 0 && from_chromosome == 0) {
            std::fprintf(plasmid_out, "@%s %s\n%s\n+\n%s\n",
                        seq->name.s, comment, seq->seq.s, seq->qual.s);
        } else if (from_chromosome > 0 && from_plasmid > 0) {
            std::fprintf(unknown_both_out, "@%s %s\n%s\n+\n%s\n",
                        seq->name.s, comment, seq->seq.s, seq->qual.s);
        } else if (params.output_chromosomal && from_chromosome > 0) {
            std::fprintf(chr_out, "@%s %s\n%s\n+\n%s\n",
                        seq->name.s, comment, seq->seq.s, seq->qual.s);
        }
    }

    kseq_destroy(seq);
    gzclose(fastq_fp);

    pclose(plasmid_out);
    pclose(unknown_neither_out);
    pclose(unknown_both_out);
    pclose(unmapped_out);
    
    if (params.output_chromosomal && chr_out) {
        pclose(chr_out);
    }

    return EXIT_SUCCESS;
}

} // namespace hyplas

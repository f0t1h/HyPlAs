/**
 * @file pipeline.cpp
 * @brief HyPlAs pipeline stage implementations
 */

#include "pipeline.hpp"
#include "file_utils.hpp"
#include "process.hpp"

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
        if (line.empty()) continue;
        
        if (line[0] == 'S') {
            // Segment line: S<tab>name<tab>sequence
            std::istringstream iss(line);
            std::string type, name, seq;
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
        
        if (empty_segments.count(from_name)) {
            outgoing[from_name].emplace_back(to_name, to_orient);
        } else if (empty_segments.count(to_name)) {
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
            segments.push_back({name, seq, tags});
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
            
            links.push_back({from_name, from_orient == "+", to_name, to_orient == "+"});
            
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
        
        forward_links[from_seg].push_back(to_seg);
        reverse_links[to_seg].push_back(from_seg);
        forward_links[-to_seg].push_back(-from_seg);
        reverse_links[-from_seg].push_back(-to_seg);
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
        if (group_1.count(edge) || group_2.count(edge)) continue;
        
        std::set<Edge> new_group_1, new_group_2;
        new_group_1.insert(edge);
        
        bool changed = true;
        while (changed) {
            changed = false;
            std::set<Edge> temp_1, temp_2;
            
            for (const auto& g1_edge : new_group_1) {
                for (const auto& match_edge : must_match[g1_edge]) {
                    if (!new_group_1.count(match_edge)) temp_1.insert(match_edge);
                }
                for (const auto& differ_edge : must_differ[g1_edge]) {
                    if (!new_group_2.count(differ_edge)) temp_2.insert(differ_edge);
                }
            }
            for (const auto& g2_edge : new_group_2) {
                for (const auto& match_edge : must_match[g2_edge]) {
                    if (!new_group_2.count(match_edge)) temp_2.insert(match_edge);
                }
                for (const auto& differ_edge : must_differ[g2_edge]) {
                    if (!new_group_1.count(differ_edge)) temp_1.insert(differ_edge);
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
    std::set<int64_t> large_trim_end, large_trim_start;
    
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
        
        int start_trim = large_trim_start.count(seg_id) ? large_half : small_half;
        int end_trim = large_trim_end.count(seg_id) ? large_half : small_half;
        
        std::string trimmed_seq = seg.sequence;
        if (static_cast<int>(trimmed_seq.length()) > start_trim + end_trim) {
            trimmed_seq = trimmed_seq.substr(start_trim, trimmed_seq.length() - start_trim - end_trim);
        }
        
        out << "S\t" << seg.name << '\t' << trimmed_seq;
        if (!seg.tags.empty()) {
            out << '\t' << seg.tags;
        }
        out << '\n';
    }
    
    // Write links with 0M overlap
    for (const auto& link : links) {
        out << "L\t" << link.from_name << '\t' << (link.from_is_forward ? '+' : '-')
            << '\t' << link.to_name << '\t' << (link.to_is_forward ? '+' : '-')
            << "\t0M\n";
    }
}

void Pipeline::extract_fasta_from_gfa(const std::filesystem::path& gfa,
                                       const std::filesystem::path& fasta,
                                       size_t min_length) {
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
        if (line.empty()) continue;
        
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
    ensure_directory(unicycler_sr_path);
    
    std::vector<std::string> cmd = {
        "unicycler_hyplas_modified",
        "-o", unicycler_sr_path.string(),
        "-t", std::to_string(config_.threads),
        "-1", config_.short_reads[0].string(),
        "--min_component_size", "10"
    };
    
    if (config_.short_reads.size() > 1) {
        cmd.push_back("-2");
        cmd.push_back(config_.short_reads[1].string());
    }
    
    run_or_die(cmd, "unicycler SR assembly");
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
        cmd.push_back("-2");
        cmd.push_back(config_.short_reads[1].string());
    }
    
    run_or_die(cmd, "SPAdes SR assembly");
    
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
    run_or_die({
        "unicycler_hyplas_modified",
        "-s", mock_fq.string(),
        "-o", unicycler_sr_path.string()
    }, "unicycler setup from existing assembly");
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
    
    run_or_die({
        "platon",
        "-c",
        "--db", config_.platon_db.string(),
        "--threads", std::to_string(config_.threads),
        "--prefix", "result",
        "--output", platon_path.string(),
        unicycler_fasta.string()
    }, "platon classification");
    
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
    int id_col = col_idx.count("ID") ? col_idx["ID"] : 0;
    int rds_col = col_idx.count("RDS") ? col_idx["RDS"] : -1;
    int circular_col = col_idx.count("Circular") ? col_idx["Circular"] : -1;
    int inc_col = col_idx.count("Inc Type(s)") ? col_idx["Inc Type(s)"] : -1;
    int rep_col = col_idx.count("# Replication") ? col_idx["# Replication"] : -1;
    int mob_col = col_idx.count("# Mobilization") ? col_idx["# Mobilization"] : -1;
    int orit_col = col_idx.count("# OriT") ? col_idx["# OriT"] : -1;
    int hits_col = col_idx.count("# Plasmid Hits") ? col_idx["# Plasmid Hits"] : -1;
    int rrna_col = col_idx.count("# rRNAs") ? col_idx["# rRNAs"] : -1;
    
    if (rds_col < 0) {
        log("ERROR", "RDS column not found in Platon output");
        std::exit(EXIT_FAILURE);
    }
    
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        
        std::vector<std::string> fields;
        std::istringstream lss(line);
        std::string field;
        while (std::getline(lss, field, '\t')) {
            fields.push_back(field);
        }
        
        if (fields.size() <= static_cast<size_t>(rds_col)) continue;
        
        std::string id = fields[id_col];
        double rds = std::stod(fields[rds_col]);
        
        bool is_chr = rds <= max_chr_rds;
        bool is_plasmid = rds >= min_plasmid_rds;
        
        // Additional plasmid indicators for ambiguous cases
        bool circular = (circular_col >= 0 && static_cast<size_t>(circular_col) < fields.size() && 
                        fields[circular_col] == "yes");
        bool has_inc = (inc_col >= 0 && static_cast<size_t>(inc_col) < fields.size() && 
                       !fields[inc_col].empty() && fields[inc_col] != "0");
        
        int rep = (rep_col >= 0 && static_cast<size_t>(rep_col) < fields.size()) ? 
                  std::stoi(fields[rep_col]) : 0;
        int mob = (mob_col >= 0 && static_cast<size_t>(mob_col) < fields.size()) ? 
                  std::stoi(fields[mob_col]) : 0;
        int orit = (orit_col >= 0 && static_cast<size_t>(orit_col) < fields.size()) ? 
                   std::stoi(fields[orit_col]) : 0;
        int hits = (hits_col >= 0 && static_cast<size_t>(hits_col) < fields.size()) ? 
                   std::stoi(fields[hits_col]) : 0;
        int rrnas = (rrna_col >= 0 && static_cast<size_t>(rrna_col) < fields.size()) ? 
                    std::stoi(fields[rrna_col]) : 0;
        
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
    
    run_or_die({
        "minigraph",
        sr_graph_fix.string(),
        config_.long_reads->string(),
        "-t", std::to_string(config_.threads),
        "-x", "lr",
        "-c"
    }, "minigraph LR to SR assembly", opts);
    
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
    
    run_or_die({
        "hyplas-utils", "split-plasmid-reads",
        "-g", graph_alignment.string(),
        "-f", config_.long_reads->string(),
        "-t", prediction_tsv.string(),
        "-p", result.plasmid_reads.string(),
        "-n", result.unknown_neither.string(),
        "-b", result.unknown_both.string(),
        "-u", result.unmapped.string()
    }, "split-plasmid-reads");
    
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
        log("ERROR", "Failed to concatenate unknown reads");
        std::exit(EXIT_FAILURE);
    }
    
    // Run innotin to filter reads
    TempFile temp_filtered(".fastq");
    {
        std::vector<std::string> innotin_cmd = {
            "hyplas-utils", "innotin",
            temp_unknown.string()
        };
        for (const auto& pf : plasmid_files) {
            innotin_cmd.push_back(pf.string());
        }
        
        RunOptions opts;
        opts.stdout_file = temp_filtered.path();
        run_or_die(innotin_cmd, "innotin filtering", opts);
    }
    
    // Clean up temp unknown file
    remove_if_exists(temp_unknown);
    
    // Run minimap2
    std::vector<std::string> minimap_cmd = {
        "minimap2",
        temp_filtered.string()
    };
    for (const auto& pf : plasmid_files) {
        minimap_cmd.push_back(pf.string());
    }
    minimap_cmd.push_back("-o");
    minimap_cmd.push_back(paf_output.string());
    minimap_cmd.push_back("-t");
    minimap_cmd.push_back(std::to_string(config_.threads));
    
    run_or_die(minimap_cmd, "minimap2 propagation round " + std::to_string(round));
    
    return paf_output;
}

std::filesystem::path Pipeline::extract_missing_long_reads(
    const std::filesystem::path& plasmid_alignment,
    const std::filesystem::path& graph_alignment,
    const std::filesystem::path& prediction_tsv,
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
        log("ERROR", "Failed to concatenate unknown reads");
        std::exit(EXIT_FAILURE);
    }
    
    run_or_die({
        "hyplas-utils", "select-missing-reads",
        "-p", plasmid_alignment.string(),
        "-g", graph_alignment.string(),
        "-f", temp_unknown.string(),
        "-o", output_path.string(),
        "-t", prediction_tsv.string()
    }, "select-missing-reads");
    
    remove_if_exists(temp_unknown);
    
    return output_path;
}

std::filesystem::path Pipeline::run_unicycler_lr_assembly(
    const std::vector<std::filesystem::path>& plasmid_files,
    int iteration) {
    
    auto unicycler_sr_path = output_dir_ / "unicycler_sr";
    auto unicycler_lr_path = output_dir_ / ("unicycler_lr_" + std::to_string(iteration));
    
    // Copy SR assembly directory as starting point
    copy_directory(unicycler_sr_path, unicycler_lr_path, true);
    
    // Remove files that will be regenerated
    remove_if_exists(unicycler_lr_path / "assembly.fasta");
    remove_if_exists(unicycler_lr_path / "assembly.gfa");
    
    // Concatenate plasmid reads
    auto temp_lr = concat_gzipped(plasmid_files);
    if (temp_lr.empty()) {
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
    TempFile empty_fq1, empty_fq2;
    if (config_.short_reads.empty()) {
        // Create empty fastq files to satisfy unicycler
        empty_fq1 = TempFile(".fq");
        empty_fq2 = TempFile(".fq");
        cmd.push_back("-1");
        cmd.push_back(empty_fq1.string());
        cmd.push_back("-2");
        cmd.push_back(empty_fq2.string());
    } else {
        cmd.push_back("-1");
        cmd.push_back(config_.short_reads[0].string());
        if (config_.short_reads.size() > 1) {
            cmd.push_back("-2");
            cmd.push_back(config_.short_reads[1].string());
        }
    }
    
    run_or_die(cmd, "unicycler LR assembly iteration " + std::to_string(iteration));
    
    remove_if_exists(temp_lr);
    
    return unicycler_lr_path / "assembly.fasta";
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
            plasmid_alignment, graph_alignment, prediction_tsv, unknown_files);
        
        if (line_count(new_reads) == 0) {
            log("INFO", "No new reads in round " + std::to_string(i + 1) + ", stopping propagation");
            symlink_remaining_iterations(i);
            break;
        }
        
        plasmid_files.push_back(new_reads);
        lr_assembly = run_unicycler_lr_assembly(plasmid_files, i + 1);
        write_circular_contigs(lr_assembly, i + 1);
    }
    
    log("INFO", "Pipeline completed successfully");
    return 0;
}

} // namespace hyplas

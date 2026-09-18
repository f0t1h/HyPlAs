/**
 * @file pipeline.cpp
 * @brief HyPlAs pipeline stage implementations
 */

#include "pipeline.hpp"
#include "error.hpp"
#include "fastx.hpp"
#include "gfa.hpp"
#include "log.hpp"
#include "contig_classification.hpp"
#include "stage.hpp"

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <sstream>
#include <string>
#include <system_error>
#include <gtl/phmap.hpp>
#include <vector>

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
    , output_dir_(config.output_directory)
    , tmp_dir_(config.output_directory / "tmp") {
    set_log_level(config.verbosity);
}

std::filesystem::path Pipeline::make_temp(const std::string& name) const {
    ensure_directory(tmp_dir_);
    auto path = tmp_dir_ / name;
    std::ofstream file(path);
    if (!file) {
        log("ERROR", "Cannot create temp file: " + path.string());
        std::exit(EXIT_FAILURE);
    }
    return path;
}

void Pipeline::log(const std::string& level, const std::string& message) const {
    hyplas::log(level, message);
}

// ============================================================================
// Final-output helpers (thin wrappers over the free gfa transforms)
// ============================================================================

void Pipeline::write_circular_plasmid_contigs(
    const std::filesystem::path& gfa_path,
    const std::filesystem::path& fasta_path,
    const std::filesystem::path& prediction_tsv,
    int iteration) const {

    auto final_path = output_dir_ / ("plasmids.final.it" + std::to_string(iteration) + ".fasta");
    std::ofstream out(final_path);
    if (!out) {
        throw HyplasError("cannot open output: " + final_path.string());
    }
    gtl::flat_hash_set<std::string> written;
    append_circular_sr_plasmids(out, gfa_path, fasta_path, prediction_tsv, written);
}

void Pipeline::write_circular_contigs(const std::filesystem::path& assembly_fasta,
                                       int iteration) const {
    auto final_path = output_dir_ / ("plasmids.final.it" + std::to_string(iteration) + ".fasta");
    std::ofstream out(final_path);
    if (!out) {
        throw HyplasError("cannot open output: " + final_path.string());
    }
    gtl::flat_hash_set<std::string> written;
    append_circular_by_header(out, assembly_fasta, written);
}

[[noreturn]] void Pipeline::soft_fail_exit() {
    log("WARNING", "Soft-fail: falling back to circular plasmid contigs from SR assembly");

    auto sr_gfa = output_dir_ / "unicycler_sr" / "assembly.gfa";
    auto sr_fasta = output_dir_ / "unicycler_sr" / "assembly.fasta";

    write_circular_plasmid_contigs(sr_gfa, sr_fasta, prediction_tsv_, 0);
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

    auto sr_stage = stage("unicycler SR assembly").expect_which("unicycler_hyplas_modified");
    for (const auto& sr : config_.short_reads) sr_stage.expect_file(sr, Expect::NON_EMPTY);
    sr_stage.proc(cmd)
        .expect_file(assembly_gfa, Expect::NON_EMPTY)
        .expect_file(assembly_fasta, Expect::NON_EMPTY)
        .expect_file(assembly_fasta, file_is_fasta, "FASTA")
        .or_die_if(!config_.soft_fail)
        .or_execute([this]{ soft_fail_exit(); });

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

    // Run SPAdes with default parameters: -k 99 --gfa11 --isolate -m 1024
    std::vector<std::string> cmd = {
        "spades.py",
        "-o", spades_path.string(),
        "-t", std::to_string(config_.threads),
        "-k", "99",
        "--gfa11",
        "--isolate",
        "-m", "1024",
        "-1", config_.short_reads[0].string()
    };

    if (config_.short_reads.size() > 1) {
        cmd.emplace_back("-2");
        cmd.emplace_back(config_.short_reads[1].string());
    }

    auto spades_stage = stage("SPAdes SR assembly").expect_which("spades.py");
    for (const auto& sr : config_.short_reads) spades_stage.expect_file(sr, Expect::NON_EMPTY);
    spades_stage.proc(cmd)
        .expect_file(spades_gfa, Expect::NON_EMPTY)
        .or_die_if(!config_.soft_fail)
        .or_execute([this]{ soft_fail_exit(); });

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
    auto mock_fq = make_temp("mock_sr.fq");
    stage("unicycler setup from existing assembly")
        .expect_which("unicycler_hyplas_modified")
        .expect_file(gfa_path, Expect::NON_EMPTY)
        .proc({
            "unicycler_hyplas_modified",
            "-s", mock_fq.string(),
            "-o", unicycler_sr_path.string()
        })
        .or_die_if(true);

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

    stage("platon classification")
        .expect_which("platon")
        .expect_file(unicycler_fasta, Expect::NON_EMPTY)
        .proc({
            "platon",
            "-c",
            "--db", config_.platon_db.string(),
            "--threads", std::to_string(config_.threads),
            "--prefix", "result",
            "--output", platon_path.string(),
            unicycler_fasta.string()
        })
        .expect_file(result_tsv, Expect::NON_EMPTY)
        .or_die_if(true);

    return platon_path;
}

std::filesystem::path Pipeline::process_platon_output(const std::filesystem::path& platon_dir) {
    auto result_tsv = platon_dir / "result.tsv";
    auto output_tsv = platon_dir / "result_p.tsv";

    std::string content = read_file(result_tsv);

    std::string body;
    try {
        body = classify_platon_tsv(content);
    } catch (const HyplasError&) {
        if (config_.soft_fail) soft_fail_exit();
        log("ERROR", "RDS column not found in Platon output");
        std::exit(EXIT_FAILURE);
    }

    std::ofstream out(output_tsv);
    if (!out) {
        log("ERROR", "Cannot open Platon output file");
        std::exit(EXIT_FAILURE);
    }
    out << body;

    return output_tsv;
}

std::filesystem::path Pipeline::run_minigraph_lr_to_sr(
    const std::filesystem::path& reads_fastq,
    const std::filesystem::path& gaf_output)
{
    auto sr_graph = output_dir_ / "unicycler_sr" / "assembly.gfa";
    auto sr_graph_fix = output_dir_ / "unicycler_sr" / "assembly_segfix.gfa";

    if (!config_.force && file_readable(gaf_output)) {
        log("WARNING", "Minigraph output exists (" + gaf_output.filename().string() +
                       "), skipping. Use --force to rerun.");
        return gaf_output;
    }

    // Fix empty segments once per pipeline (idempotent)
    if (!file_readable(sr_graph_fix)) {
        try {
            fix_gfa_empty_segments(sr_graph, sr_graph_fix);
        } catch (const HyplasError& e) {
            if (config_.soft_fail) soft_fail_exit();
            log("ERROR", e.what());
            std::exit(EXIT_FAILURE);
        }
    }

    RunOptions opts;
    opts.stdout_file = gaf_output;

    stage("minigraph LR to SR assembly")
        .expect_which("minigraph")
        .expect_file(sr_graph_fix, Expect::NON_EMPTY)
        .expect_file(reads_fastq, Expect::NON_EMPTY)
        .proc({
            "minigraph",
            sr_graph_fix.string(),
            reads_fastq.string(),
            "-t", std::to_string(config_.threads),
            "-x", "lr",
            "-c"
        }, opts)
        .expect_file(gaf_output, Expect::NON_EMPTY)
        .or_die_if(!config_.soft_fail)
        .or_execute([this]{ soft_fail_exit(); });

    return gaf_output;
}

std::filesystem::path Pipeline::run_minigraph_lr_to_sr() {
    return run_minigraph_lr_to_sr(*config_.long_reads,
                                  output_dir_ / "lr2assembly.gaf");
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

    const auto valid_cached_fastq = [](const std::filesystem::path& path) {
        return file_readable(path) && file_is_gzipped(path) && file_is_fastq(path);
    };
    if (!config_.force &&
        valid_cached_fastq(result.plasmid_reads) &&
        valid_cached_fastq(result.unknown_both) &&
        valid_cached_fastq(result.unknown_neither) &&
        valid_cached_fastq(result.unmapped)) {
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

    stage("split-plasmid-reads")
        .expect_success(split_plasmid_reads(split_params))
        .expect_file(result.plasmid_reads, file_is_gzipped_fastq, "gzipped FASTQ")
        .expect_file(result.unknown_both, file_is_gzipped_fastq, "gzipped FASTQ")
        .expect_file(result.unknown_neither, file_is_gzipped_fastq, "gzipped FASTQ")
        .expect_file(result.unmapped, file_is_gzipped_fastq, "gzipped FASTQ")
        .or_die_if(!config_.soft_fail)
        .or_execute([this]{ soft_fail_exit(); });

    return result;
}

std::filesystem::path Pipeline::find_missing_long_reads(
    const std::vector<std::filesystem::path>& plasmid_files,
    const std::vector<std::filesystem::path>& unknown_files,
    int round,
    int comp_id) {

    auto prop_dir = output_dir_ / "prop_lr";
    ensure_directory(prop_dir);

    std::string paf_name = (comp_id >= 0)
        ? ("comp_" + std::to_string(comp_id) + ".round." + std::to_string(round) + ".paf")
        : ("lr.round." + std::to_string(round) + ".paf");
    auto paf_output = prop_dir / paf_name;

    if (!config_.force && file_readable(paf_output)) {
        log("WARNING", "Propagation round " + std::to_string(round) + " output exists, skipping.");
        return paf_output;
    }

    // Tag temporary filtered output by propagation round and component.
    auto tag = (comp_id >= 0)
        ? ("r" + std::to_string(round) + "_c" + std::to_string(comp_id))
        : ("r" + std::to_string(round));

    // Run innotin to filter reads
    auto temp_filtered_path = make_temp("filtered_unknown_" + tag + ".fasta");
    {
        InnotinParams innotin_params;
        for (const auto& path : unknown_files) {
            innotin_params.main_fastqs.emplace_back(path.string());
        }
        for (const auto& pf : plasmid_files) {
            innotin_params.subset_fastqs.emplace_back(pf.string());
        }
        innotin_params.output_path = temp_filtered_path;

        stage("innotin")
            .expect_success(innotin(innotin_params))
            .expect_file(temp_filtered_path, file_is_fasta, "FASTA")
            .or_die_if(!config_.soft_fail)
            .or_execute([this]{ soft_fail_exit(); });
    }

    // Run minimap2
    std::vector<std::string> minimap_cmd = {
        "minimap2",
        temp_filtered_path.string()
    };
    for (const auto& pf : plasmid_files) {
        minimap_cmd.emplace_back(pf.string());
    }
    minimap_cmd.emplace_back("-o");
    minimap_cmd.emplace_back(paf_output.string());
    minimap_cmd.emplace_back("-t");
    minimap_cmd.emplace_back(std::to_string(config_.threads));

    auto mm_stage = stage("minimap2 propagation round " + std::to_string(round))
        .expect_which("minimap2")
        .expect_file(temp_filtered_path, Expect::NON_EMPTY);
    for (const auto& pf : plasmid_files) mm_stage.expect_file(pf);
    mm_stage.proc(minimap_cmd)
        .expect_file(paf_output)
        .or_die_if(!config_.soft_fail)
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
    SelectMissingReadsParams select_params;
    select_params.paf_path = plasmid_alignment.string();
    for (const auto& path : unknown_reads) {
        select_params.fastq_paths.emplace_back(path.string());
    }
    select_params.output_path = output_path.string();

    stage("select-missing-reads")
        .expect_success(select_missing_reads(select_params))
        .expect_file(output_path, file_is_gzipped_fastq, "gzipped FASTQ")
        .or_die_if(!config_.soft_fail)
        .or_execute([this]{ soft_fail_exit(); });

    return output_path;
}

std::filesystem::path Pipeline::run_unicycler_lr_assembly(
    const std::vector<std::filesystem::path>& plasmid_files,
    int iteration,
    int comp_id) {

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
    std::filesystem::copy(unicycler_sr_path, unicycler_lr_path,
                          std::filesystem::copy_options::recursive |
                          std::filesystem::copy_options::overwrite_existing);

    // Remove files that will be regenerated
    remove_if_exists(unicycler_lr_path / "assembly.fasta");
    remove_if_exists(unicycler_lr_path / "assembly.gfa");

    // Build tag for temp file naming (includes comp_id when in per-component context)
    auto tag = "iter" + std::to_string(iteration) +
               (comp_id >= 0 ? "_c" + std::to_string(comp_id) : "");

    // Unicycler accepts gzip input. Reuse one stream directly or concatenate
    // gzip members byte-for-byte instead of materializing plain FASTQ.
    std::filesystem::path lr_input;
    if (plasmid_files.size() == 1) {
        lr_input = plasmid_files.front();
    } else {
        lr_input = tmp_dir_ / ("lr_concat_" + tag + ".fastq.gz");
        if (!concat_files_binary(plasmid_files, lr_input)) {
            if (config_.soft_fail) soft_fail_exit();
            log("ERROR", "Failed to concatenate plasmid reads");
            std::exit(EXIT_FAILURE);
        }
    }

    std::vector<std::string> cmd = {
        "unicycler_hyplas_modified",
        "--verbosity", "1",
        "--keep", "3",
        "-o", unicycler_lr_path.string(),
        "-t", std::to_string(config_.threads),
        "-l", lr_input.string()
    };

    // Handle short reads
    if (config_.short_reads.empty()) {
        auto empty_fq1 = make_temp("empty_sr1_" + tag + ".fq");
        auto empty_fq2 = make_temp("empty_sr2_" + tag + ".fq");
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

    stage("unicycler LR assembly iteration " + std::to_string(iteration))
        .expect_which("unicycler_hyplas_modified")
        .expect_file(lr_input)
        .proc(cmd)
        .expect_file(assembly_fasta, file_is_fasta, "FASTA")
        .or_die_if(!config_.soft_fail)
        .or_execute([this]{ soft_fail_exit(); });

    return assembly_fasta;
}

// ============================================================================
// Per-component assembly helpers
// ============================================================================

std::vector<GfaComponent> Pipeline::extract_plasmid_components(
    const std::filesystem::path& gaf_path,
    const std::filesystem::path& prediction_tsv) const
{
    // 1. Parse plasmid-classified contigs
    gtl::flat_hash_set<std::string> plasmid_segs =
        plasmid_names(parse_prediction_tsv(prediction_tsv));

    // 2. Union-find over contigs co-aligned by the same long read (from GAF)
    gtl::flat_hash_map<std::string, std::string> parent;
    std::vector<std::string> all_segs;

    std::function<std::string(const std::string&)> find = [&](const std::string& x) -> std::string {
        if (parent[x] != x) parent[x] = find(parent[x]);
        return parent[x];
    };
    auto unite = [&](const std::string& a, const std::string& b) {
        parent[find(a)] = find(b);
    };

    {
        std::error_code ec;
        mio::mmap_source gaf_mmap = mio::make_mmap_source(gaf_path.string(), ec);
        if (ec) {
            log("ERROR", "Cannot mmap GAF for component extraction: " + ec.message());
            std::exit(EXIT_FAILURE);
        }
        MmapView gaf_view{gaf_mmap};

        for (AlignmentEntry entry{gaf_view}; !entry.is_terminal; entry = AlignmentEntry{gaf_view}) {
            std::vector<std::string> segs;
            segs.reserve(entry.contigs.size());
            for (const LiteView& lv : entry.contigs) {
                segs.emplace_back(static_cast<std::string>(gaf_view.sub(lv)));
            }
            for (const auto& seg : segs) {
                if (parent.find(seg) == parent.end()) {
                    parent[seg] = seg;
                    all_segs.push_back(seg);
                }
            }
            // Unioning consecutive pairs suffices; union-find handles transitivity.
            for (size_t i = 0; i + 1 < segs.size(); ++i) {
                unite(segs[i], segs[i + 1]);
            }
        }
    }

    // 3. Group segments by component root
    gtl::flat_hash_map<std::string, std::vector<std::string>> comp_map;
    for (const auto& seg : all_segs) {
        comp_map[find(seg)].push_back(seg);
    }

    // 4. Keep only components that contain at least one plasmid-classified segment
    std::vector<GfaComponent> result;
    int id = 0;
    for (auto& [root, segs] : comp_map) {
        bool has_plasmid = false;
        for (const auto& s : segs) {
            if (plasmid_segs.count(s)) { has_plasmid = true; break; }
        }
        if (has_plasmid) {
            result.push_back(GfaComponent{id++, std::move(segs)});
        }
    }

    log("INFO", "Found " + std::to_string(result.size()) + " plasmid-containing GAF component(s)");
    return result;
}

gtl::flat_hash_map<int, std::vector<std::string>> Pipeline::bin_reads_to_components(
    const std::filesystem::path& gaf_path,
    const gtl::flat_hash_map<std::string, int>& segment_to_component) const
{
    gtl::flat_hash_map<int, std::vector<std::string>> result;

    std::error_code ec;
    mio::mmap_source gaf_mmap = mio::make_mmap_source(gaf_path.string(), ec);
    if (ec) {
        log("ERROR", "Cannot mmap GAF for read binning: " + ec.message());
        std::exit(EXIT_FAILURE);
    }

    MmapView gaf_view{gaf_mmap};

    for (AlignmentEntry entry{gaf_view}; !entry.is_terminal; entry = AlignmentEntry{gaf_view}) {
        std::string read_id;
        {
            MmapView id_view = gaf_view.sub(entry.id);
            read_id = static_cast<std::string>(id_view);
        }

        gtl::flat_hash_set<int> comps;
        for (const LiteView& v : entry.contigs) {
            MmapView cv = gaf_view.sub(v);
            std::string seg = static_cast<std::string>(cv);
            auto it = segment_to_component.find(seg);
            if (it != segment_to_component.end()) {
                comps.insert(it->second);
            }
        }
        for (int cid : comps) {
            result[cid].push_back(read_id);
        }
    }

    return result;
}

gtl::flat_hash_map<int, std::vector<std::string>> Pipeline::bin_reads_via_minigraph(
    const std::filesystem::path& reads_fastq,
    const gtl::flat_hash_map<std::string, int>& segment_to_component) const
{
    auto sr_graph_fix = output_dir_ / "unicycler_sr" / "assembly_segfix.gfa";
    if (!file_readable(sr_graph_fix)) {
        log("WARNING", "assembly_segfix.gfa not found; skipping minigraph binning of new reads");
        return {};
    }

    auto temp_gaf = make_temp("binning_reads.gaf");
    RunOptions opts;
    opts.stdout_file = temp_gaf;

    (void) stage("minigraph new reads to SR assembly")
        .expect_which("minigraph")
        .expect_file(reads_fastq)
        .proc({
            "minigraph",
            sr_graph_fix.string(),
            reads_fastq.string(),
            "-t", std::to_string(config_.threads),
            "-x", "lr",
            "-c"
        }, opts);

    if (!file_readable(temp_gaf)) {
        return {};
    }

    return bin_reads_to_components(temp_gaf, segment_to_component);
}

void Pipeline::init_component_unicycler_sr(
    const std::filesystem::path& comp_gfa,
    const std::filesystem::path& comp_dir)
{
    auto unicycler_sr = comp_dir / "unicycler_sr";
    ensure_directory(unicycler_sr);

    // Place overlap-removed component GFA as the depth-filter graph.
    // No mock unicycler run — just the graph file, matching setup_from_spades_output.
    // Unicycler will process everything from scratch when given LR reads.
    remove_gfa_overlaps(comp_gfa, unicycler_sr / "002_depth_filter.gfa");
}

std::vector<std::filesystem::path> Pipeline::run_per_component_assembly(
    const std::vector<GfaComponent>& components,
    const gtl::flat_hash_map<int, std::vector<std::filesystem::path>>& comp_read_files,
    int iteration)
{
    auto sr_gfa   = output_dir_ / "unicycler_sr" / "assembly.gfa";
    auto comp_dir = output_dir_ / "components";

    std::vector<std::filesystem::path> assembly_fastas;

    for (const auto& comp : components) {
        auto cdir = comp_dir / ("comp_" + std::to_string(comp.id));

        // Get accumulated read files for this component
        auto rfIt = comp_read_files.find(comp.id);
        if (rfIt == comp_read_files.end() || rfIt->second.empty()) {
            log("INFO", "Component " + std::to_string(comp.id) + ": no LR reads, skipping");
            continue;
        }
        const auto& read_paths = rfIt->second;

        // Set up unicycler directory from component sub-GFA
        auto unicycler_comp = cdir / ("unicycler_lr_" + std::to_string(iteration));
        auto assembly_fasta = unicycler_comp / "assembly.fasta";

        if (!config_.force && file_readable(assembly_fasta)) {
            log("WARNING", "Component " + std::to_string(comp.id) + " iteration " +
                std::to_string(iteration) + " exists, skipping");
            assembly_fastas.push_back(assembly_fasta);
            continue;
        }

        // Initialise per-component SR directory (once, mirrors setup_from_existing_assembly)
        auto comp_gfa = cdir / "component.gfa";
        init_component_unicycler_sr(comp_gfa, cdir);

        // Copy SR dir as starting point for this iteration (mirrors run_unicycler_lr_assembly)
        auto unicycler_sr_comp = cdir / "unicycler_sr";
        std::filesystem::copy(unicycler_sr_comp, unicycler_comp,
                              std::filesystem::copy_options::recursive |
                              std::filesystem::copy_options::overwrite_existing);
        remove_if_exists(unicycler_comp / "assembly.fasta");
        remove_if_exists(unicycler_comp / "assembly.gfa");

        // Concatenate all accumulated gzipped reads into a persistent component-local
        // file. Gzip streams concatenate natively, so byte concat is valid.
        std::filesystem::path lr_input;
        if (read_paths.size() == 1) {
            lr_input = read_paths[0];
        } else {
            lr_input = cdir / ("lr_reads.all.it" + std::to_string(iteration) + ".fastq.gz");
            if (!concat_files_binary(read_paths, lr_input)) {
                log("WARNING", "Component " + std::to_string(comp.id) + ": failed to concat LR reads");
                continue;
            }
        }

        std::vector<std::string> cmd = {
            "unicycler_hyplas_modified",
            "--verbosity", "1",
            "--keep", "3",
            "-o", unicycler_comp.string(),
            "-t", std::to_string(config_.threads),
            "-l", lr_input.string()
        };

        std::filesystem::path empty_fq;
        if (config_.short_reads.empty()) {
            empty_fq = cdir / "empty.fq";
            std::ofstream{empty_fq};
            cmd.emplace_back("-1"); cmd.emplace_back(empty_fq.string());
            cmd.emplace_back("-2"); cmd.emplace_back(empty_fq.string());
        } else {
            cmd.emplace_back("-1"); cmd.emplace_back(config_.short_reads[0].string());
            if (config_.short_reads.size() > 1) {
                cmd.emplace_back("-2"); cmd.emplace_back(config_.short_reads[1].string());
            }
        }

        log("INFO", "Assembling component " + std::to_string(comp.id) +
            " (" + std::to_string(comp.segments.size()) + " contigs) iteration " +
            std::to_string(iteration));

        stage("unicycler component " + std::to_string(comp.id))
            .expect_which("unicycler_hyplas_modified")
            .expect_file(lr_input)
            .proc(cmd)
            .expect_file(assembly_fasta, file_is_fasta, "FASTA")
            .or_die_if(!config_.soft_fail)
            .or_execute([this, &comp](){
                log("WARNING", "Component " + std::to_string(comp.id) + " assembly failed (soft-fail)");
            });

        if (file_readable(assembly_fasta)) {
            assembly_fastas.push_back(assembly_fasta);
        }
    }

    return assembly_fastas;
}

void Pipeline::merge_circular_contigs(
    const std::vector<std::filesystem::path>& fasta_paths,
    const std::filesystem::path& sr_gfa,
    const std::filesystem::path& sr_fasta,
    const std::filesystem::path& prediction_tsv,
    int iteration) const
{
    auto final_path = output_dir_ / ("plasmids.final.it" + std::to_string(iteration) + ".fasta");
    std::ofstream out(final_path);
    if (!out) {
        throw HyplasError("cannot open merged output: " + final_path.string());
    }

    gtl::flat_hash_set<std::string> written;

    for (const auto& p : fasta_paths) {
        if (!file_readable(p)) continue;
        append_circular_by_header(out, p, written);
    }

    append_circular_sr_plasmids(out, sr_gfa, sr_fasta, prediction_tsv, written);
}

// ============================================================================
// Main pipeline execution
// ============================================================================

int Pipeline::run() {
    log("INFO", "Starting HyPlAs pipeline");
    log("INFO", "Output directory: " + output_dir_.string());

    ensure_directory(output_dir_);
    ensure_directory(tmp_dir_);

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
    prediction_tsv_ = process_platon_output(platon_dir);

    // 3. Graph alignment (if long reads provided)
    if (!config_.long_reads) {
        log("INFO", "No long reads provided, skipping propagation");
        return 0;
    }

    log("INFO", "Running minigraph alignment");
    auto graph_alignment = run_minigraph_lr_to_sr();

    // 4. Initial read selection
    log("INFO", "Running initial read selection");
    auto reads = run_long_read_selection(prediction_tsv_, graph_alignment);

    // 5. Check if any plasmid reads found
    size_t plasmid_lines = line_count(reads.plasmid_reads);
    if (plasmid_lines < 4) {
        log("INFO", "No plasmid long reads found");

        // Write circular plasmid contigs from SR assembly
        auto sr_gfa = output_dir_ / "unicycler_sr" / "assembly.gfa";
        auto sr_fasta = output_dir_ / "unicycler_sr" / "assembly.fasta";
        write_circular_plasmid_contigs(sr_gfa, sr_fasta, prediction_tsv_, 0);

        // Symlink remaining iterations
        symlink_remaining_iterations(0);
        return 0;
    }

    if (config_.per_component) {
        // ── Per-component mode ──────────────────────────────────────────────
        // Strategy: discover all plasmid long reads first (initial selection +
        // read-read overlap propagation), map only those reads to the SR graph,
        // then build components from that plasmid-only GAF. This prevents
        // chromosomal reads from collapsing the whole assembly into one
        // component via union-find over co-aligned contigs.

        auto sr_gfa   = output_dir_ / "unicycler_sr" / "assembly.gfa";
        auto sr_fasta = output_dir_ / "unicycler_sr" / "assembly.fasta";
        auto comp_dir = output_dir_ / "components";
        ensure_directory(comp_dir);

        std::vector<std::filesystem::path> unknown_files = {
            reads.unmapped, reads.unknown_both, reads.unknown_neither
        };

        // 1. Pre-grouping read discovery — read-read overlaps only, no assembly
        std::vector<std::filesystem::path> accumulated = {reads.plasmid_reads};
        for (int i = 0; i < config_.propagate_rounds; ++i) {
            log("INFO", "Pre-grouping propagation round " + std::to_string(i + 1));

            auto aln = find_missing_long_reads(accumulated, unknown_files, i);
            if (line_count(aln) == 0) {
                log("INFO", "No alignments in round " + std::to_string(i + 1) +
                            ", stopping pre-grouping discovery");
                break;
            }

            auto new_reads = extract_missing_long_reads(aln, unknown_files);
            if (line_count(new_reads) == 0) {
                log("INFO", "No new reads in round " + std::to_string(i + 1) +
                            ", stopping pre-grouping discovery");
                break;
            }

            accumulated.push_back(new_reads);
        }

        // 2. Concatenate all plasmid reads and map them to the SR graph.
        auto concat_plasmid = comp_dir / "plasmid_reads.all.fastq.gz";
        if (config_.force || !file_readable(concat_plasmid)) {
            if (!concat_files_binary(accumulated, concat_plasmid)) {
                log("ERROR", "Failed to concatenate plasmid reads for component extraction");
                std::exit(EXIT_FAILURE);
            }
        }
        auto plasmid_gaf = comp_dir / "plasmid_lr2assembly.gaf";
        run_minigraph_lr_to_sr(concat_plasmid, plasmid_gaf);

        // 3. Extract components from the plasmid-only GAF
        auto components = extract_plasmid_components(plasmid_gaf, prediction_tsv_);
        if (components.empty()) {
            log("INFO", "No plasmid components found; falling back to SR circular contigs");
            write_circular_plasmid_contigs(sr_gfa, sr_fasta, prediction_tsv_, 0);
            symlink_remaining_iterations(0);
            return 0;
        }

        gtl::flat_hash_map<std::string, int> seg_to_comp;
        for (const auto& c : components)
            for (const auto& s : c.segments) seg_to_comp[s] = c.id;

        // 4. Write component sub-GFAs
        for (const auto& comp : components) {
            auto cdir = comp_dir / ("comp_" + std::to_string(comp.id));
            ensure_directory(cdir);
            write_component_gfa(comp.segments, sr_gfa, cdir / "component.gfa");
        }

        // 5. Bin accumulated reads to components using the same plasmid GAF
        auto bins = bin_reads_to_components(plasmid_gaf, seg_to_comp);

        // 6. Write per-component read files (split the concatenated pool by bin)
        gtl::flat_hash_map<int, std::vector<std::filesystem::path>> comp_read_files;
        for (const auto& comp : components) {
            auto it = bins.find(comp.id);
            if (it == bins.end() || it->second.empty()) continue;

            auto cdir = comp_dir / ("comp_" + std::to_string(comp.id));
            auto comp_lr = cdir / "lr_reads.fastq.gz";
            if (config_.force || !file_readable(comp_lr)) {
                gtl::flat_hash_set<std::string> ids(it->second.begin(), it->second.end());
                write_reads_by_id(concat_plasmid, ids, comp_lr);
            }
            comp_read_files[comp.id] = {comp_lr};
        }

        // 7. Single assembly pass — all reads already discovered
        log("INFO", "Running per-component hybrid assembly");
        auto comp_fastas = run_per_component_assembly(components, comp_read_files, 0);
        merge_circular_contigs(comp_fastas, sr_gfa, sr_fasta, prediction_tsv_, 0);

        // 8. Mirror single-pass output to remaining iteration slots
        symlink_remaining_iterations(0);

    } else {
        // ── Standard mode ───────────────────────────────────────────────────
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

            auto new_reads = extract_missing_long_reads(plasmid_alignment, unknown_files);

            if (line_count(new_reads) == 0) {
                log("INFO", "No new reads in round " + std::to_string(i + 1) + ", stopping propagation");
                symlink_remaining_iterations(i);
                break;
            }

            plasmid_files.emplace_back(new_reads);
            lr_assembly = run_unicycler_lr_assembly(plasmid_files, i + 1);
            write_circular_contigs(lr_assembly, i + 1);
        }
    }

    if (!config_.keep_temp) {
        std::filesystem::remove_all(tmp_dir_);
    } else {
        log("INFO", "Keeping temp directory: " + tmp_dir_.string());
    }

    log("INFO", "Pipeline completed successfully");
    return 0;
}

} // namespace hyplas

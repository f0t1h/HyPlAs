/**
 * @file pipeline.hpp
 * @brief HyPlAs pipeline stage declarations
 */

#ifndef HYPLAS_PIPELINE_HPP
#define HYPLAS_PIPELINE_HPP

#include <filesystem>
#include <iosfwd>
#include <optional>
#include <string>
#include <gtl/phmap.hpp>
#include <vector>

namespace hyplas {

/**
 * @brief Pipeline configuration (from CLI arguments)
 */
struct PipelineConfig {
    // Required paths
    std::filesystem::path output_directory;
    std::filesystem::path platon_db;
    
    // Optional input paths
    std::optional<std::filesystem::path> long_reads;
    std::vector<std::filesystem::path> short_reads;
    std::optional<std::filesystem::path> sr_assembly;  // --sr-assembly shortcut
    
    // Parameters
    int propagate_rounds = 0;
    int threads = 16;
    bool force = false;
    std::string verbosity = "INFO";
    
    // Assembly options
    bool use_spades = false;  // Use SPAdes directly instead of Unicycler for SR assembly
    // SPAdes defaults: -k 53 --gfa11 --isolate -m 1024
    bool soft_fail = false;       // Fall back to SR circular contigs on post-assembly failure
    bool per_component = false;   // Bin SR contigs into GFA-connected components and assemble each independently
    bool keep_temp = false;       // Keep tmp/ directory after pipeline completes (for debugging)
};

/**
 * @brief A connected component of SR contigs in the assembly graph
 */
struct GfaComponent {
    int id;
    std::vector<std::string> segments;  // contig names in this component
};

/**
 * @brief Result of read selection stage
 */
struct ReadSelectionResult {
    std::filesystem::path plasmid_reads;
    std::filesystem::path unknown_both;
    std::filesystem::path unknown_neither;
    std::filesystem::path unmapped;
};

/**
 * @brief Pipeline execution context
 */
class Pipeline {
public:
    explicit Pipeline(const PipelineConfig& config);
    
    /// Run the complete pipeline
    int run();
    
    /// Log a message
    void log(const std::string& level, const std::string& message) const;

private:
    const PipelineConfig& config_;
    std::filesystem::path output_dir_;
    std::filesystem::path tmp_dir_;
    std::filesystem::path prediction_tsv_;

    std::filesystem::path make_temp(const std::string& name) const;
    
    // Pipeline stages
    std::filesystem::path run_unicycler_sr_assembly();
    std::filesystem::path run_spades_sr_assembly();  // SPAdes-only SR assembly
    void setup_from_existing_assembly();
    void setup_from_spades_output();  // Setup from SPAdes GFA output
    std::filesystem::path run_platon_classifier();
    std::filesystem::path process_platon_output(const std::filesystem::path& platon_dir);
    std::filesystem::path run_minigraph_lr_to_sr(
        const std::filesystem::path& reads_fastq,
        const std::filesystem::path& gaf_output);
    std::filesystem::path run_minigraph_lr_to_sr();  // uses config_.long_reads → lr2assembly.gaf
    ReadSelectionResult run_long_read_selection(
        const std::filesystem::path& prediction_tsv,
        const std::filesystem::path& graph_alignment);
    std::filesystem::path find_missing_long_reads(
        const std::vector<std::filesystem::path>& plasmid_files,
        const std::vector<std::filesystem::path>& unknown_files,
        int round,
        int comp_id = -1);
    std::filesystem::path extract_missing_long_reads(
        const std::filesystem::path& plasmid_alignment,
        const std::vector<std::filesystem::path>& unknown_reads);
    std::filesystem::path run_unicycler_lr_assembly(
        const std::vector<std::filesystem::path>& plasmid_files,
        int iteration,
        int comp_id = -1);
    
    // Per-component assembly
    std::vector<GfaComponent> extract_plasmid_components(
        const std::filesystem::path& gaf_path,
        const std::filesystem::path& prediction_tsv) const;
    gtl::flat_hash_map<int, std::vector<std::string>> bin_reads_to_components(
        const std::filesystem::path& gaf_path,
        const gtl::flat_hash_map<std::string, int>& segment_to_component) const;
    gtl::flat_hash_map<int, std::vector<std::string>> bin_reads_via_minigraph(
        const std::filesystem::path& reads_fastq,
        const gtl::flat_hash_map<std::string, int>& segment_to_component) const;
    void init_component_unicycler_sr(
        const std::filesystem::path& comp_gfa,
        const std::filesystem::path& comp_dir);
    std::vector<std::filesystem::path> run_per_component_assembly(
        const std::vector<GfaComponent>& components,
        const gtl::flat_hash_map<int, std::vector<std::filesystem::path>>& comp_read_files,
        int iteration);
    void merge_circular_contigs(
        const std::vector<std::filesystem::path>& fasta_paths,
        const std::filesystem::path& sr_gfa,
        const std::filesystem::path& sr_fasta,
        const std::filesystem::path& prediction_tsv,
        int iteration) const;

    // Final-output helpers (thin wrappers over the free gfa transforms)
    void write_circular_plasmid_contigs(
        const std::filesystem::path& gfa_path,
        const std::filesystem::path& fasta_path,
        const std::filesystem::path& prediction_tsv,
        int iteration) const;
    void write_circular_contigs(const std::filesystem::path& assembly_fasta,
                                int iteration) const;
    void symlink_remaining_iterations(int from_iteration);
    [[noreturn]] void soft_fail_exit();
};

/**
 * @brief List of required external tools
 * 
 * Includes direct dependencies plus subtools used by Unicycler:
 * - SPAdes: called by Unicycler for short-read assembly
 * - Racon: called by Unicycler for long-read polishing
 * - BLAST+ (makeblastdb, tblastn): called by Unicycler for contig rotation
 * 
 * Note: hyplas-utils subprograms (innotin, split-plasmid-reads, select-missing-reads)
 * are now baked into hyplas-pipeline and no longer require external invocation.
 */
inline const std::vector<std::string>& required_tools(bool use_spades = false) {
    static const std::vector<std::string> tools_unicycler = {
        // Direct dependencies
        "unicycler_hyplas_modified",
        "platon",
        "minigraph",
        "minimap2",
        // Unicycler subtools (called internally by Unicycler)
        "spades.py",      // SR assembly
        "racon",          // LR polishing
        "makeblastdb",    // Contig rotation (BLAST+)
        "tblastn"         // Contig rotation (BLAST+)
    };
    static const std::vector<std::string> tools_spades = {
        // Direct dependencies (no Unicycler for SR assembly)
        "spades.py",
        "unicycler_hyplas_modified",  // Still needed for LR assembly step
        "platon",
        "minigraph",
        "minimap2",
        // Unicycler subtools (still needed for hybrid assembly step)
        "racon",
        "makeblastdb",
        "tblastn"
    };
    return use_spades ? tools_spades : tools_unicycler;
}

/**
 * @brief Run tool validation (--check-deps)
 * @param use_spades If true, check for SPAdes-based dependency list
 * @return 0 if all tools found, 1 otherwise
 */
int check_dependencies(bool use_spades = false);

} // namespace hyplas

#endif // HYPLAS_PIPELINE_HPP

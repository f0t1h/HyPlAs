/**
 * @file hyplas_main.cpp
 * @brief HyPlAs pipeline orchestrator - main entry point
 */

#include "CLI11.hpp"
#include "pipeline.hpp"
#include "process.hpp"

#include <cstdio>
#include <cstdlib>

int main(int argc, char* argv[]) {
    CLI::App app{"HyPlAs - Hybrid Plasmid Assembly pipeline"};
    
    hyplas::PipelineConfig config;
    bool check_deps = false;
    
    // Utility flags (checked before required args)
    app.add_flag("--check-deps", check_deps,
        "Check that all required tools are available and exit");
    
    // Required arguments (but not when --check-deps is used)
    app.add_option("-o,--output-directory", config.output_directory,
        "Output directory for pipeline results");
    
    app.add_option("--platon-db", config.platon_db,
        "Path to Platon database")
        ->check(CLI::ExistingDirectory);
    
    // Input options
    app.add_option("-l,--long-reads", config.long_reads,
        "Long reads FASTQ file")
        ->check(CLI::ExistingFile);
    
    app.add_option("-s,--short-reads", config.short_reads,
        "Short reads FASTQ files (1 or 2 files)")
        ->check(CLI::ExistingFile);
    
    app.add_option("--sr-assembly", config.sr_assembly,
        "Pre-computed short-read assembly graph (skip SR assembly step)")
        ->check(CLI::ExistingFile);
    
    // Parameters
    app.add_option("-p,--propagate-rounds", config.propagate_rounds,
        "Number of propagation rounds")
        ->default_val(0);
    
    app.add_option("-t,--threads", config.threads,
        "Number of threads")
        ->default_val(16);
    
    app.add_option("--verbosity", config.verbosity,
        "Logging verbosity")
        ->check(CLI::IsMember({"DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"}))
        ->default_val("INFO");
    
    app.add_flag("--force", config.force,
        "Force re-run of stages even if outputs exist");
    
    // Assembly options
    app.add_flag("--use-spades", config.use_spades,
        "Use SPAdes directly for SR assembly instead of Unicycler (default: -k 53 --gfa11 --isolate -m 1024)");
    
    // Parse with allow_extras to handle --check-deps before validation
    try {
        app.parse(argc, argv);
    } catch (const CLI::ParseError& e) {
        // If --check-deps was given, run it anyway
        if (check_deps) {
            return hyplas::check_dependencies(config.use_spades);
        }
        return app.exit(e);
    }
    
    // Handle --check-deps (doesn't need other args)
    if (check_deps) {
        return hyplas::check_dependencies(config.use_spades);
    }
    
    // Now validate required arguments for pipeline run
    if (config.output_directory.empty()) {
        std::fprintf(stderr, "Error: --output-directory is required\n");
        return EXIT_FAILURE;
    }
    
    if (config.platon_db.empty()) {
        std::fprintf(stderr, "Error: --platon-db is required\n");
        return EXIT_FAILURE;
    }
    
    // Validate inputs
    if (!config.long_reads && config.short_reads.empty() && !config.sr_assembly) {
        std::fprintf(stderr, "Error: Must provide --long-reads, --short-reads, or --sr-assembly\n");
        return EXIT_FAILURE;
    }
    
    if (config.short_reads.empty() && !config.sr_assembly) {
        std::fprintf(stderr, "Warning: No short reads provided\n");
    }
    
    if (!config.long_reads) {
        std::fprintf(stderr, "Warning: No long reads provided\n");
    }
    
    // Validate tools before starting
    std::fprintf(stderr, "[INFO] Checking required tools...\n");
    if (!hyplas::validate_tools(hyplas::required_tools(config.use_spades), false)) {
        std::fprintf(stderr, "[ERROR] Missing required tools. Run with --check-deps for details.\n");
        return EXIT_FAILURE;
    }
    
    // Run pipeline
    hyplas::Pipeline pipeline(config);
    return pipeline.run();
}

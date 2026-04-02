/**
 * @file main.cpp
 * @brief HyPlAs C++ utilities - unified command-line interface
 */

#include "CLI11.hpp"
#include "subprograms.hpp"

int main(int argc, char* argv[]) {
    CLI::App app{"HyPlAs - Hybrid Plasmid Assembly utilities"};
    app.require_subcommand(1);

    // innotin
    hyplas::InnotinParams innotin_params;
    auto* innotin_cmd = app.add_subcommand(
        "innotin", 
        "Filter reads from main FASTQ that are NOT in any subset FASTQ"
    );
    innotin_cmd->add_option("main_fastq", innotin_params.main_fastq, 
        "Main FASTQ file to filter")->required()->check(CLI::ExistingFile);
    innotin_cmd->add_option("subset_fastqs", innotin_params.subset_fastqs, 
        "Subset FASTQ files to exclude")->required()->check(CLI::ExistingFile);

    // select-missing-reads
    hyplas::SelectMissingReadsParams select_params;
    auto* select_cmd = app.add_subcommand(
        "select-missing-reads",
        "Select reads that overlap with plasmid reads based on PAF alignment"
    );
    select_cmd->add_option("-p,--paf", select_params.paf_path,
        "PAF alignment file")->required()->check(CLI::ExistingFile);
    select_cmd->add_option("-g,--gaf", select_params.gaf_path,
        "GAF alignment file")->required()->check(CLI::ExistingFile);
    select_cmd->add_option("-f,--fastq", select_params.fastq_path,
        "Input FASTQ file")->required()->check(CLI::ExistingFile);
    select_cmd->add_option("-o,--output", select_params.output_path,
        "Output FASTQ file (gzipped)")->required();
    select_cmd->add_option("-t,--tsv", select_params.prediction_path,
        "Prediction TSV file")->required()->check(CLI::ExistingFile);

    // split-plasmid-reads
    hyplas::SplitPlasmidReadsParams split_params;
    auto* split_cmd = app.add_subcommand(
        "split-plasmid-reads",
        "Split reads based on alignment to plasmid/chromosome contigs"
    );
    split_cmd->add_option("-g,--gaf", split_params.gaf_path,
        "GAF alignment file")->required()->check(CLI::ExistingFile);
    split_cmd->add_option("-f,--fastq", split_params.fastq_path,
        "Input FASTQ file")->required()->check(CLI::ExistingFile);
    split_cmd->add_option("-t,--tsv", split_params.prediction_path,
        "Prediction TSV file")->required()->check(CLI::ExistingFile);
    split_cmd->add_option("-p,--plasmid-out", split_params.plasmid_out_path,
        "Output file for plasmid reads (gzipped)")->required();
    split_cmd->add_option("-c,--chr-out", split_params.chr_out_path,
        "Output file for chromosome reads (gzipped)");
    split_cmd->add_option("-n,--neither-out", split_params.unknown_neither_path,
        "Output file for reads aligning to neither (gzipped)")->required();
    split_cmd->add_option("-b,--both-out", split_params.unknown_both_path,
        "Output file for reads aligning to both (gzipped)")->required();
    split_cmd->add_option("-u,--unmapped-out", split_params.unmapped_path,
        "Output file for unmapped reads (gzipped)")->required();
    split_cmd->add_flag("--output-chromosomal", split_params.output_chromosomal,
        "Enable output of chromosome-only reads");

    CLI11_PARSE(app, argc, argv);

    if (innotin_cmd->parsed()) {
        return hyplas::run_innotin(innotin_params);
    }
    if (select_cmd->parsed()) {
        return hyplas::run_select_missing_reads(select_params);
    }
    if (split_cmd->parsed()) {
        return hyplas::run_split_plasmid_reads(split_params);
    }

    return EXIT_FAILURE;
}

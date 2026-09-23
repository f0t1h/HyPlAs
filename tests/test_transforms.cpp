/**
 * @file test_transforms.cpp
 * @brief Unit smoke tests for the pure GFA/FASTA/read transforms.
 *
 * These exercise the modules extracted from pipeline.cpp (gfa.*, fastx.*) on
 * small synthetic fixtures with hand-verified expected outputs. They require
 * no external bioinformatics tools. Build/run via `make test`.
 */

#include "fastx.hpp"
#include "gfa.hpp"
#include "contig_classification.hpp"

#include <array>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <gtl/phmap.hpp>
#include <vector>

#include <zlib.h>

namespace fs = std::filesystem;

static int g_failures = 0;
static int g_checks = 0;

#define CHECK(cond, msg)                                                    \
    do {                                                                    \
        ++g_checks;                                                         \
        if (!(cond)) {                                                      \
            ++g_failures;                                                   \
            std::fprintf(stderr, "  FAIL: %s (%s:%d)\n", msg, __FILE__, __LINE__); \
        }                                                                   \
    } while (0)

static void write_file(const fs::path& p, const std::string& content) {
    std::ofstream out(p);
    out << content;
}

static std::string read_file(const fs::path& p) {
    std::ifstream in(p);
    std::ostringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

static std::string read_gz(const fs::path& p) {
    gzFile gz = gzopen(p.string().c_str(), "rb");
    if (!gz) return {};
    std::string out;
    char buf[8192];
    int n;
    while ((n = gzread(gz, buf, sizeof(buf))) > 0) out.append(buf, static_cast<size_t>(n));
    gzclose(gz);
    return out;
}

static bool contains(const std::string& hay, const std::string& needle) {
    return hay.find(needle) != std::string::npos;
}

// ---------------------------------------------------------------------------

static void test_fix_gfa_empty_segments(const fs::path& dir) {
    std::fprintf(stderr, "test_fix_gfa_empty_segments\n");
    auto in = dir / "empty.gfa";
    auto out = dir / "empty.fixed.gfa";
    write_file(in,
        "H\tVN:Z:1.0\n"
        "S\ts1\tAAA\n"
        "S\ts2\t\n"          // empty segment, must be dropped and bypassed
        "S\ts3\tGGG\n"
        "L\ts1\t+\ts2\t+\t0M\n"
        "L\ts2\t+\ts3\t+\t0M\n");

    CHECK(hyplas::fix_gfa_empty_segments(in, out) == 0, "fix_gfa_empty_segments succeeds");
    std::string r = read_file(out);

    CHECK(contains(r, "S\ts1\tAAA"), "keeps non-empty s1");
    CHECK(contains(r, "S\ts3\tGGG"), "keeps non-empty s3");
    CHECK(!contains(r, "S\ts2\t"), "drops empty s2 segment");
    CHECK(contains(r, "L\ts1\t+\ts3\t+\t0M"), "creates bypass link s1->s3");
    CHECK(!contains(r, "\ts2\t+\t"), "no links reference removed s2");
}

static void test_remove_gfa_overlaps(const fs::path& dir) {
    std::fprintf(stderr, "test_remove_gfa_overlaps\n");

    // Odd overlap (3) so large_half(2) != small_half(1) and grouping is visible.
    auto in = dir / "ovl.gfa";
    auto out = dir / "ovl.trim.gfa";
    write_file(in,
        "H\tVN:Z:1.0\n"
        "S\ts1\tAAAGGG\n"
        "S\ts2\tGGGTTT\n"
        "L\ts1\t+\ts2\t+\t3M\n");

    CHECK(hyplas::remove_gfa_overlaps(in, out) == 0, "remove_gfa_overlaps succeeds");
    std::string r = read_file(out);

    // Deterministic golden output of the asymmetric-trim algorithm:
    // s1 trimmed 1/1 -> "AAGG"; s2 trimmed 2/1 -> "GTT"; total junction trim == overlap.
    CHECK(contains(r, "S\ts1\tAAGG"), "s1 trimmed to AAGG");
    CHECK(contains(r, "S\ts2\tGTT"), "s2 trimmed to GTT");
    CHECK(contains(r, "L\ts1\t+\ts2\t+\t0M"), "link rewritten to 0M");
    CHECK(!contains(r, "3M"), "no residual 3M overlap");

    // No-overlap passthrough: 0M input keeps sequences verbatim.
    auto in0 = dir / "noovl.gfa";
    auto out0 = dir / "noovl.out.gfa";
    write_file(in0,
        "S\ta\tAAAA\n"
        "S\tb\tCCCC\n"
        "L\ta\t+\tb\t+\t0M\n");
    CHECK(hyplas::remove_gfa_overlaps(in0, out0) == 0, "passthrough succeeds");
    std::string r0 = read_file(out0);
    CHECK(contains(r0, "S\ta\tAAAA"), "passthrough keeps a");
    CHECK(contains(r0, "S\tb\tCCCC"), "passthrough keeps b");
    CHECK(contains(r0, "L\ta\t+\tb\t+\t0M"), "passthrough keeps 0M link");
}

static void test_extract_fasta_from_gfa(const fs::path& dir) {
    std::fprintf(stderr, "test_extract_fasta_from_gfa\n");
    auto in = dir / "extract.gfa";
    auto out = dir / "extract.fasta";
    std::string longseq(250, 'C');
    write_file(in,
        "S\tshort\tAAA\n"
        "S\tlong\t" + longseq + "\n"
        "L\tshort\t+\tlong\t+\t0M\n");

    CHECK(hyplas::extract_fasta_from_gfa(in, out, 200) == 0, "extract_fasta_from_gfa succeeds");
    std::string r = read_file(out);
    CHECK(contains(r, ">long"), "emits long segment");
    CHECK(contains(r, longseq), "emits long sequence");
    CHECK(!contains(r, ">short"), "filters out short segment (< min_length)");
}

static void test_prediction_tsv(const fs::path& dir) {
    std::fprintf(stderr, "test_prediction_tsv\n");
    auto tsv = dir / "pred.tsv";
    write_file(tsv,
        "ID\tPREDICTION\n"
        "p1\tplasmid\n"
        "p2\tplasmid\n"
        "c1\tchromosome\n");

    auto parsed = hyplas::parse_prediction_tsv(tsv);
    CHECK(parsed.has_value(), "prediction TSV parsed");
    const auto& m = *parsed;
    CHECK(m.at("p1") == hyplas::ContigType::Plasmid, "p1 plasmid");
    CHECK(m.at("c1") == hyplas::ContigType::Chromosome, "c1 chromosome");
    auto names = hyplas::plasmid_names(m);
    CHECK(names.count("p1") && names.count("p2"), "plasmid_names has p1,p2");
    CHECK(!names.count("c1"), "plasmid_names excludes c1");
}

static void test_circular_by_header(const fs::path& dir) {
    std::fprintf(stderr, "test_circular_by_header\n");
    auto fasta = dir / "circ.fasta";
    write_file(fasta,
        ">c1 length=4 circular=true\nACGT\n"
        ">c2 linear\nTTTT\n"
        ">c3 circular\nGGGG\n");

    std::ostringstream out;
    gtl::flat_hash_set<std::string> written;
    CHECK(hyplas::append_circular_by_header(out, fasta, written) == 0, "append_circular_by_header succeeds");
    std::string r = out.str();
    CHECK(contains(r, ">c1 length=4 circular=true"), "keeps circular c1");
    CHECK(contains(r, ">c3 circular"), "keeps circular c3");
    CHECK(!contains(r, ">c2"), "drops linear c2");
    CHECK(written.count("c1") && written.count("c3"), "records written names");
}

static void test_circular_sr_plasmids(const fs::path& dir) {
    std::fprintf(stderr, "test_circular_sr_plasmids\n");
    auto gfa = dir / "sr.gfa";
    auto fasta = dir / "sr.fasta";
    auto tsv = dir / "sr_pred.tsv";
    write_file(gfa,
        "S\tp1\tAAAA\n"
        "S\tp2\tCCCC\n"
        "S\tc1\tGGGG\n"
        "L\tp1\t+\tp1\t+\t0M\n"   // self-loop -> circular
        "L\tp2\t+\tc1\t+\t0M\n");
    write_file(fasta,
        ">p1\nAAAA\n"
        ">p2\nCCCC\n"
        ">c1\nGGGG\n");
    write_file(tsv,
        "ID\tPREDICTION\n"
        "p1\tplasmid\n"
        "p2\tplasmid\n"
        "c1\tchromosome\n");

    std::ostringstream out;
    gtl::flat_hash_set<std::string> written;
    CHECK(hyplas::append_circular_sr_plasmids(out, gfa, fasta, tsv, written) == 0, "append_circular_sr_plasmids succeeds");
    std::string r = out.str();
    CHECK(contains(r, ">p1 circular"), "p1 circular+plasmid emitted with suffix");
    CHECK(contains(r, "AAAA"), "p1 sequence emitted");
    CHECK(!contains(r, ">p2"), "p2 not circular -> excluded");
    CHECK(!contains(r, ">c1"), "c1 not plasmid -> excluded");
}

static void test_split_plasmid_reads(const fs::path& dir) {
    std::fprintf(stderr, "test_split_plasmid_reads\n");
    auto gaf = dir / "aln.gaf";
    auto tsv = dir / "split_pred.tsv";
    auto fq = dir / "reads.fastq.gz";

    // GAF: column 6 (index 5) is the path; tokens split on '<'/'>'.
    write_file(gaf,
        "readP\t100\t0\t100\t+\t>p1>p2\t60\n"
        "readN\t100\t0\t100\t+\t>x1\t60\n"
        "readW\t100\t0\t100\t+\t>p1>c1\t60\n");
    write_file(tsv,
        "ID\tPREDICTION\n"
        "p1\tplasmid\n"
        "p2\tplasmid\n"
        "c1\tchromosome\n");

    // FASTQ, interleaved so readU1 is consumed as unmapped before readP.
    {
        hyplas::GzWriter w;
        CHECK(w.open(fq), "open gzip output");
        w.write_fastq("readU1", "", "AAAA", "IIII");
        w.write_fastq("readP", "", "CCCC", "IIII");
        w.write_fastq("readN", "", "GGGG", "IIII");
        w.write_fastq("readW", "", "TTTT", "IIII");
        w.write_fastq("readTail", "", "ACGT", "IIII");
    }

    hyplas::SplitPlasmidReadsParams p;
    p.gaf_path = gaf.string();
    p.fastq_path = fq.string();
    p.prediction_path = tsv.string();
    p.plasmid_out_path = (dir / "out.plasmid.fastq.gz").string();
    p.unknown_neither_path = (dir / "out.neither.fastq.gz").string();
    p.unknown_both_path = (dir / "out.both.fastq.gz").string();
    p.unmapped_path = (dir / "out.unmapped.fastq.gz").string();

    int rc = hyplas::split_plasmid_reads(p);
    CHECK(rc == 0, "split_plasmid_reads returns success");

    std::string plasmid = read_gz(p.plasmid_out_path);
    std::string neither = read_gz(p.unknown_neither_path);
    std::string both = read_gz(p.unknown_both_path);
    std::string unmapped = read_gz(p.unmapped_path);

    CHECK(contains(plasmid, "readP") && !contains(plasmid, "readN"),
          "plasmid stream = plasmid-only read");
    CHECK(contains(neither, "readN"), "neither stream = unclassified contig read");
    CHECK(contains(both, "readW"), "both stream = plasmid+chromosome read");
    CHECK(contains(unmapped, "readU1"), "unmapped stream = read with no GAF entry");
}

static void test_select_missing_reads(const fs::path& dir) {
    std::fprintf(stderr, "test_select_missing_reads\n");
    auto paf = dir / "prop.paf";
    auto fq1 = dir / "unknown_1.fastq.gz";
    auto fq2 = dir / "unknown_2.fastq.gz";
    // PAF: target name in column 6 (index 5), no orientation markers.
    write_file(paf, "plasmidRead\t100\t0\t100\t+\treadX\t500\tcg:Z:100M\n");
    {
        hyplas::GzWriter w;
        CHECK(w.open(fq1), "open gzip output");
        w.write_fastq("readY", "", "CCCC", "IIII");
    }
    {
        hyplas::GzWriter w;
        CHECK(w.open(fq2), "open gzip output");
        w.write_fastq("readX", "", "AAAA", "IIII");
    }

    hyplas::SelectMissingReadsParams p;
    p.paf_path = paf.string();
    p.fastq_paths = {fq1.string(), fq2.string()};
    p.output_path = (dir / "selected.fastq.gz").string();

    int rc = hyplas::select_missing_reads(p);
    CHECK(rc == 0, "select_missing_reads returns success");
    std::string out = read_gz(p.output_path);
    CHECK(contains(out, "readX"), "selects PAF target from second input");
    CHECK(!contains(out, "readY"), "excludes unreferenced read from first input");
}

static void test_innotin(const fs::path& dir) {
    std::fprintf(stderr, "test_innotin\n");
    auto main_fq1 = dir / "innotin_main_1.fastq.gz";
    auto main_fq2 = dir / "innotin_main_2.fastq.gz";
    auto subset_fq = dir / "innotin_subset.fastq.gz";
    {
        hyplas::GzWriter w;
        CHECK(w.open(main_fq1), "open gzip output");
        w.write_fastq("r1", "", "AAAA", "IIII");
        w.write_fastq("r2", "", "CCCC", "IIII");
    }
    {
        hyplas::GzWriter w;
        CHECK(w.open(main_fq2), "open gzip output");
        w.write_fastq("r3", "", "GGGG", "IIII");
    }
    {
        hyplas::GzWriter w;
        CHECK(w.open(subset_fq), "open gzip output");
        w.write_fastq("r2", "", "CCCC", "IIII");
    }

    hyplas::InnotinParams p;
    p.main_fastqs = {main_fq1.string(), main_fq2.string()};
    p.subset_fastqs = {subset_fq.string()};
    p.output_path = (dir / "innotin_out.fasta").string();

    int rc = hyplas::innotin(p);
    CHECK(rc == 0, "innotin returns success");
    std::string out = read_file(p.output_path);
    CHECK(contains(out, ">r1") && contains(out, ">r3"),
          "streams reads from every main input");
    CHECK(!contains(out, ">r2"), "drops read present in gzip subset");
}

static void test_byte_view() {
    std::fprintf(stderr, "test_byte_view\n");
    std::string s = "42\tfoo\t3.5";
    mview::byte_view v{std::string_view{s}};

    v.extend_until('\t');
    CHECK(v.to<int>().value_or(-1) == 42, "to<int> parses first token");
    CHECK(static_cast<std::string>(v) == "42", "token string is 42");

    v.skip_next('\t');
    v.extend_until('\t');
    CHECK(static_cast<std::string>(v) == "foo", "second token is foo");
    CHECK(!v.to<int>().has_value(), "to<int> fails on non-numeric token");
    CHECK(v == std::string_view("foo"), "operator== against string_view");

    v.skip_next('\t');
    v.extend_until("\t\n");
    auto d = v.to<double>();
    CHECK(d.has_value() && *d == 3.5, "to<double> parses last token");
}

static void test_split_fields() {
    std::fprintf(stderr, "test_split_fields\n");
    std::array<std::string_view, 4> f{};

    auto n = mview::split_fields("a\tb\tc\td\te", '\t', f);
    CHECK(n == 4 && f[0] == "a" && f[1] == "b" && f[2] == "c" && f[3] == "d\te",
          "last slot captures the remainder");

    std::array<std::string_view, 4> g{};
    auto m = mview::split_fields("x\ty", '\t', g);
    CHECK(m == 2 && g[0] == "x" && g[1] == "y", "fewer fields than slots");

    std::array<std::string_view, 4> h{};
    auto k = mview::split_fields("p\t\tq", '\t', h);
    CHECK(k == 3 && h[0] == "p" && h[1] == "" && h[2] == "q", "empty middle field");
}

static void test_alignment_entry() {
    std::fprintf(stderr, "test_alignment_entry\n");

    // GAF: id = col0, path = col6, tokens split on '<'/'>'.
    std::string gaf = "r1\t1\t2\t3\t4\t>s1>s2\tx\nr2\t1\t2\t3\t4\t<s3\ty\n";
    mview::byte_view gv{std::string_view{gaf}};

    hyplas::AlignmentEntry e1{gv};
    CHECK(!e1.is_terminal, "first GAF entry present");
    CHECK(static_cast<std::string>(gv.sub(e1.id)) == "r1", "entry 1 id");
    CHECK(e1.contigs.size() == 2 &&
          static_cast<std::string>(gv.sub(e1.contigs[0])) == "s1" &&
          static_cast<std::string>(gv.sub(e1.contigs[1])) == "s2",
          "entry 1 splits >s1>s2");

    hyplas::AlignmentEntry e2{gv};
    CHECK(static_cast<std::string>(gv.sub(e2.id)) == "r2", "entry 2 id");
    CHECK(e2.contigs.size() == 1 &&
          static_cast<std::string>(gv.sub(e2.contigs[0])) == "s3",
          "entry 2 splits <s3");

    hyplas::AlignmentEntry e3{gv};
    CHECK(e3.is_terminal, "GAF exhausted");

    // PAF: target name in col6 with no orientation marker.
    std::string paf = "q1\t1\t2\t3\t4\ttgt1\tcg:Z:100M\n";
    mview::byte_view pv{std::string_view{paf}};
    hyplas::AlignmentEntry p1{pv};
    CHECK(p1.contigs.size() == 1 &&
          static_cast<std::string>(pv.sub(p1.contigs[0])) == "tgt1",
          "PAF target captured without <> markers");
}

static void test_alignment_entry_short_row() {
    std::fprintf(stderr, "test_alignment_entry_short_row\n");

    // A short row (< 6 columns) between two valid rows. The fix must confine
    // tab-skipping to the row: the short row yields no contigs and must NOT
    // consume the following valid row.
    std::string gaf =
        "r1\t1\t2\t3\t4\t>s1\tx\n"
        "bad\tonly\ttwo\n"
        "r2\t1\t2\t3\t4\t>s2\ty\n";
    mview::byte_view gv{std::string_view{gaf}};

    hyplas::AlignmentEntry e1{gv};
    CHECK(static_cast<std::string>(gv.sub(e1.id)) == "r1" &&
          e1.contigs.size() == 1 &&
          static_cast<std::string>(gv.sub(e1.contigs[0])) == "s1",
          "valid row before short row unaffected");

    hyplas::AlignmentEntry e2{gv};
    CHECK(static_cast<std::string>(gv.sub(e2.id)) == "bad" && e2.contigs.empty(),
          "short row yields no contigs (no cross-row scan)");

    hyplas::AlignmentEntry e3{gv};
    CHECK(static_cast<std::string>(gv.sub(e3.id)) == "r2" &&
          e3.contigs.size() == 1 &&
          static_cast<std::string>(gv.sub(e3.contigs[0])) == "s2",
          "valid row after short row not consumed");

    hyplas::AlignmentEntry e4{gv};
    CHECK(e4.is_terminal, "GAF exhausted");
}

static void test_write_component_gfa(const fs::path& dir) {
    std::fprintf(stderr, "test_write_component_gfa\n");
    auto src = dir / "wc_src.gfa";
    auto out = dir / "wc_out.gfa";
    write_file(src,
        "H\tVN:Z:1.0\n"
        "S\ts1\tAAAA\n"
        "S\ts2\tCCCC\n"
        "S\ts3\tGGGG\n"
        "L\ts1\t+\ts2\t+\t0M\n"
        "L\ts2\t+\ts3\t+\t0M\n"
        "P\tpath1\ts1+,s2+\t*\n");
    CHECK(hyplas::write_component_gfa({"s1", "s2"}, src, out) == 0, "write_component_gfa succeeds");
    std::string r = read_file(out);
    CHECK(contains(r, "S\ts1\tAAAA") && contains(r, "S\ts2\tCCCC"), "keeps component segments");
    CHECK(!contains(r, "S\ts3"), "drops non-component segment");
    CHECK(contains(r, "L\ts1\t+\ts2\t+\t0M"), "keeps intra-component link");
    CHECK(!contains(r, "s2\t+\ts3"), "drops link leaving component");
    CHECK(!contains(r, "P\tpath1"), "drops P lines");
    CHECK(contains(r, "H\tVN:Z:1.0"), "passes through header line");
}

static void test_classify_platon() {
    std::fprintf(stderr, "test_classify_platon\n");
    std::string tsv =
        "ID\tRDS\tCircular\tInc Type(s)\t# Replication\t# Mobilization\t# OriT\t# Plasmid Hits\t# rRNAs\n"
        "c_plasmid\t2.0\tno\t\t0\t0\t0\t0\t0\n"
        "c_chr\t-9.0\tno\t\t0\t0\t0\t0\t0\n"
        "c_circ\t0.0\tyes\t\t0\t0\t0\t0\t0\n"
        "c_amb\t0.0\tno\t\t0\t0\t0\t0\t0\n"
        "c_bad\n";
    auto classified = hyplas::classify_platon_tsv(tsv);
    CHECK(classified.has_value(), "RDS column found");
    std::string out = classified.value_or("");
    CHECK(out.rfind("ID\tPREDICTION\n", 0) == 0, "output header present");
    CHECK(contains(out, "c_plasmid\tplasmid"), "high RDS -> plasmid");
    CHECK(contains(out, "c_chr\tchromosome"), "low RDS -> chromosome");
    CHECK(contains(out, "c_circ\tplasmid"), "ambiguous + circular -> plasmid");
    CHECK(!contains(out, "c_amb"), "ambiguous with no indicators omitted");
    CHECK(!contains(out, "c_bad"), "row missing RDS column skipped");
}

int main() {
    fs::path dir = fs::temp_directory_path() / "hyplas_tests";
    fs::remove_all(dir);
    fs::create_directories(dir);

    test_fix_gfa_empty_segments(dir);
    test_remove_gfa_overlaps(dir);
    test_extract_fasta_from_gfa(dir);
    test_prediction_tsv(dir);
    test_circular_by_header(dir);
    test_circular_sr_plasmids(dir);
    test_split_plasmid_reads(dir);
    test_select_missing_reads(dir);
    test_innotin(dir);
    test_byte_view();
    test_split_fields();
    test_alignment_entry();
    test_alignment_entry_short_row();
    test_write_component_gfa(dir);
    test_classify_platon();

    fs::remove_all(dir);

    std::fprintf(stderr, "\n%d checks, %d failure(s)\n", g_checks, g_failures);
    return g_failures == 0 ? 0 : 1;
}

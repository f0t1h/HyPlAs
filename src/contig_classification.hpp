/**
 * @file contig_classification.hpp
 * @brief Contig classification type and prediction-TSV parsing
 */

#ifndef HYPLAS_CONTIG_CLASSIFICATION_HPP
#define HYPLAS_CONTIG_CLASSIFICATION_HPP

#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>
#include <string_view>
#include <charconv>
#include <optional>
#include <vector>
#include <gtl/phmap.hpp>

#include "error.hpp"
#include "mview.hpp"

namespace hyplas {

/**
 * @brief Classification assigned to an assembly contig.
 */
enum class ContigType {
    Plasmid,
    Chromosome,
    Unknown
};

/**
 * @brief Parse a prediction TSV ("ID<TAB>PREDICTION", with header row).
 *
 * A row lacking a prediction column defaults to ContigType::Plasmid, matching
 * the legacy behavior of the original parser.
 *
 * @throws HyplasError if the file cannot be opened.
 */
[[nodiscard]] inline gtl::flat_hash_map<std::string, ContigType>
parse_prediction_tsv(const std::filesystem::path& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        throw HyplasError("cannot open prediction TSV: " + path.string());
    }
    std::string buf{std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>()};

    gtl::flat_hash_map<std::string, ContigType> result;
    mview::byte_view v{std::string_view{buf}};
    v.skip_next('\n'); // skip header row

    while (v.e < v.cap) {
        v.extend_until("\t\n");
        char sep = v.at_end();
        std::string name = static_cast<std::string>(v);

        ContigType type = ContigType::Plasmid;
        if (sep == '\t') {
            v.skip_next('\t');
            v.extend_until("\t\n");
            if (v == std::string_view("plasmid")) {
                type = ContigType::Plasmid;
            } else if (v == std::string_view("chromosome")) {
                type = ContigType::Chromosome;
            } else {
                type = ContigType::Unknown;
            }
        }

        if (!name.empty()) result[std::move(name)] = type;
        v.skip_next('\n');
    }

    return result;
}

/**
 * @brief Collect the names classified as plasmid from a parsed prediction map.
 */
[[nodiscard]] inline gtl::flat_hash_set<std::string>
plasmid_names(const gtl::flat_hash_map<std::string, ContigType>& predictions) {
    gtl::flat_hash_set<std::string> names;
    for (const auto& [name, type] : predictions) {
        if (type == ContigType::Plasmid) names.insert(name);
    }
    return names;
}

/**
 * @brief Thresholds for Platon's RDS-based plasmid/chromosome call.
 */
struct PlatonThresholds {
    double max_chr_rds = -7.9;
    double min_plasmid_rds = 0.7;
};

/**
 * @brief Reproduce Platon's plasmid/chromosome classification from a result.tsv
 *        buffer (header + rows), returning the "ID<TAB>PREDICTION" body for
 *        result_p.tsv.
 *
 * @throws HyplasError if the required RDS column is absent from the header.
 */
[[nodiscard]] inline std::string
classify_platon_tsv(std::string_view content, PlatonThresholds th = {}) {
    std::size_t nl = content.find('\n');
    std::string_view header = (nl == std::string_view::npos) ? content : content.substr(0, nl);
    std::string_view body   = (nl == std::string_view::npos) ? std::string_view{} : content.substr(nl + 1);

    gtl::flat_hash_map<std::string, int> col_idx;
    mview::for_each_field(header, '\t', [&](std::size_t i, std::string_view name) {
        col_idx[std::string(name)] = static_cast<int>(i);
    });
    auto col = [&](const char* n) -> int {
        auto it = col_idx.find(n);
        return it == col_idx.end() ? -1 : it->second;
    };

    int id_col       = col("ID") >= 0 ? col("ID") : 0;
    int rds_col      = col("RDS");
    int circular_col = col("Circular");
    int inc_col      = col("Inc Type(s)");
    int rep_col      = col("# Replication");
    int mob_col      = col("# Mobilization");
    int orit_col     = col("# OriT");
    int hits_col     = col("# Plasmid Hits");
    int rrna_col     = col("# rRNAs");

    if (rds_col < 0) {
        throw HyplasError("RDS column not found in Platon output");
    }

    auto to_double = [](std::string_view s) -> std::optional<double> {
        double v{};
        auto r = std::from_chars(s.data(), s.data() + s.size(), v);
        return r.ec == std::errc{} ? std::optional<double>(v) : std::nullopt;
    };
    auto to_int = [](std::string_view s, int def) -> int {
        int v{};
        auto r = std::from_chars(s.data(), s.data() + s.size(), v);
        return r.ec == std::errc{} ? v : def;
    };

    std::string out = "ID\tPREDICTION\n";
    std::vector<std::string_view> fields;

    mview::for_each_line(body, [&](std::string_view line) {
        if (line.empty()) return;
        fields.clear();
        mview::for_each_field(line, '\t', [&](std::size_t, std::string_view f) {
            fields.push_back(f);
        });

        int nfields = static_cast<int>(fields.size());
        if (nfields <= rds_col) return;
        auto at = [&](int i) -> std::string_view { return fields[static_cast<std::size_t>(i)]; };

        std::optional<double> rds_opt = to_double(at(rds_col));
        if (!rds_opt) return;  // skip rows with an unparseable RDS
        double rds = *rds_opt;

        bool is_chr     = rds <= th.max_chr_rds;
        bool is_plasmid = rds >= th.min_plasmid_rds;

        bool circular = (circular_col >= 0 && circular_col < nfields && at(circular_col) == "yes");
        bool has_inc  = (inc_col >= 0 && inc_col < nfields && !at(inc_col).empty() && at(inc_col) != "0");

        int rep   = (rep_col  >= 0 && rep_col  < nfields) ? to_int(at(rep_col), 0)  : 0;
        int mob   = (mob_col  >= 0 && mob_col  < nfields) ? to_int(at(mob_col), 0)  : 0;
        int orit  = (orit_col >= 0 && orit_col < nfields) ? to_int(at(orit_col), 0) : 0;
        int hits  = (hits_col >= 0 && hits_col < nfields) ? to_int(at(hits_col), 0) : 0;
        int rrnas = (rrna_col >= 0 && rrna_col < nfields) ? to_int(at(rrna_col), 0) : 0;

        bool repmob = (rep + mob) > 0;
        bool hit    = (rds > 0.5 && hits > 0 && rrnas == 0);

        std::string_view id = at(id_col);
        if (is_plasmid) {
            out.append(id); out += "\tplasmid\n";
        } else if (!is_chr && !is_plasmid && (circular || has_inc || repmob || orit > 0 || hit)) {
            out.append(id); out += "\tplasmid\n";
        } else if (is_chr) {
            out.append(id); out += "\tchromosome\n";
        }
    });

    return out;
}

} // namespace hyplas

#endif // HYPLAS_CONTIG_CLASSIFICATION_HPP

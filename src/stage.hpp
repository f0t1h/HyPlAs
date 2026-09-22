/**
 * @file stage.hpp
 * @brief HyPlAs glue over shrn: stage execution, output checks, file helpers.
 *
 * shrn (v0.3) reports failures as plain strings/bools; this layer applies
 * HyPlAs policy: throw HyplasError where callers cannot continue, add
 * FASTQ/FASTA checks and the "[INFO] [stage] Running:" log line.
 *
 * Gzip-aware file utilities (file_is_gzipped, file_first_byte, line_count,
 * concat_files) live here rather than in shrn so shrn stays focused on
 * subprocess execution.
 */

#ifndef HYPLAS_STAGE_HPP
#define HYPLAS_STAGE_HPP

#include "error.hpp"

#include <shrn.hpp>

#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include <zlib.h>

namespace hyplas {

using shrn::Outcome;
using shrn::RunOptions;
using shrn::RunResult;
using shrn::StageTemplate;
using shrn::file_non_empty;
using shrn::file_readable;
using shrn::many;
using shrn::slot;
using shrn::stage;
using shrn::which;

namespace detail {

/// Last-error text: errno when set, else the fallback.
inline std::string errno_string(std::string_view fallback = "I/O error") {
    int e = errno;
    return e != 0 ? std::string(std::strerror(e)) : std::string(fallback);
}

struct gz_guard {
    gzFile gz = nullptr;
    explicit gz_guard(gzFile g) noexcept : gz(g) {}
    gz_guard(const gz_guard&) = delete;
    gz_guard& operator=(const gz_guard&) = delete;
    ~gz_guard() {
        if (gz) gzclose(gz);
    }
    explicit operator bool() const noexcept { return gz != nullptr; }
};

/// Human-readable zlib failure after gzread/gzgetc.
inline std::string gz_error(gzFile gz) {
    int err = Z_OK;
    const char* msg = gzerror(gz, &err);
    if (err == Z_ERRNO) return errno_string();
    return (msg && *msg) ? std::string(msg) : std::string("zlib stream error");
}

inline std::size_t count_newlines(const char* data, std::size_t n) noexcept {
    std::size_t count = 0;
    const char* end = data + n;
    while (const char* hit = static_cast<const char*>(std::memchr(data, '\n', static_cast<std::size_t>(end - data)))) {
        ++count;
        data = hit + 1;
    }
    return count;
}

}  // namespace detail

// ---------------------------------------------------------------------------
// Filesystem helpers (owned by HyPlAs; shrn keeps only ensure_directory)
// ---------------------------------------------------------------------------

/// Whole file contents; throws HyplasError if unreadable. Sized read first,
/// then drain, so /proc-style files reporting size 0 still read correctly.
[[nodiscard]] inline std::string read_file(const std::filesystem::path& p) {
    std::ifstream in(p, std::ios::binary);
    if (!in) throw HyplasError("cannot open file: " + p.string() + ": " + detail::errno_string());

    std::string out;
    std::error_code ec;
    auto st = std::filesystem::status(p, ec);
    if (!ec && std::filesystem::is_regular_file(st)) {
        auto sz = std::filesystem::file_size(p, ec);
        if (!ec) {
            out.resize(static_cast<std::size_t>(sz));
            in.read(out.data(), static_cast<std::streamsize>(sz));
            out.resize(static_cast<std::size_t>(in.gcount()));
            if (in.bad()) throw HyplasError("cannot read file: " + p.string());
        }
    }
    char buffer[1 << 16];
    while (in.read(buffer, sizeof buffer) || in.gcount() > 0) {
        out.append(buffer, static_cast<std::size_t>(in.gcount()));
    }
    if (in.bad()) throw HyplasError("cannot read file: " + p.string());
    return out;
}

inline void ensure_directory(const std::filesystem::path& p) {
    if (!shrn::ensure_directory(p)) {
        throw HyplasError("cannot create directory: " + p.string() + ": " + detail::errno_string());
    }
}

/// Remove a file; a missing file is not an error.
inline void remove_if_exists(const std::filesystem::path& p) {
    std::error_code ec;
    std::filesystem::remove(p, ec);
    if (ec) throw HyplasError("cannot remove: " + p.string() + ": " + ec.message());
}

/// Create a symlink, replacing an existing link/file at `link`.
inline void force_symlink(const std::filesystem::path& target, const std::filesystem::path& link) {
    remove_if_exists(link);
    std::error_code ec;
    std::filesystem::create_symlink(target, link, ec);
    if (ec) throw HyplasError("cannot create symlink: " + link.string() + ": " + ec.message());
}

// ---------------------------------------------------------------------------
// Gzip-aware file utilities
// ---------------------------------------------------------------------------

/// File starts with the gzip magic bytes 1f 8b.
[[nodiscard]] inline bool file_is_gzipped(const std::filesystem::path& p) noexcept {
    std::FILE* f = std::fopen(p.c_str(), "rb");
    if (!f) return false;
    unsigned char magic[2];
    bool gz = std::fread(magic, 1, 2, f) == 2 && magic[0] == 0x1f && magic[1] == 0x8b;
    std::fclose(f);
    return gz;
}

/// First content byte of a plain or gzip file; nullopt when empty or unreadable.
[[nodiscard]] inline std::optional<int> file_first_byte(const std::filesystem::path& p) {
    detail::gz_guard gz(gzopen(p.c_str(), "rb"));
    if (!gz) return std::nullopt;
    int ch = gzgetc(gz.gz);
    return ch == -1 ? std::nullopt : std::optional<int>(ch);
}

/// Lines in a plain or gzipped file; unreadable files count as 0 (with a warning).
[[nodiscard]] inline std::size_t line_count(const std::filesystem::path& p) {
    char buffer[1 << 16];
    std::size_t count = 0;

    if (file_is_gzipped(p)) {
        detail::gz_guard gz(gzopen(p.c_str(), "rb"));
        if (!gz) {
            std::fprintf(stderr, "Warning: Cannot open %s for line count: %s\n",
                         p.c_str(), detail::errno_string().c_str());
            return 0;
        }
        int n;
        while ((n = gzread(gz.gz, buffer, sizeof buffer)) > 0) {
            count += detail::count_newlines(buffer, static_cast<std::size_t>(n));
        }
        if (n < 0) {
            std::fprintf(stderr, "Warning: Cannot count lines of %s: %s\n",
                         p.c_str(), detail::gz_error(gz.gz).c_str());
            return 0;
        }
        return count;
    }

    std::ifstream in(p, std::ios::binary);
    if (!in) {
        std::fprintf(stderr, "Warning: Cannot open %s for line count\n", p.c_str());
        return 0;
    }
    while (in.read(buffer, sizeof buffer) || in.gcount() > 0) {
        count += detail::count_newlines(buffer, static_cast<std::size_t>(in.gcount()));
    }
    if (in.bad()) {
        std::fprintf(stderr, "Warning: Cannot count lines of %s\n", p.c_str());
        return 0;
    }
    return count;
}

enum class concat_mode {
    raw,         ///< byte-for-byte; gzip members concatenate into a valid multi-member stream
    decompress,  ///< decompress gzip inputs and copy plain inputs unchanged
};

/// Concatenate `inputs` into `output` (truncated first). Logs and returns false
/// on failure.
[[nodiscard]] inline bool concat_files(const std::vector<std::filesystem::path>& inputs,
                                       const std::filesystem::path& output,
                                       concat_mode mode = concat_mode::raw) {
    auto fail = [&output](const std::string& why) {
        std::fprintf(stderr, "Error: cannot concatenate into %s: %s\n", output.c_str(), why.c_str());
        return false;
    };

    std::ofstream out(output, std::ios::binary | std::ios::trunc);
    if (!out) return fail(detail::errno_string());

    char buffer[1 << 16];
    for (const auto& input : inputs) {
        if (mode == concat_mode::decompress) {
            detail::gz_guard gz(gzopen(input.c_str(), "rb"));
            if (!gz) return fail("cannot open " + input.string() + ": " + detail::errno_string());
            int n;
            while ((n = gzread(gz.gz, buffer, sizeof buffer)) > 0) {
                out.write(buffer, n);
                if (!out) return fail("write failed");
            }
            if (n < 0) return fail(detail::gz_error(gz.gz));
            continue;
        }
        std::ifstream in(input, std::ios::binary);
        if (!in) return fail("cannot open " + input.string());
        while (in.read(buffer, sizeof buffer) || in.gcount() > 0) {
            out.write(buffer, in.gcount());
            if (!out) return fail("write failed");
        }
        if (in.bad()) return fail("read failed on " + input.string());
    }
    out.close();
    if (!out) return fail("close failed");
    return true;
}

/// Byte-for-byte concatenation (valid for gzip members).
[[nodiscard]] inline bool concat_files_binary(const std::vector<std::filesystem::path>& inputs,
                                              const std::filesystem::path& output) {
    return concat_files(inputs, output);
}

/// First byte is '@' (plain or gzipped). An empty file (0 records) is valid FASTQ.
[[nodiscard]] inline bool file_is_fastq(const std::filesystem::path& p) {
    auto b = file_first_byte(p);
    return b ? *b == '@' : true;
}

/// First byte is '>' (plain or gzipped). An empty file (0 records) is valid FASTA.
[[nodiscard]] inline bool file_is_fasta(const std::filesystem::path& p) {
    auto b = file_first_byte(p);
    return b ? *b == '>' : true;
}

/// Gzip magic bytes and a FASTQ first record (or no records).
[[nodiscard]] inline bool file_is_gzipped_fastq(const std::filesystem::path& p) {
    return file_is_gzipped(p) && file_is_fastq(p);
}

/**
 * Log every spawned command as "[INFO] [stage] Running: cmd". Call once at startup;
 * stages then read as
 *
 *   stage("minigraph alignment")
 *       .expect_which("minigraph")
 *       .expect_file(graph, file_non_empty, "non-empty")
 *       .proc({"minigraph", graph, reads}, opts)
 *       .expect_file(output, file_non_empty, "non-empty")
 *       .or_die_if(!soft_fail);
 */
inline void install_spawn_logger() {
    shrn::default_spawn_hook() = [](std::string_view stage, std::string_view cmd) {
        std::fprintf(stderr, "[INFO] [%.*s] Running: %.*s\n",
                     static_cast<int>(stage.size()), stage.data(),
                     static_cast<int>(cmd.size()), cmd.data());
    };
}

/// Check every tool resolves on PATH; prints per-tool status when verbose.
[[nodiscard]] inline bool validate_tools(const std::vector<std::string>& tools, bool verbose = true) {
    bool all_found = true;
    for (const auto& tool : tools) {
        auto path = which(tool);
        if (!path) all_found = false;
        if (verbose) {
            if (path) std::fprintf(stderr, "  [OK] %s -> %s\n", tool.c_str(), path->c_str());
            else std::fprintf(stderr, "  [MISSING] %s\n", tool.c_str());
        }
    }
    return all_found;
}

} // namespace hyplas

#endif // HYPLAS_STAGE_HPP
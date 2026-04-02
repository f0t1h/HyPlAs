/**
 * @file file_utils.hpp
 * @brief File and path utilities for HyPlAs pipeline
 * 
 * Header-only utilities for common file operations:
 * - Directory creation/copying
 * - Line counting
 * - Temp file management
 * - Gzipped file handling
 */

#ifndef HYPLAS_FILE_UTILS_HPP
#define HYPLAS_FILE_UTILS_HPP

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <random>
#include <string>
#include <vector>

#include <fcntl.h>
#include <unistd.h>
#include <zlib.h>

namespace hyplas {

/**
 * @brief Check if a file exists and is readable
 */
[[nodiscard]] inline bool file_readable(const std::filesystem::path& p) {
    return std::filesystem::exists(p) && 
           std::filesystem::is_regular_file(p) &&
           access(p.c_str(), R_OK) == 0;
}

/**
 * @brief Ensure a directory exists, creating if necessary
 * 
 * @throws std::filesystem::filesystem_error on failure
 */
inline void ensure_directory(const std::filesystem::path& p) {
    if (!std::filesystem::exists(p)) {
        std::filesystem::create_directories(p);
    }
}

/**
 * @brief Count lines in a file (like wc -l)
 * 
 * Works with both regular and gzipped files.
 */
[[nodiscard]] inline size_t line_count(const std::filesystem::path& p) {
    // Check if gzipped
    bool is_gzipped = p.extension() == ".gz";
    
    if (is_gzipped) {
        gzFile gz = gzopen(p.c_str(), "rb");
        if (!gz) {
            std::fprintf(stderr, "Warning: Cannot open %s for line count\n", p.c_str());
            return 0;
        }
        
        size_t count = 0;
        char buffer[8192];
        while (gzgets(gz, buffer, sizeof(buffer)) != nullptr) {
            ++count;
        }
        gzclose(gz);
        return count;
    }
    
    std::ifstream file(p, std::ios::binary);
    if (!file) {
        std::fprintf(stderr, "Warning: Cannot open %s for line count\n", p.c_str());
        return 0;
    }
    
    size_t count = 0;
    char buffer[8192];
    while (file.read(buffer, sizeof(buffer)) || file.gcount() > 0) {
        for (std::streamsize i = 0; i < file.gcount(); ++i) {
            if (buffer[i] == '\n') {
                ++count;
            }
        }
    }
    
    return count;
}

/**
 * @brief Copy a directory recursively (like shutil.copytree)
 * 
 * @param src Source directory
 * @param dst Destination directory
 * @param overwrite If true, overwrite existing files
 */
inline void copy_directory(const std::filesystem::path& src,
                          const std::filesystem::path& dst,
                          bool overwrite = true) {
    auto opts = std::filesystem::copy_options::recursive;
    if (overwrite) {
        opts |= std::filesystem::copy_options::overwrite_existing;
    }
    
    std::filesystem::copy(src, dst, opts);
}

/**
 * @brief Generate a unique temporary file path
 * 
 * @param suffix Optional suffix (e.g., ".fastq")
 * @param prefix Optional prefix
 * @return Path to temp file (file is NOT created)
 */
[[nodiscard]] inline std::filesystem::path make_temp_path(
    const std::string& suffix = "",
    const std::string& prefix = "hyplas_") {
    
    const char* tmpdir = std::getenv("TMPDIR");
    if (!tmpdir) {
        tmpdir = "/tmp";
    }
    
    // Generate random suffix
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_int_distribution<> dist(0, 35);
    
    const char* chars = "0123456789abcdefghijklmnopqrstuvwxyz";
    std::string random_part;
    for (int i = 0; i < 8; ++i) {
        random_part += chars[dist(gen)];
    }
    
    return std::filesystem::path(tmpdir) / (prefix + random_part + suffix);
}

/**
 * @brief Create a temporary file and return its path
 * 
 * @param suffix Optional suffix
 * @return Path to created temp file
 */
[[nodiscard]] inline std::filesystem::path make_temp_file(const std::string& suffix = "") {
    auto path = make_temp_path(suffix);
    
    // Create the file
    std::ofstream file(path);
    if (!file) {
        std::fprintf(stderr, "Error: Cannot create temp file %s\n", path.c_str());
        std::exit(EXIT_FAILURE);
    }
    
    return path;
}

/**
 * @brief RAII wrapper for temporary files
 * 
 * Automatically deletes the file when it goes out of scope.
 */
class TempFile {
public:
    explicit TempFile(const std::string& suffix = "") 
        : path_(make_temp_file(suffix)) {}
    
    ~TempFile() {
        if (std::filesystem::exists(path_)) {
            std::filesystem::remove(path_);
        }
    }
    
    // Non-copyable
    TempFile(const TempFile&) = delete;
    TempFile& operator=(const TempFile&) = delete;
    
    // Movable
    TempFile(TempFile&& other) noexcept : path_(std::move(other.path_)) {
        other.path_.clear();
    }
    TempFile& operator=(TempFile&& other) noexcept {
        if (this != &other) {
            if (std::filesystem::exists(path_)) {
                std::filesystem::remove(path_);
            }
            path_ = std::move(other.path_);
            other.path_.clear();
        }
        return *this;
    }
    
    [[nodiscard]] const std::filesystem::path& path() const { return path_; }
    [[nodiscard]] std::string string() const { return path_.string(); }
    
    /// Release ownership (file won't be deleted)
    std::filesystem::path release() {
        auto p = std::move(path_);
        path_.clear();
        return p;
    }
    
private:
    std::filesystem::path path_;
};

/**
 * @brief Concatenate multiple gzipped FASTQ files to a single uncompressed file
 * 
 * Replaces the Python pattern: subprocess.run(["zcat", *files], stdout=outfile)
 * 
 * @param inputs List of gzipped input files
 * @param output Output file path (will be uncompressed)
 * @return true on success
 */
inline bool concat_gzipped_to_file(const std::vector<std::filesystem::path>& inputs,
                                   const std::filesystem::path& output) {
    std::ofstream out(output, std::ios::binary);
    if (!out) {
        std::fprintf(stderr, "Error: Cannot open %s for writing\n", output.c_str());
        return false;
    }
    
    char buffer[65536];
    
    for (const auto& input : inputs) {
        gzFile gz = gzopen(input.c_str(), "rb");
        if (!gz) {
            std::fprintf(stderr, "Error: Cannot open %s for reading\n", input.c_str());
            return false;
        }
        
        int bytes_read;
        while ((bytes_read = gzread(gz, buffer, sizeof(buffer))) > 0) {
            out.write(buffer, bytes_read);
            if (!out) {
                std::fprintf(stderr, "Error: Write failed to %s\n", output.c_str());
                gzclose(gz);
                return false;
            }
        }
        
        if (bytes_read < 0) {
            int err;
            std::fprintf(stderr, "Error: gzread failed: %s\n", gzerror(gz, &err));
            gzclose(gz);
            return false;
        }
        
        gzclose(gz);
    }
    
    return true;
}

/**
 * @brief Concatenate gzipped files to a new temp file
 * 
 * @param inputs List of gzipped input files
 * @return Path to temp file containing concatenated data, or empty path on error
 */
[[nodiscard]] inline std::filesystem::path concat_gzipped(
    const std::vector<std::filesystem::path>& inputs) {
    
    auto temp = make_temp_path(".fastq");
    if (concat_gzipped_to_file(inputs, temp)) {
        return temp;
    }
    return {};
}

/**
 * @brief Safe file removal (ignores if doesn't exist)
 */
inline void remove_if_exists(const std::filesystem::path& p) {
    std::error_code ec;
    std::filesystem::remove(p, ec);
    // Ignore errors - we don't care if it didn't exist
}

/**
 * @brief Create a symlink, removing target if it exists
 */
inline void force_symlink(const std::filesystem::path& target,
                         const std::filesystem::path& link) {
    remove_if_exists(link);
    std::filesystem::create_symlink(target, link);
}

} // namespace hyplas

#endif // HYPLAS_FILE_UTILS_HPP

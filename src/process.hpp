/**
 * @file process.hpp
 * @brief POSIX subprocess execution with robust error handling
 * 
 * Provides a clean C++ interface for running external tools with:
 * - Full environment inheritance (conda/virtualenv compatible)
 * - Stderr capture for error reporting
 * - Signal detection (OOM kills, timeouts)
 * - Clear error messages with context
 */

#ifndef HYPLAS_PROCESS_HPP
#define HYPLAS_PROCESS_HPP

#include <filesystem>
#include <optional>
#include <string>
#include <vector>

namespace hyplas {

/**
 * @brief Result of a subprocess execution
 */
struct RunResult {
    int exit_code = -1;          ///< Exit code, -1 if not exited normally
    int signal = 0;              ///< Signal number if killed, 0 otherwise
    std::string stderr_output;   ///< Captured stderr (if requested)
    
    /// @brief Check if the process completed successfully
    [[nodiscard]] bool success() const noexcept {
        return exit_code == 0 && signal == 0;
    }
    
    /// @brief Get a human-readable error summary
    [[nodiscard]] std::string error_summary() const;
};

/**
 * @brief Options for subprocess execution
 */
struct RunOptions {
    /// Working directory (nullopt = inherit from parent)
    std::optional<std::filesystem::path> workdir;
    
    /// Capture stderr for error reporting
    bool capture_stderr = true;
    
    /// Let stdout pass through to parent's stdout
    bool inherit_stdout = true;
    
    /// Optional timeout in seconds (0 = no timeout)
    int timeout_seconds = 0;
    
    /// Optional path to redirect stdout to file
    std::optional<std::filesystem::path> stdout_file;
};

/**
 * @brief Execute a subprocess and wait for completion
 * 
 * Uses fork/execvp to run the command, inheriting the full environment.
 * This ensures compatibility with conda, virtualenvs, and module systems.
 * 
 * @param args Command and arguments (args[0] is the program name)
 * @param opts Execution options
 * @return RunResult with exit status and captured stderr
 * 
 * @note The program is found via PATH search (execvp behavior)
 * @note Environment is fully inherited - no manipulation
 */
[[nodiscard]] RunResult run(const std::vector<std::string>& args,
                            const RunOptions& opts = {});

/**
 * @brief Execute a subprocess, abort with error message on failure
 * 
 * Convenience wrapper that calls run() and exits with clear error message
 * if the command fails. Use this for pipeline stages where failure should
 * stop the entire pipeline.
 * 
 * @param args Command and arguments
 * @param stage_context Description of pipeline stage (for error messages)
 * @param opts Execution options
 * 
 * @note Calls std::exit(EXIT_FAILURE) on failure
 */
void run_or_die(const std::vector<std::string>& args,
                const std::string& stage_context,
                const RunOptions& opts = {});

/**
 * @brief Check if a tool exists in PATH
 * 
 * @param name Tool name (not a path)
 * @return true if the tool is found and executable
 */
[[nodiscard]] bool tool_exists(const std::string& name);

/**
 * @brief Get the full path to a tool in PATH
 * 
 * @param name Tool name
 * @return Full path if found, nullopt otherwise
 */
[[nodiscard]] std::optional<std::filesystem::path> which(const std::string& name);

/**
 * @brief Validate that all required tools are available
 * 
 * Checks each tool and prints status. Returns false if any are missing.
 * 
 * @param tools List of tool names to check
 * @param verbose Print status for each tool
 * @return true if all tools are available
 */
[[nodiscard]] bool validate_tools(const std::vector<std::string>& tools,
                                  bool verbose = true);

/**
 * @brief Format a command for display (logging/error messages)
 * 
 * @param args Command and arguments
 * @return Space-separated string with proper quoting
 */
[[nodiscard]] std::string format_command(const std::vector<std::string>& args);

} // namespace hyplas

#endif // HYPLAS_PROCESS_HPP

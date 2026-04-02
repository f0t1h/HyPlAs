/**
 * @file process.cpp
 * @brief POSIX subprocess execution implementation
 */

#include "process.hpp"

#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <signal.h>
#include <sys/wait.h>
#include <unistd.h>

namespace hyplas {

std::string RunResult::error_summary() const {
    if (success()) {
        return "success";
    }
    
    std::string msg;
    if (signal != 0) {
        msg = "killed by signal " + std::to_string(signal);
        // Common signals
        switch (signal) {
            case SIGKILL: msg += " (SIGKILL - possibly OOM killed)"; break;
            case SIGTERM: msg += " (SIGTERM)"; break;
            case SIGSEGV: msg += " (SIGSEGV - segmentation fault)"; break;
            case SIGABRT: msg += " (SIGABRT - aborted)"; break;
            case SIGPIPE: msg += " (SIGPIPE - broken pipe)"; break;
            default: break;
        }
    } else {
        msg = "exited with code " + std::to_string(exit_code);
    }
    
    if (!stderr_output.empty()) {
        // Truncate stderr if too long
        constexpr size_t max_stderr = 1000;
        if (stderr_output.size() > max_stderr) {
            msg += "\nstderr (truncated): " + stderr_output.substr(0, max_stderr) + "...";
        } else {
            msg += "\nstderr: " + stderr_output;
        }
    }
    
    return msg;
}

std::string format_command(const std::vector<std::string>& args) {
    std::string result;
    for (size_t i = 0; i < args.size(); ++i) {
        if (i > 0) result += ' ';
        
        // Check if quoting is needed
        bool needs_quote = false;
        for (char c : args[i]) {
            if (c == ' ' || c == '\t' || c == '"' || c == '\'' || 
                c == '\\' || c == '$' || c == '`') {
                needs_quote = true;
                break;
            }
        }
        
        if (needs_quote) {
            result += '\'';
            for (char c : args[i]) {
                if (c == '\'') {
                    result += "'\\''";
                } else {
                    result += c;
                }
            }
            result += '\'';
        } else {
            result += args[i];
        }
    }
    return result;
}

std::optional<std::filesystem::path> which(const std::string& name) {
    // Don't search if it's already a path
    if (name.find('/') != std::string::npos) {
        if (access(name.c_str(), X_OK) == 0) {
            return std::filesystem::path(name);
        }
        return std::nullopt;
    }
    
    const char* path_env = std::getenv("PATH");
    if (!path_env) {
        return std::nullopt;
    }
    
    std::string path_str(path_env);
    size_t start = 0;
    size_t end;
    
    while ((end = path_str.find(':', start)) != std::string::npos || start < path_str.size()) {
        if (end == std::string::npos) {
            end = path_str.size();
        }
        
        std::string dir = path_str.substr(start, end - start);
        if (dir.empty()) {
            dir = ".";
        }
        
        std::filesystem::path candidate = std::filesystem::path(dir) / name;
        if (access(candidate.c_str(), X_OK) == 0) {
            return candidate;
        }
        
        start = end + 1;
        if (start > path_str.size()) break;
    }
    
    return std::nullopt;
}

bool tool_exists(const std::string& name) {
    return which(name).has_value();
}

bool validate_tools(const std::vector<std::string>& tools, bool verbose) {
    bool all_found = true;
    
    for (const auto& tool : tools) {
        auto path = which(tool);
        if (path) {
            if (verbose) {
                std::fprintf(stderr, "  [OK] %s -> %s\n", tool.c_str(), path->c_str());
            }
        } else {
            all_found = false;
            if (verbose) {
                std::fprintf(stderr, "  [MISSING] %s\n", tool.c_str());
            }
        }
    }
    
    return all_found;
}

RunResult run(const std::vector<std::string>& args, const RunOptions& opts) {
    RunResult result;
    
    if (args.empty()) {
        result.stderr_output = "Empty command";
        return result;
    }
    
    // Create pipe for stderr capture if requested
    int stderr_pipe[2] = {-1, -1};
    if (opts.capture_stderr) {
        if (pipe(stderr_pipe) == -1) {
            result.stderr_output = "Failed to create pipe: " + std::string(strerror(errno));
            return result;
        }
    }
    
    // Create pipe for stdout redirection if file specified
    int stdout_fd = -1;
    if (opts.stdout_file) {
        stdout_fd = open(opts.stdout_file->c_str(), 
                         O_WRONLY | O_CREAT | O_TRUNC, 
                         0644);
        if (stdout_fd == -1) {
            result.stderr_output = "Failed to open stdout file: " + std::string(strerror(errno));
            if (opts.capture_stderr) {
                close(stderr_pipe[0]);
                close(stderr_pipe[1]);
            }
            return result;
        }
    }
    
    pid_t pid = fork();
    
    if (pid == -1) {
        // Fork failed
        result.stderr_output = "Fork failed: " + std::string(strerror(errno));
        if (opts.capture_stderr) {
            close(stderr_pipe[0]);
            close(stderr_pipe[1]);
        }
        if (stdout_fd != -1) {
            close(stdout_fd);
        }
        return result;
    }
    
    if (pid == 0) {
        // Child process
        
        // Change directory if requested
        if (opts.workdir) {
            if (chdir(opts.workdir->c_str()) == -1) {
                std::fprintf(stderr, "Failed to chdir to %s: %s\n", 
                            opts.workdir->c_str(), strerror(errno));
                _exit(127);
            }
        }
        
        // Redirect stderr if capturing
        if (opts.capture_stderr) {
            close(stderr_pipe[0]);  // Close read end
            dup2(stderr_pipe[1], STDERR_FILENO);
            close(stderr_pipe[1]);
        }
        
        // Redirect stdout if file specified
        if (stdout_fd != -1) {
            dup2(stdout_fd, STDOUT_FILENO);
            close(stdout_fd);
        } else if (!opts.inherit_stdout) {
            // Redirect to /dev/null
            int devnull = open("/dev/null", O_WRONLY);
            if (devnull != -1) {
                dup2(devnull, STDOUT_FILENO);
                close(devnull);
            }
        }
        
        // Build argv array
        std::vector<char*> argv;
        argv.reserve(args.size() + 1);
        for (const auto& arg : args) {
            argv.push_back(const_cast<char*>(arg.c_str()));
        }
        argv.push_back(nullptr);
        
        // Execute - uses PATH search
        execvp(argv[0], argv.data());
        
        // If we get here, exec failed
        std::fprintf(stderr, "Failed to execute %s: %s\n", argv[0], strerror(errno));
        _exit(127);
    }
    
    // Parent process
    
    // Close write end of stderr pipe
    if (opts.capture_stderr) {
        close(stderr_pipe[1]);
    }
    
    // Close stdout file if we opened it
    if (stdout_fd != -1) {
        close(stdout_fd);
    }
    
    // Read stderr while waiting
    if (opts.capture_stderr) {
        char buffer[4096];
        ssize_t bytes_read;
        while ((bytes_read = read(stderr_pipe[0], buffer, sizeof(buffer))) > 0) {
            result.stderr_output.append(buffer, bytes_read);
        }
        close(stderr_pipe[0]);
    }
    
    // Wait for child
    int status;
    if (waitpid(pid, &status, 0) == -1) {
        result.stderr_output += "\nwaitpid failed: " + std::string(strerror(errno));
        return result;
    }
    
    if (WIFEXITED(status)) {
        result.exit_code = WEXITSTATUS(status);
        result.signal = 0;
    } else if (WIFSIGNALED(status)) {
        result.exit_code = -1;
        result.signal = WTERMSIG(status);
    }
    
    return result;
}

void run_or_die(const std::vector<std::string>& args,
                const std::string& stage_context,
                const RunOptions& opts) {
    if (args.empty()) {
        std::fprintf(stderr, "[ERROR] [%s] Empty command\n", stage_context.c_str());
        std::exit(EXIT_FAILURE);
    }
    
    // Log the command being run
    std::fprintf(stderr, "[INFO] [%s] Running: %s\n", 
                stage_context.c_str(), 
                format_command(args).c_str());
    
    auto result = run(args, opts);
    
    if (!result.success()) {
        std::fprintf(stderr, "\n[ERROR] [%s] Command failed: %s\n",
                    stage_context.c_str(),
                    format_command(args).c_str());
        std::fprintf(stderr, "[ERROR] %s\n", result.error_summary().c_str());
        std::exit(EXIT_FAILURE);
    }
}

} // namespace hyplas

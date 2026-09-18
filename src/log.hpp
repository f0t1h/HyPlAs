/**
 * @file log.hpp
 * @brief Minimal leveled logging to stderr
 */

#ifndef HYPLAS_LOG_HPP
#define HYPLAS_LOG_HPP

#include <string>

namespace hyplas {

/**
 * @brief Set the minimum level that will be printed.
 * @param level One of "DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL".
 *              Unknown values are treated as "INFO".
 */
void set_log_level(const std::string& level);

/**
 * @brief Emit a log line to stderr if @p level meets the configured threshold.
 */
void log(const std::string& level, const std::string& message);

} // namespace hyplas

#endif // HYPLAS_LOG_HPP

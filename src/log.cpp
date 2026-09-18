/**
 * @file log.cpp
 * @brief Leveled logging implementation
 */

#include "log.hpp"

#include <cstdio>
#include <gtl/phmap.hpp>

namespace hyplas {

namespace {

int level_value(const std::string& level) {
    static const gtl::flat_hash_map<std::string, int> levels = {
        {"DEBUG", 0}, {"INFO", 1}, {"WARNING", 2}, {"ERROR", 3}, {"CRITICAL", 4}
    };
    auto it = levels.find(level);
    return it == levels.end() ? 1 /* INFO */ : it->second;
}

int g_threshold = 1; // INFO

} // namespace

void set_log_level(const std::string& level) {
    g_threshold = level_value(level);
}

void log(const std::string& level, const std::string& message) {
    if (level_value(level) >= g_threshold) {
        std::fprintf(stderr, "[%s] %s\n", level.c_str(), message.c_str());
    }
}

} // namespace hyplas

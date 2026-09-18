/**
 * @file error.hpp
 * @brief Fatal error type for HyPlAs pipeline helpers
 */

#ifndef HYPLAS_ERROR_HPP
#define HYPLAS_ERROR_HPP

#include <stdexcept>

namespace hyplas {

/**
 * @brief Fatal pipeline error.
 *
 * Thrown by pure helper functions (GFA/FASTA/read transforms) instead of
 * calling std::exit() directly. The Pipeline layer decides whether a throw
 * triggers soft-fail recovery or process termination; the main() entrypoint
 * provides a final catch-all that reports and exits non-zero.
 */
struct HyplasError : std::runtime_error {
    using std::runtime_error::runtime_error;
};

} // namespace hyplas

#endif // HYPLAS_ERROR_HPP

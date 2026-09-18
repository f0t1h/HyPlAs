/**
 * @file mview.hpp
 * @brief Zero-copy cursor over a contiguous byte buffer, plus field splitting.
 *
 * `byte_view` is a forward-scanning window into any read-only contiguous byte
 * source (std::string_view, std::span, std::vector<char>, a memory-mapped file,
 * ...). It is used to parse GAF/PAF/TSV/GFA without allocating per line/field.
 */

#ifndef HYPLAS_MVIEW_HPP
#define HYPLAS_MVIEW_HPP

#include <charconv>
#include <concepts>
#include <cstdint>
#include <optional>
#include <ostream>
#include <span>
#include <string>
#include <string_view>

namespace mview {

using u64 = uint64_t;

/// Any read-only contiguous byte source: exposes char* data() and a size().
template<class S>
concept byte_source = requires(const S& s) {
    { s.data() } -> std::convertible_to<const char*>;
    { s.size() } -> std::convertible_to<std::size_t>;
};

/**
 * @brief Forward cursor over a contiguous byte buffer.
 *
 * Invariant: [s, e) is the "current" token; `cap` is the readable extent. The
 * scanning operations (skip_next / extend_until) move `s`/`e` forward only.
 */
struct byte_view {
    /// Offset pair identifying a sub-range, independent of the base pointer.
    struct span_ref {
        u64 s;
        u64 e;
    };

    const char* base = nullptr;
    u64 s = 0;
    u64 e = 0;
    u64 cap = 0;

    byte_view() = default;

    byte_view(const char* b, u64 n) : base(b), s(0), e(0), cap(n) {}
    byte_view(const char* b, u64 cap_, u64 s_, u64 e_) : base(b), s(s_), e(e_), cap(cap_) {}
    byte_view(const char* b, u64 cap_, span_ref g) : base(b), s(g.s), e(g.e), cap(cap_) {}

    /// Construct from any contiguous byte source (string_view, span, mmap, ...).
    template<byte_source S>
    explicit byte_view(const S& src)
        : base(src.data()), s(0), e(0), cap(static_cast<u64>(src.size())) {}

    /// Materialize a sub-view over the same buffer from stored offsets.
    [[nodiscard]] byte_view sub(span_ref g) const { return byte_view{base, g.e, g.s, g.e}; }

    friend std::ostream& operator<<(std::ostream& ost, const byte_view& view) {
        for (u64 i = view.s; i < view.e; ++i) ost << view.base[i];
        return ost;
    }

    operator std::string() const { return std::string{base + s, base + e}; }
    operator span_ref() const { return {s, e}; }

    /// Current [s, e) token as a (non-owning) string_view.
    [[nodiscard]] std::string_view token_view() const {
        return std::string_view{base + s, static_cast<std::size_t>(e - s)};
    }

    /// Parse the current [s, e) token as T via std::from_chars.
    template<class T>
    [[nodiscard]] std::optional<T> to() const {
        T value{};
        auto r = std::from_chars(base + s, base + e, value);
        if (r.ec == std::errc{} && r.ptr == base + e) return value;
        return std::nullopt;
    }

    bool operator==(std::string_view other) const {
        if (other.size() != e - s) return false;
        for (u64 i = s; i < e; ++i) {
            if (base[i] != other[i - s]) return false;
        }
        return true;
    }

    char operator[](std::size_t idx) const { return base[idx + s]; }

    bool operator==(const byte_view& other) const {
        if (other.e - other.s != e - s) return false;
        for (u64 i = s; i < e; ++i) {
            if (base[i] != other.base[i - s + other.s]) return false;
        }
        return true;
    }

    bool operator==(const span_ref& other) const {
        if (other.e - other.s != e - s) return false;
        for (u64 i = s; i < e; ++i) {
            if (base[i] != base[i - s + other.s]) return false;
        }
        return true;
    }

    /// New view starting at the current position: [s, s) with cap = current e.
    [[nodiscard]] byte_view focus() const { return byte_view{base, e, s, s}; }

    /// Scanning cursor bounded to the current line: [s, next '\n'), token reset
    /// to empty. Field navigation on the result clamps at the row boundary, so a
    /// short row can never consume bytes from the following rows.
    [[nodiscard]] byte_view line() const {
        byte_view r{*this};
        r.extend_until('\n');
        return byte_view{base, r.e, s, s};
    }

    /// Advance `s` past the next occurrence of any delimiter; reset e = s.
    void skip_next(const std::string& delimiters) {
        while (s < cap) {
            for (char delimiter : delimiters) {
                if (base[s] == delimiter) {
                    s++;
                    e = s;
                    return;
                }
            }
            ++s;
        }
        s++;
        e = s;
    }

    /// Advance `s` past the next occurrence of `delimiter`; reset e = s.
    void skip_next(char delimiter) {
        while (s < cap) {
            if (base[s] == delimiter) break;
            ++s;
        }
        s++;
        e = s;
    }

    /// Extend `e` to the next occurrence of any delimiter; returns that char.
    char extend_until(const std::string& delimiters) {
        while (e < cap) {
            for (char delimiter : delimiters) {
                if (base[e] == delimiter) return delimiter;
            }
            ++e;
        }
        return 0;
    }

    /// Extend `e` to the next occurrence of `delimiter` (no return value).
    void extend_until(char delimiter) {
        while (e < cap) {
            if (base[e] == delimiter) break;
            ++e;
        }
    }

    /// Delimiter at the current end, or 0 at end-of-buffer.
    [[nodiscard]] char at_end() const { return e == cap ? char{0} : base[e]; }

    void skip_next_n(char delimiter, int n) {
        for (int i = 0; i < n; ++i) skip_next(delimiter);
    }

    void catchup() {
        ++e;
        s = e;
    }
};

/**
 * @brief Split @p s on @p delim, invoking f(index, field) for every field.
 *
 * Zero-copy: fields are std::string_view slices of @p s. An empty input yields
 * a single empty field (index 0), matching successive-delimiter semantics.
 */
template<class F>
constexpr void for_each_field(std::string_view s, char delim, F&& f) {
    std::size_t start = 0;
    std::size_t idx = 0;
    for (;;) {
        std::size_t pos = s.find(delim, start);
        if (pos == std::string_view::npos) {
            f(idx, s.substr(start));
            return;
        }
        f(idx++, s.substr(start, pos - start));
        start = pos + 1;
    }
}

/**
 * @brief Split @p s into at most out.size() fields on @p delim.
 *
 * Fields 0..n-2 are delimited tokens; the final written field holds the
 * remainder (which may itself contain @p delim), mirroring the historical
 * "read N fields, then the rest of the line" behavior. Returns the field count.
 */
inline std::size_t split_fields(std::string_view s, char delim,
                                std::span<std::string_view> out) {
    if (out.empty()) return 0;
    std::size_t n = 0;
    std::size_t start = 0;
    while (n + 1 < out.size()) {
        std::size_t pos = s.find(delim, start);
        if (pos == std::string_view::npos) break;
        out[n++] = s.substr(start, pos - start);
        start = pos + 1;
    }
    out[n++] = s.substr(start);
    return n;
}

/**
 * @brief Invoke f(line) for each '\n'-delimited line of @p buf (newline excluded).
 *
 * Zero-copy: lines are std::string_view slices of @p buf. Matches std::getline
 * semantics — a trailing newline does not yield a final empty line, and a blank
 * line yields an empty view.
 */
template<class F>
void for_each_line(std::string_view buf, F&& f) {
    byte_view v{buf};
    while (v.s < v.cap) {
        v.extend_until('\n');
        f(v.token_view());
        v.skip_next('\n');
    }
}

} // namespace mview

#endif // HYPLAS_MVIEW_HPP


#ifndef _MVIEW_H_
#define _MVIEW_H_

#include <charconv>
#include <concepts>
#include <cstdint>
#include <iterator>
#include <ostream>
#include <string>
#include <string_view>

#include "pstring.hpp"


namespace mview{
using u64 = uint64_t;
template<class MM>
concept mem_mapper_concept = requires(MM m){
  // {MM{std::string{}, u64, u64}} -> std::same_as<MM>;
  {m[int{}]} -> std::convertible_to<char>;
  {m.size()} -> std::convertible_to<size_t>;
  {m.begin()} -> std::random_access_iterator;
  {m.end()} -> std::random_access_iterator;
};

// template<class K,
//   class V, 
//   template<class, class, class, class, class> class MapIMPL=std::unordered_map>
// struct LruCache {

//   template<class ...Args>
//   constexpr auto emplace (const K &k, Args...args){
    
//   }
// };




// template<mem_mapper_concept MM>
// struct lrumm{
//   struct iterator {};
//   const std::string filepath;
//   u64 chunk_size;
//   u64 cache_size;
//   u64 elasticity;
//   LruCache<u64, MM> map_cache;
//   lrumm(const std::string &filepath, u64 chunk_size, u64 cache_size, u64 elasticity=0) : 
//     filepath  { filepath},
//     chunk_size{ chunk_size},
//     cache_size{ cache_size},
//     elasticity{ elasticity},
//     map_cache { cache_size, elasticity}
//     {}
//   constexpr u64 get_chunk(u64 idx) const {
//     return idx / chunk_size;
//   }
//   constexpr u64 get_offset(u64 idx) const {
//     return idx % chunk_size;
//   }
//   constexpr auto operator [] (u64 idx) {
//     u64 chunk = get_chunk(idx);
//     auto it = map_cache.find(chunk);
//     if(it==map_cache.end()){
//       map_cache.emplace(chunk, filepath, chunk*chunk_size, chunk_size);
//     }
//     return map_cache[chunk][get_offset(idx)];
//   }
// };

template<mem_mapper_concept MM>
struct memory_view {
  struct liteview {
    u64 s;
    u64 e;
  };
  const MM *mmap;
  u64 s;
  u64 e;
  u64 cap;

  memory_view(const MM *mmap)
      : mmap{mmap}, s{0}, e{0}, cap{static_cast<u64>(mmap->size())} {}
  memory_view(const MM *mmap, u64 s, u64 e, u64 cap)
      : mmap{mmap}, s{s}, e{e}, cap{cap} {}
  memory_view(const MM *mmap, liteview g)
      : mmap{mmap}, s{g.s}, e{g.e} {}
  friend std::ostream &operator<<(std::ostream &ost, memory_view &view) {
    for (u64 i = view.s; i < view.e; ++i) {
      ost << (char)(*view.mmap)[i];
    }
    return ost;
  }
  operator int() const {
    int value = -1;
    if (std::from_chars(&(*mmap)[s], &(*mmap)[e], value).ec == std::errc{}) {
      return value;
    }
    return -1;
  }
  template <int N> operator pistring::pstring<N>() const {
    return {&(*mmap)[s], &(*mmap)[e]};
  }
  operator std::string() const {
    return {mmap->begin() + s, mmap->begin() + e};
  }
  operator liteview() const { return {s, e}; }

  bool operator==(std::string_view other) const {
    if (other.size() != e - s) {
      return false;
    }
    for (u64 i = s; i < e; ++i) {
      u64 j = i - s;
      if ((*mmap)[i] != other[j]) {
        return false;
      }
    }
    return true;
  }

  char operator[](size_t idx) const { return (*mmap)[idx + s]; }
  bool operator==(const memory_view &other) const {
    if (other.e - other.s != e - s) {
      return false;
    }
    for (u64 i = s; i < e; ++i) {
      u64 j = i - s;
      if ((*mmap)[i] != (*other.mmap)[j + other.s]) {
        return false;
      }
    }
    return true;
  }
  bool operator==(const liteview &other) const {
    if (other.e - other.s != e - s) {
      return false;
    }
    for (u64 i = s; i < e; ++i) {
      u64 j = i - s;
      if ((*mmap)[i] != (*mmap)[j + other.s]) {
        return false;
      }
    }
    return true;
  }
  [[nodiscard]] memory_view focus() const { return memory_view{mmap, s, s, e}; }

  void skip_next(const std::string &delimiters) {
    while (s < cap) {
      for (char delimiter : delimiters) {
        if ((*mmap)[s] == delimiter) {
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
  void skip_next(char delimiter) {
    while (s < cap) {
      if ((*mmap)[s] == delimiter) {
        break;
      }
      ++s;
    }
    s++;
    e = s;
  }

  char extend_until(const std::string &delimiters) {
    while (e < cap) {
      for (char delimiter : delimiters) {
        if ((*mmap)[e] == delimiter) {
          return delimiter;
        }
      }
      ++e;
    }
    return 0;
  }

  [[nodiscard]] char at_end() const {
    if (e == cap) {
      return 0;
    }
    return (*mmap)[e];
  }
  void extend_until(char delimiter) {
    while (e < cap) {
      if ((*mmap)[e] == delimiter) {
        break;
      }
      ++e;
    }
  }
  void skip_next_n(char delimiter, int N) {
    for (int i = 0; i < N; ++i) {
      this->skip_next(delimiter);
    }
  }

  void catchup() {
    ++e;
    s = e;
  }
};
}

#endif
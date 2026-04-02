#ifndef _PSTRING_H_
#define _PSTRING_H_
#include <cassert>
#include <cstdint>
#include <cstring>
#include <string>
#include <string_view>

namespace pistring{
    namespace detail{
template<class T, size_t count>
consteval //< since C++20
T cpy_n(auto first, auto result)
{
    size_t i = 0;
    if constexpr (count > 0)
    {
        *result = *first;
        ++result;
        for (i = 1; i != count && *first; ++i, ++result)
            *result = *++first;
    }
 
    return i;
}
template<uint8_t N>
struct pstring_cap_type{
    using type = uint32_t;
    // static constexpr uint8_t minus = 3;
};

#define pscts(X,Y)   template<> struct pstring_cap_type<X>{using type = Y;}
        pscts(1,uint8_t);pscts(2,uint8_t);pscts(3,uint8_t);pscts(4,uint8_t);
        pscts(5,uint8_t);pscts(6,uint8_t);pscts(7,uint8_t);pscts(8,uint8_t);
 pscts(9,uint16_t);pscts(10,uint16_t);pscts(11,uint16_t);pscts(12,uint16_t);
pscts(13,uint16_t);pscts(14,uint16_t);pscts(15,uint16_t);pscts(16,uint16_t);
#undef pscts

inline uint32_t MurmurHash3_32(const void *key, int len, uint32_t seed) {
    const uint8_t *data = (const uint8_t *)key;
    const int nblocks = len / 4;

    uint32_t h1 = seed;

    const uint32_t c1 = 0xcc9e2d51;
    const uint32_t c2 = 0x1b873593;

    // Body: Process blocks of 4 bytes at a time
    const uint32_t *blocks = (const uint32_t *)(data + nblocks * 4);

    for (int i = -nblocks; i; i++) {
        uint32_t k1 = blocks[i];

        k1 *= c1;
        k1 = (k1 << 15) | (k1 >> (32 - 15));
        k1 *= c2;

        h1 ^= k1;
        h1 = (h1 << 15) | (h1 >> (32 - 15));
        h1 = h1 * 5 + 0xe6546b64;
    }

    // Tail: Process remaining bytes
    const uint8_t *tail = (const uint8_t *)(data + nblocks * 4);

    uint32_t k1 = 0;

    switch (len & 3) {
    case 3:
        k1 ^= tail[2] << 16;
        break;
    case 2:
        k1 ^= tail[1] << 8;
        break;
    case 1:
        k1 ^= tail[0];
        k1 *= c1;
        k1 = (k1 << 15) | (k1 >> (32 - 15));
        k1 *= c2;
        h1 ^= k1;
    }

    // Finalization: Mix the hash to ensure the last few bits are fully mixed
    h1 ^= len;

    /* fmix32 */
    h1 ^= h1 >> 16;
    h1 *= 0x85ebca6b;
    h1 ^= h1 >> 13;
    h1 *= 0xc2b2ae35;
    h1 ^= h1 >> 16;
    return h1;
}
    }
template<uint32_t N>
class pstring{
    static constexpr int cap_bits = 32 - __builtin_clz (N);
    using size_type = typename detail::pstring_cap_type<cap_bits>::type; // + (cap-1)/8 
    size_type _size;
    char _chars[N];
    uint32_t hash_value;
    public:
    constexpr pstring() : _size{0}, _chars{0} {} 
    constexpr pstring(const char *st) : _chars{0}{
        _size = detail::cpy_n<size_type, N>(st, _chars);
        if(_size < N){
            _chars[_size] = 0;
        }
        hash_value = detail::MurmurHash3_32(_chars, _size, 1453);
        assert(_size < N);
    }
    constexpr pstring(const std::string &st) : _chars{0}{
        _size = detail::cpy_n<size_type, N>(st.c_str(), _chars);
        if(_size < N){
            _chars[_size] = 0;
        }
        hash_value = detail::MurmurHash3_32(_chars, _size, 1453);
    }
    constexpr pstring(const auto *from, const auto *to){
        _size = to-from;
        assert(_size < N);
        for(auto it = from; it < to; ++it){
            _chars[it-from] = *it;
        }
        strncpy(_chars, from, to-from);
        _size = to - from;
        if(_size < N){
            _chars[_size] = 0;
        }
        hash_value = detail::MurmurHash3_32(_chars, _size, 1453);
    }
    constexpr char *data(){
        return static_cast<char *>(_chars);
    }
    constexpr const char *c_str() const{
        return static_cast<const char *>(_chars);
    }
    constexpr size_type size(){
        return _size;
    }
    constexpr operator std::string_view () const{
        return {_chars, _size};
    }

    constexpr bool operator==(const pstring<N> &other) const{
        return _size == other._size && hash_value == other.hash_value && memcmp(_chars, other._chars, _size) == 0;
    }
    struct hash{
        constexpr size_t operator ()(const pstring<N> &st) const{
            return st.hash_value;
        }
    };
    struct eq{
        constexpr size_t operator ()(const pstring<N> &st, const pstring<N> &ot) const{
            return st == ot;
        }
    };
};

}
#endif
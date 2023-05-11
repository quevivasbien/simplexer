#pragma once

#include <concepts>

namespace utils {

template <typename T>
requires std::floating_point<T>
bool approxEqual(T a, T b, T rtol, T atol) {
    return std::abs(a - b) <= (atol + rtol * std::abs(b));
}

template <typename T>
requires std::integral<T>
bool approxEqual(T a, T b, T rtol, T atol) {
    return a == b;
}

} // namespace utils
#pragma once

#include <cstdint>
#include <utility>

using u8 = std::uint8_t;
using u16 = std::uint16_t;
using u32 = std::uint32_t;
using u64 = std::uint64_t;

using i8 = std::int8_t;
using i16 = std::int16_t;
using i32 = std::int32_t;
using i64 = std::int64_t;

using byte = std::uint8_t;

template <typename T, std::size_t N>
constexpr std::size_t countof(T const (&)[N]) noexcept
{
  return N;
}

// Helper method for clamping between two values
template <typename T>
T Clamp(T val, T min, T max)
{
  return (val < min) ? min : ((val > max) ? max : val);
}

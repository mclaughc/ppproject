#pragma once
#include <cstdint>

namespace SampleConversion {
using SampleType = float;

// Helper method for clamping between two values
template<typename T>
T Clamp(T val, T min, T max)
{
  return (val < min) ? min : ((val > max) ? max : val);
}

// Conversion methods.
template<typename T>
T ConvertTo(SampleType sample);

// Signed formats.
template<>
inline int8_t ConvertTo<int8_t>(SampleType sample)
{
  return static_cast<int8_t>(Clamp(sample * 128.0f, -128.0f, 127.0f));
}
template<>
inline int16_t ConvertTo<int16_t>(SampleType sample)
{
  return static_cast<int16_t>(Clamp(sample * 32768.0f, -32768.0f, 32767.0f));
}
template<>
inline int32_t ConvertTo<int32_t>(SampleType sample)
{
  return static_cast<int32_t>(Clamp(static_cast<double>(sample * 2147483647.0), -2147483646.0, 2147483647.0));
}

// Unsigned formats.
template<>
inline uint8_t ConvertTo<uint8_t>(SampleType sample)
{
  return static_cast<uint8_t>(Clamp((sample * 0.5f + 0.5f) * 255.0f, 0.0f, 255.0f));
}
template<>
inline uint16_t ConvertTo<uint16_t>(SampleType sample)
{
  return static_cast<uint16_t>(Clamp((sample * 0.5f + 0.5f) * 65535.0f, 0.0f, 65535.0f));
}
template<>
inline uint32_t ConvertTo<uint32_t>(SampleType sample)
{
  return static_cast<uint32_t>(Clamp(static_cast<double>(sample * 0.5f + 0.5f) * 4294967295.0, 0.0, 4294967295.0));
}

// Floating-point formats.
template<>
inline float ConvertTo<float>(SampleType sample)
{
  return sample;
}
template<>
inline double ConvertTo<double>(SampleType sample)
{
  return static_cast<double>(sample);
}

// Helper methods for converting to recording samples.
template<typename T>
inline SampleType ConvertFrom(T sample);

// Signed formats.
template<>
inline SampleType ConvertFrom<int8_t>(int8_t sample)
{
  return static_cast<float>(sample) / 128.0f;
}
template<>
inline SampleType ConvertFrom<int16_t>(int16_t sample)
{
  return static_cast<float>(sample) / 32768.0f;
}
template<>
inline SampleType ConvertFrom<int32_t>(int32_t sample)
{
  return static_cast<float>(static_cast<double>(sample) / 2147483647.0);
}

// Unsigned formats.
template<>
inline SampleType ConvertFrom<uint8_t>(uint8_t sample)
{
  return ((static_cast<float>(sample) / 255.0f) - 0.5f) * 2.0f;
}
template<>
inline SampleType ConvertFrom<uint16_t>(uint16_t sample)
{
  return ((static_cast<float>(sample) / 65535.0f) - 0.5f) * 2.0f;
}
template<>
inline SampleType ConvertFrom<uint32_t>(uint32_t sample)
{
  return ((static_cast<float>(static_cast<double>(sample) / 4294967295.0)) - 0.5f * 2.0f);
}

// Floating-point formats.
template<>
inline SampleType ConvertFrom<float>(float sample)
{
  return sample;
}
template<>
inline SampleType ConvertFrom<double>(double sample)
{
  return static_cast<float>(sample);
}
} // namespace SampleConversion

#pragma once

#include <vector>

#include "shared/types.h"

// Resizable byte array type alias.
using ByteArray = std::vector<byte>;

// We use floating-point samples to match what the ogg files will store internally.
// Certainly not optimal, but at least when decoding/encoding we won't lose quality.
using RecordingSample = float;

// Sample array alias.
using RecordingSampleArray = std::vector<RecordingSample>;

// Helper methods for converting to various formats.
// These have to be defined as inline, due to the one-definition-rule.
template <typename T>
T ConvertRecordingSampleTo(RecordingSample sample);

// Signed formats.
template <>
inline i8 ConvertRecordingSampleTo<i8>(RecordingSample sample)
{
  return static_cast<i8>(Clamp(sample * 128.0f, -128.0f, 127.0f));
}
template <>
inline i16 ConvertRecordingSampleTo<i16>(RecordingSample sample)
{
  return static_cast<i16>(Clamp(sample * 32768.0f, -32768.0f, 32767.0f));
}
template <>
inline i32 ConvertRecordingSampleTo<i32>(RecordingSample sample)
{
  return static_cast<i32>(Clamp(static_cast<double>(sample * 2147483647.0), -2147483646.0, 2147483647.0));
}

// Unsigned formats.
template <>
inline u8 ConvertRecordingSampleTo<u8>(RecordingSample sample)
{
  return static_cast<u8>(Clamp((sample * 0.5f + 0.5f) * 255.0f, 0.0f, 255.0f));
}
template <>
inline u16 ConvertRecordingSampleTo<u16>(RecordingSample sample)
{
  return static_cast<u16>(Clamp((sample * 0.5f + 0.5f) * 65535.0f, 0.0f, 65535.0f));
}
template <>
inline u32 ConvertRecordingSampleTo<u32>(RecordingSample sample)
{
  return static_cast<u32>(Clamp(static_cast<double>(sample * 0.5f + 0.5f) * 4294967295.0, 0.0, 4294967295.0));
}

// Floating-point formats.
template <>
inline float ConvertRecordingSampleTo<float>(RecordingSample sample)
{
  return sample;
}
template <>
inline double ConvertRecordingSampleTo<double>(RecordingSample sample)
{
  return static_cast<double>(sample);
}

// Helper methods for converting to recording samples.
template <typename T>
RecordingSample MakeRecordingSampleFrom(T sample);

// Signed formats.
template <>
inline RecordingSample MakeRecordingSampleFrom<i8>(i8 sample)
{
  return static_cast<float>(sample) / 128.0f;
}
template <>
inline RecordingSample MakeRecordingSampleFrom<i16>(i16 sample)
{
  return static_cast<float>(sample) / 32768.0f;
}
template <>
inline RecordingSample MakeRecordingSampleFrom<i32>(i32 sample)
{
  return static_cast<float>(static_cast<double>(sample) / 2147483647.0);
}

// Unsigned formats.
template <>
inline RecordingSample MakeRecordingSampleFrom<u8>(u8 sample)
{
  return ((static_cast<float>(sample) / 255.0f) - 0.5f) * 2.0f;
}
template <>
inline RecordingSample MakeRecordingSampleFrom<u16>(u16 sample)
{
  return ((static_cast<float>(sample) / 65535.0f) - 0.5f) * 2.0f;
}
template <>
inline RecordingSample MakeRecordingSampleFrom<u32>(u32 sample)
{
  return ((static_cast<float>(static_cast<double>(sample) / 4294967295.0)) - 0.5f * 2.0f);
}

// Floating-point formats.
template <>
inline RecordingSample MakeRecordingSampleFrom<float>(float sample)
{
  return sample;
}
template <>
inline RecordingSample MakeRecordingSampleFrom<double>(double sample)
{
  return static_cast<float>(sample);
}

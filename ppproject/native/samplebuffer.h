#pragma once
#include "sampleconversion.h"
#include <Python.h>
#include <algorithm>
#include <cassert>
#include <cstdint>
#include <cstring>
#include <exception>
#include <pybind11/pybind11.h>
#include <utility>
#include <vector>

// SampleBuffer represents the storage of one or more samples, across one or
// more channels. Samples are stored in 32-bit single precision floating-point
// format, and can be added or extracted in a range of  formats, including
// signed/unsigned 8/16/32-bit integers, as well as 64-bit double precision
// floating-point. TODO: Optimize this into a circular buffer.
class SampleBuffer
{
public:
  using Sample = SampleConversion::SampleType;

  SampleBuffer(int sample_rate, int channels, int reserve_size = 0) : m_sample_rate(sample_rate), m_channels(channels)
  {
    assert(sample_rate > 0 && channels > 0);
    if (reserve_size > 0)
      m_samples.reserve(reserve_size);
  }

  ~SampleBuffer() {}

  // Accessors.
  int GetSampleRate() const { return m_sample_rate; }
  int GetChannels() const { return m_channels; }

  // Returns the "position", or how many frames have been read from the buffer.
  int GetReadPosition() const { return m_read_position; }

  // Returns the "buffer size", or how many frames have been written to the
  // buffer
  // and can be consumed by a reader.
  int GetWritePosition() const { return m_write_position; }

  // Sets the read pointer, or position to the specified frame.
  void SetReadPosition(int pos)
  {
    if (pos < 0 || pos > m_write_position)
      throw std::runtime_error("rpos is past the write position/buffer size");
    m_read_position = pos;
  }

  // Removes "junk" data at the beginning of the buffer, creating room for
  // new frames to be appended on the end.
  void Shrink()
  {
    if (m_read_position == 0)
      return;

    size_t new_size = m_write_position - m_read_position;
    if (new_size > 0)
    {
      std::memmove(m_samples.data(), m_samples.data() + m_read_position * m_channels,
                   new_size * m_channels * sizeof(Sample));
    }

    m_write_position -= m_read_position;
    m_read_position = 0;
  }

  // Returns the number of frames remaining before the end of the buffer.
  int GetSize() const { return m_write_position - m_read_position; }
  bool IsEmpty() const { return m_write_position == m_read_position; }

  // Clears or sets the "end of stream" flag.
  bool IsEndOfStream() const { return m_end_of_stream; }
  void SetEndOfStream(bool value = true) { m_end_of_stream = value; }
  void ClearEndOfStream() { m_end_of_stream = false; }

  // Removes all frames from the buffer.
  void Clear()
  {
    m_read_position = 0;
    m_write_position = 0;
    Shrink();
  }

  // Swaps with the other sample buffer.
  void Swap(SampleBuffer& rhs)
  {
    std::swap(m_samples, rhs.m_samples);
    std::swap(m_read_position, rhs.m_read_position);
    std::swap(m_write_position, rhs.m_write_position);
    std::swap(m_sample_rate, rhs.m_sample_rate);
    std::swap(m_channels, rhs.m_channels);
    std::swap(m_end_of_stream, rhs.m_end_of_stream);
  }

  // Looks ahead into the stream past the current read position, or "peeks".
  template<typename ReturnType>
  std::vector<ReturnType> Peek(int offset) const
  {
    if (offset < 0 || (m_read_position + offset) >= m_write_position)
      throw new std::runtime_error("rpos + offset is past the buffer size");
    std::vector<ReturnType> ret;
    ret.reserve(m_channels);
    size_t spos = (m_read_position + offset) * m_channels;
    for (int i = 0; i < m_channels; i++)
      ret.push_back(SampleConversion::ConvertTo<ReturnType>(m_samples[spos++]));
    return ret;
  }

  // Pushes a new value into the stream. Resizes the buffer if needed.
  // samples is assumed to be a pointer to channels samples.
  template<typename ValueType>
  void Push(const ValueType* samples)
  {
    if ((m_write_position * m_channels) == m_samples.size())
      m_samples.resize(m_samples.size() + m_channels);
    size_t spos = m_write_position * m_channels;
    for (int i = 0; i < m_channels; i++)
      m_samples[spos++] = SampleConversion::ConvertFrom<ValueType>(samples[i]);
    m_write_position++;
  }

  // Removes and returns a new value from the stream.
  template<typename ValueType>
  std::vector<ValueType> Pop()
  {
    auto buf = Peek<ValueType>(0);
    m_read_position++;
    return buf;
  }

  // Reads frames from the buffer, converting and storing them in the samples
  // pointer. samples should be an array of at least channels*num_frames
  // elements. The number of frames written to samples is returned.
  template<typename ValueType>
  int ReadFrames(ValueType* samples, int num_frames)
  {
    // TODO: This loop can be turned into num_frames * m_channels
    int written_count = 0;
    size_t in_pos = m_read_position * m_channels;
    size_t out_pos = 0;
    while (m_read_position < m_write_position)
    {
      for (int i = 0; i < m_channels; i++)
        samples[out_pos++] = SampleConversion::ConvertTo<ValueType>(m_samples[in_pos++]);
      m_read_position++;
      written_count++;
    }
    if (m_read_position >= (m_write_position / 2))
      Shrink();
    return written_count;
  }

  // Reads frames from the buffer, copying them into samples, without advancing the
  // internal read pointer.
  template<typename ValueType>
  int PeekFrames(ValueType* samples, int num_frames)
  {
    // TODO: This loop can be turned into num_frames * m_channels
    int written_count = 0;
    size_t in_pos = m_read_position * m_channels;
    size_t out_pos = 0;
    while (m_read_position < m_write_position)
    {
      for (int i = 0; i < m_channels; i++)
        samples[out_pos++] = SampleConversion::ConvertTo<ValueType>(m_samples[in_pos++]);
      written_count++;
    }
    return written_count;
  }

  // Writes frames from the buffer, converting from the samples pointer to the
  // internal
  // buffer. samples should be an array of at least channels*num_frames
  // elements.
  template<typename ValueType>
  void WriteFrames(const ValueType* samples, int num_frames)
  {
    size_t new_size = m_write_position + (size_t(num_frames) * m_channels);
    if (new_size > m_samples.size())
      m_samples.resize(new_size);

    // TODO: This loop can be turned into num_frames * m_channels
    size_t in_pos = 0;
    size_t out_pos = m_write_position * m_channels;
    for (int i = 0; i < num_frames; i++)
    {
      for (int i = 0; i < m_channels; i++)
        m_samples[out_pos++] = SampleConversion::ConvertFrom<ValueType>(samples[in_pos++]);
      m_write_position++;
    }
  }

  // Copies frames from one buffer to another.
  // TODO: Offset support.
  int CopyFrames(SampleBuffer* dest, int num_frames)
  {
    if (dest->GetChannels() != m_channels)
      throw std::invalid_argument("mismatching number of samples between buffers");

    int frames_to_copy = std::min(GetSize(), num_frames);
    size_t new_dest_size = dest->m_write_position + (size_t(frames_to_copy * m_channels));
    if (new_dest_size > dest->m_samples.size())
      dest->m_samples.resize(new_dest_size);

    size_t start_index = m_read_position * m_channels;
    size_t dest_start_index = dest->m_write_position * m_channels;
    size_t samples_to_copy = size_t(frames_to_copy) * m_channels;
    std::memcpy(&dest->m_samples[dest_start_index], &m_samples[start_index], sizeof(Sample) * samples_to_copy);
    dest->m_write_position += samples_to_copy;
    return frames_to_copy;
  }

  // Removes frames from this buffer, without reading them.
  int RemoveFrames(int num_frames)
  {
    int frames_to_remove = std::min(GetSize(), num_frames);
    m_read_position += frames_to_remove;
    if (m_read_position >= (m_write_position / 2))
      Shrink();
    return frames_to_remove;
  }

  // Returns a pointer to the internally-formatted samples. Offsetting this pointer by zero
  // results the first channel, by one the second, and so on. Valid values are only guaranteed
  // until the last channel in a single frame. Use for high-speed access to the audio data,
  // when format conversion is not necessary.
  const Sample* GetPeekPointer(int offset) const
  {
    if (offset < 0 || (m_read_position + offset) >= m_write_position)
      throw std::runtime_error("rpos + offset is past the buffer size");

    return &m_samples[size_t(offset) * m_channels];
  }

  // Returns a mutable pointer to the internally-formatted samples. Offsetting this pointer by zero
  // results the first channel, by one the second, and so on. Valid values are only guaranteed
  // until the last channel in a single frame. Use for high-speed access to the audio data,
  // when format conversion is not necessary.
  Sample* GetMutablePointer(int offset)
  {
    if (offset < 0 || (m_read_position + offset) >= m_write_position)
      throw std::runtime_error("rpos + offset is past the buffer size");

    return &m_samples[size_t(offset) * m_channels];
  }

  // Inserts n seconds of silence.
  void InsertSilence(float duration) { InsertSilenceFrames(int(duration * m_sample_rate)); }

  // Inserts n frames of silence.
  void InsertSilenceFrames(int num_frames)
  {
    if (num_frames <= 0)
      return;

    size_t new_size = m_write_position + (size_t(num_frames) * m_channels);
    if (new_size > m_samples.size())
      m_samples.resize(new_size);

    size_t out_pos = m_write_position * m_channels;
    for (int i = 0; i < num_frames; i++)
    {
      for (int i = 0; i < m_channels; i++)
        m_samples[out_pos++] = Sample(0);
      m_write_position++;
    }
  }

  // Python-wrapped versions of pop/peek/push.
  template<typename ValueType>
  void PyPush(const pybind11::tuple& samples)
  {
    if (samples.size() != m_channels)
      throw std::invalid_argument("samples contains an incorrect number of channels");

    if ((m_write_position * m_channels) == m_samples.size())
      m_samples.resize(m_samples.size() + m_channels);
    size_t out_pos = m_write_position * m_channels;
    for (int i = 0; i < m_channels; i++)
      m_samples[out_pos++] = SampleConversion::ConvertFrom<ValueType>(samples[i].cast<ValueType>());
    m_write_position++;
  }
  template<typename ValueType>
  pybind11::tuple PyPeek(int offset)
  {
    if (offset < 0 || (m_read_position + offset) >= m_write_position)
      throw std::invalid_argument("rpos + offset is past the buffer size");

    pybind11::tuple res(m_channels);
    size_t spos = (m_read_position + offset) * m_channels;
    for (int i = 0; i < m_channels; i++)
      res[i] = pybind11::cast(SampleConversion::ConvertTo<ValueType>(m_samples[spos++]));
    return res;
  }
  template<typename ValueType>
  pybind11::tuple PyPop()
  {
    auto res = PyPeek<ValueType>(0);
    m_read_position++;
    return res;
  }

private:
  std::vector<Sample> m_samples;
  size_t m_read_position = 0;
  size_t m_write_position = 0;
  int m_sample_rate;
  int m_channels;
  bool m_end_of_stream = false;
};

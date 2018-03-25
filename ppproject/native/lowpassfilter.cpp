#include "samplebuffer.h"
#include "sampleconversion.h"
#include <cmath>
#include <cstdlib>
#include <exception>
#include <pybind11/pybind11.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

class LowPassFilter
{
public:
  LowPassFilter(float cutoff_freq, int num_taps);
  ~LowPassFilter();

  float GetCutoffFrequency() const { return m_cutoff_freq; }
  int GetNumTaps() const { return m_num_taps; }

  void SetCutoffFrequency(float freq);
  void SetNumTaps(int taps);
  void Reset();

  void Run(SampleBuffer* in_buf, SampleBuffer* out_buf);

private:
  void CalculateTaps();

  int m_num_channels = 0;
  int m_last_sample_rate = 44100;
  int m_num_taps = 32;
  float m_cutoff_freq = 22050.0f;
  float m_lambda = 0.0f;

  std::vector<SampleBuffer::Sample> m_taps;
  std::vector<SampleBuffer::Sample> m_outputs;
  std::vector<std::vector<SampleBuffer::Sample>> m_prev;
};

LowPassFilter::LowPassFilter(float cutoff_freq, int num_taps) : m_num_taps(num_taps), m_cutoff_freq(cutoff_freq)
{
  if (SampleConversion::NearEqual(m_cutoff_freq, 0.0f))
    throw std::invalid_argument("cutoff_frequency should be a positive number");
  if (m_num_taps < 1)
    throw std::invalid_argument("taps must be positive");

  m_cutoff_freq = cutoff_freq;
  CalculateTaps();
}

LowPassFilter::~LowPassFilter() {}

void LowPassFilter::SetCutoffFrequency(float freq)
{
  if (SampleConversion::NearEqual(freq, 0.0f))
    throw std::invalid_argument("freq should be a positive number");

  m_cutoff_freq = freq;
  CalculateTaps();
  Reset();
}

void LowPassFilter::SetNumTaps(int taps)
{
  if (taps < 1)
    throw std::invalid_argument("taps must be positive");

  m_num_taps = taps;
  CalculateTaps();
  Reset();
}

void LowPassFilter::Reset()
{
  m_outputs.resize(m_num_channels);
  std::fill(m_outputs.begin(), m_outputs.end(), 0.0f);
  m_prev.resize(m_num_channels);
  for (auto& prev : m_prev)
  {
    prev.resize(m_num_taps);
    std::fill(prev.begin(), prev.end(), 0.0f);
  }
}

void LowPassFilter::CalculateTaps()
{
  m_lambda = float(M_PI * m_cutoff_freq / (m_last_sample_rate / 2.0f));
  m_taps.resize(m_num_taps);
  for (int i = 0; i < m_num_taps; i++)
  {
    float mm = float(i) - ((m_num_taps - 1) / 2.0f);
    if (SampleConversion::NearEqual(mm, 0.0f))
      m_taps[i] = SampleBuffer::Sample(m_lambda / M_PI);
    else
      m_taps[i] = SampleBuffer::Sample(std::sin(mm * m_lambda) / (mm * M_PI));
  }
}

void LowPassFilter::Run(SampleBuffer* in_buf, SampleBuffer* out_buf)
{
  if (in_buf->GetSampleRate() != out_buf->GetSampleRate() || in_buf->GetChannels() != out_buf->GetChannels())
    throw std::runtime_error("input and output buffers should have the same sample rates and channels");

  if (in_buf->GetSampleRate() != m_last_sample_rate || in_buf->GetChannels() != m_num_channels)
  {
    m_last_sample_rate = in_buf->GetSampleRate();
    m_num_channels = in_buf->GetChannels();
    CalculateTaps();
    Reset();
  }

  const int num_frames = in_buf->GetSize();
  for (int frame = 0; frame < num_frames; frame++)
  {
    const SampleBuffer::Sample* in_samples = in_buf->GetPeekPointer(frame);
    for (int channel = 0; channel < m_num_channels; channel++)
    {
      // TODO: This can be optimized, we can peek ahead and remove fewer items.
      auto& prev = m_prev[channel];
      std::memmove(&prev[1], &prev[0], sizeof(SampleBuffer::Sample) * (m_num_taps - 1));
      prev[0] = in_samples[channel];

      SampleBuffer::Sample val = SampleBuffer::Sample(0);
      for (int i = 0; i < m_num_taps; i++)
        val += prev[i] * m_taps[i];
      m_outputs[channel] = val;
    }

    out_buf->Push(m_outputs.data());
  }
  in_buf->RemoveFrames(num_frames);
}

namespace py = pybind11;
PYBIND11_MODULE(lowpassfilter, m)
{
  py::class_<LowPassFilter> lowpassfilter(m, "LowPassFilter");
  lowpassfilter.def(py::init<float, int>(), py::arg("cutoff_freq") = 22050.0f, py::arg("taps") = 32)
    .def("get_cutoff_frequency", &LowPassFilter::GetCutoffFrequency)
    .def("set_cutoff_frequency", &LowPassFilter::SetCutoffFrequency)
    .def("get_num_taps", &LowPassFilter::GetNumTaps)
    .def("set_num_taps", &LowPassFilter::SetNumTaps)
    .def("reset", &LowPassFilter::Reset)
    .def("run", &LowPassFilter::Run);
}

// Implements the "Noise removal from the waveform" algoritm as described in "The calculation of acoustic
// indices derived from long-duration recordings" (2017, Towsey, Michael W.).

#include "samplebuffer.h"
#include "sampleconversion.h"
#include "shared/types.h"
#include <array>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <pybind11/pybind11.h>
#include <tuple>
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

using Sample = double;
using SampleVector = std::vector<Sample>;

// Computes min/max of a sample array.
static std::tuple<Sample, Sample> ComputeMinMax(const SampleVector& samples)
{
  Sample minval = samples[0];
  Sample maxval = minval;
  for (size_t i = 1; i < samples.size(); i++)
  {
    minval = std::min(minval, samples[i]);
    maxval = std::max(maxval, samples[i]);
  }

  return std::tie(minval, maxval);
}

// Computes a 100-bin histogram of sample values.
constexpr int NUM_HISTOGRAM_BINS = 100;
struct Histogram
{
  struct Bin
  {
    Sample bin_min;
    Sample bin_max;
    int num_values;

    Sample GetIntensity() const
    {
      // Take the mid-way point between min/max.
      return bin_min + (bin_max - bin_min) / 2.0;
    }
  };

  Sample all_min;
  Sample all_max;
  Sample range_per_bin;
  std::array<Bin, NUM_HISTOGRAM_BINS> bins;
};

Histogram ComputeHistogram(const SampleVector& samples, Sample min_val, Sample max_val)
{
  Histogram histogram;
  {
    // Distribute the min/max over 100 samples.
    histogram.range_per_bin = (max_val - min_val) / NUM_HISTOGRAM_BINS;
    Sample current_min = min_val;
    for (auto& bin : histogram.bins)
    {
      bin.bin_min = current_min;
      bin.bin_max = current_min + histogram.range_per_bin;
      bin.num_values = 0;
      current_min += histogram.range_per_bin;
    }
  }

  // Fit the samples to the histogram.
  for (const Sample& sample : samples)
  {
    int bin_index = int((sample - histogram.all_min) / histogram.range_per_bin);
    bin_index = Clamp(bin_index, 0, NUM_HISTOGRAM_BINS - 1);
    histogram.bins[bin_index].num_values++;
  }

  return histogram;
}

// Smooths the histogram using a moving average filter.
// TODO: How does this fit into counts?
template<typename T>
static std::vector<T> SmoothHistogram(const std::vector<T>& input_values, int window_size)
{
  std::vector<T> smoothed_histogram;
  smoothed_histogram.reserve(input_values.size());

  const int mid_window = window_size / 2;
  const int value_count = static_cast<int>(input_values.size());
  for (int i = 0; i < NUM_HISTOGRAM_BINS; i++)
  {
    T sum = T(0);
    for (int j = 0; j < window_size; j++)
    {
      // i=3, window_size=5, indices=1,2,3,4,5
      int index = i - mid_window + j;
      sum += input_values[Clamp(index, 0, value_count - 1)];
    }
    smoothed_histogram.push_back(sum / T(window_size));
  }

  return smoothed_histogram;
}

static void SmoothHistogram(const Histogram& histogram, int window_size = 5)
{
  // TODO: Work out how this works with integer counts.
}

// Finds the mode sample intensity.
static Sample ComputeModeIntensity(const Histogram& histogram)
{
  int pos = 0;
  int count = 0;
  for (int i = 0; i < NUM_HISTOGRAM_BINS; i++)
  {
    const Histogram::Bin& bin = histogram.bins[i];
    if (bin.num_values > count)
    {
      pos = i;
      count = bin.num_values;
    }
  }

  return histogram.bins[pos].GetIntensity();
}

static double ComputeStandardDeviation(const Histogram& histogram)
{
  Sample mode_intensity = ComputeModeIntensity(histogram);

  // Work out 68% of the total counts below the mode.
  int total_below_mode = 0;
  for (int i = 0; i < NUM_HISTOGRAM_BINS; i++)
  {
    const Histogram::Bin& bin = histogram.bins[i];
    if (bin.GetIntensity() < mode_intensity)
      total_below_mode += bin.num_values;
  }

  // TODO: This does not seem right...
  return int(double(total_below_mode) * 0.68);
}

// Computes the background noise intensity.
static Sample ComputeBackgroundNoise(const SampleVector& samples, double N)
{
  // (1) Find minimum and maximum values of waveform.
  Sample min_val, max_val;
  std::tie(min_val, max_val) = ComputeMinMax(samples);
  std::fprintf(stderr, "Min: %f, Max: %f\n", min_val, max_val);

  // (2) Compute 100-bin histogram.
  Histogram histogram = ComputeHistogram(samples, min_val, max_val);

  // (3) Smooth histogram with moving average filter.
  SmoothHistogram(histogram);

  // (4) Calculate mode of histogram.
  Sample mode_intensity = ComputeModeIntensity(histogram);
  std::fprintf(stderr, "Mode Intensity: %f\n", noise);

  // (5) Calculate standard deviation.
  Sample stddev = ComputeStandardDeviation(histogram);
  std::fprintf(stderr, "stddev: %f\n", noise);

  // (6) Calculate background noise.
  Sample noise = mode_intensity + (stddev * N);
  std::fprintf(stderr, "Background noise: %f\n", noise);
  return noise;
}

static Sample ScaleInputRange(SampleBuffer::Sample val)
{
  // Scale from -1..1 to 0..1.
  //const double linear = Clamp(double(val) * 0.5 + 0.5, 0.0, 1.0);
  const double linear = Clamp(std::abs(val), 0.0, 1.0);

  // Convert to logarithmic/db scale
  const double db = std::log10(20.0 * std::max(linear, std::numeric_limits<double>::epsilon()));
  return std::abs(std::max(db, -80.0));
}

// Re-scales from 0..1 to -1..1.
// TODO: This should be decibels.
static Sample ScaleOutputRange(SampleBuffer::Sample orig_sample, Sample sample)
{
  // Convert logarithmic back to linear scale.
  double linear = std::pow(10.0, sample / 20.0);

  // Preserve the original sign.
  return (orig_sample < 0.0) ? -linear : linear;
}

// Performs noise removal on a single channel.
static void ChannelNoiseRemoval(SampleBuffer* buffer, int channel, double N, int start_pos, int num_frames)
{
  if (num_frames == 0)
  {
    if (buffer->IsEmpty())
      throw std::invalid_argument("Input sample buffer must not be empty.");
    num_frames = buffer->GetSize();
  }

  // Pull samples from the buffer, scaling to the range of the noise removal algorithm.
  SampleVector samples;
  samples.reserve(num_frames);
  for (int i = 0; i < num_frames; i++)
    samples.push_back(ScaleInputRange(buffer->GetPeekPointer(i)[channel]));

  // (6) Compute background noise, (7) subtract it from all samples, truncating negative values to zero.
  Sample noise = ComputeBackgroundNoise(samples, N);
  for (Sample& sample : samples)
    sample = std::max(sample - noise, Sample(0.0));

  // Write samples back to the buffer, scaling to the range of the sample buffer.
  for (int i = 0; i < num_frames; i++)
    buffer->GetMutablePointer(i)[channel] = ScaleOutputRange(samples[i], buffer->GetMutablePointer(i)[channel]);
}

// Main entry point from python.
static void WaveformNoiseRemoval(SampleBuffer* input_samples, SampleBuffer* output_samples, double N)
{
  if (input_samples->IsEmpty())
    throw std::invalid_argument("Input sample buffer must not be empty.");

  // Copy the samples from the input to the output buffer.
  // These are merely here as placeholders, we overwrite them as we go.
  // This is necessary for multi-channel streams.
  int start_pos = output_samples->GetWritePosition();
  int num_frames = input_samples->GetSize();
  output_samples->CopyFrames(input_samples, num_frames);

  // Perform noise removal on each channel in sequence.
  for (int i = 0; i < output_samples->GetChannels(); i++)
    ChannelNoiseRemoval(output_samples, i, N, start_pos, num_frames);
}

namespace py = pybind11;
PYBIND11_MODULE(noiseremoval, m)
{
  m.def("waveform_noise_removal", WaveformNoiseRemoval,
        "Removes noise across all channels from the input sample buffer (waveform), and writes to the output sample "
        "buffer.",
        py::arg("input_samples"), py::arg("output_samples"), py::arg("N") = 0.0);
  m.def("channel_noise_removal", ChannelNoiseRemoval,
        "Removes noise across a single channel from the input sample buffer (waveform).", py::arg("samples"),
        py::arg("channel"), py::arg("N") = 0.0, py::arg("start_pos") = 0, py::arg("num_frames") = 0);
}

// Implements the "Noise removal from the waveform" algoritm as described in "The calculation of acoustic
// indices derived from long-duration recordings" (2017, Towsey, Michael W.).

#include "samplebuffer.h"
#include "sampleconversion.h"
#include "shared/log.h"
#include "shared/types.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <numeric>
#include <pybind11/pybind11.h>
#include <tuple>
#include <vector>

#define DEBUG
#ifdef DEBUG
#define DBG_LOG(...)                                                                                                   \
  do                                                                                                                   \
  {                                                                                                                    \
    std::fprintf(stderr, __VA_ARGS__);                                                                                 \
  } while (0)
#else
#define DBG_LOG(...)                                                                                                   \
  do                                                                                                                   \
  {                                                                                                                    \
  } while (0)
#endif

// Computes min/max of a sample array.
template<typename T>
static std::tuple<T, T> ComputeMinMax(const std::vector<T>& samples)
{
  T minval = samples[0];
  T maxval = minval;
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
    double bin_min;
    double bin_max;
    int num_values;

    double GetIntensity() const
    {
      // Take the mid-way point between min/max.
      return bin_min + (bin_max - bin_min) / 2.0;
    }
  };

  double all_min;
  double all_max;
  double range_per_bin;
  std::array<Bin, NUM_HISTOGRAM_BINS> bins;
};
struct SmoothedHistogram
{
  double all_min;
  double all_max;
  double range_per_bin;
  std::array<double, NUM_HISTOGRAM_BINS> bins;
};

Histogram ComputeHistogram(const std::vector<double>& samples, double min_val, double max_val)
{
  Histogram histogram;
  {
    // Distribute the min/max over 100 samples.
    histogram.all_min = min_val;
    histogram.all_max = max_val;
    histogram.range_per_bin = (max_val - min_val) / NUM_HISTOGRAM_BINS;
    double current_min = min_val;
    for (auto& bin : histogram.bins)
    {
      bin.bin_min = current_min;
      bin.bin_max = current_min + histogram.range_per_bin;
      bin.num_values = 0;
      current_min += histogram.range_per_bin;
    }
  }

  // Fit the samples to the histogram.
  for (const double& sample : samples)
  {
    int bin_index = int((sample - histogram.all_min) / histogram.range_per_bin);
    bin_index = Clamp(bin_index, 0, NUM_HISTOGRAM_BINS - 1);
    histogram.bins[bin_index].num_values++;
  }

#ifdef DEBUG
  DBG_LOG("histogram: %f %f %f\n", histogram.all_min, histogram.all_max, histogram.range_per_bin);
  for (int i = 0; i < NUM_HISTOGRAM_BINS; i++)
    DBG_LOG("%d ", histogram.bins[i].num_values);
  DBG_LOG("\n");
#endif

  return histogram;
}

// Smooths the histogram using a moving average filter.
// TODO: How does this fit into counts?
static SmoothedHistogram SmoothHistogram(const Histogram& input_histogram, int window_size)
{
  // Convert int array -> double array for counts.
  SmoothedHistogram smoothed_histogram;
  smoothed_histogram.all_min = input_histogram.all_min;
  smoothed_histogram.all_max = input_histogram.all_max;
  smoothed_histogram.range_per_bin = input_histogram.range_per_bin;
  for (int i = 0; i < NUM_HISTOGRAM_BINS; i++)
    smoothed_histogram.bins[i] = double(input_histogram.bins[i].num_values);

  const int mid = window_size / 2;
  if (window_size >= NUM_HISTOGRAM_BINS)
    return smoothed_histogram;

  for (int i = 0; i < mid; i++)
  {
    smoothed_histogram.bins[i] =
      std::accumulate(&smoothed_histogram.bins[0], &smoothed_histogram.bins[i + mid + 1], 0.0) / double(i + mid + 1);
  }

  for (int i = mid; i < (NUM_HISTOGRAM_BINS - mid); i++)
  {
    smoothed_histogram.bins[i] =
      std::accumulate(&smoothed_histogram.bins[i - mid], &smoothed_histogram.bins[i + mid], 0.0) / double(window_size);
  }

  for (int i = NUM_HISTOGRAM_BINS - mid; i < NUM_HISTOGRAM_BINS; i++)
  {
    smoothed_histogram.bins[i] =
      std::accumulate(&smoothed_histogram.bins[i], &smoothed_histogram.bins[NUM_HISTOGRAM_BINS], 0.0) /
      double(NUM_HISTOGRAM_BINS - i);
  }

#ifdef DEBUG
  DBG_LOG("smoothed histogram: \n");
  for (int i = 0; i < NUM_HISTOGRAM_BINS; i++)
    DBG_LOG("%f ", smoothed_histogram.bins[i]);
  DBG_LOG("\n");
#endif

  return smoothed_histogram;
}

// Finds the mode sample intensity.
static std::tuple<double, double> ComputeModeAndStdDev(const SmoothedHistogram& hist)
{
  double max_count = hist.bins[0];
  int index_of_mode = 0;
  for (size_t i = 1; i < hist.bins.size(); i++)
  {
    if (hist.bins[i] > max_count)
    {
      max_count = hist.bins[i];
      index_of_mode = int(i);
    }
  }

  index_of_mode = std::min(index_of_mode, int(double(hist.bins.size()) * 0.95));
  double area_under_curve = std::accumulate(&hist.bins[0], &hist.bins[index_of_mode + 1], 0.0);
  DBG_LOG("index_of_mode=%d\narea_under_curve=%f\n", index_of_mode, area_under_curve);

  int index_of_one_sd = index_of_mode;
  double sum = 0.0;
  double threshold = area_under_curve * 0.68;
  for (int i = index_of_mode; i > 0; i--)
  {
    sum += hist.bins[i];
    index_of_one_sd = i;
    if (sum > threshold)
      break;
  }

  double mode = hist.all_min + (double(index_of_mode + 1) * hist.range_per_bin);
  double stddev = double(std::max(index_of_mode - index_of_one_sd, 1)) * hist.range_per_bin;
  DBG_LOG("index_of_one_sd=%d,mode=%f,stddev=%f\n", index_of_one_sd, mode, stddev);

  return std::tie(mode, stddev);
}

// Computes the background noise intensity.
static double ComputeBackgroundNoise(const std::vector<double>& samples, double N, int window_size)
{
  // (1) Find minimum and maximum values of waveform.
  double min_val, max_val;
  std::tie(min_val, max_val) = ComputeMinMax(samples);
  DBG_LOG("Min: %f, Max: %f\n", min_val, max_val);

  // (2) Compute 100-bin histogram.
  Histogram histogram = ComputeHistogram(samples, min_val, max_val);

  // (3) Smooth histogram with moving average filter.
  SmoothedHistogram smoothed_histogram = SmoothHistogram(histogram, window_size);

  // (4) Calculate mode of histogram.
  double mode, stddev;
  std::tie(mode, stddev) = ComputeModeAndStdDev(smoothed_histogram);
  DBG_LOG("Mode Intensity: %f\n", mode);
  DBG_LOG("stddev: %f\n", stddev);

  // (6) Calculate background noise.
  double noise = mode + (stddev * N);
  DBG_LOG("Background noise: %f\n", noise);
  return noise;
}

static double ConvertToDecibels(SampleBuffer::Sample val)
{
  const double linear = Clamp(std::abs(double(val)), std::numeric_limits<double>::epsilon(), 1.0);

  // Convert to logarithmic/db scale
  const double db = 20.0 * std::log10(linear);
  return std::max(db, -80.0);
}

// Re-scales from 0..1 to -1..1.
// TODO: This should be decibels.
static SampleBuffer::Sample ConvertToLinear(SampleBuffer::Sample orig_sample, double sample)
{
  // Convert logarithmic back to linear scale.
  double linear = std::pow(10.0, sample / 20.0);

  // Since we take the absolute linear value to convert to decibels, preserve the original sign.
  double linear_signed = static_cast<SampleBuffer::Sample>((orig_sample < 0.0) ? -linear : linear);
  // printf("%f -> %f -> %f\n", orig_sample, sample, linear_signed);
  return linear_signed;
}

// Performs noise removal on a single channel.
static void ChannelNoiseRemoval(SampleBuffer* buffer, int channel, double N, int window_size, int start_pos,
                                int num_frames)
{
  if (num_frames == 0)
  {
    if (buffer->IsEmpty())
      throw std::invalid_argument("Input sample buffer must not be empty.");
    num_frames = buffer->GetSize();
  }

  // Pull samples from the buffer, scaling to the range of the noise removal algorithm.
  std::vector<double> samples;
  samples.reserve(num_frames);
  for (int i = 0; i < num_frames; i++)
    samples.push_back(ConvertToDecibels(buffer->GetPeekPointer(i)[channel]));

  // (6) Compute background noise, (7) subtract it from all samples, truncating negative values to zero.
  double noise = ComputeBackgroundNoise(samples, N, window_size);
  for (double& sample : samples)
  {
    // sample = std::max(sample - noise, Sample(0.0));
    sample = std::max(-80.0, sample + noise);
  }

  // Convert from decimals back to linear volume.
  for (int i = 0; i < num_frames; i++)
    buffer->GetMutablePointer(i)[channel] = ConvertToLinear(buffer->GetPeekPointer(i)[channel], samples[i]);
}

// Main entry point from python.
static void WaveformNoiseRemoval(SampleBuffer* input_samples, SampleBuffer* output_samples, double N, int window_size)
{
  if (input_samples->IsEmpty())
    throw std::invalid_argument("Input sample buffer must not be empty.");

  // Copy the samples from the input to the output buffer.
  // These are merely here as placeholders, we overwrite them as we go.
  // This is necessary for multi-channel streams.
  int start_pos = output_samples->GetWritePosition();
  int num_frames = input_samples->GetSize();
  input_samples->CopyFrames(output_samples, num_frames);

  // Perform noise removal on each channel in sequence.
  for (int i = 0; i < output_samples->GetChannels(); i++)
    ChannelNoiseRemoval(output_samples, i, N, window_size, start_pos, num_frames);
}

namespace py = pybind11;
PYBIND11_MODULE(noiseremoval, m)
{
  m.def("waveform_noise_removal", WaveformNoiseRemoval,
        "Removes noise across all channels from the input sample buffer (waveform), and writes to the output sample "
        "buffer.",
        py::arg("input_samples"), py::arg("output_samples"), py::arg("N") = 0.0, py::arg("window_size") = 3);
  m.def("channel_noise_removal", ChannelNoiseRemoval,
        "Removes noise across a single channel from the input sample buffer (waveform).", py::arg("samples"),
        py::arg("channel"), py::arg("N") = 0.0, py::arg("window_size") = 3, py::arg("start_pos") = 0,
        py::arg("num_frames") = 0);
}

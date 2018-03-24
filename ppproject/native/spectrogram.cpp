/*
 * Originally based off spectrogram.c from sndfile-tools
 * https://raw.githubusercontent.com/erikd/sndfile-tools/master/src/spectrogram.c
 * Copyright (C) 2007-2016 Erik de Castro Lopo <erikd@mega-nerd.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 or version 3 of the
 * License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "samplebuffer.h"
#include "shared/string_helpers.h"
#include "shared/types.h"
#include <Python.h>
#include <algorithm>
#include <array>
#include <cairo.h>
#include <cassert>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <fftw3.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <string>
#include <vector>

enum WINDOW_FUNCTION
{
  RECTANGULAR = 0,
  KAISER = 1,
  NUTTALL = 2,
  HANN = 3
};

static void calc_kaiser_window(double* data, int datalen, double beta);
static void calc_nuttall_window(double* data, int datalen);
static void calc_hann_window(double* data, int datalen);

static double besseli0(double x);
static double factorial(int k);

void calc_kaiser_window(double* data, int datalen, double beta)
{
  /*
  **          besseli0 (beta * sqrt (1- (2*x/N).^2))
  ** w (x) =  --------------------------------------,  -N/2 <= x <= N/2
  **                 besseli0 (beta)
  */

  double denom = besseli0(beta);
  if (!std::isfinite(denom))
  {
    printf("besseli0 (%f) : %f\nExiting\n", beta, denom);
    exit(1);
  }

  for (int k = 0; k < datalen; k++)
  {
    double n = k + 0.5 - 0.5 * datalen;
    double two_n_on_N = (2.0 * n) / datalen;
    data[k] = besseli0(beta * std::sqrt(1.0 - two_n_on_N * two_n_on_N)) / denom;
  }
}

void calc_nuttall_window(double* data, int datalen)
{
  // Nuttall window function from http://en.wikipedia.org/wiki/Window_function
  static const std::array<double, 4> a = {{0.355768, 0.487396, 0.144232, 0.012604}};
  for (int k = 0; k < datalen; k++)
  {
    double scale = M_PI * k / (datalen - 1);
    data[k] = a[0] - a[1] * std::cos(2.0 * scale) + a[2] * std::cos(4.0 * scale) - a[3] * std::cos(6.0 * scale);
  }
}

void calc_hann_window(double* data, int datalen)
{
  // Hann window function from http://en.wikipedia.org/wiki/Window_function
  for (int k = 0; k < datalen; k++)
    data[k] = 0.5 * (1.0 - std::cos(2.0 * M_PI * k / (datalen - 1)));
}

/*==============================================================================
 */

double besseli0(double x)
{
  double result = 0.0;
  for (int k = 1; k < 25; k++)
  {
    double temp;

    temp = pow(0.5 * x, k) / factorial(k);
    result += temp * temp;
  }

  return 1.0 + result;
}

double factorial(int val)
{
  static double memory[64] = {1.0};
  static int have_entry = 0;

  assert(val > 0 && val <= (sizeof(memory) / sizeof(memory[0])));
  if (val < have_entry)
    return memory[val];

  for (int k = have_entry + 1; k <= val; k++)
    memory[k] = k * memory[k - 1];

  have_entry = val;
  return memory[val];
}

class FrequencyDomain
{
public:
  FrequencyDomain(int speclen_, WINDOW_FUNCTION wfunc_);
  ~FrequencyDomain();

  int GetLength() const { return speclen; }
  WINDOW_FUNCTION GetWindowFunction() const { return wfunc; }

  int GetInputCount() const { return speclen * 2; }
  double* GetInputData() const { return time_domain; }
  double* GetOutputMagnitudes() const { return mag_spec; }
  double GetMaxMagnitude() const { return max; }

  void Calculate();

private:
  int speclen;
  WINDOW_FUNCTION wfunc;
  fftw_plan plan;

  double* time_domain;
  double* window;
  double* freq_domain;
  double* mag_spec;

  double max;
};

FrequencyDomain::FrequencyDomain(int speclen_, WINDOW_FUNCTION wfunc_) : speclen(speclen_), wfunc(wfunc_)
{
  /* mag_spec has values from [0..speclen] inclusive for 0Hz to Nyquist.
  ** time_domain has an extra element to be able to interpolate between
  ** samples for better time precision, hoping to eliminate artifacts.
  */
  time_domain = new double[2 * speclen + 1];
  window = new double[2 * speclen];
  freq_domain = new double[2 * speclen];
  mag_spec = new double[speclen + 1];

  plan = fftw_plan_r2r_1d(2 * speclen, time_domain, freq_domain, FFTW_R2HC, FFTW_MEASURE | FFTW_PRESERVE_INPUT);
  if (!plan)
    throw std::runtime_error("failed to create fftw plan");

  switch (wfunc)
  {
    case RECTANGULAR:
      break;
    case KAISER:
      calc_kaiser_window(window, 2 * speclen, 20.0);
      break;
    case NUTTALL:
      calc_nuttall_window(window, 2 * speclen);
      break;
    case HANN:
      calc_hann_window(window, 2 * speclen);
      break;
    default:
      throw std::invalid_argument("unknown window function");
  }
}

FrequencyDomain::~FrequencyDomain()
{
  fftw_destroy_plan(plan);
  delete[] time_domain;
  delete[] window;
  delete[] freq_domain;
  delete[] mag_spec;
}

void FrequencyDomain::Calculate()
{
  int freqlen = 2 * speclen;
  if (wfunc != RECTANGULAR)
  {
    for (int k = 0; k < freqlen; k++)
      time_domain[k] *= window[k];
  }

  fftw_execute(plan);

  /* Convert from FFTW's "half complex" format to an array of magnitudes.
  ** In HC format, the values are stored:
  ** r0, r1, r2 ... r(n/2), i(n+1)/2-1 .. i2, i1
  **/
  max = mag_spec[0] = fabs(freq_domain[0]);

  for (int k = 1; k < speclen; k++)
  {
    double re = freq_domain[k];
    double im = freq_domain[freqlen - k];
    mag_spec[k] = sqrt(re * re + im * im);
    max = std::max(max, mag_spec[k]);
  }
  /* Lastly add the point for the Nyquist frequency */
  mag_spec[speclen] = fabs(freq_domain[speclen]);
}

struct RENDER
{
  cairo_surface_t* surface = nullptr;
  int width, height;
  bool log_freq, gray_scale;
  double min_freq, max_freq, fft_freq;
  WINDOW_FUNCTION window_function;
  double spec_floor_db;

  ~RENDER()
  {
    if (surface)
      cairo_surface_destroy(surface);
  }
};

static void get_colour_map_value(float value, double spec_floor_db, unsigned char colour[3], bool gray_scale)
{
  static unsigned char map[][3] = {
    /* These values were originally calculated for a dynamic range of 180dB.
     */
    {255, 255, 255}, /* -0dB */
    {240, 254, 216}, /* -10dB */
    {242, 251, 185}, /* -20dB */
    {253, 245, 143}, /* -30dB */
    {253, 200, 102}, /* -40dB */
    {252, 144, 66},  /* -50dB */
    {252, 75, 32},   /* -60dB */
    {237, 28, 41},   /* -70dB */
    {214, 3, 64},    /* -80dB */
    {183, 3, 101},   /* -90dB */
    {157, 3, 122},   /* -100dB */
    {122, 3, 126},   /* -110dB */
    {80, 2, 110},    /* -120dB */
    {45, 2, 89},     /* -130dB */
    {19, 2, 70},     /* -140dB */
    {1, 3, 53},      /* -150dB */
    {1, 3, 37},      /* -160dB */
    {1, 2, 19},      /* -170dB */
    {0, 0, 0},       /* -180dB */
  };

  float rem;
  int indx;

  if (gray_scale)
  {           /* "value" is a negative value in decibels.
               * black (0,0,0) is for <= -180.0, and the other 255 values
               * should cover the range from -180 to 0 evenly.
               * (value/spec_floor_db) is >=0.0  and <1.0
               * because both value and spec_floor_db are negative.
               * (v/s) * 255.0 goes from 0.0 to 254.9999999 and
               * floor((v/s) * 255) gives us 0 to 254
               * converted to 255 to 1 by subtracting it from 255. */
    int gray; /* The pixel value */

    if (value <= spec_floor_db)
      gray = 0;
    else
    {
      gray = 255 - lrint(floor((value / spec_floor_db) * 255.0));
      assert(gray >= 1 && gray <= 255);
    }
    colour[0] = colour[1] = colour[2] = gray;
    return;
  }

  if (value >= 0.0)
  {
    colour[0] = colour[1] = colour[2] = 255;
    return;
  }

  value = fabs(value * (-180.0 / spec_floor_db) * 0.1);

  indx = lrintf(floor(value));

  if (indx < 0)
  {
    printf("\nError : colour map array index is %d\n\n", indx);
    exit(1);
  }

  if (indx >= (sizeof(map) / sizeof(map[0])))
  {
    colour[0] = colour[1] = colour[2] = 0;
    return;
  }

  rem = fmod(value, 1.0);

  colour[0] = lrintf((1.0 - rem) * map[indx][0] + rem * map[indx + 1][0]);
  colour[1] = lrintf((1.0 - rem) * map[indx][1] + rem * map[indx + 1][1]);
  colour[2] = lrintf((1.0 - rem) * map[indx][2] + rem * map[indx + 1][2]);
}

static void render_spectrogram(cairo_surface_t* surface, double spec_floor_db, float** mag2d, double maxval,
                               double left, double top, double width, double height, bool gray_scale)
{
  unsigned char colour[3] = {0, 0, 0};
  unsigned char* data;
  double linear_spec_floor;
  int w, h, stride;

  stride = cairo_image_surface_get_stride(surface);

  data = cairo_image_surface_get_data(surface);
  memset(data, 0, stride * cairo_image_surface_get_height(surface));

  linear_spec_floor = pow(10.0, spec_floor_db / 20.0);

  for (w = 0; w < width; w++)
    for (h = 0; h < height; h++)
    {
      int x, y;

      mag2d[w][h] = mag2d[w][h] / maxval;
      mag2d[w][h] = (mag2d[w][h] < linear_spec_floor) ? spec_floor_db : 20.0 * log10(mag2d[w][h]);

      get_colour_map_value(mag2d[w][h], spec_floor_db, colour, gray_scale);

      y = height + top - 1 - h;
      x = (w + left) * 4;
      data[y * stride + x + 0] = colour[2];
      data[y * stride + x + 1] = colour[1];
      data[y * stride + x + 2] = colour[0];
      data[y * stride + x + 3] = 0;
    }

  cairo_surface_mark_dirty(surface);
}

/* Helper function:
** Map the index for an output pixel in a column to an index into the
** FFT result representing the same frequency.
** magindex is from 0 to maglen-1, representing min_freq to max_freq Hz.
** Return values from are from 0 to speclen representing frequencies from
** 0 to the Nyquist frequency.
** The result is a floating point number as it may fall between elements,
** allowing the caller to interpolate onto the input array.
*/
static double magindex_to_specindex(int speclen, int maglen, int magindex, double min_freq, double max_freq,
                                    int samplerate, bool log_freq)
{
  double freq; /* The frequency that this output value represents */

  if (!log_freq)
    freq = min_freq + (max_freq - min_freq) * magindex / (maglen - 1);
  else
    freq = min_freq * pow(max_freq / min_freq, (double)magindex / (maglen - 1));

  return (freq * speclen / (samplerate / 2));
}

/* Map values from the spectrogram onto an array of magnitudes, the values
** for display. Reads spec[0..speclen], writes mag[0..maglen-1].
*/
static void interp_spec(float* mag, int maglen, const double* spec, int speclen, const RENDER* render, int samplerate)
{
  int k;

  /* Map each output coordinate to where it depends on in the input array.
  ** If there are more input values than output values, we need to average
  ** a range of inputs.
  ** If there are more output values than input values we do linear
  ** interpolation between the two inputs values that a reverse-mapped
  ** output value's coordinate falls between.
  **
  ** spec points to an array with elements [0..speclen] inclusive
  ** representing frequencies from 0 to samplerate/2 Hz. Map these to the
  ** scale values min_freq to max_freq so that the bottom and top pixels
  ** in the output represent the energy in the sound at min_ and max_freq Hz.
  */

  for (k = 0; k < maglen; k++)
  {
    // Average the pixels in the range it comes from
    double this_val =
      magindex_to_specindex(speclen, maglen, k, render->min_freq, render->max_freq, samplerate, render->log_freq);
    double next =
      magindex_to_specindex(speclen, maglen, k + 1, render->min_freq, render->max_freq, samplerate, render->log_freq);

    // Range check: can happen if --max-freq > samplerate / 2
    if (this_val > speclen)
    {
      mag[k] = 0.0;
      return;
    }

    if (next > this_val + 1)
    {
      /* The output indices are more sparse than the input
       *indices
       ** so average the range of input indices that map to
       *this output,
       ** making sure not to exceed the input array
       *(0..speclen inclusive)
       */
      /* Take a proportional part of the first sample */
      double count = 1.0 - (this_val - floor(this_val));
      double sum = spec[(int)this_val] * count;

      while ((this_val += 1.0) < next && (int)this_val <= speclen)
      {
        sum += spec[(int)this_val];
        count += 1.0;
      }
      /* and part of the last one */
      if ((int)next <= speclen)
      {
        sum += spec[(int)next] * (next - floor(next));
        count += next - floor(next);
      }

      mag[k] = sum / count;
    }
    else
    /* The output indices are more densely packed than the input indices
    ** so interpolate between input values to generate more output values.
    */
    { /* Take a weighted average of the nearest values */
      mag[k] = spec[(int)this_val] * (1.0 - (this_val - floor(this_val))) +
               spec[(int)this_val + 1] * (this_val - floor(this_val));
    }
  }
}

/* Pick the best FFT length good for FFTW?
**
** We use fftw_plan_r2r_1d() for which the documantation
** http://fftw.org/fftw3_doc/Real_002dto_002dReal-Transforms.html says:
**
** "FFTW is generally best at handling sizes of the form
** 2^a 3^b 5^c 7^d 11^e 13^f
** where e+f is either 0 or 1, and the other exponents are arbitrary."
*/

/* Helper function: does N have only 2, 3, 5 and 7 as its factors? */
static bool is_2357(int n)
{
  /* Just eliminate all factors os 2, 3, 5 and 7 and see if 1 remains */
  while (n % 2 == 0)
    n /= 2;
  while (n % 3 == 0)
    n /= 3;
  while (n % 5 == 0)
    n /= 5;
  while (n % 7 == 0)
    n /= 7;
  return (n == 1);
}

/* Helper function: is N a "fast" value for the FFT size? */
static bool is_good_speclen(int n)
{
  /* It wants n, 11*n, 13*n but not (11*13*n)
  ** where n only has as factors 2, 3, 5 and 7
  */
  if (n % (11 * 13) == 0)
    return 0; /* No good */

  return is_2357(n) || ((n % 11 == 0) && is_2357(n / 11)) || ((n % 13 == 0) && is_2357(n / 13));
}

static void render_to_surface(const RENDER* render, const SampleBuffer* inbuf)
{
  const int samplerate = inbuf->GetSampleRate();
  const int filelen = inbuf->GetSize();
  float** mag_spec = NULL; // Indexed by [w][h]

  if (render->width < 1 || render->height < 1)
    throw std::invalid_argument("width and height must be >= 1");

  /*
  ** Choose a speclen value, the spectrum length.
  ** The FFT window size is twice this.
  */
  int speclen;
  if (render->fft_freq != 0.0)
    /* Choose an FFT window size of 1/fft_freq seconds of audio */
    speclen = (samplerate / render->fft_freq + 1) / 2;
  else
    /* Long enough to represent frequencies down to 20Hz. */
    speclen = render->height * (samplerate / 20 / render->height + 1);

  /* Find the nearest fast value for the FFT size. */
  for (int d = 0; /* Will terminate */; d++)
  { /* Logarithmically, the integer above is closer than
     ** the integer below, so prefer it to the one below.
     */
    if (is_good_speclen(speclen + d))
    {
      speclen += d;
      break;
    }
    /* FFT length must also be >= the output height,
    ** otherwise repeated pixel rows occur in the output.
    */
    if (speclen - d >= render->height && is_good_speclen(speclen - d))
    {
      speclen -= d;
      break;
    }
  }

  mag_spec = new float*[render->width];
  for (int w = 0; w < render->width; w++)
    mag_spec[w] = new float[render->height];

  FrequencyDomain spec(speclen, render->window_function);

  double max_mag = 0.0;
  for (int current_x = 0; current_x < render->width; current_x++)
  {
    double* data = spec.GetInputData();
    int datalen = spec.GetInputCount();
    std::memset(data, 0, datalen * sizeof(data[0]));

    // Watch out integer overflow here, as indx * filelen can produce a number greater than 2^31.
    int start = int((u64(current_x) * u64(filelen)) / render->width) - datalen / 2;
    if (start < 0)
    {
      // Fill negative indices with zeros
      data += -start;
      datalen -= -start;
      start = 0;
    }
    if ((start + datalen) > filelen)
      datalen -= (start + datalen) - filelen;
    for (int i = 0; i < datalen; i++)
      data[i] = SampleConversion::ConvertTo<double>(*inbuf->GetPeekPointer(start + i));

    spec.Calculate();
    max_mag = std::max(max_mag, spec.GetMaxMagnitude());

    interp_spec(mag_spec[current_x], render->height, spec.GetOutputMagnitudes(), speclen, render, samplerate);
  }

  render_spectrogram(render->surface, render->spec_floor_db, mag_spec, max_mag, 0, 0, render->width, render->height,
                     render->gray_scale);

  for (int w = 0; w < render->width; w++)
    delete[] mag_spec[w];
  delete[] mag_spec;
}

static RENDER setup_render(SampleBuffer* buf, int width, int height, bool log_freq, bool grayscale, float min_freq,
                           float max_freq, float fft_freq, float dyn_range, const std::string& window_func)
{
  if (width < 1 || height < 1)
    throw std::invalid_argument("width/height must be positive");
  if (min_freq < 0.0f)
    throw std::invalid_argument("min_freq cannot be negative");
  if (fft_freq < 0.0f)
    throw std::invalid_argument("fft_freq cannot be negative");

  RENDER render = {};
  render.width = width;
  render.height = height;
  render.log_freq = log_freq;
  render.gray_scale = grayscale;
  render.min_freq = min_freq;
  render.max_freq = max_freq;
  render.fft_freq = fft_freq;
  render.spec_floor_db = -std::abs(dyn_range);

  if (window_func == "rectangular")
    render.window_function = RECTANGULAR;
  else if (window_func == "kaiser")
    render.window_function = KAISER;
  else if (window_func == "nuttall")
    render.window_function = NUTTALL;
  else if (window_func == "hann")
    render.window_function = HANN;
  else
    throw std::invalid_argument("unknown window function");

  if (render.max_freq == 0.0)
    render.max_freq = static_cast<double>(buf->GetSampleRate()) / 2.0;
  if (render.min_freq == 0.0 && render.log_freq)
    render.min_freq = 20.0;

  /* Do this sanity check here, as soon as max_freq has its default value */
  if (render.min_freq >= render.max_freq)
  {
    throw std::invalid_argument(
      StringFromFormat("min_freq (%g) must be less than max_freq (%g)", render.min_freq, render.max_freq));
  }

  render.surface = cairo_image_surface_create(CAIRO_FORMAT_RGB24, width, height);
  if (!render.surface || cairo_surface_status(render.surface) != CAIRO_STATUS_SUCCESS)
  {
    throw std::runtime_error(StringFromFormat("Error while creating surface : %s",
                                              cairo_status_to_string(cairo_surface_status(render.surface))));
  }

  cairo_surface_flush(render.surface);
  return render;
}

namespace py = pybind11;

static py::array_t<float> Py_render_to_array(SampleBuffer* buf, int width = 640, int height = 480,
                                             bool log_freq = false, bool grayscale = false, float min_freq = 0.0f,
                                             float max_freq = 0.0f, float fft_freq = 0.0f, float dyn_range = 180.0f,
                                             const std::string& window_func = "kaiser")
{
  RENDER render =
    setup_render(buf, width, height, log_freq, grayscale, min_freq, max_freq, fft_freq, dyn_range, window_func);
  render_to_surface(&render, buf);

  const int stride = cairo_image_surface_get_stride(render.surface);
  const u8* data = cairo_image_surface_get_data(render.surface);

  // Allocate the output array. If we did a round-trip to png and back, imread() would give us a WxHx3 float array.
  py::array_t<float> out_array({height, width, 3});
  auto r = out_array.mutable_unchecked<3>();
  for (int y = 0; y < r.shape(0); y++)
  {
    for (int x = 0; x < r.shape(1); x++)
    {
      u32 rgba;
      std::memcpy(&rgba, data + stride * y + x * sizeof(rgba), sizeof(rgba));
      r(y, x, 0) = float((rgba >> 16) & 0xFF) * (1.0f / 255.0f);
      r(y, x, 1) = float((rgba >> 8) & 0xFF) * (1.0f / 255.0f);
      r(y, x, 2) = float(rgba & 0xFF) * (1.0f / 255.0f);
    }
  }

  return out_array;
}

PYBIND11_MODULE(spectrogram, m)
{
  m.def("render_to_array", Py_render_to_array,
        "Render a spectrogram image to a numpy array, from the given audio buffer", py::arg("buf"),
        py::arg("width") = 640, py::arg("height") = 480, py::arg("log_freq") = false, py::arg("grayscale") = false,
        py::arg("min_freq") = 0.0f, py::arg("max_freq") = 0.0f, py::arg("fft_freq") = 0.0f,
        py::arg("dyn_range") = 180.0f, py::arg("window_func") = "kaiser");
}
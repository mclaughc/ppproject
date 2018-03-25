from enum import Enum
from matplotlib import pyplot
from ppproject.native.spectrogram import render_to_array as native_render_to_array

# enumeration of window functions
class WindowFunction(Enum):
  RECTANGULAR = 0
  KAISER = 1
  NUTTALL = 2
  HANN = 3

def render_to_array(buf, width = 640, height = 480, log_freq = False, grayscale = False,
                    min_freq = 0.0, max_freq = 0.0, fft_freq = 0.0, dyn_range = 180.0,
                    window_func = WindowFunction.KAISER):
  """Renders a spectrogram from buf to a numpy array."""
  return native_render_to_array(buf, width = width, height = height, log_freq = log_freq, grayscale = grayscale,
                                min_freq = min_freq, max_freq = max_freq, fft_freq = fft_freq, dyn_range = dyn_range,
                                window_func = window_func.value)

def render_to_file(buf, filename, width = 640, height = 480, log_freq = False, grayscale = False,
                   min_freq = 0.0, max_freq = 0.0, fft_freq = 0.0, dyn_range = 180.0,
                   window_func = WindowFunction.KAISER):
  """Renders a spectrogram from buf to an external png file."""
  spec_array = render_to_array(buf, width = width, height = height, log_freq = log_freq, grayscale = grayscale,
                               min_freq = min_freq, max_freq = max_freq, fft_freq = fft_freq, dyn_range = dyn_range,
                               window_func = window_func)
  pyplot.imsave(filename, spec_array, vmin = 0.0, vmax = 1.0, cmap = pyplot.cm.gray)

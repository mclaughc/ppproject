from ppproject.channels.audiochannel import AudioChannel
from ppproject.native.samplebuffer import SampleBuffer
from ppproject.metadata import Metadata
from ppproject.pipeline import PipelineStage
from ppproject import spectrogram

import logging
logger = logging.getLogger("outputs.spectrogramimage")

class SpectrogramImage(PipelineStage):
  input_channel = None
  output_filename_template = ""
  output_metadata = Metadata()
  next_output_counter = 1
  image_width = 640
  image_height = 480
  border = False
  grayscale = False
  log_freq = False
  min_freq = 0.0
  max_freq = 0.0
  fft_freq = 0.0
  dyn_range = 180.0
  window_func = spectrogram.WindowFunction.KAISER
  ready = False
  done = False

  def __init__(self,
               output_filename_template = "spec-%output-counter%.png",
               image_width = 640,
               image_height = 480,
               border = False,
               grayscale = False,
               log_freq = False,
               min_freq = 0.0,
               max_freq = 0.0,
               fft_freq = 0.0,
               dyn_range = 180.0,
               window_func = spectrogram.WindowFunction.KAISER,
               name = ""):
    super().__init__(name)
    self.output_filename_template = output_filename_template
    self.image_width = image_width
    self.image_height = image_height
    self.border = border
    self.grayscale = grayscale
    self.log_freq = log_freq
    self.min_freq = min_freq
    self.max_freq = max_freq
    self.fft_freq = fft_freq
    self.dyn_range = dyn_range
    self.window_func = window_func

  def set_input_channel(self, name, channel):
    if (name == "input"):
      if (not isinstance(channel, AudioChannel)):
        raise ValueError("channel must be an audio channel")
      self.input_channel = channel
      self.output_buffer = SampleBuffer(channel.get_sample_rate(), channel.get_channels())
    else:
      raise ValueError("unknown input %s" % (name))

  def update_metadata(self):
    self.output_metadata.copy_properties(self.input_channel.get_metadata(), clear_properties=True)
    self.output_metadata.set_property("output-counter", self.next_output_counter)

  def render_spectrogram(self):
    filename = self.output_metadata.format_string(self.output_filename_template)
    logger.debug("rendering spectrogram of %d frames to %s", self.output_buffer.get_size(), filename)
    spectrogram.render_to_file(buf = self.output_buffer, filename = filename, width = self.image_width,
                               height = self.image_height, grayscale = self.grayscale, log_freq = self.log_freq,
                               min_freq = self.min_freq, max_freq = self.max_freq, fft_freq = self.fft_freq,
                               dyn_range = self.dyn_range, window_func = self.window_func)
    
    # spectrogram.render_to_file does not clear the buffer, so clear it here
    # otherwise we end up creating a spectrogram of the entire input
    self.output_buffer.clear()

  ######################################################################################################################
  # Pipeline interface
  ######################################################################################################################
  def make_ready(self):
    if (self.input_channel is None):
      raise RuntimeError("configured without input channel")
    self.ready = True

  def is_ready(self):
    return self.ready

  def is_done(self):
    return self.done

  def run(self):
    # Copy samples into the output buffer.
    frames_to_copy = self.input_channel.get_buffer().get_size()
    if (frames_to_copy > 0):
      self.input_channel.get_buffer().copy_frames(self.output_buffer, frames_to_copy)
      self.input_channel.get_buffer().remove_frames(frames_to_copy)

    # If it's the end of the stream, render the spectrogram and increment the counter.
    if (self.input_channel.end_of_stream()):
      self.update_metadata()
      self.render_spectrogram()
      self.next_output_counter += 1
      self.done = True
    else:
      self.done = False

from ppproject.filter import AudioFilter
from ppproject.channels import AudioChannel
from ppproject.native.lowpassfilter import LowPassFilter as NativeLPF

class LowPassFilter(AudioFilter):
  input_channel = None
  output_channel = None
  native_filter = None
  ready = False

  def __init__(self, cutoff_freq = 1000.0, num_taps = 32, name = ""):
    super().__init__(name)
    self.native_filter = NativeLPF(cutoff_freq, num_taps)

  def get_cutoff_frequency(self):
    return self.native_filter.get_cutoff_frequency()

  def set_cutoff_frequency(self, freq):
    self.native_filter.set_cutoff_frequency(freq)

  def set_input_channel(self, name, channel):
    super().set_input_channel(name, channel)
    if (name == "input"):
      self.input_channel = channel
    else:
      raise ValueError("unknown input %s" % (name))

  def get_output_channel(self, name):
    self.make_ready()
    if (name == "output"):
      return self.output_channel
    else:
      raise ValueError("unknown output %s" % (name))

  def is_ready(self):
    return self.ready

  def make_ready(self):
    if (self.ready):
      return
    
    # Create an output channel with the same format as the input channel.
    self.output_channel = AudioChannel(self.input_channel.get_sample_rate(), self.input_channel.get_channels())
    self.ready = True

  def is_done(self):
    # We can't use is_empty_and_end_of_stream due to overlaps.
    return self.input_channel.end_of_stream()

  def run(self):
    # Run the filter.
    self.native_filter.run(self.input_channel.get_buffer(), self.output_channel.get_buffer())

    # Propogate the end-of-stream through.
    if (self.input_channel.end_of_stream()):
      self.output_channel.set_end_of_stream()
    else:
      self.output_channel.clear_end_of_stream()

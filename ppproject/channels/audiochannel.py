from ppproject.native.samplebuffer import SampleBuffer
from ppproject.metadata import Metadata

class AudioChannel:
  def __init__(self, sample_rate, channels, buffer_size = 0, name = ""):
    if (sample_rate <= 0):
      raise ValueError("sample_rate must be a positive integer")
    if (channels <= 0):
      raise ValueError("channels must be a positive integer")
    if (buffer_size < 0):
      raise ValueError("buffer_size must not be a negative integer")
    self.sample_rate = sample_rate
    self.channels = channels
    self.buffer = SampleBuffer(sample_rate, channels)
    self.metadata = Metadata()
    self.name = name
    self.end_of_stream_flag = False

  def get_sample_rate(self):
    """Returns the sample rate of the audio streamed through this channel."""
    return self.sample_rate

  def get_channels(self):
    """Returns the number of channels for audio streamed through this channel."""
    return self.channels

  def get_buffer(self):
    """Returns the buffer containing the samples streamed through this channel."""
    return self.buffer

  def get_buffer_size(self):
    """Returns the number of frames in the audio buffer."""
    return self.buffer.get_size()

  def is_buffer_empty(self):
    """Returns true if the sample buffer is empty."""
    return self.buffer.is_empty()

  def create_matching_buffer(self):
    """Returns a new, empty sample buffer, with the same format as this channel, suitable for swapping."""
    return SampleBuffer(self.sample_rate, self.channels)

  def detach_buffer(self):
    """Resets the sample buffer for this channel to an empty buffer, and returns the old sample buffer.
    The caller now owns the returned sample buffer."""
    buffer = self.create_matching_buffer()
    buffer.swap(self.buffer)
    return buffer

  def get_metadata(self):
    """Returns the metadata container for audio streamed through this channel."""
    return self.metadata
  
  def get_name(self):
    """Returns the (optional) name of this audio channel."""
    return self.name

  def end_of_stream(self):
    """Returns True if the end-of-stream flag has been set."""
    return self.end_of_stream_flag

  def set_end_of_stream(self):
    """Sets the end-of-stream flag, instructing future streams to complete processing."""
    self.end_of_stream_flag = True

  def clear_end_of_stream(self):
    """Clears the end-of-stream flag, instructing future streams that a new file is being processed."""
    self.end_of_stream_flag = False

  def is_empty_and_end_of_stream(self):
    """Returns true when there are no frames in the buffer, and the end-of-stream flag is set."""
    return (self.end_of_stream_flag and self.buffer.get_size() == 0)

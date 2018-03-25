from ppproject.channels.audiochannel import AudioChannel
from ppproject.native.audiofile import open_file_writer
from ppproject.native.samplebuffer import SampleBuffer
from ppproject.metadata import Metadata
from ppproject.pipeline import PipelineStage

import logging
logger = logging.getLogger("outputs.iterator")

######################################################################################################################
# Value returned by the iterator
######################################################################################################################
class IteratorValue:
  metadata = None
  sample_buffer = None
  def __init__(self, metadata, sample_buffer):
    self.metadata = metadata
    self.sample_buffer = sample_buffer

  def get_metadata(self):
    return self.metadata

  def get_sample_buffer(self):
    return self.sample_buffer

######################################################################################################################
# Main logic of iterator
######################################################################################################################
class IteratorOutput(PipelineStage):
  input_channel = None
  output_filename_template = ""
  next_output_counter = 1
  ready = False
  done = False

  def __init__(self, output_filename_template = "%output-counter%.wav", name = ""):
    super().__init__(name)
    self.output_filename_template = output_filename_template

  def set_input_channel(self, name, channel):
    if (name == "input"):
      if (not isinstance(channel, AudioChannel)):
        raise ValueError("channel must be an audio channel")
      self.input_channel = channel
    else:
      raise ValueError("unknown input %s" % (name))

  def get_iterator(self):
    return IteratorType(self)

  def next_iterator_value(self):
    pipeline_done = False
    while (not self.done):
      # Keep running a single until we hit the end of a stream.
      pipeline_done = self.get_pipeline().run(1)
      if (pipeline_done):
        break

    # If we're done, and have no frames, this is the end.
    if (pipeline_done and self.input_channel.is_buffer_empty()):
      return None

    # Not ideal, but clear the done flag here.
    # This way, the next iteration will actually run.
    self.done = False

    # Construct metadata for the iterator value.
    metadata = Metadata()
    metadata.copy_properties(self.input_channel.get_metadata(), clear_properties=True)
    metadata.set_property("output-counter", self.next_output_counter)
    self.next_output_counter += 1

    # Create a new sample buffer with the same format as the current sample buffer.
    # This is modifying the audio channel, which isn't ideal, but it saves creating a copy.
    # Creating a copy would be bad because this buffer could be rather large.
    buffer = self.input_channel.detach_buffer()
    logger.debug("Returning iterator segment with size of %d", buffer.get_size())
    return IteratorValue(metadata, buffer)

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
    # Hold onto samples until we're ready to return them.
    self.done = self.input_channel.end_of_stream()

######################################################################################################################
# Iterator class called by Python
######################################################################################################################
class IteratorType:
  output = None
  def __init__(self, output):
    self.output = output

  def __iter__(self):
    # We are the iterator class.
    return self

  def __next__(self):
    # The main meat of the iterator.
    val = self.output.next_iterator_value()
    if (val is None):
      raise StopIteration
    else:
      return val

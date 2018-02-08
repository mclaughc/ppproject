from ppproject.channels.audiochannel import AudioChannel
from ppproject.native.audiofile import open_file_writer
from ppproject.metadata import Metadata
from ppproject.pipeline import PipelineStage

import logging
logger = logging.getLogger("outputs.fileoutput")

class FileOutput(PipelineStage):
  input_channel = None
  output_filename_template = ""
  output_metadata = Metadata()
  next_output_counter = 1
  output_file = None
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

  def open_file(self):
    if (self.output_file is not None):
      return
    output_filename = self.output_metadata.format_string(self.output_filename_template)
    logger.info("opening new output file %s", output_filename)
    self.output_file = open_file_writer(filename = output_filename,
                                        sample_rate = self.input_channel.get_sample_rate(),
                                        channels = self.input_channel.get_channels())

  def close_file(self):
    if (self.output_file is None):
      return
    logger.debug("closing output file")
    self.output_file.close()
    self.output_file = None

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
    # Update metadata
    self.output_metadata.copy_properties(self.input_channel.get_metadata(), clear_properties=True)
    self.output_metadata.set_property("output-counter", self.next_output_counter)

    # Open output file if it isn't already
    self.open_file()

    # Push samples to the output file
    buf = self.input_channel.get_buffer()
    num_frames = buf.get_size()
    if (num_frames > 0):
      logger.debug("writing %d frames", num_frames)
      self.output_file.write_frames(buf, num_frames)

    # If it's the end of the stream, close the output file and increment the counter.
    # This way the next write will go to a new file.
    if (self.input_channel.end_of_stream()):
      self.close_file()
      self.next_output_counter += 1
      self.done = True
    else:
      self.done = False

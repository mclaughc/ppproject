from ppproject.native.audiofile import open_file_reader
from ppproject.channels.audiochannel import AudioChannel
from ppproject.metadata import Metadata
from ppproject.pipeline import PipelineStage
import math
import logging

MINIMUM_CHUNK_LENGTH = 1.0
MAXIMUM_CHUNK_LENGTH = 3600.0
logger = logging.getLogger("FileSource")

class FileSource(PipelineStage):
  def __init__(self, filename = "", chunk_size = 10.0, name = ""):
    super().__init__(name)
    if (chunk_size < MINIMUM_CHUNK_LENGTH or chunk_size > MAXIMUM_CHUNK_LENGTH):
      raise ValueError("chunk_size must be between %f and %f" % (MINIMUM_CHUNK_LENGTH, MAXIMUM_CHUNK_LENGTH))

    self.filename = filename
    self.chunk_size = chunk_size
    self.ready = False
    self.reader = None
    self.output_channel = None
    self.channels = 0
    self.sample_rate = 0
    self.chunk_size_frames = 0

  def set_filename(self, filename):
    if (len(filename) == 0):
      raise ValueError("filename must not be an empty string")
    self.filename = filename

  def set_chunk_size(self, chunk_size):
    if (chunk_size < MINIMUM_CHUNK_LENGTH or chunk_size > MAXIMUM_CHUNK_LENGTH):
      raise ValueError("chunk_size must be between %f and %f" % (MINIMUM_CHUNK_LENGTH, MAXIMUM_CHUNK_LENGTH))
    self.chunk_size = chunk_size

  def init(self):
    # Ensure a file was specified.
    if (len(self.filename) == 0):
      raise ValueError("filename must be set before creation")

    # Create the FileReader.
    self.reader = open_file_reader(self.filename)

    # Read properties from the input file.
    self.channels = self.reader.get_channels()
    self.sample_rate = self.reader.get_sample_rate()
    logger.info("%s: opened file %s: %d hz, %d channels, %d frames", self.name, self.filename,
                self.channels, self.sample_rate, self.reader.get_total_frames())

    # Create the channel
    self.chunk_size_frames = int(math.ceil(self.chunk_size * self.sample_rate))
    logger.info("%s: using a chunk size of %d frames (%f seconds)",
                self.name, self.chunk_size_frames, self.chunk_size)
    self.output_channel = AudioChannel(self.sample_rate, self.channels, self.chunk_size_frames,
                                       self.name)
    self.ready = True

  def set_initial_metadata(self):
    pass

  def update_metadata(self):
    pass

  def push_frames(self):
    # Push as many frames as possible to the output buffer, while respecting the chunk size.
    frames_to_push = min(self.reader.get_remaining_frames(),
                         self.chunk_size_frames - self.output_channel.get_buffer().get_size())
    if (frames_to_push > 0):
      logger.debug("push %d frames to output", frames_to_push)
      actual_frames = self.reader.read_frames(self.output_channel.get_buffer(), frames_to_push)
      logger.debug("read %d of %d frames", actual_frames, frames_to_push)

  ######################################################################################################################
  # Pipeline interface
  ######################################################################################################################
  def get_output_channel(self, name):
    if (name == "output"):
      self.make_ready()
      return self.output_channel

    raise ValueError("this node does not have an output named %s" % (name))

  def make_ready(self):
    if (not self.ready):
      self.init()

  def is_ready(self):
    return self.ready

  def is_done(self):
    return self.output_channel.end_of_stream()

  def run(self):
    # Push frames, and update metadata with new positional information, etc.
    if (not self.output_channel.end_of_stream()):
      self.push_frames()
      self.update_metadata()

      # Set end-of-stream flag if there are no remaining samples.
      if (self.reader.get_remaining_frames() == 0):
        logger.debug("setting end-of-stream")
        self.output_channel.set_end_of_stream()

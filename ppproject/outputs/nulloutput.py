from ppproject.channels.audiochannel import AudioChannel
from ppproject.pipeline import PipelineStage

import logging
logger = logging.getLogger("NullOutput")

class NullOutput(PipelineStage):
  is_audio_input = False
  input_channel = None
  ready = False
  done = False

  def __init__(self, name = ""):
    super().__init__(name)

  def set_input_channel(self, name, input_channel):
    if (name == "input"):
      self.input_channel = input_channel
      self.is_audio_input = isinstance(input_channel, AudioChannel)
      logger.debug("input channel is '%s' (%s)", input_channel.get_name(), type(input_channel))
    else:
      raise ValueError("unknown input %s" % (name))

  def pull_audio(self):
    buf = self.input_channel.get_buffer()
    num_frames = buf.get_size()
    logger.debug("dropping %d frames", num_frames)
    buf.remove_frames(num_frames)
    return self.input_channel.end_of_stream()

  ######################################################################################################################
  # Pipeline interface
  ######################################################################################################################
  def get_name(self):
    return self.name

  def make_ready(self):
    if (self.input_channel is None):
      raise RuntimeError("configured without input channel")
    self.ready = True

  def is_ready(self):
    return self.ready

  def is_done(self):
    return self.done

  def run(self):
    if (self.is_audio_input):
      self.done = self.pull_audio()
    else:
      pass

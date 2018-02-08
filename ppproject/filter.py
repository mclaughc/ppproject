from ppproject.channels.audiochannel import AudioChannel
from ppproject.pipeline import PipelineStage

########################################################################################################################
# Filter Base Class
########################################################################################################################
class Filter(PipelineStage):
  def __init__(self, name = ""):
    super().__init__(name)

  def set_input_channel(self, name, channel):
    raise RuntimeError("set_input_channel not implemented by subclass")

  def get_output_channel(self, name, channel):
    raise RuntimeError("get_output_channel not implemented by subclass")

########################################################################################################################
# Audio Filter Base Class
########################################################################################################################
class AudioFilter(Filter):
  def __init__(self, name = ""):
    super().__init__(name)

  def set_input_channel(self, name, channel):
    if (not isinstance(channel, AudioChannel)):
      raise ValueError("input channels for an audio filter must be an AudioChannel")
    else:
      pass

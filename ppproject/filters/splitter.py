from ppproject.filter import AudioFilter
from ppproject.channels import AudioChannel
import math
import logging

logger = logging.getLogger("filters.Splitter")

class Splitter(AudioFilter):
  input_channel = None
  output_channel = None
  split_length = 0.0
  overlap_duration = 0.0
  overlap_frames = 0
  chunk_frames = 0
  chunk_frames_remaining = 0
  next_segment_number = 1
  pad_short_segments = True
  ready = False

  def __init__(self, split_length = 10.0, overlap_duration = 0.0, pad_short_segments = True, name = ""):
    super().__init__(name)
    if (split_length <= 0.0):
      raise ValueError("split_length must be positive")
    if (overlap_duration < 0.0):
      raise ValueError("overlap_duration must not be negative")
    if (overlap_duration >= split_length):
      raise ValueError("overlap duration cannot be greater than the split length")
    self.split_length = split_length
    self.overlap_duration = overlap_duration
    self.pad_short_segments = pad_short_segments

  def set_split_length(self, split_length):
    if (split_length <= 0.0):
      raise ValueError("split_length must be positive")
    if (self.overlap_duration >= split_length):
      raise ValueError("overlap duration cannot be greater than the split length")
    self.split_length = split_length

  def set_overlap_duration(self, overlap_duration):
    if (self.ready):
      raise RuntimeError("cannot change after ready")
    if (overlap_duration < 0.0):
      raise ValueError("overlap_duration must not be negative")
    if (overlap_duration >= self.split_length):
      raise ValueError("overlap duration cannot be greater than the split length")
    self.overlap_duration = overlap_duration

  def set_pad_short_segments(self, enabled):
    self.pad_short_segments = enabled

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
    self.output_channel.get_metadata().set_property("segment-number", self.next_segment_number)
    self.ready = True

    # Determine chunk size in frames
    self.chunk_frames = math.ceil(self.split_length * self.input_channel.get_sample_rate())
    self.chunk_frames_remaining = self.chunk_frames
    logger.debug("%d frames per chunk for %f seconds", self.chunk_frames, self.split_length)

    # Determine how many frames of overlap we need.
    self.overlap_frames = math.ceil(self.overlap_duration * self.input_channel.get_sample_rate())
    logger.debug("%d frames of overlap for %f seconds", self.overlap_frames, self.overlap_duration)

  def is_done(self):
    # We can't use is_empty_and_end_of_stream due to overlaps.
    return self.input_channel.end_of_stream() and self.input_channel.get_buffer().get_size() <= self.overlap_frames

  def run(self):
    input_buf = self.input_channel.get_buffer()
    end_of_input = self.input_channel.end_of_stream()
    output_buf = self.output_channel.get_buffer()

    # We need at least overlap_frames remaining in the buffer, in case we push a chunk.
    if (not end_of_input and input_buf.get_size() < self.overlap_frames):
      return

    # Figure out how many frames to write to the output.
    frames_to_copy = min(self.chunk_frames_remaining, input_buf.get_size())
    frames_copied = input_buf.copy_frames(output_buf, frames_to_copy)
    assert(frames_copied == frames_to_copy)
    self.chunk_frames_remaining -= frames_to_copy
    logger.debug("copied %d frames, %d frames remaining in chunk, %d in buffer",
                 frames_copied, self.chunk_frames_remaining, input_buf.get_size())

    # How many frames do we pop, to preserve the overlap?
    frames_to_remove = min(self.overlap_frames, input_buf.get_size()) if self.overlap_frames > 0 else frames_to_copy
    frames_removed = input_buf.remove_frames(frames_to_remove)
    assert(frames_removed == frames_to_remove)
    logger.debug("removed %d frames", frames_removed)

    # Do we have a complete chunk?
    if (self.chunk_frames_remaining == 0 or end_of_input):
      if (self.pad_short_segments and self.chunk_frames_remaining > 0 and frames_copied > 0):
        logger.debug("inserting %d frames of silence", self.chunk_frames_remaining)
        output_buf.insert_silence_frames(self.chunk_frames_remaining)

      # Set the end-of-stream flag. We'll remove this next time around.
      logger.debug("pushing chunk")
      self.output_channel.set_end_of_stream()
      if (self.chunk_frames_remaining == 0):
        self.chunk_frames_remaining = self.chunk_frames
        self.output_channel.get_metadata().set_property("segment-number", self.next_segment_number)
        self.next_segment_number += 1
    else:
      self.output_channel.clear_end_of_stream()

class segment:
  def __init__(self, filename, num_channels, sample_rate, num_frames):
    self.filename = filename
    self.num_channels = num_channels
    self.sample_rate = sample_rate
    self.num_frames = num_samples

  def get_channels(self):
    return self.num_channels

  def get_sample_rate(self):
    return self.sample_rate

  def get_num_frames(self):
    return self.num_frames

  def get_length(self):
    assert(self.sample_rate > 0)
    return self.num_frames / self.sample_rate

  def get_length_float(self):
    assert(self.sample_rate > 0)
    return float(self.num_frames) / float(self.sample_rate)

from ppproject.native.audiofile import open_file_reader, open_file_writer

def write_buffer_to_file(buffer, filename):
  """Encodes all samples in buffer to the specified filename.
  The samples will not be consumed by this function, so the buffer must be manually
  cleared afterwards, if more data is to be appended."""
  if (buffer.is_empty()):
    raise RuntimeError("buffer contains no samples")

  file = open_file_writer(filename, buffer.get_sample_rate(), buffer.get_channels())
  file.write_frames(buffer, buffer.get_size())
  file.close()

def read_file_to_buffer(filename):
  """Reads all samples in a file to the returned buffer."""
  file = open_file_reader(filename)
  return file.read_all()

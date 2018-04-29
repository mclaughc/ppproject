import ppproject
import logging
from matplotlib import image as mpimage

INPUT_FILE = "tests/wavs/ff.wav"
SEGMENT_LENGTH = 60.0
SEGMENT_OVERLAP_LENGTH = 0.0
SPECTROGRAM_WIDTH = 640
SPECTROGRAM_HEIGHT = 480
SPECTROGRAM_GRAYSCALE = False

# disable debug logging, it's quite noisy.
#logging.basicConfig(format="[%(asctime)s] %(levelname)s:%(name)s: %(message)s", level = logging.DEBUG)
logging.basicConfig(format="[%(asctime)s] %(levelname)s:%(name)s: %(message)s", level = logging.INFO)

# Create a pipeline for this whole process. Everything has to go through a pipeline.
pipeline = ppproject.Pipeline()

# Set up our input audio stream, and add it as the first stage of the pipeline.
source = ppproject.sources.FileSource(INPUT_FILE, name = "input")
pipeline.add_stage(source)

# Split the audio stream into segments, with optional overlap.
# This is added after the source in the pipeline, so the splitter splits the source directly.
# We have to connect the output of the file source to the splitter, so it knows where to pull its
# data from. The pipeline just determines the ordering, not the connections.
splitter = ppproject.filters.Splitter(SEGMENT_LENGTH, SEGMENT_OVERLAP_LENGTH, "splitter")
splitter.set_input_channel("input", source.get_output_channel("output"))
pipeline.add_stage(splitter)

# Construct our iterator output. This enables us to loop through the segments below with high efficiency,
# processing the data directly, instead of writing it out to a file each time.
output = ppproject.outputs.IteratorOutput()
output.set_input_channel("input", splitter.get_output_channel("output"))
pipeline.add_stage(output)

# Loop through the generated segments. This iterator ends when we reach the end of the source file.
for segment in output.get_iterator():
  segid = segment.get_metadata().get_property("output-counter")
  buf = segment.get_sample_buffer()
  print("got a segment %d of %d frames" % (segid, buf.get_size()))

  # Create spectrogram array.
  spec_array = ppproject.spectrogram.render_to_array(buf = buf, width = SPECTROGRAM_WIDTH,
                                                     height = SPECTROGRAM_HEIGHT,
                                                     grayscale = SPECTROGRAM_GRAYSCALE)

  # spec_array is numpy array of 480x680x3. We could feed this to tensorflow, or whatever.
  # For now, let's just write it out to a file.
  filename = "segment-%d-spec.png" % (segid)
  print("Saving spectrogram to %s" % (filename))
  mpimage.imsave(filename, spec_array)

  # Also write the wav out to a file
  filename = "segment-%d-audio.wav" % (segid)
  print("Saving audio to %s" % (filename))
  ppproject.util.bufferio.write_buffer_to_file(buf, filename)

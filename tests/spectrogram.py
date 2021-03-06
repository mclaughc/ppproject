from ppproject import open_file_reader
from ppproject import SampleBuffer, spectrogram
from numpy import array
from matplotlib import pyplot

SAMPLERATE = 44100
CHANNELS = 1

import math
def make_sine_wave(frequency, buffer, count):
  for i in range(0, count):
    val = math.sin(frequency * (2.0 * math.pi) * i / SAMPLERATE)
    buffer.push_float((val,))

#print("Making sine wave of 10KHz for five seconds")
#buf = SampleBuffer(SAMPLERATE, CHANNELS)
#make_sine_wave(10000, buf, SAMPLERATE * 5)
#print("Popping")
#while buf.get_remaining_frames() > 0:
  #print(buf.pop_float())

#reader = open_file_reader("tests/wavs/ibis.wav")
#reader = open_file_reader("tests/wavs/audiocheck.net_sweep_10Hz_100Hz_-3dBFS_3s.wav")
reader = open_file_reader("tests/wavs/audiocheck.net_sweep_500Hz_1000Hz_-3dBFS_5s.wav")
#reader = open_file_reader("tests/wavs/sweep.wav")
#reader = open_file_reader("tests/wavs/long.wav")
#reader = open_file_reader("tests/wavs/comb.wav")

print("Length in frames", reader.get_total_frames())
print("Length in time", reader.get_duration())

print("Reading to buffer...")
buf = reader.read_all()
print("Buffer size", buf.get_size())

print("Rendering spectrogram...")
#spectrogram.render_to_file(buf=buf, filename="spec.png")

spec_array = spectrogram.render_to_array(buf=buf, grayscale = True, width = 640, height = 480, dyn_range = 180.0)
#print(spec_array)

print("Saving file...")
pyplot.imsave("test.png", spec_array, vmin = 0.0, vmax = 1.0, cmap = pyplot.cm.gray)

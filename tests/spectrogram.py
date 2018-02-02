import sys
sys.path.append("..")

from ppproject import SampleBuffer, render_spectrogram_to_file

SAMPLERATE = 44100
CHANNELS = 1

import math
def make_sine_wave(frequency, buffer, count):
  for i in range(0, count):
    val = math.sin(frequency * (2.0 * math.pi) * i / SAMPLERATE)
    buffer.push_float((val,))

print("Making sine wave of 10KHz for five seconds")
buf = SampleBuffer(SAMPLERATE, CHANNELS)
make_sine_wave(10000, buf, SAMPLERATE * 5)
#print("Popping")
#while buf.get_remaining_frames() > 0:
  #print(buf.pop_float())

print("Rendering spectrogram...")
render_spectrogram_to_file(buf=buf, filename="spec.png")


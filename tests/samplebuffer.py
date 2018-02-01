import sys
sys.path.append("..")

from ppproject import SampleBuffer

SAMPLERATE = 44100
CHANNELS = 2

import math
def make_sine_wave(frequency, buffer, count):
  for i in range(0, count):
    val = math.sin(frequency * (2.0 * math.pi) * i / SAMPLERATE)
    buffer.push_float((val, val))

buf = SampleBuffer(SAMPLERATE, CHANNELS)
print("Position", buf.get_position());
print("Size", buf.get_size());

print("Pushing")
buf.push_float((0.0, 0.0))
buf.push_float((0.5, -0.5))
buf.push_float((1.0, -1.0))
buf.push_float((-0.5, 0.5))
buf.push_float((0.0, 0.0))
print("Position", buf.get_position());
print("Size", buf.get_size());
print("Remaining", buf.get_remaining_frames());

print("Popping")
for i in range(0, 5):
  print(buf.pop_s16())

print("Making sine wave of 10KHz for five seconds")
make_sine_wave(10000, buf, SAMPLERATE * 5)
print("Popping")
while buf.get_remaining_frames() > 0:
  print(buf.pop_float())

import sys
sys.path.append("..")

import ppproject

print("Opening file...")
reader = ppproject.make_file_reader("tests/wavs/ibis.wav")
#reader = ppproject.make_file_reader("tests/wavs/audiocheck.net_sweep_10Hz_100Hz_-3dBFS_3s.wav")
#reader = ppproject.make_file_reader("tests/wavs/audiocheck.net_sweep_500Hz_1000Hz_-3dBFS_5s.wav")

print("Length in frames", reader.get_total_frames())
print("Length in time", reader.get_duration())

print("Reading to buffer...")
buf = reader.read_all()
print("Buffer size", buf.get_size())

#print("Popping")
#while buf.get_remaining_frames() > 0:
  #print(buf.pop_float())

print("Rendering spectrogram...")
ppproject.render_spectrogram_to_file(buf=buf, filename="spec.png", window_func="hann")


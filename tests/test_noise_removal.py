from ppproject import open_file_reader, open_file_writer
from ppproject import SampleBuffer, spectrogram
from ppproject.native.lowpassfilter import LowPassFilter as NativeLowPassFilter
from ppproject.native.noiseremoval import waveform_noise_removal
from numpy import array
from matplotlib import pyplot

reader = open_file_reader("segment-1-wav.wav")

print("Length in frames", reader.get_total_frames())
print("Length in time", reader.get_duration())

print("Reading to buffer...")
buf = reader.read_all()
print("Buffer size", buf.get_size())

print("Noise removal...")
new_buf = SampleBuffer(buf.get_sample_rate(), buf.get_channels())
waveform_noise_removal(buf, new_buf, N=-0.5, window_size=3)
#waveform_noise_removal(buf, new_buf, N=0.0)

print("Rendering spectrogram...")
#spectrogram.render_to_file(buf=buf, filename="spec.png")
#spec_array = spectrogram.render_to_array(buf=new_buf, grayscale = True, width = 640, height = 480, dyn_range = 180.0)
spec_array = spectrogram.render_to_array(buf=new_buf, grayscale = False, width = 640, height = 480, dyn_range = 180.0)
#print(spec_array)

print("Saving file...")
#pyplot.imsave("test.png", spec_array, vmin = 0.0, vmax = 1.0, cmap = pyplot.cm.gray)
pyplot.imsave("test.png", spec_array, vmin = 0.0, vmax = 1.0)
of = open_file_writer("test.wav", sample_rate = new_buf.get_sample_rate(), channels = new_buf.get_channels())
of.write_frames(new_buf, new_buf.get_size())
of.close()

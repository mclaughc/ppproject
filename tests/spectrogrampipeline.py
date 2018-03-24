import ppproject
import logging

logging.basicConfig(format="[%(asctime)s] %(levelname)s:%(name)s: %(message)s", level = logging.DEBUG)

pipeline = ppproject.Pipeline()

current = ppproject.sources.FileSource("tests/wavs/ibis.wav", name = "input")
pipeline.add_stage(current)

next = ppproject.filters.Splitter(2.0, 0.0, "splitter")
next.set_input_channel("input", current.get_output_channel("output"))
pipeline.add_stage(next)
current = next

#next = ppproject.filters.LowPassFilter(cutoff_freq = 1000.0, num_taps = 128)
#next.set_input_channel("input", current.get_output_channel("output"))
#pipeline.add_stage(next)
#current = next

#next = ppproject.outputs.NullOutput(name = "dummy sink")
#next = ppproject.outputs.FileOutput()
next = ppproject.outputs.SpectrogramImage()
next.set_input_channel("input", current.get_output_channel("output"))
pipeline.add_stage(next)
current = next

#pipeline.run(5)
pipeline.run_until_complete()



import ppproject
import logging

logging.basicConfig(format="[%(asctime)s] %(levelname)s:%(name)s: %(message)s", level = logging.DEBUG)

pipeline = ppproject.Pipeline()

source = ppproject.sources.FileSource("tests/wavs/ibis.wav", name = "input")
pipeline.add_stage(source)

splitter = ppproject.filters.Splitter(2.0, 1.0, "splitter")
splitter.set_input_channel("input", source.get_output_channel("output"))
pipeline.add_stage(splitter)

sink = ppproject.outputs.NullOutput(name = "dummy sink")
sink.set_input_channel(splitter.get_output_channel("output"))
pipeline.add_stage(sink)

#pipeline.run(5)
pipeline.run_until_complete()

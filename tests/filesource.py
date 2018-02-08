import ppproject
import logging

logging.basicConfig(format="[%(asctime)s] %(levelname)s:%(name)s: %(message)s", level = logging.DEBUG)

pipeline = ppproject.Pipeline()

source = ppproject.sources.FileSource("tests/wavs/ibis.wav", name = "input")
pipeline.add_stage(source)

sink = ppproject.outputs.NullOutput(name = "dummy sink")
sink.set_input_channel(source.get_output_channel("output"))
pipeline.add_stage(sink)

print("running until completion...")
pipeline.run_until_complete()

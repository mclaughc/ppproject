# Pipeline System Design Overview
The framework represents a data processing system as a series of "nodes", known as *stages*, linked together in a graph to form a *pipeline*.

A stage represents a node in the pipeline graph. Each node in the graph can perform one or more roles, including producing data, writing data to an output file, or processing/manipulating the data in some way. Usually, this data will be in the form of audio samples. These audio samples are passed between nodes through *channels*.

In addition to these channels, metadata can also be passed between filter nodes. In this instance, metadata is a simple key-value pair store. By providing such a store, stages which execute later in the pipeline can be passed additional information from earlier stages, altering behavior. For example, the output filename after a Splitter filter could be based on current the segment number, or offset in the original file.

Channels also provide an "end-of-stream" flag. When this flag is set by a preceding filter, it means that any data produced in the next pipeline iteration is not related to previously-seen data, and stages should reset their internal state. The end-of-stream flag can also be used by outputs to signal that the current file should be closed, and a new output file opened.

Nodes which perform data processing/manipulating should use *Filter* as a base class. The filter class provides a common interface for assigning input and output channels to the node. A filter can provide more than one input or output channel, however in most cases a single channel will be sufficient. An example where this is required could be a filter which splits the incoming audio stream into two output streams, one above and below a specified cutoff frequency.

Pipeline stages have two flags; "ready" and "done". The ready flag indicates that the pipeline has been configured correctly, all required properties have been set, and no errors were detected which would prevent processing. The done flag is set when a node has completed processing all available input. When all stages within a pipeline have completed and set the corresponding done flag, the program will exit. It is possible for the done flag to revert to False if additional input becomes available after previously being set to True.

# Getting Started
TODO - Write me!

# Creating a Pipeline
TODO - Write me!

# Creating a New Filter
TODO - Write me!

# Built-in Sources
The following section describes the 
## pproject.sources.FileSource
FileSource reads audio samples from both uncompressed and compressed file sources, including .wav.  
*Outputs:* _output_
- *source.set_filename(filename)*
- *source.set_chunk_size(chunk_size)*

# Built-in Outputs
## ppproject.outputs.FileOutput
FileOutput writes incoming audio samples to an output file, based on a template filename.  
*Inputs:* _input_
- *output.set_filename_template(filename_template)*

## ppproject.outputs.NullOutput
NullOutput acceps incoming audio samples, but immediately discards them.  
*Inputs:* _input_

## ppproject.outputs.SpectrogramImage
SpectrogramImage generates a spectrogram in .png format of the incoming audio stream, and saves images based on a template filename. Each spectrogram is generated when the preceding filter is "done", or the end-of-stream flag is set.
*Inputs:* _input_
- *output.set_filename_template(filename_template)*
- *output.set_image_size(width, height)*
- *output.set_border_enabled(enable)*
- *output.set_greyscale_enabled(enable)*
- *output.set_log_freq_enabled(enable)*
- *output.set_freq_range(min_freq, max_freq)*
- *output.set_fft_freq(fft_freq)*
- *output.set_dynamic_range(range_db)*
- *output.set_window_func(func)*

# Built-in Filters
## pproject.filters.Splitter
Splitter breaks an incoming audio stream into multiple output streams by manipulating the end-of-stream flag. The splitter filter can optionally copy overlapping samples of a specified duration between streams, providing a "sliding window".
*Inputs:* _input_
*Outputs:* _output_
- *output.set_split_length(split_length)*
- *output.set_overlap_duration(overlap_duration)*
 
# Interface Specification
The following lists define the methods available as part of the pipeline interface.
## ppproject.Pipeline
 - *pipeline.add_stage(stage)*  
   Adds the specified stage to the processing pipeline.
 - *pipeline.run(num_iterations = 1)*  
   Executes the processing pipeline for the specified number of iterations. If the pipeline finishes before _num_iterations_ is executed, this method will return early.
 - *pipeline.run_until_complete()*  
   Executes the processing pipeline until all filters return a "done" flag. This ensures that all available input has been processed, and all output has been produced.

## ppproject.PipelineStage
 - *stage.get_name()*  
   Returns the name of the stage object in pipeline. This name must be unique within the pipeline graph the stage is added to.
 - *stage.is_ready()*  
   Virtual method which returns true if the stage is considered "ready". If the ready state of a stage is False when it is added to a the pipeline graph, *make_ready()* will be called.
 - *stage.make_ready()*  
   Virtual method which is called when this pipeline stage is added to the pipeline graph. This enables derived classes to perform "late initialization", for example, opening a file from the disk. This would not be possible to do at object construction time, as the required attributes would not have been configured yet. If an error occurs during late initialization, an appropriate exception should be raised.
 - *stage.is_done()*  
   Virtual method which is called by the pipeline executor to determine whether the pipeline has completed. As aforementioned, the done state can return to False after being set to True.
 - *stage.run()*  
   Called once per pipeline iteration. Stages should consume as much input is available, unless they wish to insert an end-of-stream marker.

## ppproject.Filter
 - *filter.set_input_channel(name, channel)*  
   Connects the specified channel to the named input of a filter. Usually a filter will only have one input, named "input". Only a single channel can be connected to each input.
 - *filter.get_output_channel(name, channel)*  
   Returns the specified output channel of a filter. Usually a filter will only have one output, named "output".

## ppproject.channels.AudioChannel
 - *channel.get_sample_rate()*   
   Returns the sample rate of the audio streamed through this channel.  
 - *channel.get_channels()*  
   Returns the number of channels for audio streamed through this channel.
 - *channel.get_buffer()*  
   Returns the buffer containing the samples streamed through this channel.
 - *channel.get_metadata()*  
   Returns the metadata container for audio streamed through this channel.
 - *channel.get_name()*  
   Returns the (optional) name of this audio channel.
 - *channel.end_of_stream()*  
   Returns True if the end-of-stream flag has been set.
 - *channel.set_end_of_stream()*  
   Sets the end-of-stream flag, instructing future streams to complete processing.
 - *channel.clear_end_of_stream()*  
   Clears the end-of-stream flag, instructing future streams that a new file is being processed.
 - *channel.is_empty_and_end_of_stream()*  
   Returns true when there are no frames in the buffer, and the end-of-stream flag is set.

## ppproject.metadata
The Metadata class represents a series of key/value properties associated with an object of some description. These objects could include audio chunks, images, or binary data.
 - *metadata.get_property(name)*  
  Returns the property named name from within this object. If a property with this name does not exist, KeyError will be raised.
 - *metadata.get_property_default(name, default_value)*  
Returns the property named name from within this object. If a property with this name does not exist, the value default will be returned.
 - *metadata.set_property(name, value)*  
Sets property name to value within this object.
 - *metadata.clear_properties()*  
Removes all properties from this object.
 - *metadata.copy_properties(from_metadata, clear_properties = False)*  
Copies the properties from from_metadata to this metadata object. All properties are "deep copied", i.e. no references are maintained to the original object. If clear_properties is set, any existing properties set in this object will be removed.
 - *metadata.remove_property(name)*  
Removes the property name from this object, if it exists.

- *metadata.format string(name)*  
Returns a string with fields in the template string substituted with property names. Property names are represented by percent symbol in the template string. For example, template_for_%name%. %name% will be replaced with the property 'name' from this object. If the property does not exist in this object, in empty string will be substituted instead. If the caller wishes to insert a percent symbol in the template string, use %%.
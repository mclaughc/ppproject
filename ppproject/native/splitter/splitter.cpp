#include <Python.h>
#include <algorithm>
#include <cassert>
#include <cstdio>
#include <sndfile.h>
#include <vector>

#include "shared/format.h"
#include "shared/log.h"
#include "shared/string_helpers.h"
#include "shared/types.h"
Log_SetChannel(_splitter);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// _splitter.split_file
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
static SNDFILE* OpenInputWaveFile(const char* filename, SF_INFO* info)
{
  info->frames = 0;
  info->samplerate = 0;
  info->channels = 0;
  info->format = 0;
  info->sections = 0;
  info->seekable = 1;

  SNDFILE* f = sf_open(filename, SFM_READ, info);
  if (!f)
  {
    Log_ErrorPrintf("Failed to open %s as wave input", filename);
    return nullptr;
  }

  return f;
}

static SNDFILE* OpenOutputWaveFile(const char* filename, u32 channels, u32 sample_rate, SF_INFO* info)
{
  info->frames = 0;
  info->samplerate = sample_rate;
  info->channels = channels;
  info->format = SF_FORMAT_WAV | SF_FORMAT_PCM_16;
  info->sections = 0;
  info->seekable = 0;

  SNDFILE* f = sf_open(filename, SFM_WRITE, info);
  if (!f)
  {
    Log_ErrorPrintf("Failed to open %s as wave output", filename);
    return nullptr;
  }

  return f;
}

static bool WriteSplitFile(const SF_INFO* input_info, SNDFILE* input_file, sf_count_t start_frame,
                           sf_count_t num_frames, const char* output_prefix)
{
  u32 input_offset_time = u32(start_frame) / u32(input_info->samplerate);
  u32 input_offset_mins = input_offset_time / 60;
  u32 input_offset_secs = input_offset_time % 60;

  std::string output_filename =
    StringFromFormat("%s_split_%04u_%02u.wav", output_prefix, input_offset_mins, input_offset_secs);
  SF_INFO output_info = {};
  SNDFILE* output_file =
    OpenOutputWaveFile(output_filename.c_str(), input_info->channels, input_info->samplerate, &output_info);
  if (!output_file)
    return false;

  Log_DevPrintf("Writing %u frames starting from frame %u (%04u:%02u) to %s", u32(num_frames), u32(start_frame),
                input_offset_mins, input_offset_secs, output_filename.c_str());

  if (sf_seek(input_file, start_frame, SEEK_SET) < 0)
  {
    Log_ErrorPrintf("sf_seek() failed: %s", sf_strerror(input_file));
    return false;
  }

  constexpr sf_count_t CHUNK_SIZE = 1024;
  std::vector<float> chunk_data(CHUNK_SIZE * input_info->channels);

  for (sf_count_t current_frame = 0; current_frame < num_frames; current_frame += CHUNK_SIZE)
  {
    sf_count_t this_chunk = std::min(CHUNK_SIZE, num_frames - current_frame);

    if (sf_readf_float(input_file, chunk_data.data(), this_chunk) != this_chunk)
    {
      Log_ErrorPrintf("sf_readf_float() failed: %s", sf_strerror(input_file));
      sf_close(output_file);
      return false;
    }

    if (sf_writef_float(output_file, chunk_data.data(), this_chunk) != this_chunk)
    {
      Log_ErrorPrintf("sf_writef_float() failed: %s", sf_strerror(input_file));
      sf_close(output_file);
      return false;
    }
  }

  // sync before closing, needed?
  sf_write_sync(output_file);
  sf_close(output_file);
  return true;
}

static PyObject* Py_splitter_split_file(PyObject* self, PyObject* args)
{
  char* input_filename;
  char* output_filename_prefix;
  float segment_length;
  float overlap_length;

  if (!PyArg_ParseTuple(args, "ssff", &input_filename, &output_filename_prefix, &segment_length, &overlap_length))
    return nullptr;

  assert(segment_length > 0.0f && overlap_length >= 0.0f);

  SF_INFO input_info = {};
  SNDFILE* input_file = OpenInputWaveFile(input_filename, &input_info);
  if (!input_file)
    return nullptr;

  // Determine how many frames per segment based on the input sample rate.
  sf_count_t frames_per_segment = static_cast<sf_count_t>(double(segment_length) * input_info.samplerate);
  sf_count_t overlap_frames = static_cast<sf_count_t>(double(overlap_length) * input_info.samplerate);
  assert(frames_per_segment > 0 && overlap_frames < frames_per_segment);

  // Process frames in the input.
  const sf_count_t increment_frame_count = frames_per_segment - overlap_frames;
  const sf_count_t num_frames = input_info.frames;
  u32 num_segments = 0;

  for (sf_count_t start_frame = 0; start_frame < num_frames; start_frame += increment_frame_count)
  {
    sf_count_t this_segment = std::min(frames_per_segment, num_frames - start_frame);
    if (!WriteSplitFile(&input_info, input_file, start_frame, this_segment, output_filename_prefix))
    {
      sf_close(input_file);
      return nullptr;
    }

    num_segments++;
  }

  sf_close(input_file);

  PyObject* ret = Py_BuildValue("i", num_segments);
  return ret;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// _splitter
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
static const char module_docstring[] = "Module documentation string";

static PyMethodDef module_methods[] = {
  {"split_file", Py_splitter_split_file, METH_VARARGS, "docs"}, {},
};

static PyModuleDef module_def = {
  PyModuleDef_HEAD_INIT,
  "_splitter",      // name
  module_docstring, // docstring
  0,                // size
  module_methods,   // methods
  nullptr,          // slots
  nullptr,          // gc traverse
  nullptr,          // gc clear
  nullptr           // free
};

PyMODINIT_FUNC PyInit__splitter()
{
  PyObject* mod = PyModule_Create(&module_def);
  if (!mod)
    return nullptr;

  Log::SetFilterLevel(LOGLEVEL_TRACE);
  Log::SetConsoleOutputParams(true);

  return mod;
}

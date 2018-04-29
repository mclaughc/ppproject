#include "samplebuffer.h"
#include "shared/string_helpers.h"
#include <algorithm>
#include <cassert>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <memory>
#include <pybind11/pybind11.h>
#include <sndfile.h>
#include <vector>

// TODO: Make this class have a header interface, accessible to other modules.
class FileReader
{
public:
  ~FileReader();

  // Opens an existing file for reading.
  static FileReader* OpenFile(const std::string& filename);

  // Closes the file, preventing further operations.
  bool IsOpen() const { return m_file != nullptr; }
  void Close();

  // Accessors.
  const std::string& GetFileName() const { return m_filename; }
  int GetSampleRate() const { return info.samplerate; }
  int GetChannels() const { return info.channels; }
  int64_t GetTotalFrames() const { return info.frames; }
  int64_t GetCurrentFrame() const { return m_current_frame; }
  int64_t GetRemainingFrames() const { return info.frames - m_current_frame; }
  double GetDuration() const { return info.frames / static_cast<double>(info.samplerate); }
  double GetCurrentTime() const { return m_current_frame / static_cast<double>(info.samplerate); }
  double GetRemainingTime() const { return GetRemainingFrames() / static_cast<double>(info.samplerate); }

  // Mutators
  void Seek(int64_t frame_number);
  void SeekTime(double time);

  // Reads the specified number of frames to the provided buffer.
  // Returns the actual number of frames read.
  int64_t ReadFrames(SampleBuffer* buf, int64_t count);

  // Reads the entire file to a single buffer.
  // Restores the previous position afterwards.
  SampleBuffer* ReadAll() const;

private:
  FileReader(const std::string& filename);
  void CheckIsOpen();

  std::string m_filename;
  SNDFILE* m_file = nullptr;
  SF_INFO info = {};
  int64_t m_current_frame = 0;
  int m_chunk_size = 0;
};

FileReader::FileReader(const std::string& filename) : m_filename(filename) {}

FileReader::~FileReader()
{
  if (m_file)
    sf_close(m_file);
}

void FileReader::CheckIsOpen()
{
  if (!IsOpen())
    throw std::runtime_error("Attempt to manipulate a FileReader which is already closed.");
}

FileReader* FileReader::OpenFile(const std::string& filename)
{
  std::unique_ptr<FileReader> fr = std::unique_ptr<FileReader>(new FileReader(filename));
  fr->info.seekable = 1;
  fr->m_file = sf_open(filename.c_str(), SFM_READ, &fr->info);
  if (!fr->m_file)
    throw std::runtime_error(StringFromFormat("Failed to open %s: %s", filename.c_str(), sf_strerror(nullptr)));

  return fr.release();
}

void FileReader::Close()
{
  CheckIsOpen();
  sf_close(m_file);
  m_file = nullptr;
}

void FileReader::Seek(int64_t frame_number)
{
  if (frame_number < 0)
    throw std::invalid_argument("attempt to seek to a negative frame number");
  if (frame_number > info.frames)
    throw std::invalid_argument("attempt to seek past the end of the file");

  if (sf_seek(m_file, frame_number, SEEK_SET) != frame_number)
    throw std::runtime_error("seek failed");
  m_current_frame = frame_number;
}

void FileReader::SeekTime(double time)
{
  int64_t frame_number = std::floor(time * static_cast<double>(info.samplerate));
  Seek(frame_number);
}

int64_t FileReader::ReadFrames(SampleBuffer* buf, int64_t count)
{
  // TODO: Optimize this. Remove the extra copy.
  std::unique_ptr<double[]> tmp = std::make_unique<double[]>(count * buf->GetChannels());
  int64_t actual_count = sf_readf_double(m_file, tmp.get(), count);
  if (actual_count > 0)
    buf->WriteFrames<double>(tmp.get(), static_cast<int>(actual_count));
  m_current_frame += actual_count;
  return actual_count;
}

SampleBuffer* FileReader::ReadAll() const
{
  // seek back to start, remembering original position
  int64_t old_position = m_current_frame;
  if (sf_seek(m_file, 0, SEEK_SET) < 0)
    throw std::runtime_error("seek to start failed");

  // read the samples into a temporary buffer
  int64_t frame_count = info.frames;
  int64_t sample_count = frame_count * GetChannels();
  std::unique_ptr<double[]> tmp = std::make_unique<double[]>(sample_count);
  if (sample_count > 0 && sf_readf_double(m_file, tmp.get(), frame_count) != frame_count)
    throw std::runtime_error("reading samples from file failed");

  // return to the original position
  if (sf_seek(m_file, 0, SEEK_SET) != m_current_frame)
    throw std::runtime_error("seek to original position failed");

  // copy into a sample buffer
  std::unique_ptr<SampleBuffer> buf = std::make_unique<SampleBuffer>(GetSampleRate(), GetChannels());
  if (sample_count > 0)
    buf->WriteFrames<double>(tmp.get(), static_cast<int>(frame_count));
  return buf.release();
}

// TODO: Output formats.
class FileWriter
{
public:
  ~FileWriter();

  // Opens an existing file for writing.
  static FileWriter* OpenFile(const std::string& filename, int sample_rate, int channels);

  // Closes the file, preventing further operations.
  bool IsOpen() const { return m_file != nullptr; }
  void Close();

  // Accessors.
  const std::string& GetFileName() const { return m_filename; }
  int GetSampleRate() const { return info.samplerate; }
  int GetChannels() const { return info.channels; }
  int64_t GetTotalFrames() const { return info.frames; }
  double GetDuration() const { return info.frames / static_cast<double>(info.samplerate); }

  // Writes the specified number of frames to the provided buffer.
  // Returns the actual number of frames written.
  int64_t WriteFrames(SampleBuffer* buf, int64_t count);

  // Writes the entire buffer to the file, and advances its read pointer.
  int64_t WriteAll(SampleBuffer* buf);

private:
  FileWriter(const std::string& filename);
  void CheckIsOpen();

  std::string m_filename;
  SNDFILE* m_file = nullptr;
  SF_INFO info = {};
  int64_t m_current_frame = 0;
  int m_chunk_size = 0;
};

FileWriter::FileWriter(const std::string& filename) : m_filename(filename) {}

FileWriter::~FileWriter()
{
  if (m_file)
    sf_close(m_file);
}

void FileWriter::CheckIsOpen()
{
  if (!IsOpen())
    throw std::runtime_error("Attempt to manipulate a FileWriter which is already closed.");
}

FileWriter* FileWriter::OpenFile(const std::string& filename, int sample_rate, int channels)
{
  std::unique_ptr<FileWriter> fr = std::unique_ptr<FileWriter>(new FileWriter(filename));
  fr->info.format = SF_FORMAT_WAV | SF_FORMAT_PCM_16;
  fr->info.samplerate = sample_rate;
  fr->info.channels = channels;
  fr->m_file = sf_open(filename.c_str(), SFM_WRITE, &fr->info);
  if (!fr->m_file)
    throw std::runtime_error(StringFromFormat("Failed to open %s: %s", filename.c_str(), sf_strerror(nullptr)));

  return fr.release();
}

void FileWriter::Close()
{
  CheckIsOpen();
  sf_close(m_file);
  m_file = nullptr;
}

int64_t FileWriter::WriteFrames(SampleBuffer* buf, int64_t count)
{
  if (buf->IsEmpty())
    return 0;

  return sf_writef_float(m_file, buf->GetPeekPointer(0), buf->GetSize());
}

int64_t FileWriter::WriteAll(SampleBuffer* buf)
{
  if (buf->GetSize() == 0)
    return 0;

  int64_t actual_count = WriteFrames(buf, buf->GetSize());
  if (actual_count > 0)
    buf->RemoveFrames(static_cast<int>(actual_count));

  return actual_count;
}

namespace py = pybind11;
PYBIND11_MODULE(audiofile, m)
{
  m.def("open_file_reader", &FileReader::OpenFile, py::arg("filename"));
  m.def("open_file_writer", &FileWriter::OpenFile, py::arg("filename"), py::arg("sample_rate"), py::arg("channels"));

  py::class_<FileReader> filereader(m, "FileReader");
  filereader.def("is_open", &FileReader::IsOpen)
    .def("close", &FileReader::Close)
    .def("get_file_name", &FileReader::GetFileName)
    .def("get_sample_rate", &FileReader::GetSampleRate)
    .def("get_channels", &FileReader::GetChannels)
    .def("get_total_frames", &FileReader::GetTotalFrames)
    .def("get_remaining_frames", &FileReader::GetRemainingFrames)
    .def("get_duration", &FileReader::GetDuration)
    .def("get_current_time", &FileReader::GetCurrentTime)
    .def("get_remaining_time", &FileReader::GetRemainingTime)
    .def("seek", &FileReader::Seek)
    .def("seek_time", &FileReader::SeekTime)
    .def("read_frames", &FileReader::ReadFrames)
    .def("read_all", &FileReader::ReadAll);

  py::class_<FileWriter> filewriter(m, "FileWriter");
  filewriter.def("is_open", &FileWriter::IsOpen)
    .def("close", &FileWriter::Close)
    .def("get_file_name", &FileWriter::GetFileName)
    .def("get_sample_rate", &FileWriter::GetSampleRate)
    .def("get_channels", &FileWriter::GetChannels)
    .def("get_total_frames", &FileWriter::GetTotalFrames)
    .def("get_duration", &FileWriter::GetDuration)
    .def("write_frames", &FileWriter::WriteFrames)
    .def("write_all", &FileWriter::WriteAll);
}
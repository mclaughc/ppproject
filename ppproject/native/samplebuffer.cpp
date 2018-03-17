#include "samplebuffer.h"
#include <pybind11/pybind11.h>

namespace py = pybind11;
PYBIND11_MODULE(samplebuffer, m)
{
  py::class_<SampleBuffer> samplebuffer(m, "SampleBuffer");
  samplebuffer.def(py::init<int, int>())
    .def("get_read_position", &SampleBuffer::GetReadPosition)
    .def("get_write_position", &SampleBuffer::GetWritePosition)
    .def("set_read_position", &SampleBuffer::SetReadPosition)
    .def("get_size", &SampleBuffer::GetSize)
    .def("is_empty", &SampleBuffer::IsEmpty)
    .def("copy_frames", &SampleBuffer::CopyFrames)
    .def("remove_frames", &SampleBuffer::RemoveFrames)
    .def("clear", &SampleBuffer::Clear)
    .def("swap", &SampleBuffer::Swap)
    .def("insert_silence", &SampleBuffer::InsertSilence)
    .def("insert_silence_frames", &SampleBuffer::InsertSilenceFrames)
    .def("push_float", &SampleBuffer::PyPush<float>)
    .def("peek_float", &SampleBuffer::PyPeek<float>)
    .def("pop_float", &SampleBuffer::PyPop<float>)
    .def("push_double", &SampleBuffer::PyPush<double>)
    .def("peek_double", &SampleBuffer::PyPeek<double>)
    .def("pop_double", &SampleBuffer::PyPop<double>)
    .def("push_s8", &SampleBuffer::PyPush<int8_t>)
    .def("peek_s8", &SampleBuffer::PyPeek<int8_t>)
    .def("pop_s8", &SampleBuffer::PyPop<int8_t>)
    .def("push_s16", &SampleBuffer::PyPush<int16_t>)
    .def("peek_s16", &SampleBuffer::PyPeek<int16_t>)
    .def("pop_s16", &SampleBuffer::PyPop<int16_t>)
    .def("push_s32", &SampleBuffer::PyPush<int32_t>)
    .def("peek_s32", &SampleBuffer::PyPeek<int32_t>)
    .def("pop_s32", &SampleBuffer::PyPop<int32_t>)
    .def("push_u8", &SampleBuffer::PyPush<uint8_t>)
    .def("peek_u8", &SampleBuffer::PyPeek<uint8_t>)
    .def("pop_u8", &SampleBuffer::PyPop<uint8_t>)
    .def("push_u16", &SampleBuffer::PyPush<uint16_t>)
    .def("peek_u16", &SampleBuffer::PyPeek<uint16_t>)
    .def("pop_u16", &SampleBuffer::PyPop<uint16_t>)
    .def("push_u32", &SampleBuffer::PyPush<uint32_t>)
    .def("peek_u32", &SampleBuffer::PyPeek<uint32_t>)
    .def("pop_u32", &SampleBuffer::PyPop<uint32_t>);
}

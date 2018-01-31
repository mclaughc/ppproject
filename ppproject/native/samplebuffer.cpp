#include "samplebuffer.h"
#include <pybind11/pybind11.h>

namespace py = pybind11;
PYBIND11_MODULE(samplebuffer, m)
{
  py::class_<SampleBuffer> samplebuffer(m, "SampleBuffer");
  samplebuffer.def(py::init<int, int>())
    .def("get_position", &SampleBuffer::GetPosition)
    .def("set_position", &SampleBuffer::SetPosition)
    .def("get_size", &SampleBuffer::GetSize)
    .def("get_remaining_frames", &SampleBuffer::GetRemainingFrames);
}

#include <nanobind/nanobind.h>
#include <parameterization.hpp>
#include <nanobind/eigen/dense.h>

namespace nb = nanobind;

void init_harmonic(nb::module_ &m) {
  m.def("harmonic", &harmonic);
}

#include <nanobind/nanobind.h>
#include <parameterization.hpp>
#include <nanobind/eigen/dense.h>

namespace nb = nanobind;

void init_tutte(nb::module_ &m) {
  m.def("tutte", &tutte_ext);
}

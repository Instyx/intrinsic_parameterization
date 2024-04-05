#include <nanobind/nanobind.h>
#include <parameterization.hpp>
#include <nanobind/eigen/dense.h>

namespace nb = nanobind;

void init_conformal(nb::module_ &m) {
  m.def("conformal", &LSCM);
}

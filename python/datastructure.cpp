#include <nanobind/nanobind.h>
#include <datageo.hpp>
#include <nanobind/eigen/dense.h>

namespace nb = nanobind;

void init_datastructure(nb::module_ &m) {
  nb::class_<DataGeo>(m, "DataGeo")
            .def(nb::init<>())
            .def_ro("V", &DataGeo::V)
            .def_ro("F", &DataGeo::F)
            ;
}

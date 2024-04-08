#include <nanobind/nanobind.h>
#include <nanobind/eigen/dense.h>
#include <nanobind/eigen/sparse.h>

namespace nb = nanobind;

void init_datastructure(nb::module_ &m);
void init_create_datastructure(nb::module_ &m);
void init_harmonic(nb::module_ &m);
void init_conformal(nb::module_ &m);
void init_tutte(nb::module_ &m);
void init_adjacency(nb::module_ &m);

NB_MODULE(pyiparam, m) {
    init_datastructure(m);
    init_create_datastructure(m);
    init_harmonic(m);
    init_conformal(m);
    init_tutte(m);
    init_adjacency(m);
}

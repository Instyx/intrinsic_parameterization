#include <nanobind/nanobind.h>
#include <harmonic.hpp>
#include <nanobind/eigen/dense.h>
#include <nanobind/eigen/sparse.h>

namespace nb = nanobind;

Eigen::Matrix<double,-1,2> wrap_harmonic(
           const Eigen::SparseMatrix<double>& L,
           const nb::DRef<const Eigen::VectorXi>& B){
   Eigen::Matrix<double,-1,2> X;
   harmonic<Eigen::Matrix<double,-1,2>>(L, B, X, 0);
   return X;
}

void init_harmonic(nb::module_ &m) {
  m.def("harmonic", &wrap_harmonic);
}

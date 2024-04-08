#include <nanobind/nanobind.h>
#include <nanobind/eigen/dense.h>
#include <nanobind/eigen/sparse.h>
#include <nanobind/stl/tuple.h>

#include <adjacency.hpp>

namespace nb = nanobind;

Eigen::VectorXi wrap_bdy_loop(
    const nb::DRef<Eigen::Matrix<int,-1,3>>& F,
    const nb::DRef<Eigen::Matrix<int,-1,3>>& TT,
    const nb::DRef<Eigen::VectorXi>& VT,
    const nb::DRef<Eigen::VectorXi>& VTi){
   Eigen::VectorXi B;
   bdy_loop(F, TT, VT, VTi, B);
   return B;
}

Eigen::Matrix<int,-1,3> wrap_tt_adjacency(const nb::DRef<Eigen::Matrix<int,-1,3>>& F,
                  const nb::DRef<Eigen::VectorXi>& VT,
                  const nb::DRef<Eigen::VectorXi>& VTi){
   Eigen::Matrix<int,-1,3> TT;
   tt_adjacency(F, VT, VTi, TT);
   return TT;
}

std::tuple<Eigen::VectorXi,Eigen::VectorXi> wrap_vv_adjacency(
                  const nb::DRef<Eigen::Matrix<int,-1,3>>& F,
                  const nb::DRef<Eigen::Matrix<double,-1,3>>& V,
                  const nb::DRef<Eigen::Matrix<int,-1,3>>& TT){
    Eigen::VectorXi VV, VVi;
   vv_adjacency(F, V, TT, VV, VVi);
   return std::make_tuple(VV,VVi);
}

std::tuple<Eigen::VectorXi,Eigen::VectorXi> wrap_vt_adjacency(
      const nb::DRef<Eigen::Matrix<int,-1,-1>>& F,
      const nb::DRef<Eigen::Matrix<double,-1,-1>>& V){
   Eigen::VectorXi VT, VTi;
   vt_adjacency(F, V, VT, VTi);
   return std::make_tuple(VT,VTi);
}

void init_adjacency(nb::module_ &m) {
  m.def("tt_adjacency", &wrap_tt_adjacency);
  m.def("vv_adjacency", &wrap_vv_adjacency);
  m.def("vt_adjacency", &wrap_vt_adjacency);
  m.def("bdy_loop", &wrap_bdy_loop);
}

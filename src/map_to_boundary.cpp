#include "map_to_boundary.hpp"

#include <type_traits>
#include <iostream>
#include <math.h>


template <typename DerivedIn,typename DerivedOut>
void circle_boundary_proportional(
  const Eigen::MatrixBase<DerivedIn>& V,
  const Eigen::Ref<const Eigen::VectorXi> B,
  Eigen::MatrixBase<DerivedOut>& X)
{
  typedef typename DerivedOut::Scalar Scalar;
  // Map boundary to unit circle
  std::vector<double> len(B.size());
  len[0] = 0.;

  for (int i = 1; i < B.size(); i++)
  {
    len[i] = len[i-1] + (V.row(B[i-1]) - V.row(B[i])).norm();
  }
  double total_len = len[len.size()-1] + (V.row(B[0]) - V.row(B[B.size()-1])).norm();
  X.derived().resize(V.rows(),2);
  for (int i = 0; i < B.size(); i++)
  {
    double frac = -len[i] * 2. * M_PI / total_len;
    X.row(B[i]) << (Scalar) cos(frac), (Scalar) sin(frac);
  }
}


template <typename Derived>
void circle_boundary(
  const Eigen::Ref<const Eigen::VectorXi> B,
  Eigen::MatrixBase<Derived>& X)
{
  typedef typename Derived::Scalar Scalar;
  int nv = X.rows();
  int bc = B.size();
  for (int i = 0; i < B.rows(); i++){
    // boundary construction:
    double phi = (2.0*M_PI * -i) / (double)B.rows();
    X.row(B[i]) << (Scalar) cos(phi), (Scalar) sin(phi);
  }
}


template <typename Derived>
void square_boundary(
  const Eigen::Ref<const Eigen::VectorXi> B,
  Eigen::MatrixBase<Derived>& X)
{
  typedef typename Derived::Scalar Scalar;
  int nv = X.rows();
  int bc = B.size();
  X.row(B(0)) << 0, 0;
  X.row(B(bc / 4)) << 0, 1;
  X.row(B(bc*2 / 4)) << 1, 1;
  X.row(B(bc*3 / 4)) << 1, 0;
  for (int i = 0; i < B.rows(); i++){
    // boundary construction:
    int cornerA = bc*3 / 4;
    int cornerB = bc;
    if (i < bc / 4)
    {
        cornerA = 0;
        cornerB = bc / 4;
    }
    else if (i < bc*2 / 4)
    {
        cornerA = bc / 4;
        cornerB = bc*2 / 4;
    }
    else if (i < bc*3 / 4)
    {
        cornerA = bc*2 / 4;
        cornerB = bc*3 / 4;
    }
    double alpha = ((double) i-cornerA)/(cornerB-cornerA);
    X.row(B(i)) << (1-alpha) * X.row(B(cornerA)) + alpha * X.row(B(cornerB%bc));
  }
}


template void circle_boundary_proportional(
  const Eigen::MatrixBase<Eigen::Matrix<double,-1,3>>& V,
  const Eigen::Ref<const Eigen::VectorXi> B,
  Eigen::MatrixBase<Eigen::Matrix<double,-1,2>>& X);

template void circle_boundary_proportional(
  const Eigen::MatrixBase<Eigen::Matrix<double,-1,3>>& V,
  const Eigen::Ref<const Eigen::VectorXi> B,
  Eigen::MatrixBase<Eigen::Matrix<double,-1,-1>>& X);

template void circle_boundary_proportional(
  const Eigen::MatrixBase<Eigen::Matrix<double,-1,-1>>& V,
  const Eigen::Ref<const Eigen::VectorXi> B,
  Eigen::MatrixBase<Eigen::Matrix<double,-1,-1>>& X);

template void circle_boundary_proportional(
  const Eigen::MatrixBase<Eigen::Matrix<double,-1,-1>>& V,
  const Eigen::Ref<const Eigen::VectorXi> B,
  Eigen::MatrixBase<Eigen::Matrix<double,-1,2>>& X);


template void circle_boundary(
  const Eigen::Ref<const Eigen::VectorXi> B,
  Eigen::MatrixBase<Eigen::Matrix<double,-1,2>>& X);

template void circle_boundary(
  const Eigen::Ref<const Eigen::VectorXi> B,
  Eigen::MatrixBase<Eigen::Matrix<double,-1,-1>>& X);


template void square_boundary(
  const Eigen::Ref<const Eigen::VectorXi> B,
  Eigen::MatrixBase<Eigen::Matrix<double,-1,2>>& X);

template void square_boundary(
  const Eigen::Ref<const Eigen::VectorXi> B,
  Eigen::MatrixBase<Eigen::Matrix<double,-1,-1>>& X);


template void circle_boundary_proportional(
  const Eigen::MatrixBase<Eigen::Matrix<float,-1,3>>& V,
  const Eigen::Ref<const Eigen::VectorXi> B,
  Eigen::MatrixBase<Eigen::Matrix<float,-1,2>>& X);

template void circle_boundary_proportional(
  const Eigen::MatrixBase<Eigen::Matrix<float,-1,3>>& V,
  const Eigen::Ref<const Eigen::VectorXi> B,
  Eigen::MatrixBase<Eigen::Matrix<float,-1,-1>>& X);

template void circle_boundary_proportional(
  const Eigen::MatrixBase<Eigen::Matrix<float,-1,-1>>& V,
  const Eigen::Ref<const Eigen::VectorXi> B,
  Eigen::MatrixBase<Eigen::Matrix<float,-1,-1>>& X);

template void circle_boundary_proportional(
  const Eigen::MatrixBase<Eigen::Matrix<float,-1,-1>>& V,
  const Eigen::Ref<const Eigen::VectorXi> B,
  Eigen::MatrixBase<Eigen::Matrix<float,-1,2>>& X);


template void circle_boundary(
  const Eigen::Ref<const Eigen::VectorXi> B,
  Eigen::MatrixBase<Eigen::Matrix<float,-1,2>>& X);

template void circle_boundary(
  const Eigen::Ref<const Eigen::VectorXi> B,
  Eigen::MatrixBase<Eigen::Matrix<float,-1,-1>>& X);


template void square_boundary(
  const Eigen::Ref<const Eigen::VectorXi> B,
  Eigen::MatrixBase<Eigen::Matrix<float,-1,2>>& X);

template void square_boundary(
  const Eigen::Ref<const Eigen::VectorXi> B,
  Eigen::MatrixBase<Eigen::Matrix<float,-1,-1>>& X);

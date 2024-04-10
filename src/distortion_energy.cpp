#include "distortion_energy.hpp"
#include <math.h>
#include <svd.hpp>
#include <iostream>
#include <igl/polar_svd.h>

double dirichlet(const Eigen::Matrix2d &J){
  return J.norm()*J.norm();
}

double symmetric_dirichlet(const Eigen::Matrix2d &J){
  // check also with polar svd instead of SSVD2x2
  Eigen::Matrix2d U, S, VV;
  SSVD2x2(J, U, S, VV);
  return S(0,0)*S(0,0) + 1/(S(0,0)*S(0,0)) + S(1,1)*S(1,1)+1/(S(1,1)*S(1,1));
}

double symmetric_dirichlet_alt(const Eigen::Matrix2d &J){
  Eigen::Matrix2d ri, ti, ui, vi;
  Eigen::Vector2d sing;
  igl::polar_svd(J, ri, ti, ui, sing, vi);
  double s1 = sing(0);
  double s2 = sing(1);
  return pow(s1, 2) + pow(s1, -2) + pow(s2, 2) + pow(s2, -2);
}

double arap(const Eigen::Matrix2d &J){
  Eigen::Matrix2d U, S, VV;
  SSVD2x2(J, U, S, VV);
  Eigen::Matrix2d R = U*VV.transpose();
  return (J-R).norm() * (J-R).norm();
}

double asap(const Eigen::Matrix2d &J){
  /*
  Eigen::Matrix2d ri, ti, ui, vi;
  Eigen::Vector2d sing;
  igl::polar_svd(J, ri, ti, ui, sing, vi);
  double s1 = sing(0);
  double s2 = sing(1);
  return ((pow(s1, 2) + pow(s2, 2)) / (2 * s1 * s2));
  */
  return (J(0,0)-J(1,1)) *  (J(0,0)-J(1,1)) +  (J(1,0)+J(0,1)) * (J(1,0)+J(0,1));
}

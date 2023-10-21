#include "distortion_energy.hpp"
#include <math.h>
#include <svd.hpp>
#include <iostream>
double dirichlet(Eigen::Matrix2d J){
  return J.norm()*J.norm();
}

double symmetric_dirichlet(Eigen::Matrix2d J){
  if(J.determinant()>=0){
    Eigen::Matrix2d U, S, VV;
    SSVD2x2(J, U, S, VV);
    return S(0,0)*S(0,0) + 1/(S(0,0)*S(0,0)) + S(1,1)*S(1,1)+1/(S(1,1)*S(1,1));
  }
  else{
    std::cout << " negative dirichlet infinity " << std::endl;
    return 1e8; // should be inf
  }
}

double arap(Eigen::Matrix2d J){
  Eigen::Matrix2d U, S, VV;
  SSVD2x2(J, U, S, VV);
  Eigen::Matrix2d R = U*VV.transpose();
  return (J-R).norm() * (J-R).norm();
   
}

double asap(Eigen::Matrix2d J){
  return (J(0,0)-J(1,1)) *  (J(0,0)-J(1,1)) +  (J(1,0)+J(0,1)) * (J(1,0)+J(0,1)); 

}



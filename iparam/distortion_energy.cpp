#include "distortion_energy.hpp"
#include <math.h>

void SSVD2x2(const Eigen::Matrix2d& J, Eigen::Matrix2d& U, Eigen::Matrix2d& S, Eigen::Matrix2d& V) {
	double e = (J(0) + J(3))*0.5;
	double f = (J(0) - J(3))*0.5;
	double g = (J(1) + J(2))*0.5;
	double h = (J(1) - J(2))*0.5;
	double q = sqrt((e*e) + (h*h));
	double r = sqrt((f*f) + (g*g));
	double a1 = atan2(g, f);
	double a2 = atan2(h, e);
	double rho = (a2 - a1)*0.5;
	double phi = (a2 + a1)*0.5;

	S(0) = q + r;
	S(1) = 0;
	S(2) = 0;
	S(3) = q - r;

	double c = cos(phi);
	double s = sin(phi);
	U(0) = c;
	U(1) = s;
	U(2) = -s;
	U(3) = c;

	c = cos(rho);
	s = sin(rho);
	V(0) = c;
	V(1) = -s;
	V(2) = s;
	V(3) = c;
}


double dirichlet(Eigen::Matrix2d J){
  return J.norm()*J.norm();
}

double symmetricDirichlet(Eigen::Matrix2d J){
  if(J.determinant()>=0){
    return J.norm() * J.norm() + J.inverse().norm()*J.inverse().norm();
  }
  else{
    return 1e7; // should be inf
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



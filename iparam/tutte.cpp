#include "tutte.hpp"

#include <Eigen/CholmodSupport>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>

#include <type_traits>
#include <iostream>
#include <math.h>


template <typename Derived>
void tutte_circle_boundary(
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
void tutte_square_boundary(
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



template <typename Derived>
void tutte(const Eigen::Ref<const Eigen::VectorXi> VV,
           const Eigen::Ref<const Eigen::VectorXi> VVi,
           const Eigen::Ref<const Eigen::VectorXi> B,
           //Eigen::Ref<Eigen::MatrixXd> X)
           Eigen::MatrixBase<Derived> & X,
           int boundary_type)
{
    typedef typename Derived::Scalar Scalar;

    const int nv = VV.maxCoeff() + 1;
    // return positions
    // Eigen::MatrixBase<Derived>& X = const_cast< Eigen::MatrixBase<Derived>& >(X_);
    X.derived().resize(nv,2);
    switch (boundary_type) {
      case 2:
        break;
      case 1:
        tutte_square_boundary(B, X);
        break;
      case 0:
      default:
        tutte_circle_boundary(B, X);
    }
    std::vector<int> idx(nv,0);
    int bc = B.size();
    if (bc == nv){
      return;
    }
    for (int i = 0; i < bc; i++){
      idx[B[i]] = -i-1;
    }
    int row = 0;
    for (int i = 0; i < nv; i++) if (idx[i]==0) idx[i] = row++;

    std::vector<Eigen::Triplet<Scalar> > coef;
    coef.reserve( (nv-B.rows())*6 );
    Eigen::Matrix<Scalar,Eigen::Dynamic,2> b = Eigen::Matrix<Scalar,Eigen::Dynamic,2>::Zero(nv-B.rows(),2);
    row = 0;
    for (int i = 0; i < nv; i++)
        if (idx[i] >= 0)
        {
            coef.push_back( Eigen::Triplet<Scalar>(row,row,(Scalar)(VVi[i+1]-VVi[i])) );
            for (int j = VVi[i]; j < VVi[i+1]; j++)
            {
                if (idx[VV[j]] >= 0)
                    coef.push_back(Eigen::Triplet<Scalar>(row,idx[VV[j]],-1.0));
                else
                {
                    // boundary
                    b(row,0) += X(VV[j],0);
                    b(row,1) += X(VV[j],1);
                }
            }
            row++;
        }
    Eigen::SparseMatrix<Scalar> A(row,row);
    A.setFromTriplets(coef.begin(), coef.end());

    Eigen::Matrix<Scalar,Eigen::Dynamic,2> x;
    if constexpr(std::is_same_v<Scalar, double>)
    {
      Eigen::CholmodSupernodalLLT<Eigen::SparseMatrix<Scalar> > solver;
      solver.compute(A);
      x = solver.solve(b);
    } else {
      Eigen::SimplicialLLT<Eigen::SparseMatrix<Scalar> > solver;
      solver.compute(A);
      x = solver.solve(b);
    }

    row = 0;
    for (int i = 0; i < nv; i++)
        if (idx[i] >= 0)
            X.row(i) = x.row(row++);
}


template void tutte<Eigen::Matrix<double,-1,2>>(const Eigen::Ref<const Eigen::VectorXi> VV,
           const Eigen::Ref<const Eigen::VectorXi> VVi,
           const Eigen::Ref<const Eigen::VectorXi> B,
           Eigen::MatrixBase<Eigen::Matrix<double,-1,2>> & X, int boundary_type=0);

template void tutte<Eigen::Matrix<float,-1,2>>(const Eigen::Ref<const Eigen::VectorXi> VV,
          const Eigen::Ref<const Eigen::VectorXi> VVi,
          const Eigen::Ref<const Eigen::VectorXi> B,
          Eigen::MatrixBase<Eigen::Matrix<float,-1,2>> & X, int boundary_type=0);

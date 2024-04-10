#include "gauss_seidel.hpp"
#include <iostream>

// Code taken from Derek adopeted from Hendrik
// https://github.com/HTDerekLiu/surface_multigrid_code/
// https://github.com/HTDerekLiu/surface_multigrid_code/blob/adc49abe10a97200a3a5735bd370e98560091a37/src/mg_VCycle.cpp#L114


template <typename T, typename DerivedB, typename DerivedX>
bool gauss_seidel_compute(const Eigen::SparseMatrix<T> &A, const Eigen::MatrixBase<DerivedB> &B, const int maxIter,
                      Eigen::PlainObjectBase<DerivedX> &X) {
    Eigen::Matrix<typename DerivedX::Scalar, 1, DerivedX::ColsAtCompileTime> sum(B.cols());
    for (int iter = 0; iter < maxIter; ++iter) {
        // currently works only for symmetric matrices
        for (int oIdx = 0; oIdx < A.outerSize(); ++oIdx) {
            sum.setZero();
            for (typename Eigen::SparseMatrix<T>::InnerIterator it(A, oIdx); it; ++it) {
                if (it.index() != it.outer()) {
                    sum += it.value() * X.row(it.index());
                }
            }
            X.row(oIdx) = (B.row(oIdx) - sum) / A.coeff(oIdx, oIdx);
        }
    }

    return true;
}


template <typename T, typename DerivedB, typename DerivedX>
bool gauss_seidel_tillconverges(const Eigen::SparseMatrix<T> &A, const Eigen::MatrixBase<DerivedB> &B,
                      Eigen::PlainObjectBase<DerivedX> &X) {
    Eigen::Matrix<typename DerivedX::Scalar, 1, DerivedX::ColsAtCompileTime> sum(B.cols());
    Eigen::Matrix<typename DerivedX::Scalar, DerivedX::RowsAtCompileTime, DerivedX::ColsAtCompileTime> X_old = X;
    int iter = 0;
    double tol = 1e-6;
    double rel_tol = 1e-4;
    do{
        X_old = X;
        // currently works only for symmetric matrices
        for (int oIdx = 0; oIdx < A.outerSize(); ++oIdx) {
            sum.setZero();
            for (typename Eigen::SparseMatrix<T>::InnerIterator it(A, oIdx); it; ++it) {
                if (it.index() != it.outer()) {
                    sum += it.value() * X.row(it.index());
                }
            }
            X.row(oIdx) = (B.row(oIdx) - sum) / A.coeff(oIdx, oIdx);
        }
        ++iter;
    }while(iter < 100 && (X-X_old).norm() > tol && (X-X_old).norm()/X_old.norm() > rel_tol);
    std::cout << std::endl;
    std::cout << "gs converged in " <<iter << " iterations" << std::endl;
    return true;
}


template <typename T, typename DerivedB, typename DerivedX, typename TMC>
bool gauss_seidel_mc(const Eigen::SparseMatrix<T> &A, const Eigen::MatrixBase<DerivedB> &B, const int maxIter,
                     const Eigen::SparseMatrix<TMC> &MC, Eigen::PlainObjectBase<DerivedX> &X) {
    const typename Eigen::SparseMatrix<TMC>::StorageIndex *outer = MC.outerIndexPtr();
    const typename Eigen::SparseMatrix<TMC>::StorageIndex *inner = MC.innerIndexPtr();

    const typename Eigen::SparseMatrix<TMC>::StorageIndex *innerStart;
    int innerSize;
    Eigen::Index idx;

    Eigen::Matrix<typename DerivedX::Scalar, 1, DerivedX::ColsAtCompileTime> sum(B.cols());
    for (int iter = 0; iter < maxIter; ++iter) {
        // currently works only for symmetrix matrices
        for (int mcIdx = 0; mcIdx < MC.outerSize(); ++mcIdx) {
            innerStart = inner + *(outer + mcIdx);
            innerSize = (mcIdx < MC.outerSize() - 1 ? *(outer + mcIdx + 1) : MC.nonZeros()) - *(outer + mcIdx);
#pragma omp parallel for shared(innerSize, innerStart, A, X, B) firstprivate(sum) private(idx) default(none)
            for (int i = 0; i < innerSize; ++i) {
                idx = *(innerStart + i);
                sum.setZero();
                for (typename Eigen::SparseMatrix<T>::InnerIterator it(A, idx); it; ++it) {
                    if (it.index() != it.outer()) {
                        sum += it.value() * X.row(it.index());
                    }
                }
                X.row(idx) = (B.row(idx) - sum) / A.coeff(idx, idx);
            }
        }
    }
    return true;
}


template <typename T, typename TMC>
int coloring(Eigen::SparseMatrix<T>& A, Eigen::SparseMatrix<TMC>& MC) {
    const auto m = A.rows();
    assert(A.cols() == A.rows() && "A should be square");

    std::vector<Eigen::Triplet<TMC>> IJV;
    IJV.reserve(m);

    // m  means not yet visited
    Eigen::Vector<Eigen::Index, Eigen::Dynamic> C = Eigen::Vector<Eigen::Index, Eigen::Dynamic>::Constant(m, m);

    std::vector<bool> avail(0);

    for (Eigen::Index f = 0; f < m; f++) {
        for (typename Eigen::SparseMatrix<T>::InnerIterator it(A, f); it; ++it) {
            Eigen::Index c = C(it.index());
            if (c < m) avail[c] = false;
        }

        C(f) = std::distance(avail.begin(), std::find(avail.begin(), avail.end(), true));
        IJV.emplace_back(f, C(f), 1);

        if (C(f) == avail.size()) {
            avail.resize(C(f) + 1);
        }

        std::fill(avail.begin(), avail.end(), true);
    }

    MC.resize(m, avail.size());
    MC.setFromTriplets(IJV.begin(), IJV.end());
    MC.makeCompressed();

    return avail.size();
}

void gauss_seidel(const Eigen::SparseMatrix<double> &A, const Eigen::Matrix<double, -1, -1> &B, Eigen::Matrix<double, -1, -1> &X){
  gauss_seidel_tillconverges(A, B, X);
}

void gauss_seidel(const Eigen::SparseMatrix<double> &A, const Eigen::Matrix<double, -1, -1> &B, const int maxIter, Eigen::Matrix<double, -1, -1> &X){
  gauss_seidel_compute(A, B, maxIter, X);
}


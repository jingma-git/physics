#pragma once
#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace egl
{
    template <typename DerivedF, typename T>
    inline void adjacency_matrix(const Eigen::MatrixBase<DerivedF> &F,
                                 Eigen::SparseMatrix<T> &A)
    {
        using namespace std;
        using namespace Eigen;
        typedef typename DerivedF::Scalar Index;

        vector<Triplet<T>> ijv;
        ijv.reserve(F.size() * 2);
        for (size_t i = 0; i < F.cols(); ++i) // element
        {
            for (size_t j = 0; j < F.rows(); ++j)
                for (size_t k = j + 1; k < F.rows(); ++k)
                {
                    Index s = F(j, i);
                    Index d = F(k, i);
                    ijv.emplace_back(s, d, 1);
                    ijv.emplace_back(d, s, 1);
                }
        }

        const Index n = F.maxCoeff() + 1;
        A.resize(n, n);
        const Index capacity = (F.rows() == 3 ? 6 * n : 26 * n);
        A.reserve(capacity);
        A.setFromTriplets(ijv.begin(), ijv.end());

        for (size_t k = 0; k < A.outerSize(); ++k)
            for (typename SparseMatrix<T>::InnerIterator it(A, k); it; ++it)
                A.coeffRef(it.row(), it.col()) = 1;
    }
}
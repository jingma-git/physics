#pragma once
#include "adjacency_matrix.h"

namespace egl
{
    template <typename DerivedF, typename DerivedE>
    void extract_edges(Eigen::MatrixBase<DerivedF> &F,
                       Eigen::PlainObjectBase<DerivedE> &E)
    {
        typedef typename DerivedF::Scalar Index;
        using namespace Eigen;

        SparseMatrix<Index> A;
        egl::adjacency_matrix(F, A);
        assert(A.nonZeros() % 2 == 0);

        E.resize(A.nonZeros() / 2, 2);
        typename DerivedE::Scalar i = 0;
        for (size_t k = 0; k < A.outerSize(); ++k)
            for (typename SparseMatrix<Index>::InnerIterator it(A, k); it; ++it)
            {
                if (it.row() < it.col())
                {
                    E(i, 0) = it.row();
                    E(i, 1) = it.col();
                    ++i;
                }
            }
    }
}
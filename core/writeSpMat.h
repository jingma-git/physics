#pragma once
#include <Eigen/Sparse>
#include <fstream>

namespace egl
{
    template <typename T>
    int writeSpMat(const char *path, const Eigen::SparseMatrix<T> &mat)
    {
        //-> only for col major sparse matrix
        Eigen::SparseMatrix<T> mat_tm = mat;
        if (!mat_tm.isCompressed())
            mat_tm.makeCompressed();

        std::ofstream ofs(path, std::ios::binary);
        if (ofs.fail())
            return __LINE__;

        const Eigen::Index rows = mat_tm.rows(), cols = mat_tm.cols(), nnz = mat_tm.nonZeros();
        ofs.write((const char *)&rows, sizeof(Eigen::Index));
        ofs.write((const char *)&cols, sizeof(Eigen::Index));
        ofs.write((const char *)&nnz, sizeof(Eigen::Index));
        //-> the default SparseMatrix::StorageIndex is int
        ofs.write((const char *)mat_tm.outerIndexPtr(), (cols + 1) * sizeof(int));
        ofs.write((const char *)mat_tm.innerIndexPtr(), nnz * sizeof(int));
        ofs.write((const char *)mat_tm.valuePtr(), nnz * sizeof(T));
        ofs.close();

        return 0;
    }
}

#include "util.h"

namespace egl
{
    void polar(const Eigen::Matrix3d &F, Eigen::Matrix3d &Q)
    {
        Eigen::JacobiSVD<Eigen::Matrix3d> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
        Eigen::Matrix3d U = svd.matrixU();
        Eigen::Matrix3d V = svd.matrixV();
        Q = U * V.transpose();
    }
}
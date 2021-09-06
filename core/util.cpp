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

    Eigen::Vector3d cross(const Eigen::Vector3d &v1, const Eigen::Vector3d &v2)
    {
        return v1.cross(v2);
    }
}
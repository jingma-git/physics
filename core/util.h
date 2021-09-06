#pragma once

#include <Eigen/Eigen>

namespace egl
{
    void polar(const Eigen::Matrix3d &F, Eigen::Matrix3d &Q);
    Eigen::Vector3d cross(const Eigen::Vector3d &v1, const Eigen::Vector3d &v2);

    template <typename FT, typename IT>
    void ele_verts(const Eigen::MatrixBase<FT> &nods,
                   const Eigen::MatrixBase<IT> &tet,
                   Eigen::MatrixBase<FT> &eleV)
    {
        eleV.resize(3, 4);
        eleV.col(0) = nods.col(tet(0));
        eleV.col(1) = nods.col(tet(1));
        eleV.col(2) = nods.col(tet(2));
        eleV.col(3) = nods.col(tet(3));
    }

    template <typename FT, typename IT>
    void ele_coord(const Eigen::MatrixBase<FT> &nods,
                   const Eigen::MatrixBase<IT> &tet,
                   Eigen::MatrixBase<FT> &dm)
    {
        dm.resize(3, 3);
        dm.col(0) = nods.col(tet(1)) - nods.col(tet(0));
        dm.col(1) = nods.col(tet(2)) - nods.col(tet(0));
        dm.col(2) = nods.col(tet(3)) - nods.col(tet(0));
    }

}
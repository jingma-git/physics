#pragma once

#include <Eigen/Eigen>

namespace egl
{
    void polar(const Eigen::Matrix3d &F, Eigen::Matrix3d &Q);
    inline Eigen::Vector3d cross(const Eigen::Vector3d &v1, const Eigen::Vector3d &v2) { return v1.cross(v2); }

    template <typename FT, typename IT>
    void tri_verts(const Eigen::MatrixBase<FT> &nods,
                   const Eigen::MatrixBase<IT> &tri,
                   Eigen::MatrixBase<FT> &eleV)
    {
        eleV.col(0) = nods.col(tri(0));
        eleV.col(1) = nods.col(tri(1));
        eleV.col(2) = nods.col(tri(2));
    }

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

    template <typename DerivedV, typename DerivedF, typename DerivedN>
    void face_normal(const Eigen::MatrixBase<DerivedV> &nods,
                     const Eigen::MatrixBase<DerivedF> &tris,
                     Eigen::PlainObjectBase<DerivedN> &fn,
                     bool is_normalize = true)
    {
        using namespace Eigen;
        size_t m = tris.cols();
        fn.resize(3, m);
        for (size_t i = 0; i < m; ++i)
        {
            fn.col(i) = cross(nods.col(tris(1, i)) - nods.col(tris(0, i)), nods.col(tris(2, i)) - nods.col(tris(0, i)));
            if (is_normalize)
                fn.col(i).normalize();
        }
    }

    template <typename DerivedV>
    double calc_tri_area(const Eigen::MatrixBase<DerivedV> &vert)
    {
        return 0.5 * (cross(vert.col(1) - vert.col(0), vert.col(2) - vert.col(1))).norm();
    }
}
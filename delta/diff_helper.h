#pragma once

#include <Eigen/Eigen>

template <class T, int d>
struct diff_helper
{
    typedef Eigen::Matrix<T, d, d> MatT;
    typedef Eigen::Matrix<T, d, 1> VecT;

    MatT H;
    VecT A;
    VecT B;
    bool flip;

    // T_iikk = H_ik
    // T_ikik = (A_ik + B_ik)/2
    // T_ikki = (A_ik - B_ik)/2
    // if flip, B_ik=-B_ki; otherwise, B_ik=B_ki

    MatT apply(const MatT &M) const
    {
    }
};
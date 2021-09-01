#pragma once

#include <Eigen/Eigen>
#include <iostream>
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
        std::cout << "M\n"
                  << M << std::endl;
        MatT R;
        R.setZero();
        for (int i = 0; i < 3; i++)
            for (int k = 0; k < 3; k++)
                R(i, i) += H(i, k) * M(k, k);

        for (int i = 0; i < 3; i++)
        {
            int j = (i + 1) % 3;
            int k = (i + 2) % 3;
            T a = A[i] / 2, b = B[i] / 2, c = flip ? -b : b;
            R(j, k) += (a + b) * M(j, k);
            R(j, k) += (a - b) * M(k, j);
            R(k, j) += (a + c) * M(k, j);
            R(k, j) += (a - c) * M(j, k);
        }

        return R;
    }
};
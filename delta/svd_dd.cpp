#include "diff_helper.h"
#include <iostream>
#include <random>
#include <chrono>
using namespace std;
using namespace Eigen;

typedef mt19937 RandomT;
const double eps = 1e-6;

template <class T, int m, int n>
void fill_random(RandomT &gen,
                 Matrix<T, m, n> &mat,
                 T lo, T hi)
{
    std::uniform_real_distribution<> dist(lo, hi);
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            mat(i, j) = dist(gen);
}

Matrix3d cpm(const Vector3d &v)
{
    Matrix3d M;
    M(0, 0) = 0;
    M(0, 1) = -v[2];
    M(0, 2) = v[1];
    M(1, 0) = v[2];
    M(1, 1) = 0;
    M(1, 2) = -v[0];
    M(2, 0) = -v[1];
    M(2, 1) = v[0];
    M(2, 2) = 0;
    return M;
}

Matrix3d random_rotation(RandomT &gen, double max_angle)
{
    Vector3d u;
    fill_random(gen, u, -max_angle, max_angle);
    double t = u.norm();
    Matrix3d R = cpm(u / t);
    return (1 - cos(t)) * R * R + sin(t) * R + Matrix3d::Identity();
}

template <class T, int m, int n>
Matrix<T, m, m> outer(const Matrix<T, m, n> &u, const Matrix<T, m, n> &v)
{
    assert(m > 1 && n == 1 && "Must be a vector!");
    Matrix<T, m, m> A;
    assert(n == u.rows() && u.cols() == 1 && v.cols() == 1);
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < m; ++j)
            A(i, j) = u(i) * v(j);
    return A;
}

double compute_diag(const Vector3d &sig,
                    Vector3d &dE,
                    Matrix3d &ddE,
                    Vector3d &div)
{
    double a = sig.sum();
    double b = sig.dot(sig);
    double E = cos(a) + sin(b);
    cout << "E:" << E << endl;

    Vector3d da = Vector3d::Ones();
    Vector3d db = 2 * sig;
    dE = -sin(a) * da + cos(b) * db;
    cout << "db:" << db.transpose() << " dE:" << dE.transpose() << endl;

    Matrix3d ddb = 2 * Matrix3d::Identity();
    ddE = -cos(a) * outer(da, da) - sin(b) * outer(db, db) + cos(b) * ddb;
    cout << "ddE(sigma):\n"
         << ddE << endl;
    for (int i = 0; i < 3; ++i)
    {
        int j = (i + 1) % 3;
        int k = (i + 2) % 3;
        div(i) = (dE(j) - dE(k)) / (sig(j) - sig(k));
    }
    cout << "div=" << div.transpose() << endl;
    return E;
}

// ddE = Hessian : df
double compute(const Vector3d &sig,
               const Matrix3d &UU,
               const Matrix3d &VV,
               const Matrix3d &dF,
               Matrix3d &dE,
               Matrix3d &ddE)
{
    diff_helper<double, 3> h;
    Vector3d diag_dE;
    double E = compute_diag(sig, diag_dE, h.H, h.A);
    for (int i = 0; i < 3; ++i)
    {
        int j = (i + 1) % 3;
        int k = (i + 2) % 3;
        h.B(i) = (diag_dE(j) + diag_dE(k)) / (sig(j) + sig(k));
    }
    cout << "h.A=" << h.A.transpose() << endl;
    cout << "h.B=" << h.B.transpose() << endl;

    h.flip = false;
    dE = UU * diag_dE.asDiagonal() * VV.transpose(); // PK1 stress tensor
    cout << "dE:\n"
         << dE << endl;
    ddE = UU * h.apply(UU.transpose() * dF * VV) * VV.transpose();
    cout << "ddE\n"
         << ddE << endl;
    return E;
}

int main()
{
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    RandomT gen(seed);
    Matrix3d UU = random_rotation(gen, M_PI);
    Matrix3d VV = random_rotation(gen, M_PI);

    Vector3d Sig;
    fill_random(gen, Sig, -1., 1.);

    Matrix3d dUU = random_rotation(gen, eps);
    Matrix3d dVV = random_rotation(gen, eps);
    Matrix3d UU1 = UU * dUU;
    Matrix3d VV1 = VV * dUU;
    Vector3d dSig;
    fill_random(gen, dSig, -eps, eps);
    Vector3d Sig1 = Sig + dSig;

    Matrix3d F = UU * Sig.asDiagonal() * VV.transpose();
    Matrix3d F1 = UU1 * Sig1.asDiagonal() * VV1.transpose();
    Matrix3d dF = F1 - F;

    {
        cout << "UU" << endl;
        cout << UU << endl;
        cout << "Sig:" << Sig.transpose() << endl;
        cout << "F\n"
             << F << endl;

        cout << "UU1" << endl;
        cout << UU1 << endl;
        cout << "Sig1:" << Sig1.transpose() << endl;
        cout << "F1\n"
             << F1 << endl;
    }

    Matrix3d dE0, dE1;
    Matrix3d ddE0, ddE1;
    double E0 = compute(Sig, UU, VV, dF, dE0, ddE0);
    double E1 = compute(Sig1, UU1, VV1, dF, dE1, ddE1);
    double d0 = (E1 - E0) / eps;
    double d1 = ((dE0 + dE1).array() * dF.array()).sum() / (2 * eps);
    printf("%g %g %g\n", d0, d1, fabs(d0 - d1) / std::max(std::max(fabs(d0), fabs(d1)), 1e-30));

    Matrix3d u = dE1 - dE0;
    Matrix3d v = (ddE0 + ddE1) / 2;
    double e0 = u.norm() / eps;
    double e1 = v.norm() / eps;
    double e2 = (v - u).norm() / eps;
    printf("%g %g %g\n", e0, e1, e2 / std::max(std::max(e0, e1), 1e-30));
    return 0;
}
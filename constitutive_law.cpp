#include <Eigen/Eigen>
#include <iostream>
using namespace Eigen;
using namespace std;

extern "C"
{
    void tet_linear_(double *val, const double *x, const double *Dm, const double *vol, const double *lam, const double *miu);
    void tet_linear_jac_(double *jac, const double *x, const double *Dm, const double *vol, const double *lam, const double *miu);
    void tet_linear_hes_(double *hes, const double *x, const double *Dm, const double *vol, const double *lam, const double *miu);
}

int main()
{
    typedef Matrix<double, 3, 4> MatT;

    MatT Xini = MatT::Random();
    Matrix3d dm;
    dm.col(0) = Xini.col(1) - Xini.col(0);
    dm.col(1) = Xini.col(2) - Xini.col(0);
    dm.col(2) = Xini.col(3) - Xini.col(0);
    dm = dm.inverse().eval();

    double vol = 1.0, lam = 1.0, miu = 1.0;
    MatT x = MatT::Random();
    MatrixXd jac(12, 1);
    MatrixXd hes(12, 12);
    double val;
    tet_linear_(&val, x.data(), dm.data(), &vol, &lam, &miu);
    tet_linear_jac_(jac.data(), x.data(), dm.data(), &vol, &lam, &miu);
    tet_linear_hes_(hes.data(), x.data(), dm.data(), &vol, &lam, &miu);

    VectorXd X = Map<VectorXd>(x.data(), 12, 1);
    VectorXd X0 = Map<VectorXd>(Xini.data(), 12, 1);
    double err = (hes * (X - X0) - jac).norm();
    cout << "energy:" << val << endl;
    cout << "err:" << err << endl;
    return 0;
}
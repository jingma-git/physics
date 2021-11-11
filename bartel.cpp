#include <linear_tet_dphi_dX.h>

using namespace Eigen;
using namespace std;

int main(int argc, char *argv[])
{
    MatrixXd V(4, 3), X;
    V << 0, 0, 0,
        0.5, 0, 0,
        0, 0.5, 0,
        0, 0, 0.5;
    MatrixXi T(1, 4); // Tet
    T << 0, 1, 2, 3;

    Matrix43d dphi;
    sim::linear_tet_dphi_dX(dphi, V, T.row(0), X);
    cout << dphi << endl;
    return 0;
}
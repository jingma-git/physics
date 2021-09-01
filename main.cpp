#include <iostream>
#include <Eigen/Eigen>
using namespace std;
using namespace Eigen;

int main()
{
    Matrix3d F;
    F << 0, 1, 2,
        3, 4, 5,
        6, 7, 8;
    // cout << F.array() * F.array() << endl;
    // cout << F.transpose() * F << endl
    //      << endl;
    // cout << (F.array() * F.array()).sum() << endl;
    // cout << (F.transpose() * F).trace() << endl;
    // cout << sqrt((F.array() * F.array()).sum()) << endl;
    // cout << F.norm() << endl;
    return 0;
}
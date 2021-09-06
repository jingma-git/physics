#include <iostream>
#include <Eigen/Eigen>
using namespace std;
using namespace Eigen;

int main(int argc, char **argv)
{
    Matrix3d F;
    F << 0, 1, 2,
        3, 4, 5,
        6, 7, 8;
    int flag = atoi(argv[1]);
    cout << "flag=" << flag << endl;
    int a, b, c;
    switch (flag)
    {

    case 1:
        a = 1;
        b = 2;
        c = 3;
        break;
    case 2:
        a = 3;
        b = 4;
        c = 5;
        break;
    case 3:
        a = 7;
        b = 8;
        c = 9;
        break;
    default:
        cout << "no" << endl;
    }
    cout << a << " " << b << " " << c << endl;

    // cout << F.array() * F.array() << endl;
    // cout << F.transpose() * F << endl
    //      << endl;
    // cout << (F.array() * F.array()).sum() << endl;
    // cout << (F.transpose() * F).trace() << endl;
    // cout << sqrt((F.array() * F.array()).sum()) << endl;
    // cout << F.norm() << endl;
    return 0;
}
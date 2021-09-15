#include <iostream>
#include <Eigen/Eigen>
#include <bitset>
using namespace std;
using namespace Eigen;

int main(int argc, char **argv)
{
     // Matrix3d F;
     // F << 0, 1, 2,
     //     3, 4, 5,
     //     6, 7, 8;
     // int flag = atoi(argv[1]);
     // cout << "flag=" << flag << endl;
     // int a, b, c;
     // switch (flag)
     // {

     // case 1:
     //     a = 1;
     //     b = 2;
     //     c = 3;
     //     break;
     // case 2:
     //     a = 3;
     //     b = 4;
     //     c = 5;
     //     break;
     // case 3:
     //     a = 7;
     //     b = 8;
     //     c = 9;
     //     break;
     // default:
     //     cout << "no" << endl;
     // }
     // cout << a << " " << b << " " << c << endl;

     // VectorXd X(3);
     // X << 0, 1, 2;
     // std::swap(X(0), X(1));
     // cout << X << endl;
     // cout << X.cwiseSqrt() << endl;

     // cout << F.array() * F.array() << endl;
     // cout << F.transpose() * F << endl
     //      << endl;
     // cout << (F.array() * F.array()).sum() << endl;
     // cout << (F.transpose() * F).trace() << endl;
     // cout << sqrt((F.array() * F.array()).sum()) << endl;
     // cout << F.norm() << endl;

     // Eigen::PermutationMatrix<-1, -1> P;
     // P.resize(10);
     // P.setIdentity();
     // Eigen::MatrixXd p = P;
     // cout << "p\n"
     //      << p << endl;
     // P.indices()[1] = 0;
     // P.indices()[0] = 1;
     // Eigen::MatrixXd pp = P;
     // cout << "pp\n"
     // << pp << endl;
     // P.indices().reverseInPlace();
     // cout << P.indices() << endl;
     // Eigen::MatrixXd pi = P.inverse();
     // cout << "pi\n"
     //      << pi << endl;

     // typedef Eigen::SparseMatrix<int, Eigen::ColMajor, int> PatternType;
     // int nnz = 8;
     // int a[6] = {0, 2, 4, 5, 6, 8};
     // int b[8] = {1, 2, 0, 2, 4, 2, 1, 4};
     // int val[8] = {22, 7, 3, 5, 14, 1, 17, 8};
     // vector<int> outer_indices(a, a + 6);
     // vector<int> inner_indices(b, b + 8);
     // vector<int> data(val, val + 8);
     // PatternType s = Map<PatternType>(5, 5, nnz, &outer_indices[0], &inner_indices[0], &data[0]);
     // cout << s.toDense() << endl;
     // cout << "nonZero\n"
     //      << s.nonZeros() << endl;
     // cout << "iteration\n";
     // for (int k = 0; k < s.outerSize(); ++k)
     // {
     //      for (PatternType::InnerIterator it(s, k); it; ++it)
     //      {
     //           cout << it.row() << ", " << it.col() << ": " << it.value() << endl;
     //      }
     // }
     // cout << "---- selfAdjoint" << endl;
     // PatternType s_adj = s.selfadjointView<Eigen::Upper>();
     // cout << s_adj.toDense() << endl;
     // cout << "iteration2\n"
     //      << endl;
     // for (int j = 0; j < s.cols(); ++j)
     // {
     //      cout << "col" << j << endl;
     //      for (int iter = s.outerIndexPtr()[j]; iter < s.outerIndexPtr()[j + 1]; ++iter)
     //      {
     //           int i = s.innerIndexPtr()[iter];
     //           cout << iter << ":" << i << endl;
     //      }
     // }

     // SparseMatrix<int> tp;
     // vector<Triplet<int>> trips;
     // for (int i = 0; i < 10; ++i)
     // {
     //      trips.emplace_back(i, i, i);
     // }
     // tp.setFromTriplets(trips.begin(), trips.end());
     // cout << "after xP" << endl;
     // SparseMatrix<int> ttp = tp.twistedBy(P);
     // cout << ttp.toDense() << endl;

     std::bitset<64> val("0010");
     cout << val.none() << endl;
     // for (int i = 0; i < 64; ++i)
     //      val[i] = 1;
     cout << val << endl;
     return 0;
}
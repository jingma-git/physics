// LAPACK test code
//compile with: g++ main.cpp -llapack -lblas -o testprog

#include <iostream>
#include <vector>
#include <cblas.h>
using namespace std;

extern "C" void dgetrf_(int *dim1, int *dim2, double *a, int *lda, int *ipiv, int *info);
extern "C" void dgetrs_(char *TRANS, int *N, int *NRHS, double *A, int *LDA, int *IPIV, double *B, int *LDB, int *INFO);

void test_LU()
{
    char trans = 'N';
    int dim = 2;
    int nrhs = 1;
    int LDA = dim;
    int LDB = dim;
    int info;

    vector<double> a, b;

    a.push_back(1);
    a.push_back(1);
    a.push_back(1);
    a.push_back(-1);

    b.push_back(2);
    b.push_back(0);

    int ipiv[3];

    dgetrf_(&dim, &dim, &*a.begin(), &LDA, ipiv, &info);
    dgetrs_(&trans, &dim, &nrhs, &*a.begin(), &LDA, ipiv, &*b.begin(), &LDB, &info);

    std::cout << "solution is:";
    std::cout << "[" << b[0] << ", " << b[1] << ", "
              << "]" << std::endl;
    std::cout << "Info = " << info << std::endl;
}

void test_mm()
{
    double P[3][3] = {{0, 1, 0},
                      {1, 0, 0},
                      {0, 0, 1}};
    double A[3][3] = {{1, 2, 3},
                      {4, 5, 6},
                      {7, 8, 9}};
    double C[3][3];
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 3, 3, 3, 1, &P[0][0], 3, &A[0][0], 3, 0, &C[0][0], 3);
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
            cout << C[i][j] << " ";
        cout << endl;
    }
}

int main()
{
    // test_LU();
    test_mm();
    return (0);
}
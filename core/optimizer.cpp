#include "optimizer.h"
#include <Eigen/Sparse>
#include <iostream>

namespace egl
{
    void newton_solve(double *x, const size_t dof, const pfunc &f, const opt_args &args)
    {
        using namespace Eigen;
        using namespace std;
        assert(dof == f->Dof() && "[Error] dof does not match\n");

        SimplicialCholesky<SparseMatrix<double>> sol;
        sol.setMode(SimplicialCholeskyLDLT);

        Map<VectorXd> X(x, dof);
        VectorXd xstar = X;
        VectorXd dx(dof);

        for (size_t iter = 0; iter < args.max_iter; ++iter)
        {
            // cout << "........................iter" << iter << endl;
            double value = 0;
            {
                f->Val(&xstar[0], &value);
                if (iter % 100 == 0)
                    cout << "iter" << iter << " energy value:" << value << endl;
            }

            VectorXd grad = VectorXd::Zero(dof);
            {
                f->Gra(&xstar[0], &grad[0]);
                if (grad.norm() <= args.eps)
                    cout << "Gradient converged" << endl;
            }

            SparseMatrix<double> H(dof, dof);
            {
                vector<Triplet<double>> trips;
                f->Hes(&xstar[0], &trips);
                H.reserve(trips.size());
                H.setFromTriplets(trips.begin(), trips.end());
            }
            sol.compute(H);
            assert(sol.info() == Success && "Fail to factorize the system!");
            dx = -sol.solve(grad);
            assert(sol.info() == Success && "Fail to solve the system!");
            double xstar_norm = xstar.norm();
            xstar += dx;
            if (dx.norm() <= args.eps * xstar_norm)
            {
                cout << "Converged after " << iter << " iterations!" << endl;
                break;
            }
        }

        X = xstar;
    }
}
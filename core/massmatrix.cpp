#include "massmatrix.h"
#include "util.h"
#include <iostream>
#include <vector>
#include <unsupported/Eigen/KroneckerProduct>

namespace egl
{
    int calc_mass_matrix(const matd_t &nods,
                         const mati_t &cell,
                         const double rho,
                         const size_t dim,
                         spmat_t &M,
                         bool lumped)
    {
        using namespace Eigen;
        using namespace std;

        const size_t SIMPLEX = cell.rows();
        std::vector<Triplet<double>> trips;
        for (size_t i = 0; i < cell.cols(); ++i)
        {
            double coeff = 0;
            switch (SIMPLEX)
            {
            case 2:
            {
                matd_t edge = nods.col(cell(0, i)) - nods.col(cell(1, i));
                double length = edge.norm();
                coeff = rho * length / 6.0;
                break;
            }
            case 3:
            {
                matd_t edge(nods.rows(), 2);
                edge.col(0) = nods.col(cell(1, i)) - nods.col(cell(0, i));
                edge.col(1) = nods.col(cell(2, i)) - nods.col(cell(0, i));
                double area = 0.5 * edge.block(0, 0, 2, 2).determinant();
                coeff = rho * area / 12.0;
                break;
            }
            case 4:
            {
                matd_t dm(3, 3);
                ele_coord(nods, cell.col(i), dm);
                double volume = fabs(dm.determinant()) / 6.0;
                coeff = rho * volume / 20.0;
                break;
            }
            default:
            {
                std::cerr << "[info] mesh type not supported\n";
                exit(EXIT_FAILURE);
            }
            }
            for (size_t p = 0; p < cell.rows(); ++p)
            {
                for (size_t q = p; q < cell.rows(); ++q)
                {
                    trips.push_back(Triplet<double>(cell(p, i), cell(q, i), coeff));
                    trips.push_back(Triplet<double>(cell(q, i), cell(p, i), coeff));
                }
            }
        }
        if (lumped)
        {
            for (size_t i = 0; i < trips.size(); ++i)
                trips[i] = Triplet<double>(trips[i].row(), trips[i].row(), trips[i].value());
        }
        M.resize(nods.cols(), nods.cols());
        M.reserve(trips.size());
        M.setFromTriplets(trips.begin(), trips.end());

        if (dim > 1)
        {
            SparseMatrix<double> I(dim, dim);
            I.setIdentity();
            SparseMatrix<double> Mtemp = kroneckerProduct(M, I);
            M = Mtemp;
        }
        return 0;
    }
}

#include <iostream>

#include <Eigen/Cholesky>
#include <boost/filesystem.hpp>

#include "core/def.h"
#include "core/io.h"
#include "core/energy.h"
#include "core/util.h"
#include "core/massmatrix.h"
#include "json/json.h"

using namespace std;
using namespace Eigen;
using namespace egl;

struct SimulationParameters
{
    int total_frame;
    size_t max_iter;
    double eps;
    std::string output_folder;

    double dt, density, ym, pr;
    double we, wg, wp; // weight for elastic, grativational, positional potential
};

void newton_solve(double *x, const size_t dof, const pfunc &f, const SimulationParameters &params);
bool readInputFile(const char *fname,
                   SimulationParameters &params,
                   matd_t &nods, mati_t &tets, vector<size_t> &fixed);

int main()
{
    SimulationParameters params;
    matd_t nods;
    mati_t tets;
    std::vector<size_t> fixed;

    readInputFile("inputs/simplefem.json", params, nods, tets, fixed);
    if (!boost::filesystem::exists(params.output_folder))
        boost::filesystem::create_directories(params.output_folder);

    vector<pfunc> ebf(5);
    pfunc energy;
    ebf[0] = std::make_shared<momentum_potential_imp_euler>(nods, tets, params.density, params.dt, 1e0);
    ebf[1] = std::make_shared<elastic_potential>(nods, tets, elastic_potential::COROTATIONAL, params.ym, params.pr, params.we);
    ebf[2] = std::make_shared<gravitational_potential>(nods, tets, params.density, params.wg);
    ebf[3] = std::make_shared<positional_potential>(nods, params.wp);
    ebf[4] = std::make_shared<ext_force_energy>(nods, 1e0);
    energy = std::make_shared<energy_t<double>>(ebf);

    for (size_t ele : fixed)
        dynamic_pointer_cast<positional_potential>(ebf[3])->Pin(ele, &nods(0, ele));

    const vector<size_t> driver{149, 150, 151, 152,
                                153, 154, 155, 156, 157, 158, 159, 160,
                                161, 162, 163, 164, 165, 166, 167};
    const double intensity = 35;

    char outfile[256];
    for (int i = 0; i < params.total_frame; ++i)
    {
        cout << "[info] frame " << i << endl;
        sprintf(outfile, (params.output_folder + "simplefem%03d.vtk").c_str(), i);
        ofstream os(outfile);
        tet2vtk(os, nods.data(), nods.cols(), tets.data(), tets.cols());
        os.close();

        Vector3d n = cross((nods.col(167) - nods.col(166)), (nods.col(161) - nods.col(167)));
        Vector3d o = (nods.col(167) + nods.col(166) + nods.col(161)) / 3.0;

        if (i < 80)
        {
            for (size_t pi : driver)
            {
                Vector3d force = cross(n, nods.col(pi) - o);
                force = intensity * force / force.norm();
                dynamic_pointer_cast<ext_force_energy>(ebf[4])->ApplyForce(pi, force.data());
            }
        }

        if (i == 80)
        {
            for (size_t pi : driver)
            {
                dynamic_pointer_cast<ext_force_energy>(ebf[4])->RemoveForce(pi);
            }
        }

        newton_solve(nods.data(), nods.size(), energy, params);
        dynamic_pointer_cast<momentum_potential_imp_euler>(ebf[0])->Update(nods.data());
    }

    return 0;
}

void newton_solve(double *x, const size_t dof, const pfunc &f, const SimulationParameters &params)
{
    assert(dof == f->Dof() && "[Error] dof does not match\n");

    SimplicialCholesky<SparseMatrix<double>> sol;
    sol.setMode(SimplicialCholeskyLDLT);

    Map<VectorXd> X(x, dof);
    VectorXd xstar = X;
    VectorXd dx(dof);

    for (size_t iter = 0; iter < params.max_iter; ++iter)
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
            if (grad.norm() <= params.eps)
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
        if (dx.norm() <= params.eps * xstar_norm)
        {
            cout << "Converged after " << iter << " iterations!" << endl;
            break;
        }
    }

    X = xstar;
}

bool readInputFile(const char *fname,
                   SimulationParameters &params,
                   matd_t &nods, mati_t &tets, vector<size_t> &fixed)
{
    std::ifstream in(fname, std::ios::in);

    Json::Value root;
    Json::Reader jReader;

    if (!jReader.parse(in, root))
    {
        std::cout << "couldn't read input file: " << fname << '\n'
                  << jReader.getFormattedErrorMessages() << std::endl;
        exit(1);
    }

    params.dt = root.get("dt", 0.01).asDouble();
    params.total_frame = root.get("total_frame", 200).asInt();
    params.max_iter = root.get("max_iter", 10000).asInt();
    params.eps = root.get("eps", 1e-8).asDouble();
    params.density = root.get("density", 1.0).asDouble();
    params.ym = root.get("yong_modulus", 1e4).asDouble();
    params.pr = root.get("poisson_ratio", 0.45).asDouble();
    params.we = root.get("we", 1.0).asDouble();
    params.wg = root.get("wg", 1.0).asDouble();
    params.wp = root.get("wp", 1e3).asDouble();
    params.output_folder = root.get("output_folder", std::string("output/simplefem/")).asString();

    bool flag = tet_mesh_read_from_vtk(root.get("tets", "").asString().c_str(), nods, tets);
    flag = read_fixed_verts(root.get("fixed", "").asString().c_str(), fixed);
    return flag;
}
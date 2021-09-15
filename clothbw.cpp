#include <iostream>
#include <boost/filesystem.hpp>
#include "json/json.h"

#include "core/def.h"
#include "core/io.h"
#include "core/energy.h"
#include "core/util.h"
#include "core/readOBJ.h"
#include "core/edge_topology.h"
#include "core/optimizer.h"

using namespace std;
using namespace Eigen;
using namespace egl;

static opt_args optparam = {200, 1e-8, false};
struct SimulationParameters
{
    int total_frame;
    size_t max_iter;
    double eps;
    std::string output_folder;

    double dt, density;
    double wp, wg, ws, wb; // weight positional_potential, gravitational_potential, stretch, bending
};

bool readInputFile(const char *fname,
                   SimulationParameters &params,
                   matd_t &nods, mati_t &tris, vector<size_t> &fixed);

int main()
{
    SimulationParameters params;
    matd_t nods;
    mati_t tris;
    std::vector<size_t> fixed;

    readInputFile("inputs/clothsim.json", params, nods, tris, fixed);
    cout << "nods:" << nods.cols() << " tris:" << tris.cols() << endl;
    if (!boost::filesystem::exists(params.output_folder))
        boost::filesystem::create_directories(params.output_folder);

    mati_t edges, diams;
    {
        MeshTopology<matd_t, mati_t> mesh(nods, tris);
        mesh.get_edges(edges);
        mesh.get_diams(diams);
    }

    vector<pfunc> ebf(7);
    ebf[0] = std::make_shared<momentum_potential_imp_euler>(nods, tris, params.density, params.dt);
    ebf[1] = std::make_shared<positional_potential>(nods, params.wp);
    ebf[2] = std::make_shared<gravitational_potential>(nods, tris, params.density, params.wg);
    ebf[3] = std::make_shared<ext_force_energy>(nods, 1e0);

    ebf[4] = std::make_shared<bw98_stretch_energy>(nods, tris, params.ws);
    ebf[5] = std::make_shared<bw98_shear_energy>(nods, tris, 0.75 * params.ws);
    ebf[6] = std::make_shared<surf_bending_potential>(diams, nods, params.wb);
    pfunc energy = std::make_shared<energy_t<double>>(ebf);

    for (size_t ele : fixed)
        dynamic_pointer_cast<positional_potential>(ebf[1])->Pin(ele, &nods(0, ele));

    const double intense = 200;
    double f[3] = {-200, 0, -200};
    char outfile[256];
    for (size_t i = 0; i < params.total_frame; ++i)
    {
        cout << "[info] frame " << i << endl;
        sprintf(outfile, "%s/clothbw_%zu.vtk", params.output_folder.c_str(), i);
        ofstream os(outfile);
        tri2vtk(os, nods.data(), nods.cols(), tris.data(), tris.cols());
        os.close();

        //Apply force
        {
            if (i == 0)
                dynamic_pointer_cast<ext_force_energy>(ebf[3])->ApplyForce(3, f);
            if (i == 40)
                dynamic_pointer_cast<ext_force_energy>(ebf[3])->RemoveForce(3);
            // if (i == 160)
            //     dynamic_pointer_cast<positional_potential>(ebf[1])->Release(2);
        }

        newton_solve(nods.data(), nods.size(), energy, optparam);
        dynamic_pointer_cast<momentum_potential_imp_euler>(ebf[0])->Update(nods.data());
    }
}

bool readInputFile(const char *fname,
                   SimulationParameters &params,
                   matd_t &nods, mati_t &tris, vector<size_t> &fixed)
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
    params.ws = root.get("ws", 1.2e3).asDouble();
    params.wb = root.get("wb", 1e-3).asDouble();
    params.wp = root.get("wp", 1e3).asDouble();
    params.wg = root.get("wg", 1.0).asDouble();
    params.density = root.get("density", 1.0).asDouble();
    params.output_folder = root.get("output_folder", std::string("output/clothbw/")).asString();

    bool flag;
    flag = readOBJ(root.get("mesh", "").asString().c_str(), nods, tris);
    flag = read_fixed_verts(root.get("fixed", "").asString().c_str(), fixed);
    return flag;
}
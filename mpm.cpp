#include <Eigen/Eigen>
#include <random>
#include <iostream>
#include <fstream>
using namespace std;
using namespace Eigen;

const int n_particles = 8192;
const int n_grid = 128;
const double dx = 1.0 / double(n_grid);
const double inv_dx = 1 / dx;
const double dt = 2e-4;
const double p_vol = (dx * 0.5) * (dx * 0.5);
const double p_rho = 1.;
const double p_mass = p_vol * p_rho;
const double Ym = 400; // Yong's modulus

struct Particle
{
    Vector2d pos, vel;
    Matrix2d C; // ?
    double J = 1;
};
Vector2d grid_v[n_grid + 1][n_grid + 1]; // velocity
double grid_m[n_grid + 1][n_grid + 1];   // mass
vector<Particle> particles(n_particles);

void point2vtk(std::ofstream &os, const std::vector<Particle> &particles)
{
    os << "# vtk DataFile Version 2.0\nTRI\nASCII\n\nDATASET UNSTRUCTURED_GRID\n";

    os << "POINTS " << particles.size() << " float\n";
    for (size_t i = 0; i < particles.size(); ++i)
    {
        os << particles[i].pos.x() << " " << particles[i].pos.y() << " " << 0 << "\n";
    }

    os << "CELLS " << particles.size() << " " << particles.size() * 2 << "\n";
    for (size_t i = 0; i < particles.size(); ++i)
        os << 1 << " " << i << "\n";

    os << "CELL_TYPES " << particles.size() << "\n";
    for (size_t i = 0; i < particles.size(); ++i)
        os << 1 << "\n";
    os.close();
}

int main()
{
    std::mt19937 gen;
    std::uniform_real_distribution<double> dis(0.0, 1.0);

    for (int i = 0; i < n_particles; ++i)
    {
        double x = dis(gen) * 0.4 + 0.2, y = dis(gen) * 0.4 + 0.2;
        particles[i].pos = Vector2d(x, y);
        particles[i].vel = Vector2d(0, -1);
        particles[i].C.setZero();
    }

    ofstream os("output/mpm/rest.vtk");
    point2vtk(os, particles);

    char output_fname[256];
    for (int frame = 0; frame < 20000; ++frame)
    {
        sprintf(output_fname, "output/mpm/mpm-%04d.vtk", frame);
        ofstream os(output_fname);
        point2vtk(os, particles);
        cout << "frame" << frame << endl;

        for (int s = 0; s < 50; ++s)
        {
            memset(grid_v, 0, sizeof(grid_v));
            memset(grid_m, 0, sizeof(grid_m));
            // particle to grid
            for (Particle &p : particles)
            {
                Vector2d base = (((p.pos * inv_dx).array() - 0.5).cast<int>()).cast<double>();
                Vector2d fx = (p.pos * inv_dx).array() - base.array();
                int bi = base.x();
                int bj = base.y();
                // cout << "base:" << base.transpose() << " fx:" << fx.transpose()
                //      << " float:" << (p.pos * inv_dx).array().transpose() << endl;
                Matrix<double, 2, 3> w;
                w.col(0) = 0.5 * (1.5 - fx.array()) * (1.5 - fx.array());
                w.col(1) = 0.75 - (fx.array() - 1) * (fx.array() - 1);
                w.col(2) = 0.5 * (fx.array() - 0.5) * (fx.array() - 0.5);
                double stress = -dt * p_vol * (p.J - 1) * 4 * inv_dx * inv_dx * Ym; // why?
                Matrix2d affine;
                affine << stress + p_mass * p.C(0, 0), 0 + p_mass * p.C(0, 1),
                    0 + p_mass * p.C(1, 0), stress + p_mass * p.C(1, 1);
                for (int i = 0; i < 3; ++i)
                {
                    for (int j = 0; j < 3; ++j)
                    {
                        Vector2d dpos((i - fx(0)) * dx, (j - fx(1)) * dx);
                        double weight = w(0, i) * w(1, j);
                        grid_v[bi + i][bj + j].array() += weight * ((p_mass * p.vel).array() + (affine * dpos).array());
                        grid_m[bi + i][bj + j] += weight * p_mass;
                    }
                }
            }
            // constraint
            for (int i = 0; i < n_grid + 1; ++i)
                for (int j = 0; j < n_grid + 1; ++j)
                {
                    if (grid_m[i][j] > 0)
                    {
                        int bound = 3;
                        double inv_m = 1 / grid_m[i][j];
                        grid_v[i][j] *= inv_m;
                        grid_v[i][j].y() -= dt * 9.8;
                        if (i < bound && grid_v[i][j].x() < 0)
                            grid_v[i][j].x() = 0;
                        if (i > n_grid - bound && grid_v[i][j].x() > 0)
                            grid_v[i][j].x() = 0;
                        if (j < bound && grid_v[i][j].y() < 0)
                            grid_v[i][j].y() = 0;
                        if (i > n_grid - bound && grid_v[i][j].y() > 0)
                            grid_v[i][j].y() = 0;
                    }
                }
            // grid to particle
            for (Particle &p : particles)
            {
                Vector2d base = (((p.pos * inv_dx).array() - 0.5).cast<int>()).cast<double>();
                Vector2d fx = (p.pos * inv_dx).array() - base.array();
                int bi = base.x();
                int bj = base.y();
                Matrix<double, 2, 3> w;
                w.col(0) = 0.5 * (1.5 - fx.array()) * (1.5 - fx.array());
                w.col(1) = 0.75 - (fx.array() - 1) * (fx.array() - 1);
                w.col(2) = 0.5 * (fx.array() - 0.5) * (fx.array() - 0.5);

                Vector2d newV = Vector2d::Zero();
                Matrix2d newC = Matrix2d::Zero();
                for (int i = 0; i < 3; ++i)
                    for (int j = 0; j < 3; ++j)
                    {
                        Vector2d dpos((i - fx(0)) * dx, (j - fx(1)) * dx);
                        double weight = w(0, i) * w(1, j);
                        Vector2d g_v = grid_v[bi + i][bj + j];
                        newV.array() += weight * g_v.array();
                        newC.array() += (4 * weight * (g_v * dpos.transpose()) * inv_dx).array(); // why?
                    }

                p.vel = newV;
                p.pos.array() += dt * p.vel.array();
                p.J *= 1 + dt * newC.trace();
                p.C = newC;
            }
        }
    }
    cout << "done!" << endl;
    return 0;
}

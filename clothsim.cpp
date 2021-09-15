#include <core/meshgen.h>
#include <core/extract_edges.h>
#include <core/scalar_op.h>
#include <core/io.h>
#include <Eigen/Geometry>

using namespace std;
using namespace Eigen;
using namespace egl;

typedef double Float;
typedef size_t Int;
typedef Matrix<Float, -1, -1> matd_t;
typedef Matrix<Int, -1, -1> mati_t;

const Int N = 101; // number vertices in one row
const double dt = 8e-3;
const double stiff = 20;
const double damp = 2.6;
const double gravity = 0.8;

matd_t nods, v, a; // position, velocity, accelaration
mati_t eles;
mati_t edges;
vector<double> restl;

Vector3d ballBoundReflect(const Vector3d &pos, const Vector3d &vel, const Vector3d &center, const double radius)
{
    Vector3d ret = vel;
    double dis2surf = (pos - center).norm() - radius;
    if (dis2surf <= 0)
    {
        Vector3d normal = pos - center;
        double projn = vel.dot(normal) - 6 * egl::hermite_interp(dis2surf, 0, -0.1);
        if (projn < 0)
            ret -= projn * normal;
    }
    return vel;
}

void collide(matd_t &x, matd_t &vel)
{
    for (int i = 0; i < vel.cols(); ++i)
        vel(1, i) -= gravity * dt;
    for (int i = 0; i < x.cols(); ++i)
        vel.col(i) = ballBoundReflect(x.col(i), v.col(i), Vector3d(0, 0.2, -0.0), 0.4);
}

void update_pos(matd_t &x, matd_t &vel, matd_t &acc)
{
    for (size_t i = 0; i < x.cols(); ++i)
    {
        v.col(i) *= std::exp(dt * -damp);
        x.col(i) += dt * vel.col(i);
    }
}

void acc(const matd_t &x, matd_t &vel, matd_t &acc)
{
    acc.setZero();
    for (size_t e = 0; e < edges.rows(); ++e)
    {
        Int i = edges(e, 0);
        Int j = edges(e, 1);
        Vector3d dx = x.col(j) - x.col(i);
        double s = (dx.norm() - restl[e]) / (restl[e] * restl[e]); // why: E = k*x^2, f=-dE/(dx)
        acc.col(i) += s * dx;
        acc.col(j) -= s * dx; // why
    }
    for (size_t i = 0; i < nods.cols(); ++i)
    {
        vel.col(i) += stiff * acc.col(i) * dt;
        // vel.col(i) *= std::exp(-damp * dt);
    }
}

void explicit_euler()
{
    acc(nods, v, a);
    collide(nods, v);
    update_pos(nods, v, a);
}

int main()
{
    v = matd_t::Zero(3, N * N);
    a = matd_t::Zero(3, N * N);
    {
        generate_square_mesh("regular", N * N, 3, nods, eles);
        matd_t rot = AngleAxis<Float>(-M_PI / 2.0, Vector3d(1, 0, 0)).matrix();
        nods = (rot * nods).colwise() + Vector3d(0, 0.8, 0);
        std::ofstream os("output/clothsim/rest.vtk");
        egl::tri2vtk(os, nods.data(), nods.cols(), eles.data(), eles.cols());

        egl::extract_edges(eles, edges);
        restl.resize(edges.rows());
        for (int i = 0; i < edges.rows(); ++i)
            restl[i] = (nods.col(edges(i, 1)) - nods.col(edges(i, 0))).norm();
    }

    char outf[256];
    for (int frame = 0; frame < 2000; ++frame)
    {
        cout << "frame" << frame << endl;
        sprintf(outf, "output/clothsim/cloth-%d.vtk", frame);
        std::ofstream os(outf);
        egl::tri2vtk(os, nods.data(), nods.cols(), eles.data(), eles.cols());
        explicit_euler();
    }

    return 0;
}
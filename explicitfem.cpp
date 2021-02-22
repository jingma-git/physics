#include <iostream>
#include "explicitfem.h"

using namespace std;
using namespace Eigen;

void polar_decomp(const Eigen::Matrix3d &F, Eigen::Matrix3d &Q)
{
    Eigen::JacobiSVD<Matrix3d> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Q = svd.matrixU() * svd.matrixV().transpose();
}

void initialize(Object &obj, double density)
{
    for (std::vector<Particle>::iterator p = obj.particles.begin(); p != obj.particles.end(); p++)
    {
        p->mass = 0.0;
    }

    for (std::vector<Element>::iterator e = obj.elements.begin(); e != obj.elements.end(); e++)
    {
        Eigen::Vector3d &x0 = obj.particles[(*e)[0]].pos;
        Eigen::Vector3d &x1 = obj.particles[(*e)[1]].pos;
        Eigen::Vector3d &x2 = obj.particles[(*e)[2]].pos;
        Eigen::Vector3d &x3 = obj.particles[(*e)[3]].pos;

        e->normals[0] = (x2 - x1).cross(x3 - x1);
        e->normals[1] = (x3 - x0).cross(x2 - x0);
        e->normals[2] = (x1 - x0).cross(x3 - x0);
        e->normals[3] = (x2 - x0).cross(x1 - x0);

        e->basis << x1 - x0, x2 - x0, x3 - x0;
        double det = e->basis.determinant();
        double mass = density * det / 24.0; // density * volume
        e->basis = e->basis.inverse().eval();

        obj.particles[(*e)[0]].mass += mass;
        obj.particles[(*e)[1]].mass += mass;
        obj.particles[(*e)[2]].mass += mass;
        obj.particles[(*e)[3]].mass += mass;
    }
}

int main(int argc, char *argv[])
{
    char fname[80];
    SimulationParameters params;
    std::vector<Object> objects;
    readInputFile(argv[1], params, objects);

    for (std::vector<Object>::iterator obj = objects.begin(); obj != objects.end(); obj++)
    {
        initialize(*obj, params.density);
    }

    double time = 0;
    int frame = 0;
    double frameTime = -1.0;
    double dt = params.dt;
    while (time < params.total_time)
    {
        // write frame
        if (frameTime < 0)
        {
            for (unsigned int o = 0; o < objects.size(); o++)
            {
                sprintf(fname, params.output_fname.c_str(), o, frame);
                writeObj(fname, objects[o].particles, objects[o].triangles);
                frameTime = 1.0 / 30.0 - 0.0001;
                frame++;
            }
        }

        // take a time step
        for (vector<Object>::iterator obj = objects.begin(); obj != objects.end(); obj++)
        {
            // apply gravity
            for (std::vector<Particle>::iterator p = obj->particles.begin(); p != obj->particles.end(); p++)
            {
                p->frc.setZero();
                p->frc[2] = -9.8 * p->mass;
            }

            for (std::vector<Element>::iterator e = obj->elements.begin(); e != obj->elements.end(); e++)
            {
                Matrix3d X, F, Q, Ftide, strain, stress;
                Matrix3d Xdot, Fdot, Fdottide, strainrate, stressrate, I;
                I = Matrix3d::Identity();

                // deformation gradient, stress, strain
                Vector3d &x0 = obj->particles[(*e)[0]].pos;
                Vector3d &x1 = obj->particles[(*e)[1]].pos;
                Vector3d &x2 = obj->particles[(*e)[2]].pos;
                Vector3d &x3 = obj->particles[(*e)[3]].pos;

                X << x1 - x0, x2 - x0, x3 - x0;
                F = X * e->basis;
                polar_decomp(F, Q);
                Ftide = Q.transpose() * F;
                strain = 0.5 * (Ftide + Ftide.transpose()) - I;
                stress = params.lambda * strain.trace() * I + 2 * params.mu * strain;
                // compute Elatic Force
                // damping
                Vector3d v0 = obj->particles[(*e)[0]].vel;
                Vector3d v1 = obj->particles[(*e)[1]].vel;
                Vector3d v2 = obj->particles[(*e)[2]].vel;
                Vector3d v3 = obj->particles[(*e)[3]].vel;
                Xdot << v1 - v0, v2 - v0, v3 - v0;
                Fdot = Xdot * e->basis;
                Fdottide = Q.transpose() * F;
                // strainrate = 0.5 * (Fdot + Fdot.transpose()) - I;
                strainrate = 0.5 * (Fdottide + Fdottide.transpose()) - I;
                stressrate = params.damp * (params.lambda * strainrate.trace() * I + 2 * params.mu * strainrate);
                // assert(abs((Q - Q.transpose()).sum()) < 1e-8);
                obj->particles[(*e)[0]].frc += Q * (stress + stressrate) * e->normals[0] / 6.0;
                obj->particles[(*e)[1]].frc += Q * (stress + stressrate) * e->normals[1] / 6.0;
                obj->particles[(*e)[2]].frc += Q * (stress + stressrate) * e->normals[2] / 6.0;
                obj->particles[(*e)[3]].frc += Q * (stress + stressrate) * e->normals[3] / 6.0;
            }

            // time integration
            for (std::vector<Particle>::iterator p = obj->particles.begin(); p != obj->particles.end(); p++)
            {
                // update velocity and position
                p->vel += dt * (p->frc / p->mass);
                p->pos += dt * p->vel;
                // collision
                if (p->pos[2] < 0.0)
                {
                    p->pos[2] = 0.0;
                    p->vel[2] = 0.0;
                }
            }
        }

        time += params.dt;
        frameTime -= params.dt;
    }

    return 0;
}

#include "json/json.h"
#include <fstream>

bool readInputFile(const char *fname,
                   SimulationParameters &params,
                   std::vector<Object> &objects)
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

    params.dt = root.get("dt", 1.0 / 300.0).asDouble();
    params.total_time = root.get("total_time", 1.0).asDouble();
    params.density = root.get("density", 1e4).asDouble();
    params.lambda = root.get("lambda", 1e4).asDouble();
    params.mu = root.get("mu", 1e4).asDouble();
    params.damp = root.get("damping", 1e-2).asDouble();
    params.output_fname = root.get("output_fname", std::string("output-%02d.%04d.obj")).asString();

    Json::Value objectsIn = root["objects"];
    objects.resize(objectsIn.size());
    for (unsigned int i = 0; i < objectsIn.size(); i++)
    {
        readObject((objectsIn[i]["filename"]).asString().c_str(), objects[i]);
    }
    return true;
}

bool readObject(const char *fname, Object &object)
{
    char ch;
    Particle p;
    Element e;
    Tri t;

    std::ifstream in(fname, std::ios::in);
    while (in >> ch)
    {
        if (ch == 'p')
        {
            in >> p.pos[0] >> p.pos[1] >> p.pos[2] >> p.vel[0] >> p.vel[1] >> p.vel[2];
            object.particles.push_back(p);
            continue;
        }
        if (ch == 'e')
        {
            in >> e[0] >> e[1] >> e[2] >> e[3];
            object.elements.push_back(e);
            continue;
        }
        if (ch == 't')
        {
            in >> t[0] >> t[1] >> t[2];
            object.triangles.push_back(t);
            continue;
        }
    }
    in.close();
    std::cout << "inputfile " << fname << " read" << std::endl;
    return true;
}

void writeObj(char *fname, const std::vector<Particle> &meshPts, const std::vector<Tri> &triangles)
{
    std::cout << "writing " << fname << std::endl;
    std::ofstream out;
    std::vector<Particle>::const_iterator p;
    std::vector<Tri>::const_iterator t;

    out.open(fname);

    for (p = meshPts.begin(); p != meshPts.end(); p++)
        out << "v " << p->pos[0] << " " << p->pos[1] << " " << p->pos[2] << std::endl;

    for (t = triangles.begin(); t != triangles.end(); t++)
        out << "f " << (*t)[0] + 1 << " " << (*t)[1] + 1 << " " << (*t)[2] + 1 << std::endl;

    out.close();
}

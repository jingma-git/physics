#pragma once
// https://www.toptal.com/game/video-game-physics-part-iii-constrained-rigid-body-simulation
// http://www.cs.unc.edu/~lin/COMP768-F07/LEC/rbd2.pdf

#include <iostream>
#include <vector>
#include <string>
#include <Eigen/Dense>

struct Particle
{
    Eigen::Vector3d pos;
};

// record state and force
struct RigidBody
{
    double mass;
    Eigen::Matrix3d I0;          // Inertia
    Eigen::Vector3d pos, P, frc; // position, Linear Momentum, force
    Eigen::Quaterniond q;        // orientation
    Eigen::Vector3d L, trq;      // Angular Momentum, torque
};

class Tri
{
public:
    int indices[3];
    inline int &operator[](const unsigned int &i) { return indices[i]; };
    inline int operator[](const unsigned int &i) const { return indices[i]; };
};

struct SimulationParameters
{
    double dt, total_time, cor; // https://en.wikipedia.org/wiki/Coefficient_of_restitution
    std::string output_fname;
};

struct Object
{
    RigidBody rb;
    std::vector<Particle> particles;
    std::vector<Tri> triangles;
};

// I/O
#include "json/json.h"
#include <fstream>

bool readObject(const char *fname, Object &object);
bool readInputFile(const char *fname, SimulationParameters &params, std::vector<Object> &objects);
void writeObj(char *fname, const Object &obj);
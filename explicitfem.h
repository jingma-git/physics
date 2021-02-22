#pragma once
#include <iostream>
#include <vector>
#include <string>
#include <Eigen/Dense>

struct Particle
{
    double mass;
    Eigen::Vector3d pos, vel, frc;
};

class Element
{
    unsigned int indices[4];

public:
    Eigen::Matrix3d basis;
    Eigen::Vector3d normals[4]; // normals of opposite faces, for computing forces

    inline unsigned int &operator[](const unsigned int &i) { return indices[i]; }
    inline unsigned int operator[](const unsigned int &i) const { return indices[i]; }
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
    double dt, total_time, density, lambda, mu, damp;
    std::string output_fname;
};

struct Object
{
    std::vector<Particle> particles;
    std::vector<Element> elements;
    std::vector<Tri> triangles;
};

bool readObject(const char *fname, Object &object);
bool readInputFile(const char *fname, SimulationParameters &params, std::vector<Object> &objects);
void writeObj(char *fname, const std::vector<Particle> &meshPts, const std::vector<Tri> &triangles);
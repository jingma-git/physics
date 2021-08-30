#include <cmath>
#include <iostream>
using namespace std;

void compute(double x, double &z, double &dz, double &ddz)
{
    double a = cos(x);
    double b = 1. / a;
    double c = x * b;
    double d = 3. / 2. * c;
    double e = d * d - c;
    double f = sqrt(e);
    z = d + f;

    double da = -sin(x);
    double g = b * b;
    double db = -g * da;
    double dc = b + x * db;
    double dd = 3. / 2. * dc;
    double de = 2 * d * dd - dc;
    double h = 1. / f;
    double df = de * h / 2;
    dz = dd + df;

    double dda = -a;
    double dg = 2 * b * db;
    double ddb = -dg * da - g * dda;
    double ddc = 2 * db + x * ddb;
    double ddd = 3. / 2. * ddc;
    double dde = 2 * (dd * dd + d * ddd) - ddc;
    double dh = -h * h * df;
    double ddf = (dde * h + de * dh) / 2.;
    ddz = ddd + ddf;
}

int main()
{
    double x = 1.341, dx = 1e-6;

    double z0, dz0, ddz0;
    double z1, dz1, ddz1;

    compute(x - dx, z0, dz0, ddz0);
    compute(x + dx, z1, dz1, ddz1);

    double a = (z1 - z0) / dx;
    double b = dz0 + dz1;
    double err = fabs(a - b) / max(max(fabs(a), fabs(b)), 1e-30);
    printf("%g %g %g", a, b, err);
    double c = (dz1 - dz0) / dx;
    double d = ddz0 + ddz1;
    double err0 = fabs(c - d) / max(max(fabs(c), fabs(d)), 1e-30);
    return 0;
}
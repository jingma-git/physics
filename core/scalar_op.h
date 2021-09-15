#pragma once
namespace egl
{
    inline double clamp(double x, double xmin = 0, double xmax = 1)
    {
        return std::min(xmax, std::max(xmin, x));
    }

    inline double hermite_interp(double x, double a = 0, double b = 1)
    {
        double t = (x - a) / (b - a);
        return t * t * (3 - 2 * t);
    }
}
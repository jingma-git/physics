#pragma once
#include "def.h"

namespace egl
{
    struct opt_args
    {
        size_t max_iter;
        double eps;
        bool lineseach;
    };

    void newton_solve(double *x, const size_t dim, const pfunc &f, const opt_args &args);
}
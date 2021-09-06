#pragma once

#include "def.h"

namespace egl
{
    int calc_mass_matrix(const matd_t &nods,
                         const mati_t &cell,
                         const double rho,
                         const size_t dim,
                         spmat_t &M,
                         bool lumped);
}

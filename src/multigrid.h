#ifndef MULTIGRID_H
#define MULTIGRID_H

#include <vector>

std::vector<double>
restrictResidual(
    const std::vector<double>& fine,
    int nxf,
    int nyf);

#endif
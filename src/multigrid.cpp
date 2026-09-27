#include "multigrid.h"

#include <vector>

using namespace std;

int coarseIdx(int i,int j,int nx)
{
    return j*nx + i;
}

vector<double>
restrictResidual(
    const vector<double>& fine,
    int nxf,
    int nyf)
{
    int nxc = nxf/2;
    int nyc = nyf/2;

    vector<double> coarse(
        nxc*nyc,
        0.0);

    for(int jc=0;jc<nyc;jc++)
    {
        for(int ic=0;ic<nxc;ic++)
        {
            int ifine = 2*ic;
            int jfine = 2*jc;

            coarse[
                coarseIdx(ic,jc,nxc)
            ]
            =
            fine[
                coarseIdx(ifine,jfine,nxf)
            ];
        }
    }

    return coarse;
}
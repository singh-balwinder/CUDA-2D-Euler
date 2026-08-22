#include "definitions.h"
#include <cuda.h>

__global__ void timestep(device_edgeSoA deviceedge, device_cellSoA devicecell, int numEdges, int numIntcells, double CFL)
{
    int tid = blockIdx.x*blockDim.x + threadIdx.x;
    while(tid < numIntcells)
    {
        double const gamma = 1.4;

        double den = devicecell.den[tid];
		double xMom = devicecell.xMom[tid];
		double yMom = devicecell.yMom[tid];
		double energy = devicecell.energy[tid];

        double l1_nx = devicecell.l1_nx[tid];
        double l1_ny = devicecell.l1_ny[tid];
        double l1_len = devicecell.l1_len[tid];

        double l2_nx = devicecell.l2_nx[tid];
        double l2_ny = devicecell.l2_ny[tid];
        double l2_len = devicecell.l2_len[tid];

        double l3_nx = devicecell.l3_nx[tid];
        double l3_ny = devicecell.l3_ny[tid];
        double l3_len = devicecell.l3_len[tid];

        double l4_nx = devicecell.l4_nx[tid];
        double l4_ny = devicecell.l4_ny[tid];
        double l4_len = devicecell.l4_len[tid];

        double xVel = xMom/den;
		double yVel = yMom/den;
		double pressure = (gamma - 1.0)*(energy - ((xMom*xMom)/(2.0*den)) - ((yMom*yMom)/(2.0*den)));

        double soundspeed = sqrt(gamma*pressure/den);

        double I_nx = (l1_nx - l3_nx)/2.0;
        double I_ny = (l1_ny - l3_ny)/2.0;
        double aa = fabs(xVel*I_nx + yVel*I_ny) + soundspeed;

        double J_nx = (l2_nx - l4_nx)/2.0;
        double J_ny = (l2_ny - l4_ny)/2.0;
        double bb = fabs(xVel*J_nx + yVel*J_ny) + soundspeed;

        double cc = (l1_len + l3_len)/2.0;
        double dd = (l2_len + l4_len)/2.0;

        double cflden = (aa*cc) + (bb*dd);

        devicecell.dt[tid] = (CFL*devicecell.area[tid])/cflden;

        tid += gridDim.x*blockDim.x;
    }
}
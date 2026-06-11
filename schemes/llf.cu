#include "definitions.h"
#include "cuda.h"

__global__ void llf(device_edgeSoA deviceedge, device_cellSoA devicecell, int numEdges)
{
	int tid = blockIdx.x*blockDim.x + threadIdx.x;
	while (tid < numEdges)
	{
		double const gamma = 1.4;

		int parentcell = deviceedge.pCell[tid];
		int neighbourcell = deviceedge.nCell[tid];

		double nx = deviceedge.nx[tid];
		double ny = deviceedge.ny[tid];

		double length = deviceedge.length[tid];

		double parent_den = devicecell.den[parentcell];
		double neighbour_den = devicecell.den[neighbourcell];
		double parent_xMom = devicecell.xMom[parentcell];
		double neighbour_xMom = devicecell.xMom[neighbourcell];
		double parent_yMom = devicecell.yMom[parentcell];
		double neighbour_yMom = devicecell.yMom[neighbourcell];
		double parent_energy = devicecell.energy[parentcell];
		double neighbour_energy = devicecell.energy[neighbourcell];

		double parent_xVel = parent_xMom/parent_den;
		double parent_yVel = parent_yMom/parent_den;
		double parent_pressure = (gamma - 1.0)*(parent_energy - ((parent_xMom*parent_xMom)/(2.0*parent_den)) - ((parent_yMom*parent_yMom)/(2.0*parent_den)));

		double neighbour_xVel = neighbour_xMom/neighbour_den;
		double neighbour_yVel = neighbour_yMom/neighbour_den;
		double neighbour_pressure = (gamma - 1.0)*(neighbour_energy - (neighbour_xMom*neighbour_xMom)/(2.0*neighbour_den) - (neighbour_yMom*neighbour_yMom)/(2.0*neighbour_den));

		double norm_parent_den = parent_den;
		double norm_parent_xVel = (parent_xVel*nx) + (parent_yVel*ny);
		double norm_parent_yVel = (-parent_xVel*ny) + (parent_yVel*nx);
		double norm_parent_pressure = parent_pressure;

		double norm_neighbour_den = neighbour_den;
		double norm_neighbour_xVel = (neighbour_xVel*nx) + (neighbour_yVel*ny);
		double norm_neighbour_yVel = (-neighbour_xVel*ny) + (neighbour_yVel*nx);
		double norm_neighbour_pressure = neighbour_pressure;

		double norm_parent_xMom = norm_parent_den*norm_parent_xVel;
		double norm_parent_yMom = norm_parent_den*norm_parent_yVel;
		double norm_parent_energy = (norm_parent_pressure/(gamma - 1.0)) + (0.5*norm_parent_den*norm_parent_xVel*norm_parent_xVel) + (0.5*norm_parent_den*norm_parent_yVel*norm_parent_yVel);

		double norm_neighbour_xMom = norm_neighbour_den*norm_neighbour_xVel;
		double norm_neighbour_yMom = norm_neighbour_den*norm_neighbour_yVel;
		double norm_neighbour_energy = (norm_neighbour_pressure/(gamma - 1.0)) + (0.5*norm_neighbour_den*norm_neighbour_xVel*norm_neighbour_xVel) + (0.5*norm_neighbour_den*norm_neighbour_yVel*norm_neighbour_yVel);

		double norm_parent_ana_denFlux = norm_parent_den*norm_parent_xVel;
		double norm_parent_ana_xMomFlux = (norm_parent_den*norm_parent_xVel*norm_parent_xVel) + norm_parent_pressure;
		double norm_parent_ana_yMomFlux = (norm_parent_den*norm_parent_xVel*norm_parent_yVel);
		double norm_parent_ana_energyFlux = (norm_parent_pressure + norm_parent_energy)*norm_parent_xVel;

		double norm_neighbour_ana_denFlux = norm_neighbour_den*norm_neighbour_xVel;
		double norm_neighbour_ana_xMomFlux = (norm_neighbour_den*norm_neighbour_xVel*norm_neighbour_xVel) + norm_neighbour_pressure;
		double norm_neighbour_ana_yMomFlux = (norm_neighbour_den*norm_neighbour_xVel*norm_neighbour_yVel);
		double norm_neighbour_ana_energyFlux = (norm_neighbour_pressure + norm_neighbour_energy)*norm_neighbour_xVel;

		double norm_parent_soundspeed = sqrt(gamma*norm_parent_pressure/norm_parent_den);
		double norm_neighbour_soundspeed = sqrt(gamma*norm_neighbour_pressure/norm_neighbour_den);

		double norm_parent_eigenvalue = fabs(norm_parent_xVel) + norm_parent_soundspeed;
		double norm_neighbour_eigenvalue = fabs(norm_neighbour_xVel) + norm_neighbour_soundspeed;

		double alpha;

		if(norm_parent_eigenvalue >= norm_neighbour_eigenvalue)
		{
			alpha = norm_parent_eigenvalue;
		}
		else
		{
			alpha = norm_neighbour_eigenvalue;
		}

		double flux0 = ((0.5*(norm_parent_ana_denFlux + norm_neighbour_ana_denFlux)) - (0.5*alpha*(norm_neighbour_den - norm_parent_den)))*length;
		double flux1 = ((0.5*(norm_parent_ana_xMomFlux + norm_neighbour_ana_xMomFlux)) - (0.5*alpha*(norm_neighbour_xMom - norm_parent_xMom)))*length;
		double flux2 = ((0.5*(norm_parent_ana_yMomFlux + norm_neighbour_ana_yMomFlux)) - (0.5*alpha*(norm_neighbour_yMom - norm_parent_yMom)))*length;
		double flux3 = ((0.5*(norm_parent_ana_energyFlux + norm_neighbour_ana_energyFlux)) - (0.5*alpha*(norm_neighbour_energy - norm_parent_energy)))*length;

		deviceedge.denFlux[tid] = flux0;
		deviceedge.xMomFlux[tid] = (flux1*nx) - (flux2*ny);
		deviceedge.yMomFlux[tid] = (flux1*ny) + (flux2*nx);
		deviceedge.energyFlux[tid] = flux3;

		tid += gridDim.x*blockDim.x;
	}
}

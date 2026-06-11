#include "definitions.h"

void device_SoA(vector <edge> &E, vector <cell> &convol, device_edgeSoA &deviceedge, device_cellSoA &devicecell, host_edgeSoA &hostedge, host_cellSoA &hostcell)
{
    	int numEdges = E.size();

	cudaMalloc(&deviceedge.pCell, numEdges*sizeof(int));
	cudaMalloc(&deviceedge.nCell, numEdges*sizeof(int));

	cudaMalloc(&deviceedge.nx, numEdges*sizeof(double));
	cudaMalloc(&deviceedge.ny, numEdges*sizeof(double));

	cudaMalloc(&deviceedge.length, numEdges*sizeof(double));

	cudaMalloc(&deviceedge.denFlux, numEdges*sizeof(double));
	cudaMalloc(&deviceedge.xMomFlux, numEdges*sizeof(double));
	cudaMalloc(&deviceedge.yMomFlux, numEdges*sizeof(double));
	cudaMalloc(&deviceedge.energyFlux, numEdges*sizeof(double));

	int numCells = convol.size();

	cudaMalloc(&devicecell.area, numCells*sizeof(double));
	cudaMalloc(&devicecell.den, numCells*sizeof(double));
	cudaMalloc(&devicecell.xMom, numCells*sizeof(double));
	cudaMalloc(&devicecell.yMom, numCells*sizeof(double));
	cudaMalloc(&devicecell.energy, numCells*sizeof(double));

	cudaMemcpy(deviceedge.pCell, hostedge.pCell.data(), numEdges*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(deviceedge.nCell, hostedge.nCell.data(), numEdges*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(deviceedge.nx, hostedge.nx.data(), numEdges*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(deviceedge.ny, hostedge.ny.data(), numEdges*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(deviceedge.length, hostedge.length.data(), numEdges*sizeof(double), cudaMemcpyHostToDevice);

	cudaMemcpy(devicecell.area, hostcell.area.data(), numCells*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(devicecell.den, hostcell.den.data(), numCells*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(devicecell.xMom, hostcell.xMom.data(), numCells*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(devicecell.yMom, hostcell.yMom.data(), numCells*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(devicecell.energy, hostcell.energy.data(), numCells*sizeof(double), cudaMemcpyHostToDevice);
}

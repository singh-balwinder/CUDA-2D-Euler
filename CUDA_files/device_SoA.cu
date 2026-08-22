#include "definitions.h"

void device_SoA(vector <edge> &E, vector <cell> &convol, device_edgeSoA &deviceedge, device_cellSoA &devicecell, host_edgeSoA &hostedge, host_cellSoA &hostcell, int numIntcells)
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
	cudaMalloc(&devicecell.l1_nx, numIntcells*sizeof(double));
	cudaMalloc(&devicecell.l1_ny, numIntcells*sizeof(double));
	cudaMalloc(&devicecell.l1_len, numIntcells*sizeof(double));
	cudaMalloc(&devicecell.l2_nx, numIntcells*sizeof(double));
	cudaMalloc(&devicecell.l2_ny, numIntcells*sizeof(double));
	cudaMalloc(&devicecell.l2_len, numIntcells*sizeof(double));
	cudaMalloc(&devicecell.l3_nx, numIntcells*sizeof(double));
	cudaMalloc(&devicecell.l3_ny, numIntcells*sizeof(double));
	cudaMalloc(&devicecell.l3_len, numIntcells*sizeof(double));
	cudaMalloc(&devicecell.l4_nx, numIntcells*sizeof(double));
	cudaMalloc(&devicecell.l4_ny, numIntcells*sizeof(double));
	cudaMalloc(&devicecell.l4_len, numIntcells*sizeof(double));
	cudaMalloc(&devicecell.dt, numIntcells*sizeof(double));

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
	cudaMemcpy(devicecell.l1_nx, hostcell.l1_nx.data(), numIntcells*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(devicecell.l1_ny, hostcell.l1_ny.data(), numIntcells*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(devicecell.l1_len, hostcell.l1_len.data(), numIntcells*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(devicecell.l2_nx, hostcell.l2_nx.data(), numIntcells*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(devicecell.l2_ny, hostcell.l2_ny.data(), numIntcells*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(devicecell.l2_len, hostcell.l2_len.data(), numIntcells*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(devicecell.l3_nx, hostcell.l3_nx.data(), numIntcells*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(devicecell.l3_ny, hostcell.l3_ny.data(), numIntcells*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(devicecell.l3_len, hostcell.l3_len.data(), numIntcells*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(devicecell.l4_nx, hostcell.l4_nx.data(), numIntcells*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(devicecell.l4_ny, hostcell.l4_ny.data(), numIntcells*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(devicecell.l4_len, hostcell.l4_len.data(), numIntcells*sizeof(double), cudaMemcpyHostToDevice);
}
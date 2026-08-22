#include "definitions.h"

void cuda_free(device_edgeSoA &deviceedge, device_cellSoA &devicecell, host_edgeSoA &hostedge, host_cellSoA &hostcell)
{

	cudaFree(deviceedge.pCell);
	cudaFree(deviceedge.nCell);

	cudaFree(deviceedge.nx);
	cudaFree(deviceedge.ny);

	cudaFree(deviceedge.length);

	cudaFree(deviceedge.denFlux);
	cudaFree(deviceedge.xMomFlux);
	cudaFree(deviceedge.yMomFlux);
	cudaFree(deviceedge.energyFlux);

	cudaFree(devicecell.area);
	cudaFree(devicecell.den);
	cudaFree(devicecell.xMom);
	cudaFree(devicecell.yMom);
	cudaFree(devicecell.energy);
	cudaFree(devicecell.l1_nx);
	cudaFree(devicecell.l1_ny);
	cudaFree(devicecell.l1_len);
	cudaFree(devicecell.l2_nx);
	cudaFree(devicecell.l2_ny);
	cudaFree(devicecell.l2_len);
	cudaFree(devicecell.l3_nx);
	cudaFree(devicecell.l3_ny);
	cudaFree(devicecell.l3_len);
	cudaFree(devicecell.l4_nx);
	cudaFree(devicecell.l4_ny);
	cudaFree(devicecell.l4_len);
	cudaFree(devicecell.dt);
}
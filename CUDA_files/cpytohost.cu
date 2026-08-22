#include "definitions.h"

void cpytohost(device_edgeSoA &deviceedge, device_cellSoA &devicecell, host_edgeSoA &hostedge, host_cellSoA &hostcell, vector <edge> &E, int numIntcells, double &del_t)
{
	int numEdges = E.size();
	cudaMemcpy(hostedge.cart_flux0.data(), deviceedge.denFlux, numEdges*sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(hostedge.cart_flux1.data(), deviceedge.xMomFlux, numEdges*sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(hostedge.cart_flux2.data(), deviceedge.yMomFlux, numEdges*sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(hostedge.cart_flux3.data(), deviceedge.energyFlux, numEdges*sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(hostcell.dt.data(), devicecell.dt, numIntcells*sizeof(double), cudaMemcpyDeviceToHost);

	for(int i=0; i<numEdges; i++)
	{
		E[i].cartesian_interface_flux[0] = hostedge.cart_flux0[i];
		E[i].cartesian_interface_flux[1] = hostedge.cart_flux1[i];
		E[i].cartesian_interface_flux[2] = hostedge.cart_flux2[i];
		E[i].cartesian_interface_flux[3] = hostedge.cart_flux3[i];
	}

	del_t = hostcell.dt[0];
	for (int i=1; i<numIntcells; i++) {
    	if (hostcell.dt[i] <= del_t)
        	del_t = hostcell.dt[i];
	}
	
}

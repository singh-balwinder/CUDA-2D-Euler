#include "definitions.h"

void cpytohost(device_edgeSoA &deviceedge, host_edgeSoA &hostedge, vector <edge> &E)
{
	int numEdges = E.size();
	cudaMemcpy(hostedge.cart_flux0.data(), deviceedge.denFlux, numEdges*sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(hostedge.cart_flux1.data(), deviceedge.xMomFlux, numEdges*sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(hostedge.cart_flux2.data(), deviceedge.yMomFlux, numEdges*sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(hostedge.cart_flux3.data(), deviceedge.energyFlux, numEdges*sizeof(double), cudaMemcpyDeviceToHost);

	for(int i=0; i<numEdges; i++)
	{
		E[i].cartesian_interface_flux[0] = hostedge.cart_flux0[i];
		E[i].cartesian_interface_flux[1] = hostedge.cart_flux1[i];
		E[i].cartesian_interface_flux[2] = hostedge.cart_flux2[i];
		E[i].cartesian_interface_flux[3] = hostedge.cart_flux3[i];
	}
}

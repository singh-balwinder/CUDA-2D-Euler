#include "definitions.h"

int main()
{
    // variable declaration
    vector <point> p;
    vector <cell> convol;
    vector <edge> E;

    int numEdges, numIntcells, numCells;

    bool converged = false;
    double del_t, total_time = 0.0;
    double CFL = 0.4;
    int loopIter = 0;
    // double h = 0.0;

    // Pre-processing
    read_vtk(p, E, convol, numEdges, numIntcells);
    ghostcells(E, convol, numCells);

    // Initialization
    testcase1(convol);
    
    // update loop
    while (!converged)
    {
        loopIter = loopIter + 1;
        boundary(E, convol);

        host_edgeSoA hostedge(numEdges);
        host_cellSoA hostcell(numCells);
        host_SoA(E, convol, hostedge, hostcell);

        device_edgeSoA deviceedge;
        device_cellSoA devicecell;
        device_SoA(E, convol, deviceedge, devicecell, hostedge, hostcell);

        llf<<<1000, 256>>> (deviceedge, devicecell, numEdges);
	cudaDeviceSynchronize();
	cpytohost(deviceedge, hostedge, E);
        timestep(convol, E, del_t, CFL);
        // h = h + ((20/sqrt(3.0))*del_t);
        tempupdate(convol, E, del_t, converged);
        total_time=total_time+del_t;
        cout << total_time << endl;
        
        if (total_time >= 0.1561)
        {
            break;
        }
        
        /*
        if (loopIter >= 50000)
        {
            break;
        }
        */
    }

    string vtkFilename = "finalSolution.vtk";
    writevtk(vtkFilename, p, E, convol);
    cout << loopIter << endl;

    return 0;
}

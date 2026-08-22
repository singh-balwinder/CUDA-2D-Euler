#include <iostream>
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include "cuda.h"
#include "geometry.h"
#include "cell.h"

using namespace std;

// functions
void read_vtk(vector <point> &p, vector <edge> &E, vector <cell> &convol, int &numEdges, int &numIntcells);
void ghostcells(vector <edge> &E, vector <cell> &convol, int &numCells);
void testcase1(vector <cell> &convol);
void testcase2(vector <cell> &convol);
void boundary(vector <edge> &E, vector <cell> &convol);
// void boundary(vector <edge> &E, vector <cell> &convol, double &h);
void tempupdate(vector <cell> &convol, vector <edge> &E, double &del_t, bool &converged);
void writevtk(string &vtkFilename, vector <point> &p, vector <edge> &E, vector <cell> &convol);
void host_SoA(vector <edge> &E, vector <cell> &convol, host_edgeSoA &hostedge, host_cellSoA &hostcell, int numIntcells);
void device_SoA(vector <edge> &E, vector <cell> &convol, device_edgeSoA &deviceedge, device_cellSoA &devicecell, host_edgeSoA &hostedge, host_cellSoA &hostcell, int numIntcells);
void cpytohost(device_edgeSoA &deviceedge, device_cellSoA &devicecell, host_edgeSoA &hostedge, host_cellSoA &hostcell, vector <edge> &E, int numIntcells, double &del_t);
void cuda_free(device_edgeSoA &deviceedge, device_cellSoA &devicecell, host_edgeSoA &hostedge, host_cellSoA &hostcell);

#ifdef __CUDACC__
__global__ void llf(device_edgeSoA deviceedge, device_cellSoA devicecell, int numEdges);
__global__ void timestep(device_edgeSoA deviceedge, device_cellSoA devicecell, int numEdges, int numIntcells, double CFL);
#endif
void kfvs(vector <edge> &E, vector <cell> &convol);
void pvu(vector <edge> &E, vector <cell> &convol);
void pvuc(vector <edge> &E, vector <cell> &convol);
void mpvu(vector <edge> &E, vector <cell> &convol);
void tkfvs(vector <edge> &E, vector <cell> &convol);
void perthame(vector <edge> &E, vector <cell> &convol);
void pvuc2(vector <edge> &E, vector <cell> &convol);
void aravind(vector <edge> &E, vector <cell> &convol);
void blendedPVU(vector <edge> &E, vector <cell> &convol);
void roe(vector <edge> &E, vector <cell> &convol);
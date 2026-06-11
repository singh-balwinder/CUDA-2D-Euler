#include "definitions.h"

void host_SoA(vector <edge> &E, vector <cell> &convol, host_edgeSoA &hostedge, host_cellSoA &hostcell)
{
    int numEdges = E.size();

    for (int i=0; i<numEdges; i++)
    {
        hostedge.pCell[i]   = E[i].parentid;
        hostedge.nCell[i]   = E[i].neighbourid;
        hostedge.nx[i]      = E[i].nx;
        hostedge.ny[i]      = E[i].ny;
        hostedge.length[i]  = E[i].L.length;
    }

    int numCells = convol.size();

    for (int i=0; i<numCells; i++)
    {
        hostcell.area[i] = convol[i].area;
        hostcell.den[i] = convol[i].cons_var[0];
        hostcell.xMom[i] = convol[i].cons_var[1];
        hostcell.yMom[i] = convol[i].cons_var[2];
        hostcell.energy[i] = convol[i].cons_var[3];
    }
}
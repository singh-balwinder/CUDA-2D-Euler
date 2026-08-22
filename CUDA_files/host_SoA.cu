#include "definitions.h"

void host_SoA(vector <edge> &E, vector <cell> &convol, host_edgeSoA &hostedge, host_cellSoA &hostcell, int numIntcells)
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

    for (int i=0; i<numIntcells; i++)
    {
        hostcell.l1_nx[i] = convol[i].quad.l1.nx;
        hostcell.l1_ny[i] = convol[i].quad.l1.ny;
        hostcell.l1_len[i] = convol[i].quad.l1.getLength();

        hostcell.l2_nx[i] = convol[i].quad.l2.nx;
        hostcell.l2_ny[i] = convol[i].quad.l2.ny;
        hostcell.l2_len[i] = convol[i].quad.l2.getLength();

        hostcell.l3_nx[i] = convol[i].quad.l3.nx;
        hostcell.l3_ny[i] = convol[i].quad.l3.ny;
        hostcell.l3_len[i] = convol[i].quad.l3.getLength();

        hostcell.l4_nx[i] = convol[i].quad.l4.nx;
        hostcell.l4_ny[i] = convol[i].quad.l4.ny;
        hostcell.l4_len[i] = convol[i].quad.l4.getLength();
    }
}
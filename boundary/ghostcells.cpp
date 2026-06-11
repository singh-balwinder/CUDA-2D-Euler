#include "definitions.h"

void ghostcells(vector <edge> &E, vector <cell> &convol, int &numCells)
{
    int i;
    int numEdges = E.size();
    numCells = convol.size();

    for (i=0; i<numEdges; i++)
    {
        if (E[i].edgeType == "BOUNDARY")
        {
            convol.push_back(cell(numCells, "GHOST"));
            convol[numCells].area = convol[E[i].parentid].getArea();
            E[i].neighbourid = numCells;
            numCells = convol.size();
        }
    }

    // print statements
    /*
    cout << "total number of cells" << "  " << numCell << endl;

    for (i=0;i<numEdge;i++)
    {
        cout << i << "  " << E[i].parentid << "  " << E[i].neighbourid << endl;
    }
    */
    
    /*
    cout << convol[1].quad.l1.nx << "  " << convol[1].quad.l1.ny << endl;
    cout << convol[1].quad.l2.nx << "  " << convol[1].quad.l2.ny << endl;
    cout << convol[1].quad.l3.nx << "  " << convol[1].quad.l3.ny << endl;
    cout << convol[1].quad.l4.nx << "  " << convol[1].quad.l4.ny << endl;
    */
}
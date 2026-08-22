#include "definitions.h"

void writevtk(string &vtkFilename, vector <point> &p, vector <edge> &E, vector <cell> &convol)
{
    int i;
    int numPoints = p.size();
    int numEdge = E.size();
    int numGhostcell = 0;
    
    for (i=0;i<numEdge;i++)
    {
        if (E[i].edgeType == "BOUNDARY")
        {
            numGhostcell = numGhostcell + 1;
        }
        else
        {
            continue;
        }
    }

    // cout << numGhostcell << endl;

    int numCells = convol.size() - numGhostcell;

    ofstream vtkFile(vtkFilename);

    
    /* check to see .vtk file is open or not */
    if (!vtkFile.is_open())
    {
        cerr << "Failed to open the .vtk file" << endl;
        // return 1;
    }

    /* VTK file header */
    vtkFile << "# vtk DataFile Version 2.0" << endl;
    vtkFile << "oblique_shock, Created by Solver" << endl;
    vtkFile << "ASCII" << endl;
    vtkFile << "DATASET UNSTRUCTURED_GRID" << endl;

    vtkFile << endl;

    /* POINTS data */
    vtkFile << "POINTS" << " " << numPoints << " " << "double" << endl;

    for (i=0;i<numPoints;i++)
    {
        vtkFile << p[i].x << " " << p[i].y << " " << p[i].z << endl;
    }

    vtkFile << endl;

    vtkFile << "CELLS" << " " << numCells << " " << (numCells*5) << endl;
    for (i=0;i<numCells;i++)
    {
        vtkFile << convol[i].numPolyVertices << " " << convol[i].quad.p1.pointid << " " << convol[i].quad.p2.pointid << " " << convol[i].quad.p3.pointid << " " << convol[i].quad.p4.pointid << endl;
    }

    vtkFile << endl;

    vtkFile << "CELL_TYPES" << " " << numCells << endl;
    for (i=0;i<numCells;i++)
    {
        vtkFile << "9" << endl;
    }

    vtkFile << endl;

    vtkFile << "CELL_DATA" << " " << numCells << endl;
    vtkFile << "SCALARS" << " " << "DENSITY" << " " << "double" << " " << "1" << endl;
    vtkFile << "LOOKUP_TABLE" << " " << "default" << endl;

    for (i=0;i<numCells;i++)
    {
        vtkFile << convol[i].prim_var[0] << endl;
    }

    vtkFile << endl;

    vtkFile << "SCALARS" << " " << "PRESSURE" << " " << "double" << " " << "1" << endl;
    vtkFile << "LOOKUP_TABLE" << " " << "default" << endl;

    for (i=0;i<numCells;i++)
    {
        vtkFile << convol[i].prim_var[3] << endl;
    }

    vtkFile << endl;

    vtkFile << "VECTORS" << " " << "VELOCITY" << " " << "double" << endl;
    // vtkFile << "LOOKUP_TABLE" << " " << "default" << endl;

    for (i=0;i<numCells;i++)
    {
        vtkFile << convol[i].prim_var[1] << " " << convol[i].prim_var[2] << " " << 0 << endl;
    }

    vtkFile.close();
}
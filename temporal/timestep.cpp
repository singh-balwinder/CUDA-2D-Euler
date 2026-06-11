#include "definitions.h"

// FOTU: First Order Temporal Update
void timestep(vector <cell> &convol, vector <edge> &E, double &del_t, double &CFL)
{
    // timestep calculation
    int i;
    int numCell = convol.size();
    int numEdge = E.size();
    int numGhostcell=0;
    double aa, cc;
    double bb, dd;
    double cflDen;

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

    int dt_size = numCell - numGhostcell;
    vector <double> dt(dt_size);
    
    for (i=0;i<numCell;i++)
    {
        if (convol[i].cellType != "GHOST")
        {
            const double gamma = 1.4;
            double soundspeed = sqrt(gamma*convol[i].prim_var[3]/convol[i].prim_var[0]);
            double I_nx = (convol[i].quad.l1.nx-convol[i].quad.l3.nx)/2.0;
            double I_ny = (convol[i].quad.l1.ny-convol[i].quad.l3.ny)/2.0;
            aa = fabs((convol[i].prim_var[1]*I_nx)+(convol[i].prim_var[2]*I_ny))+(soundspeed);

            double J_nx = (convol[i].quad.l2.nx - convol[i].quad.l4.nx)/2.0;
            double J_ny = (convol[i].quad.l2.ny - convol[i].quad.l4.ny)/2.0;
            bb = fabs((convol[i].prim_var[1]*J_nx) + (convol[i].prim_var[2]*J_ny))+(soundspeed);

            cc = (convol[i].quad.l1.length + convol[i].quad.l3.length)/2.0;
            dd = (convol[i].quad.l2.length + convol[i].quad.l4.length)/2.0;

            cflDen = (aa*cc) + (bb*dd);

            dt[i] = (CFL*convol[i].getArea())/cflDen;
        }
        else
        {
            continue;
        }
    }

    del_t = *min_element(dt.begin(), dt.end());

    // cout << del_t << endl;

    // cout << numGhostcell << endl;
}

/*
void timestep(vector <cell> &convol, vector <edge> &E, double &del_t, double &CFL)
{
    // kinetic timestep calculation (for perthame scheme)
    int i;
    int numCell = convol.size();
    int numEdge = E.size();
    int numGhostcell=0;
    double del_x;
    double theta;
    double dt1, dt2;
    double CFLDen1, CFLDen2;

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

    int dt_size = numCell - numGhostcell;
    vector <double> dt(dt_size);
    
    for (i=0;i<numCell;i++)
    {
        if (convol[i].cellType != "GHOST")
        {
            del_x = convol[i].quad.l1.length;
            theta = (convol[i].prim_var[3]/convol[i].prim_var[0]);
            CFLDen1 = (fabs(convol[i].prim_var[1]) + sqrt(3.0*theta));
            CFLDen2 = (fabs(convol[i].prim_var[2]) + sqrt(3.0*theta));

            dt1 = (CFL*del_x)/(CFLDen1);

            dt2 = (CFL*del_x)/(CFLDen2);

            dt[i] = min(dt1, dt2);
        }
        else
        {
            continue;
        }
    }

    del_t = *min_element(dt.begin(), dt.end());

    // cout << del_t << endl;

    // cout << numGhostcell << endl;
}
*/
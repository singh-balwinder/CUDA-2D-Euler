#include "definitions.h"

/*
// oblique shock test case
void boundary(vector <edge> &E, vector <cell> &convol)
{
    int i, j;
    int numEdge=E.size();
    // int numCell=convol.size();

    for (i=0;i<numEdge;i++)
    {
        if (E[i].edgeType == "BOUNDARY")
        {
            if (E[i].boundaryType == "INLET1")
            {
                // Top inlet
                convol[E[i].neighbourid].set(1.69997, 2.61934, -0.50633, 1.52819);
            }
            if (E[i].boundaryType == "INLET2")
            {
                // Left inlet
                convol[E[i].neighbourid].set(1.0, 2.9, 0.0, (1/1.4));
            }
            if (E[i].boundaryType == "OUTLET")
            {
                // Outlet
                for (j=0;j<4;j++)
                {
                    convol[E[i].neighbourid].prim_var[j]=convol[E[i].parentid].prim_var[j];
                }
                convol[E[i].neighbourid].convert_p2c();
            }
            if (E[i].boundaryType == "WALL")
            {
                // Wall
                double temp_den, temp_u, temp_v, temp_pr;
                double term1, term2, term3;
                double term4, term5, term6;
                term1=2.0*((E[i].tx)*(E[i].ty)*(convol[E[i].parentid].prim_var[2]));
                term2=-(E[i].ty)*(E[i].ty)*(convol[E[i].parentid].prim_var[1]);
                term3=(E[i].tx)*(E[i].tx)*(convol[E[i].parentid].prim_var[1]);
                term4=-(E[i].tx)*(E[i].tx)*(convol[E[i].parentid].prim_var[2]);
                term5=2.0*((E[i].tx)*(E[i].ty)*(convol[E[i].parentid].prim_var[1]));
                term6=(E[i].ty)*(E[i].ty)*(convol[E[i].parentid].prim_var[2]);
                temp_den=convol[E[i].parentid].prim_var[0];
                temp_u=term1 + term2 + term3;
                temp_v=term4 + term5 + term6;
                temp_pr=convol[E[i].parentid].prim_var[3];
                convol[E[i].neighbourid].set(temp_den, temp_u, temp_v, temp_pr);
            }
        }
    }
}
*/

/*
// An oblique shock in the shock wave-boundary layer interaction
void boundary(vector <edge> &E, vector <cell> &convol)
{
    int i, j;
    int numEdge=E.size();
    // int numCell=convol.size();

    for (i=0;i<numEdge;i++)
    {
        if (E[i].edgeType == "BOUNDARY")
        {
            // Top inlet
            if (E[i].boundaryType == "INLET1")
            {
                convol[E[i].neighbourid].set(1.1706, 0.9617, -0.064, 0.1916);
            }

            // Left inlet
            if (E[i].boundaryType == "INLET2")
            {
                if (convol[E[i].parentid].quad.centroidy <= 0.765)
                {
                    convol[E[i].neighbourid].set(1.0, 1.0, 0.0, 0.154);
                }
                else
                {
                    convol[E[i].neighbourid].set(1.1706, 0.9617, -0.064, 0.1916);
                }
            }

            // Outlet
            if (E[i].boundaryType == "OUTLET")
            {
                for (j=0;j<4;j++)
                {
                    convol[E[i].neighbourid].prim_var[j]=convol[E[i].parentid].prim_var[j];
                }
                convol[E[i].neighbourid].convert_p2c();
            }

            // Wall
            if (E[i].boundaryType == "WALL")
            {
                double temp_den, temp_u, temp_v, temp_pr;
                double term1, term2, term3;
                double term4, term5, term6;
                term1=2.0*((E[i].tx)*(E[i].ty)*(convol[E[i].parentid].prim_var[2]));
                term2=-(E[i].ty)*(E[i].ty)*(convol[E[i].parentid].prim_var[1]);
                term3=(E[i].tx)*(E[i].tx)*(convol[E[i].parentid].prim_var[1]);
                term4=-(E[i].tx)*(E[i].tx)*(convol[E[i].parentid].prim_var[2]);
                term5=2.0*((E[i].tx)*(E[i].ty)*(convol[E[i].parentid].prim_var[1]));
                term6=(E[i].ty)*(E[i].ty)*(convol[E[i].parentid].prim_var[2]);
                temp_den=convol[E[i].parentid].prim_var[0];
                temp_u=term1 + term2 + term3;
                temp_v=term4 + term5 + term6;
                temp_pr=convol[E[i].parentid].prim_var[3];
                convol[E[i].neighbourid].set(temp_den, temp_u, temp_v, temp_pr);
            }
        }
    }
}
*/

/*
// horizontal flow test case
void boundary(vector <edge> &E, vector <cell> &convol)
{
    int i, j;
    int numEdge=E.size();
    // int numCell=convol.size();

    for (i=0;i<numEdge;i++)
    {
        if (E[i].edgeType == "BOUNDARY")
        {
            if (E[i].boundaryType == "INLET1")
            {
                // Top inlet
                convol[E[i].neighbourid].set(1.4, 2.0, 0.0, 1.0);
            }
            if (E[i].boundaryType == "INLET2")
            {
                // Left inlet
                convol[E[i].neighbourid].set(1.4, 3.0, 0.0, 1.0);
            }
            if (E[i].boundaryType == "OUTLET")
            {
                // Outlet
                for (j=0;j<4;j++)
                {
                    convol[E[i].neighbourid].prim_var[j]=convol[E[i].parentid].prim_var[j];
                }
                convol[E[i].neighbourid].convert_p2c();
            }
            if (E[i].boundaryType == "WALL")
            {
                // Wall
                double temp_den, temp_u, temp_v, temp_pr;
                double term1, term2, term3;
                double term4, term5, term6;
                term1=2.0*((E[i].tx)*(E[i].ty)*(convol[E[i].parentid].prim_var[2]));
                term2=-(E[i].ty)*(E[i].ty)*(convol[E[i].parentid].prim_var[1]);
                term3=(E[i].tx)*(E[i].tx)*(convol[E[i].parentid].prim_var[1]);
                term4=-(E[i].tx)*(E[i].tx)*(convol[E[i].parentid].prim_var[2]);
                term5=2.0*((E[i].tx)*(E[i].ty)*(convol[E[i].parentid].prim_var[1]));
                term6=(E[i].ty)*(E[i].ty)*(convol[E[i].parentid].prim_var[2]);
                temp_den=convol[E[i].parentid].prim_var[0];
                temp_u=term1 + term2 + term3;
                temp_v=term4 + term5 + term6;
                temp_pr=convol[E[i].parentid].prim_var[3];
                convol[E[i].neighbourid].set(temp_den, temp_u, temp_v, temp_pr);
            }
        }
    }
}
*/

/*
// odd-even decoupling test case
void boundary(vector <edge> &E, vector <cell> &convol)
{
    int i, j;
    int numEdge=E.size();
    // int numCell=convol.size();

    for (i=0;i<numEdge;i++)
    {
        if (E[i].edgeType == "BOUNDARY")
        {
            if (E[i].boundaryType == "INLET")
            {
                // Top inlet
                convol[E[i].neighbourid].set(7.364, 4.8602, 0.0, 41.833);
            }

            if (E[i].boundaryType == "OUTLET")
            {
                // Outlet
                for (j=0;j<4;j++)
                {
                    convol[E[i].neighbourid].prim_var[j]=convol[E[i].parentid].prim_var[j];
                }
                convol[E[i].neighbourid].convert_p2c();
            }

            if (E[i].boundaryType == "WALL")
            {
                // Wall
                double temp_den, temp_u, temp_v, temp_pr;
                double term1, term2, term3;
                double term4, term5, term6;
                term1=2.0*((E[i].tx)*(E[i].ty)*(convol[E[i].parentid].prim_var[2]));
                term2=-(E[i].ty)*(E[i].ty)*(convol[E[i].parentid].prim_var[1]);
                term3=(E[i].tx)*(E[i].tx)*(convol[E[i].parentid].prim_var[1]);
                term4=-(E[i].tx)*(E[i].tx)*(convol[E[i].parentid].prim_var[2]);
                term5=2.0*((E[i].tx)*(E[i].ty)*(convol[E[i].parentid].prim_var[1]));
                term6=(E[i].ty)*(E[i].ty)*(convol[E[i].parentid].prim_var[2]);
                temp_den=convol[E[i].parentid].prim_var[0];
                temp_u=term1 + term2 + term3;
                temp_v=term4 + term5 + term6;
                temp_pr=convol[E[i].parentid].prim_var[3];
                convol[E[i].neighbourid].set(temp_den, temp_u, temp_v, temp_pr);
            }
        }
    }
}
*/

// general B.C.
void boundary(vector <edge> &E, vector <cell> &convol)
{
    int i, j;
    int numEdge=E.size();

    for (i=0;i<numEdge;i++)
    {
        if (E[i].edgeType == "BOUNDARY")
        {
            if (E[i].boundaryType == "INLET")
            {
                // Inlet
                convol[E[i].neighbourid].set(7.041132, 4.07795, 0.0, 30.05945);
            }

            if (E[i].boundaryType == "OUTLET")
            {
                // Outlet
                for (j=0;j<4;j++)
                {
                    convol[E[i].neighbourid].prim_var[j]=convol[E[i].parentid].prim_var[j];
                }
                convol[E[i].neighbourid].convert_p2c();
            }

            if (E[i].boundaryType == "WALL")
            {
                // Wall
                double temp_den, temp_u, temp_v, temp_pr;
                double term1, term2, term3;
                double term4, term5, term6;
                term1=2.0*((E[i].tx)*(E[i].ty)*(convol[E[i].parentid].prim_var[2]));
                term2=-(E[i].ty)*(E[i].ty)*(convol[E[i].parentid].prim_var[1]);
                term3=(E[i].tx)*(E[i].tx)*(convol[E[i].parentid].prim_var[1]);
                term4=-(E[i].tx)*(E[i].tx)*(convol[E[i].parentid].prim_var[2]);
                term5=2.0*((E[i].tx)*(E[i].ty)*(convol[E[i].parentid].prim_var[1]));
                term6=(E[i].ty)*(E[i].ty)*(convol[E[i].parentid].prim_var[2]);
                temp_den=convol[E[i].parentid].prim_var[0];
                temp_u=term1 + term2 + term3;
                temp_v=term4 + term5 + term6;
                temp_pr=convol[E[i].parentid].prim_var[3];
                convol[E[i].neighbourid].set(temp_den, temp_u, temp_v, temp_pr);
            }
        }
    }
}

/*
// kinked mach stem test case
void boundary(vector <edge> &E, vector <cell> &convol, double &h)
{
    int i, j;
    int numEdge=E.size();

    for (i=0;i<numEdge;i++)
    {
        if (E[i].edgeType == "BOUNDARY")
        {
            if (E[i].boundaryType == "INLET")
            {
                // Inlet
                convol[E[i].neighbourid].set(7.999, 7.14461, -4.12494, 116.5);
            }

            if (E[i].boundaryType == "OUTLET")
            {
                // Outlet
                for (j=0;j<4;j++)
                {
                    convol[E[i].neighbourid].prim_var[j]=convol[E[i].parentid].prim_var[j];
                }
                convol[E[i].neighbourid].convert_p2c();
            }

            if (E[i].boundaryType == "WALL")
            {
                // Wall
                double temp_den, temp_u, temp_v, temp_pr;
                double term1, term2, term3;
                double term4, term5, term6;
                term1=2.0*((E[i].tx)*(E[i].ty)*(convol[E[i].parentid].prim_var[2]));
                term2=-(E[i].ty)*(E[i].ty)*(convol[E[i].parentid].prim_var[1]);
                term3=(E[i].tx)*(E[i].tx)*(convol[E[i].parentid].prim_var[1]);
                term4=-(E[i].tx)*(E[i].tx)*(convol[E[i].parentid].prim_var[2]);
                term5=2.0*((E[i].tx)*(E[i].ty)*(convol[E[i].parentid].prim_var[1]));
                term6=(E[i].ty)*(E[i].ty)*(convol[E[i].parentid].prim_var[2]);
                temp_den=convol[E[i].parentid].prim_var[0];
                temp_u=term1 + term2 + term3;
                temp_v=term4 + term5 + term6;
                temp_pr=convol[E[i].parentid].prim_var[3];
                convol[E[i].neighbourid].set(temp_den, temp_u, temp_v, temp_pr);
            }

            if (E[i].boundaryType == "DYNAMIC")
            {
                double pcell = E[i].parentid;
                double ncell = E[i].neighbourid;
                double init_dist = (1.0/6.0) + (1.0/sqrt(3.0));
                if (convol[pcell].quad.centroidx < (init_dist+h))
                {
                    convol[ncell].set(7.999, 7.14461, -4.12494, 116.5);
                }
                else
                {
                    double temp_den1, temp_u1, temp_v1, temp_pr1;
                    double term11, term22, term33;
                    double term44, term55, term66;
                    term11=2.0*((E[i].tx)*(E[i].ty)*(convol[E[i].parentid].prim_var[2]));
                    term22=-(E[i].ty)*(E[i].ty)*(convol[E[i].parentid].prim_var[1]);
                    term33=(E[i].tx)*(E[i].tx)*(convol[E[i].parentid].prim_var[1]);
                    term44=-(E[i].tx)*(E[i].tx)*(convol[E[i].parentid].prim_var[2]);
                    term55=2.0*((E[i].tx)*(E[i].ty)*(convol[E[i].parentid].prim_var[1]));
                    term66=(E[i].ty)*(E[i].ty)*(convol[E[i].parentid].prim_var[2]);
                    temp_den1=convol[E[i].parentid].prim_var[0];
                    temp_u1=term11 + term22 + term33;
                    temp_v1=term44 + term55 + term66;
                    temp_pr1=convol[E[i].parentid].prim_var[3];
                    convol[ncell].set(temp_den1, temp_u1, temp_v1, temp_pr1);
                }
            }
        }
    }
}
*/

// print statements
/*
for (i=0;i<numCell;i++)
{
    cout << convol[i].prim_var[2] << endl;
}

for (i=0;i<numEdge;i++)
{
    cout << i << "  " << "nx = " << E[i].nx << "  " << "ny = " << E[i].ny << endl;
}
*/
    
/*
cout << convol[1].quad.l1.nx << "  " << convol[1].quad.l1.ny << endl;
cout << convol[1].quad.l2.nx << "  " << convol[1].quad.l2.ny << endl;
cout << convol[1].quad.l3.nx << "  " << convol[1].quad.l3.ny << endl;
cout << convol[1].quad.l4.nx << "  " << convol[1].quad.l4.ny << endl;
*/

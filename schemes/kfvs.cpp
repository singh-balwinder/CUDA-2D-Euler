#include "definitions.h"

// kinetic flux vector splitting scheme (Deshpande et. al.)
void kfvs(vector <edge> &E, vector <cell> &convol)
{
    int i, j;
    int numEdge=E.size();
    double const pi = 3.141592;
    double const gamma = 1.4;

    double parent_s, neighbour_s;
    double parent_X, neighbour_X;
    double parent_Y, neighbour_Y;
    vector <double> positive_flux(4);
    vector <double> negative_flux(4);
    double temp1, temp2, temp3, temp4;
    double temp5, temp6, temp7, temp8;

    for (i=0;i<numEdge;i++)
    {
        int parentcellid=E[i].parentid;
        int neighbourcellid=E[i].neighbourid;

        convol[parentcellid].ncv(E[i].nx, E[i].ny);
        convol[neighbourcellid].ncv(E[i].nx, E[i].ny);

        parent_s = convol[parentcellid].norm_prim_var[1]*sqrt(convol[parentcellid].norm_prim_var[0]/(2.0*convol[parentcellid].norm_prim_var[3]));
        neighbour_s = convol[neighbourcellid].norm_prim_var[1]*sqrt(convol[neighbourcellid].norm_prim_var[0]/(2.0*convol[neighbourcellid].norm_prim_var[3]));
        parent_X = (1.0 + erf(parent_s))/2.0;
        neighbour_X = (1.0 - erf(neighbour_s))/2.0;
        parent_Y = exp(-powf(parent_s,2))*sqrt(convol[parentcellid].norm_prim_var[3]/(2.0*pi*convol[parentcellid].norm_prim_var[0]));
        neighbour_Y = exp(-powf(neighbour_s,2))*sqrt(convol[neighbourcellid].norm_prim_var[3]/(2.0*pi*convol[neighbourcellid].norm_prim_var[0]));

        positive_flux[0] = (convol[parentcellid].norm_prim_var[0]*convol[parentcellid].norm_prim_var[1]*parent_X) + (convol[parentcellid].norm_prim_var[0]*parent_Y);
        negative_flux[0] = (convol[neighbourcellid].norm_prim_var[0]*convol[neighbourcellid].norm_prim_var[1]*neighbour_X) - (convol[neighbourcellid].norm_prim_var[0]*neighbour_Y);
        positive_flux[1] = ((convol[parentcellid].norm_prim_var[3] + (convol[parentcellid].norm_prim_var[0]*powf(convol[parentcellid].norm_prim_var[1],2)))*parent_X) + (convol[parentcellid].norm_prim_var[0]*convol[parentcellid].norm_prim_var[1]*parent_Y);
        negative_flux[1] = ((convol[neighbourcellid].norm_prim_var[3] + (convol[neighbourcellid].norm_prim_var[0]*powf(convol[neighbourcellid].norm_prim_var[1],2)))*neighbour_X) - (convol[neighbourcellid].norm_prim_var[0]*convol[neighbourcellid].norm_prim_var[1]*neighbour_Y);
        positive_flux[2] = (convol[parentcellid].norm_prim_var[0]*convol[parentcellid].norm_prim_var[1]*convol[parentcellid].norm_prim_var[2]*parent_X) + (convol[parentcellid].norm_prim_var[0]*convol[parentcellid].norm_prim_var[2]*parent_Y);
        negative_flux[2] = (convol[neighbourcellid].norm_prim_var[0]*convol[neighbourcellid].norm_prim_var[1]*convol[neighbourcellid].norm_prim_var[2]*neighbour_X) - (convol[neighbourcellid].norm_prim_var[0]*convol[neighbourcellid].norm_prim_var[2]*neighbour_Y);

        temp1 = (gamma/(gamma-1.0))*convol[parentcellid].norm_prim_var[3]*convol[parentcellid].norm_prim_var[1]*parent_X;
        temp2 = (0.5*convol[parentcellid].norm_prim_var[0]*convol[parentcellid].norm_prim_var[1]*parent_X)*(powf(convol[parentcellid].norm_prim_var[1],2)+powf(convol[parentcellid].norm_prim_var[2],2));
        temp3 = ((gamma+1.0)/(2.0*(gamma-1.0)))*convol[parentcellid].norm_prim_var[3]*parent_Y;
        temp4 = (0.5*convol[parentcellid].norm_prim_var[0]*parent_Y)*(powf(convol[parentcellid].norm_prim_var[1],2)+powf(convol[parentcellid].norm_prim_var[2],2));
        positive_flux[3] = temp1 + temp2 + temp3 +temp4;

        temp5 = (gamma/(gamma-1.0))*convol[neighbourcellid].norm_prim_var[3]*convol[neighbourcellid].norm_prim_var[1]*neighbour_X;
        temp6 = (0.5*convol[neighbourcellid].norm_prim_var[0]*convol[neighbourcellid].norm_prim_var[1]*neighbour_X)*(powf(convol[neighbourcellid].norm_prim_var[1],2)+powf(convol[neighbourcellid].norm_prim_var[2],2));
        temp7 = ((gamma+1.0)/(2.0*(gamma-1.0)))*convol[neighbourcellid].norm_prim_var[3]*neighbour_Y;
        temp8 = (0.5*convol[neighbourcellid].norm_prim_var[0]*neighbour_Y)*(powf(convol[neighbourcellid].norm_prim_var[1],2)+powf(convol[neighbourcellid].norm_prim_var[2],2));
        negative_flux[3] = temp5 + temp6 - temp7 -temp8;

        for (j=0;j<4;j++)
        {
            E[i].normal_interface_flux[j] = (positive_flux[j] + negative_flux[j])*(E[i].L.length);
        }

        E[i].cart_interface_flux();
    }
}
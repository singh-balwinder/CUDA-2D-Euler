#include "definitions.h"

// peculiar velocity-based upwind scheme (Raghurama Rao et. al.)
void pvu(vector <edge> &E, vector <cell> &convol)
{
    int i, j;
    int numEdge=E.size();
    double const pi = 3.141592;
    // double const gamma = 1.4;
    
    // c: convective
    vector <double> c_positive_flux(4);
    vector <double> c_negative_flux(4);
    vector <double> c_flux(4);

    // p : pressure
    vector <double> p_positive_flux(4);
    vector <double> p_negative_flux(4);
    vector <double> p_flux(4);

    // double pos_speed, neg_speed;
    double temp1, temp2, temp3;
    double temp4, temp5, temp6;

    for (i=0;i<numEdge;i++)
    {
        int parentcellid=E[i].parentid;
        int neighbourcellid=E[i].neighbourid;

        convol[parentcellid].ncv(E[i].nx, E[i].ny);
        convol[neighbourcellid].ncv(E[i].nx, E[i].ny);
        
        // convective flux
        
        /*
        // first type
        double avgVelocity = (convol[parentcellid].norm_prim_var[1] + convol[neighbourcellid].norm_prim_var[1])/2.0;
        pos_speed = (avgVelocity + fabs(avgVelocity))/2.0;
        neg_speed = (avgVelocity - fabs(avgVelocity))/2.0;
        for (j=0;j<4;j++)
        {
            c_positive_flux[j] = pos_speed*(convol[parentcellid].norm_cons_var[j]);
            c_negative_flux[j] = neg_speed*(convol[neighbourcellid].norm_cons_var[j]);
            c_flux[j] = c_positive_flux[j] + c_negative_flux[j];
        }
        */

        /*
        // second type
        pos_speed = (convol[parentcellid].norm_prim_var[1] + fabs(convol[parentcellid].norm_prim_var[1]))/2.0;
        neg_speed = (convol[neighbourcellid].norm_prim_var[1] - fabs(convol[neighbourcellid].norm_prim_var[1]))/2.0;
        for (j=0;j<4;j++)
        {
            c_positive_flux[j] = pos_speed*(convol[parentcellid].norm_cons_var[j]);
            c_negative_flux[j] = neg_speed*(convol[neighbourcellid].norm_cons_var[j]);
            c_flux[j] = c_positive_flux[j] + c_negative_flux[j];
        }
        */

        // third type
        double alpha = (fabs(convol[parentcellid].norm_prim_var[1]) + fabs(convol[neighbourcellid].norm_prim_var[1]))/2.0;
        for (j=0;j<4;j++)
        {
            c_positive_flux[j] = ((convol[parentcellid].norm_cons_var[j]*convol[parentcellid].norm_prim_var[1])/2.0) + (alpha*convol[parentcellid].norm_cons_var[j])/2.0;
            c_negative_flux[j] = ((convol[neighbourcellid].norm_cons_var[j]*convol[neighbourcellid].norm_prim_var[1])/2.0) - (alpha*convol[neighbourcellid].norm_cons_var[j])/2.0;
            c_flux[j] = c_positive_flux[j] + c_negative_flux[j];
        }

        // Pressure flux
        double parent_beta, neighbour_beta;
        parent_beta = (0.5*convol[parentcellid].norm_prim_var[0])/convol[parentcellid].norm_prim_var[3];
        neighbour_beta = (0.5*convol[neighbourcellid].norm_prim_var[0])/convol[neighbourcellid].norm_prim_var[3];

        p_positive_flux[0] = convol[parentcellid].norm_prim_var[0]/(2.0*sqrt(pi*parent_beta));
        p_negative_flux[0] = -convol[neighbourcellid].norm_prim_var[0]/(2.0*sqrt(pi*neighbour_beta));
        p_positive_flux[1] = (convol[parentcellid].norm_prim_var[3]/2.0) + ((convol[parentcellid].norm_prim_var[0]/(2.0*sqrt(pi*parent_beta)))*(convol[parentcellid].norm_prim_var[1]));
        p_negative_flux[1] = (convol[neighbourcellid].norm_prim_var[3]/2.0) - ((convol[neighbourcellid].norm_prim_var[0]/(2.0*sqrt(pi*neighbour_beta)))*(convol[neighbourcellid].norm_prim_var[1]));
        p_positive_flux[2] = ((convol[parentcellid].norm_prim_var[0]/(2.0*sqrt(pi*parent_beta)))*(convol[parentcellid].norm_prim_var[2]));
        p_negative_flux[2] = -((convol[neighbourcellid].norm_prim_var[0]/(2.0*sqrt(pi*neighbour_beta)))*(convol[neighbourcellid].norm_prim_var[2]));

        temp1 = 0.5*convol[parentcellid].norm_prim_var[3]*convol[parentcellid].norm_prim_var[1];
        temp2 = convol[parentcellid].norm_prim_var[3]/(4.0*sqrt(pi*parent_beta));
        temp3 = convol[parentcellid].norm_cons_var[3]/(2.0*sqrt(pi*parent_beta));
        p_positive_flux[3] = temp1 + temp2 + temp3;

        temp4 = 0.5*convol[neighbourcellid].norm_prim_var[3]*convol[neighbourcellid].norm_prim_var[1];
        temp5 = -convol[neighbourcellid].norm_prim_var[3]/(4.0*sqrt(pi*neighbour_beta));
        temp6 = -convol[neighbourcellid].norm_cons_var[3]/(2.0*sqrt(pi*neighbour_beta));
        p_negative_flux[3] = temp4 + temp5 + temp6;

        for (j=0;j<4;j++)
        {
            p_flux[j] = p_positive_flux[j] + p_negative_flux[j];
        }

        for (j=0;j<4;j++)
        {
            E[i].normal_interface_flux[j] = (c_flux[j] + p_flux[j])*(E[i].L.length);
        }

        E[i].cart_interface_flux();
    }
}
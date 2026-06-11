#include "definitions.h"

// modified peculiar velocity-based upwind scheme
void blendedPVU(vector <edge> &E, vector <cell> &convol)
{
    int i, j;
    int numEdge=E.size();
    double const pi = 3.141592;
    double const gamma = 1.4;
    
    // c: convective
    vector <double> c_positive_flux(4);
    vector <double> c_negative_flux(4);
    vector <double> c_flux(4);

    // p1: pressure flux: PVU+ variant 3
    vector <double> p1_positive_flux(4);
    vector <double> p1_negative_flux(4);
    vector <double> p1_flux(4);

    // p2: pressure flux: PVU
    vector <double> p2_positive_flux(4);
    vector <double> p2_negative_flux(4);
    vector <double> p2_flux(4);

    double temp1, temp2, temp3, temp4;
    double temp5, temp6, temp7, temp8;

    for (i=0;i<numEdge;i++)
    {
        int parentcellid=E[i].parentid;
        int neighbourcellid=E[i].neighbourid;

        convol[parentcellid].ncv(E[i].nx, E[i].ny);
        convol[neighbourcellid].ncv(E[i].nx, E[i].ny);

        // convective flux
        /*
        // first type
        double pos_speed, neg_speed;
        pos_speed = (convol[parentcellid].norm_prim_var[1] + abs(convol[parentcellid].norm_prim_var[1]))/2.0;
        neg_speed = (convol[neighbourcellid].norm_prim_var[1] - abs(convol[neighbourcellid].norm_prim_var[1]))/2.0;

        for (j=0;j<4;j++)
        {
            c_positive_flux[j] = pos_speed*(convol[parentcellid].norm_cons_var[j]);
            c_negative_flux[j] = neg_speed*(convol[neighbourcellid].norm_cons_var[j]);
            c_flux[j] = c_positive_flux[j] + c_negative_flux[j];
        }
        */

        /*
        // second type
        double avg_speed, pos_speed, neg_speed;
        avg_speed = (convol[parentcellid].norm_prim_var[1] + convol[neighbourcellid].norm_prim_var[1])/2.0;
        pos_speed = (avg_speed + fabs(avg_speed))/2.0;
        neg_speed = (avg_speed - fabs(avg_speed))/2.0;

        for (j=0;j<4;j++)
        {
            c_positive_flux[j] = pos_speed*(convol[parentcellid].norm_cons_var[j]);
            c_negative_flux[j] = neg_speed*(convol[neighbourcellid].norm_cons_var[j]);
            c_flux[j] = c_positive_flux[j] + c_negative_flux[j];
        }
        */

        // third type
        double alpha = max(fabs(convol[parentcellid].norm_prim_var[1]), fabs(convol[neighbourcellid].norm_prim_var[1]));
        for (j=0;j<4;j++)
        {
            c_positive_flux[j] = ((convol[parentcellid].norm_cons_var[j]*convol[parentcellid].norm_prim_var[1])/2.0) + (alpha*convol[parentcellid].norm_cons_var[j])/2.0;
            c_negative_flux[j] = ((convol[neighbourcellid].norm_cons_var[j]*convol[neighbourcellid].norm_prim_var[1])/2.0) - (alpha*convol[neighbourcellid].norm_cons_var[j])/2.0;
            c_flux[j] = c_positive_flux[j] + c_negative_flux[j];
        }


        // pressure flux: PVU+ variant 3
        double parent_beta, neighbour_beta;
        double pcell_n_mach, ncell_n_mach;
        double parent_phi, neighbour_phi;

        parent_beta = (0.5*convol[parentcellid].norm_prim_var[0])/convol[parentcellid].norm_prim_var[3];
        neighbour_beta = (0.5*convol[neighbourcellid].norm_prim_var[0])/convol[neighbourcellid].norm_prim_var[3];
        pcell_n_mach = convol[parentcellid].norm_prim_var[1]*sqrt(2.0*parent_beta/gamma);
        ncell_n_mach = convol[neighbourcellid].norm_prim_var[1]*sqrt(2.0*neighbour_beta/gamma);
        parent_phi = ((gamma/2.0)*powf(pcell_n_mach,2))/(1.0+((gamma/2.0)*powf(pcell_n_mach,2)));
        neighbour_phi = ((gamma/2.0)*powf(ncell_n_mach,2))/(1.0+((gamma/2.0)*powf(ncell_n_mach,2)));

        p1_positive_flux[0] = (convol[parentcellid].norm_prim_var[0]/(2.0*sqrt(pi*parent_beta)))*parent_phi;
        p1_negative_flux[0] = (-convol[neighbourcellid].norm_prim_var[0]/(2.0*sqrt(pi*neighbour_beta)))*neighbour_phi;
        p1_positive_flux[1] = (convol[parentcellid].norm_prim_var[3]/2.0) + (((convol[parentcellid].norm_prim_var[0]/(2.0*sqrt(pi*parent_beta)))*(convol[parentcellid].norm_prim_var[1]))*parent_phi);
        p1_negative_flux[1] = (convol[neighbourcellid].norm_prim_var[3]/2.0) - (((convol[neighbourcellid].norm_prim_var[0]/(2.0*sqrt(pi*neighbour_beta)))*(convol[neighbourcellid].norm_prim_var[1]))*neighbour_phi);
        p1_positive_flux[2] = ((convol[parentcellid].norm_prim_var[0]/(2.0*sqrt(pi*parent_beta)))*(convol[parentcellid].norm_prim_var[2]))*parent_phi;
        p1_negative_flux[2] = -((convol[neighbourcellid].norm_prim_var[0]/(2.0*sqrt(pi*neighbour_beta)))*(convol[neighbourcellid].norm_prim_var[2]))*neighbour_phi;

        temp1 = 0.5*convol[parentcellid].norm_prim_var[3]*convol[parentcellid].norm_prim_var[1];
        temp2 = (convol[parentcellid].norm_prim_var[3]/(4.0*sqrt(pi*parent_beta)))*((3.0-gamma)/(gamma-1.0))*parent_phi;
        temp3 = ((convol[parentcellid].norm_prim_var[0]*powf(convol[parentcellid].norm_prim_var[1],2))/(4.0*sqrt(pi*parent_beta)))*parent_phi;
        temp4 = (convol[parentcellid].norm_prim_var[3]/(2.0*sqrt(pi*parent_beta)))*powf(parent_phi,2);
        p1_positive_flux[3] = temp1 + temp2 + temp3 +temp4;

        temp5 = 0.5*convol[neighbourcellid].norm_prim_var[3]*convol[neighbourcellid].norm_prim_var[1];
        temp6 = -(convol[neighbourcellid].norm_prim_var[3]/(4.0*sqrt(pi*neighbour_beta)))*((3.0-gamma)/(gamma-1.0))*neighbour_phi;
        temp7 = -((convol[neighbourcellid].norm_prim_var[0]*powf(convol[neighbourcellid].norm_prim_var[1],2))/(4.0*sqrt(pi*neighbour_beta)))*neighbour_phi;
        temp8 = -(convol[neighbourcellid].norm_prim_var[3]/(2.0*sqrt(pi*neighbour_beta)))*powf(neighbour_phi,2);
        p1_negative_flux[3] = temp5 + temp6 + temp7 + temp8;

        for (j=0;j<4;j++)
        {
            p1_flux[j] = p1_positive_flux[j] + p1_negative_flux[j];
        }

        // pressure flux: PVU
        parent_beta = (0.5*convol[parentcellid].norm_prim_var[0])/convol[parentcellid].norm_prim_var[3];
        neighbour_beta = (0.5*convol[neighbourcellid].norm_prim_var[0])/convol[neighbourcellid].norm_prim_var[3];

        p2_positive_flux[0] = convol[parentcellid].norm_prim_var[0]/(2.0*sqrt(pi*parent_beta));
        p2_negative_flux[0] = -convol[neighbourcellid].norm_prim_var[0]/(2.0*sqrt(pi*neighbour_beta));
        p2_positive_flux[1] = (convol[parentcellid].norm_prim_var[3]/2.0) + ((convol[parentcellid].norm_prim_var[0]/(2.0*sqrt(pi*parent_beta)))*(convol[parentcellid].norm_prim_var[1]));
        p2_negative_flux[1] = (convol[neighbourcellid].norm_prim_var[3]/2.0) - ((convol[neighbourcellid].norm_prim_var[0]/(2.0*sqrt(pi*neighbour_beta)))*(convol[neighbourcellid].norm_prim_var[1]));
        p2_positive_flux[2] = ((convol[parentcellid].norm_prim_var[0]/(2.0*sqrt(pi*parent_beta)))*(convol[parentcellid].norm_prim_var[2]));
        p2_negative_flux[2] = -((convol[neighbourcellid].norm_prim_var[0]/(2.0*sqrt(pi*neighbour_beta)))*(convol[neighbourcellid].norm_prim_var[2]));

        temp1 = 0.5*convol[parentcellid].norm_prim_var[3]*convol[parentcellid].norm_prim_var[1];
        temp2 = convol[parentcellid].norm_prim_var[3]/(4.0*sqrt(pi*parent_beta));
        temp3 = convol[parentcellid].norm_cons_var[3]/(2.0*sqrt(pi*parent_beta));
        p2_positive_flux[3] = temp1 + temp2 + temp3;

        temp4 = 0.5*convol[neighbourcellid].norm_prim_var[3]*convol[neighbourcellid].norm_prim_var[1];
        temp5 = -convol[neighbourcellid].norm_prim_var[3]/(4.0*sqrt(pi*neighbour_beta));
        temp6 = -convol[neighbourcellid].norm_cons_var[3]/(2.0*sqrt(pi*neighbour_beta));
        p2_negative_flux[3] = temp4 + temp5 + temp6;

        for (j=0;j<4;j++)
        {
            p2_flux[j] = p2_positive_flux[j] + p2_negative_flux[j];
        }

        /*
        double shockSwitch = 1.0 - exp(-50.0*fabs(convol[parentcellid].norm_prim_var[3] - convol[neighbourcellid].norm_prim_var[3])/(convol[parentcellid].norm_prim_var[3] + convol[neighbourcellid].norm_prim_var[3]));
        */

        double shockSwitch;
        if (fabs(convol[parentcellid].norm_prim_var[1] - convol[neighbourcellid].norm_prim_var[1]) <= powf(10,-4) && fabs(convol[parentcellid].norm_prim_var[3] - convol[neighbourcellid].norm_prim_var[3]) <= powf(10,-4))
        {
            shockSwitch = 0.0;
        }
        else
        {
            shockSwitch = 1.0;
        }

        // total flux
        for (j=0;j<4;j++)
        {
            E[i].normal_interface_flux[j] = (c_flux[j] + ((1.0 - shockSwitch)*p1_flux[j]) + (shockSwitch*p2_flux[j]))*(E[i].L.length);
        }

        E[i].cart_interface_flux();
    }
}
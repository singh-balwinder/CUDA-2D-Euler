#include "definitions.h"

// peculiar velocity-based upwind scheme to capture contacts
// variant 2
void pvuc2(vector <edge> &E, vector <cell> &convol)
{
    int i, j;
    double const pi = 3.141592;
    double const gamma = 1.4;

    double avg_sound_speed;
    int numEdge=E.size();
    
    /*
     int numCell = convol.size();
    int numGhostcell = 0;
    for(j=0;j<numEdge;j++)
    {
        if (E[j].edgeType == "BOUNDARY")
        {
            numGhostcell = numGhostcell + 1;
        }
        else
        {
            continue;
        }
    }
    numCell = numCell - numGhostcell;

    vector <double> soundspeed(numCell);
    for (j=0;j<numCell;j++)
    {
        soundspeed[j] = sqrt((gamma*convol[j].prim_var[3])/convol[j].prim_var[0]);
    }
    avg_sound_speed = (*max_element(soundspeed.begin(),soundspeed.end()))/2.0;
    */
    
    // c: convective
    vector <double> c_positive_flux(4);
    vector <double> c_negative_flux(4);
    vector <double> c_flux(4);

    // p : pressure
    vector <double> p_positive_flux(4);
    vector <double> p_negative_flux(4);
    vector <double> p_flux(4);
    
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
        double u_plus = (convol[parentcellid].norm_prim_var[1] + fabs(convol[parentcellid].norm_prim_var[1]))/2.0;
        double u_minus = (convol[neighbourcellid].norm_prim_var[1] - fabs(convol[neighbourcellid].norm_prim_var[1]))/2.0;

        for (j=0;j<4;j++)
        {
            c_positive_flux[j] = u_plus*convol[parentcellid].norm_cons_var[j];
            c_negative_flux[j] = u_minus*convol[neighbourcellid].norm_cons_var[j];
            c_flux[j] = c_positive_flux[j] + c_negative_flux[j];
        }
        */
        
        // second type
        double u_average = (convol[parentcellid].norm_prim_var[1] + convol[neighbourcellid].norm_prim_var[1])/2.0;
        double u_plus = (u_average + fabs(u_average))/2.0;
        double u_minus = (u_average - fabs(u_average))/2.0;

        for (j=0;j<4;j++)
        {
            c_positive_flux[j] = u_plus*convol[parentcellid].norm_cons_var[j];
            c_negative_flux[j] = u_minus*convol[neighbourcellid].norm_cons_var[j];
            c_flux[j] = c_positive_flux[j] + c_negative_flux[j];
        }
        
        /*
        // third type
        double alpha = max(fabs(convol[parentcellid].norm_prim_var[1]), fabs(convol[neighbourcellid].norm_prim_var[1]));
        for (j=0;j<4;j++)
        {
            c_positive_flux[j] = ((convol[parentcellid].norm_cons_var[j]*convol[parentcellid].norm_prim_var[1])/2.0) + (alpha*convol[parentcellid].norm_cons_var[j])/2.0;
            c_negative_flux[j] = ((convol[neighbourcellid].norm_cons_var[j]*convol[neighbourcellid].norm_prim_var[1])/2.0) - (alpha*convol[neighbourcellid].norm_cons_var[j])/2.0;
            c_flux[j] = c_positive_flux[j] + c_negative_flux[j];
        }
        */

        // pressure flux
        double pcell_sound_speed, ncell_sound_speed;
        pcell_sound_speed = sqrt((gamma*convol[parentcellid].norm_prim_var[3])/convol[parentcellid].norm_prim_var[0]);
        ncell_sound_speed = sqrt((gamma*convol[neighbourcellid].norm_prim_var[3])/convol[neighbourcellid].norm_prim_var[0]);
        avg_sound_speed = (pcell_sound_speed + ncell_sound_speed)/2.0;

        p_positive_flux[0] = sqrt(gamma/(2.0*pi))*(convol[parentcellid].norm_prim_var[3]/(2.0*avg_sound_speed));
        p_negative_flux[0] = -sqrt(gamma/(2.0*pi))*(convol[neighbourcellid].norm_prim_var[3]/(2.0*avg_sound_speed));
        p_positive_flux[1] = (convol[parentcellid].norm_prim_var[3]/2.0) + (sqrt(gamma/(2.0*pi))*((convol[parentcellid].norm_prim_var[3]*convol[parentcellid].norm_prim_var[1])/(2.0*avg_sound_speed)));
        p_negative_flux[1] = (convol[neighbourcellid].norm_prim_var[3]/2.0) - (sqrt(gamma/(2.0*pi))*((convol[neighbourcellid].norm_prim_var[3]*convol[neighbourcellid].norm_prim_var[1])/(2.0*avg_sound_speed)));
        p_positive_flux[2] = (sqrt(gamma/(2.0*pi))*((convol[parentcellid].norm_prim_var[3]*convol[parentcellid].norm_prim_var[2])/(2.0*avg_sound_speed)));
        p_negative_flux[2] = -(sqrt(gamma/(2.0*pi))*((convol[neighbourcellid].norm_prim_var[3]*convol[neighbourcellid].norm_prim_var[2])/(2.0*avg_sound_speed)));

        temp1 = (convol[parentcellid].norm_prim_var[3]*convol[parentcellid].norm_prim_var[1])/2.0;
        temp2 = ((gamma + 1.0)/(4.0*(gamma - 1.0)))*sqrt(gamma/(2.0*pi))*(convol[parentcellid].norm_prim_var[3]*avg_sound_speed);
        temp3 = (sqrt(gamma/(2.0*pi))*((convol[parentcellid].norm_prim_var[3]*powf(convol[parentcellid].norm_prim_var[1],2))/(4.0*avg_sound_speed))) + (sqrt(gamma/(2.0*pi))*((convol[parentcellid].norm_prim_var[3]*powf(convol[parentcellid].norm_prim_var[2],2))/(4.0*avg_sound_speed)));

        temp4 = (convol[neighbourcellid].norm_prim_var[3]*convol[neighbourcellid].norm_prim_var[1])/2.0;
        temp5 = -((gamma + 1.0)/(4.0*(gamma - 1.0)))*sqrt(gamma/(2.0*pi))*(convol[neighbourcellid].norm_prim_var[3]*avg_sound_speed);
        temp6 = -(sqrt(gamma/(2.0*pi))*((convol[neighbourcellid].norm_prim_var[3]*powf(convol[neighbourcellid].norm_prim_var[1],2))/(4.0*avg_sound_speed))) - (sqrt(gamma/(2.0*pi))*((convol[neighbourcellid].norm_prim_var[3]*powf(convol[neighbourcellid].norm_prim_var[2],2))/(4.0*avg_sound_speed)));

        p_positive_flux[3] = temp1 + temp2 + temp3;
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
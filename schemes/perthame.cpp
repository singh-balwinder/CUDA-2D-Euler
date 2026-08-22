#include "definitions.h"

// perthame scheme
void perthame(vector <edge> &E, vector <cell> &convol)
{
    int i, j;
    int numEdge=E.size();
    double const gamma = 1.4;

    vector <double> positive_flux(4);
    vector <double> negative_flux(4);
    double pcell_norm_mach, ncell_norm_mach;
    double pcell_theta, ncell_theta;
    double temp1, temp2, temp3, temp4;

    for (i=0;i<numEdge;i++)
    {
        int parentcellid=E[i].parentid;
        int neighbourcellid=E[i].neighbourid;

        convol[parentcellid].ncv(E[i].nx, E[i].ny);
        convol[neighbourcellid].ncv(E[i].nx, E[i].ny);

        pcell_theta = (convol[parentcellid].norm_prim_var[3]/convol[parentcellid].norm_prim_var[0]);
        ncell_theta = (convol[neighbourcellid].norm_prim_var[3]/convol[neighbourcellid].norm_prim_var[0]);

        pcell_norm_mach = convol[parentcellid].norm_prim_var[1]/sqrt(3.0*pcell_theta);
        ncell_norm_mach = convol[neighbourcellid].norm_prim_var[1]/sqrt(3.0*ncell_theta);

        // positive flux
        if (pcell_norm_mach > 1.0)
        {
            for (j=0;j<4;j++)
            {
                positive_flux[j] = convol[parentcellid].norm_analytic_flux[j];
            }
        }

        if (fabs(pcell_norm_mach) <= 1.0)
        {
            positive_flux[0] = (convol[parentcellid].norm_prim_var[0]/(4.0*sqrt(3.0*pcell_theta)))*powf((convol[parentcellid].norm_prim_var[1]+sqrt(3.0*pcell_theta)),2);
            positive_flux[1] = (convol[parentcellid].norm_prim_var[0]/(6.0*sqrt(3.0*pcell_theta)))*powf((convol[parentcellid].norm_prim_var[1]+sqrt(3.0*pcell_theta)),3);
            positive_flux[2] = (convol[parentcellid].norm_prim_var[0]/(4.0*sqrt(3.0*pcell_theta)))*(convol[parentcellid].norm_prim_var[2])*(powf((convol[parentcellid].norm_prim_var[1]+sqrt(3.0*pcell_theta)),2));
            
            temp1 = (convol[parentcellid].norm_prim_var[0]/(16.0*sqrt(3.0*pcell_theta)))*powf((convol[parentcellid].norm_prim_var[1]+sqrt(3.0*pcell_theta)),4);
            temp2 = (1.0/(8.0*sqrt(3.0*pcell_theta)))*((((3.0-gamma)/(gamma-1.0))*convol[parentcellid].norm_prim_var[3])+(convol[parentcellid].norm_prim_var[0]*powf(convol[parentcellid].norm_prim_var[2],2)))*powf((convol[parentcellid].norm_prim_var[1]+sqrt(3.0*pcell_theta)),2);
            positive_flux[3] = temp1 + temp2;
        }

        if (pcell_norm_mach < -1.0)
        {
            for (j=0;j<4;j++)
            {
                positive_flux[j] = 0.0;
            }
        }
        
        // negative flux
        if (ncell_norm_mach > 1.0)
        {
            for (j=0;j<4;j++)
            {
                negative_flux[j] = 0.0;
            }
        }

        if (fabs(ncell_norm_mach) <= 1.0)
        {
            negative_flux[0] = -(convol[neighbourcellid].norm_prim_var[0]/(4.0*sqrt(3.0*ncell_theta)))*powf((convol[neighbourcellid].norm_prim_var[1]-sqrt(3.0*ncell_theta)),2);
            negative_flux[1] = -(convol[neighbourcellid].norm_prim_var[0]/(6.0*sqrt(3.0*ncell_theta)))*powf((convol[neighbourcellid].norm_prim_var[1]-sqrt(3.0*ncell_theta)),3);
            negative_flux[2] = -(convol[neighbourcellid].norm_prim_var[0]/(4.0*sqrt(3.0*ncell_theta)))*(convol[neighbourcellid].norm_prim_var[2])*(powf((convol[neighbourcellid].norm_prim_var[1]-sqrt(3.0*ncell_theta)),2));
            
            temp3 = -(convol[neighbourcellid].norm_prim_var[0]/(16.0*sqrt(3.0*ncell_theta)))*powf((convol[neighbourcellid].norm_prim_var[1]-sqrt(3.0*ncell_theta)),4);
            temp4 = -(1.0/(8.0*sqrt(3.0*ncell_theta)))*((((3.0-gamma)/(gamma-1.0))*convol[neighbourcellid].norm_prim_var[3])+(convol[neighbourcellid].norm_prim_var[0]*powf(convol[neighbourcellid].norm_prim_var[2],2)))*powf((convol[neighbourcellid].norm_prim_var[1]-sqrt(3.0*ncell_theta)),2);
            negative_flux[3] = temp3 + temp4;
        }

        if (ncell_norm_mach <= -1.0)
        {
            for (j=0;j<4;j++)
            {
                negative_flux[j] = convol[neighbourcellid].norm_analytic_flux[j];
            }
        }

        for (j=0;j<4;j++)
        {
            E[i].normal_interface_flux[j] = (positive_flux[j] + negative_flux[j])*(E[i].L.length);
        }

        E[i].cart_interface_flux();
    }
}
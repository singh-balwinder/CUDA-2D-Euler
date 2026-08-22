#include "definitions.h"

// Roe's approximate Riemann solver without entropy fix
void roe(vector <edge> &E, vector <cell> &convol)
{
    int i, j;
    int numEdge = E.size();
    double const gamma = 1.4;

    double totalEnthalpyP, totalEnthalpyN;
    double roe_averaged_nvelocity, roe_averaged_tvelocity, roe_averaged_sound_speed, roe_averaged_total_enthalpy;
    double eigenvalue[4];
    double K1[4], K2[4], K3[4], K4[4];
    double wave_strength[4];

    for (i=0;i<numEdge;i++)
    {
        int parentcellid = E[i].parentid;
        int neighbourcellid = E[i].neighbourid;

        convol[parentcellid].ncv(E[i].nx, E[i].ny);
        convol[neighbourcellid].ncv(E[i].nx, E[i].ny);

        totalEnthalpyP = (0.5*(powf(convol[parentcellid].norm_prim_var[1], 2) + powf(convol[parentcellid].norm_prim_var[2], 2))) + ((gamma/(gamma - 1.0))*(convol[parentcellid].norm_prim_var[3]/convol[parentcellid].norm_prim_var[0]));
        totalEnthalpyN = (0.5*(powf(convol[neighbourcellid].norm_prim_var[1], 2) + powf(convol[neighbourcellid].norm_prim_var[2], 2))) + ((gamma/(gamma - 1.0))*(convol[neighbourcellid].norm_prim_var[3]/convol[neighbourcellid].norm_prim_var[0]));

        roe_averaged_nvelocity = ((sqrt(convol[parentcellid].norm_prim_var[0])*convol[parentcellid].norm_prim_var[1]) + (sqrt(convol[neighbourcellid].norm_prim_var[0])*convol[neighbourcellid].norm_prim_var[1]))/(sqrt(convol[parentcellid].norm_prim_var[0]) + sqrt(convol[neighbourcellid].norm_prim_var[0]));
        roe_averaged_tvelocity = ((sqrt(convol[parentcellid].norm_prim_var[0])*convol[parentcellid].norm_prim_var[2]) + (sqrt(convol[neighbourcellid].norm_prim_var[0])*convol[neighbourcellid].norm_prim_var[2]))/(sqrt(convol[parentcellid].norm_prim_var[0]) + sqrt(convol[neighbourcellid].norm_prim_var[0]));
        roe_averaged_total_enthalpy = ((sqrt(convol[parentcellid].norm_prim_var[0])*totalEnthalpyP) + (sqrt(convol[neighbourcellid].norm_prim_var[0])*totalEnthalpyN))/(sqrt(convol[parentcellid].norm_prim_var[0]) + sqrt(convol[neighbourcellid].norm_prim_var[0]));
        roe_averaged_sound_speed = sqrt((gamma - 1.0)*(roe_averaged_total_enthalpy - (0.5*(powf(roe_averaged_nvelocity, 2) + powf(roe_averaged_tvelocity, 2)))));

        eigenvalue[0] = roe_averaged_nvelocity - roe_averaged_sound_speed;
        eigenvalue[1] = roe_averaged_nvelocity;
        eigenvalue[2] = roe_averaged_nvelocity;
        eigenvalue[3] = roe_averaged_nvelocity + roe_averaged_sound_speed;

        K1[0] = 1.0;
        K1[1] = eigenvalue[0];
        K1[2] = roe_averaged_tvelocity;
        K1[3] = (roe_averaged_total_enthalpy - (roe_averaged_nvelocity*roe_averaged_sound_speed));

        K2[0] = 1.0;
        K2[1] = eigenvalue[1];
        K2[2] = roe_averaged_tvelocity;
        K2[3] = 0.5*(powf(roe_averaged_nvelocity, 2) + powf(roe_averaged_tvelocity, 2));

        K3[0] = 0.0;
        K3[1] = 0.0;
        K3[2] = 1.0;
        K3[3] = roe_averaged_tvelocity;

        K4[0] = 1.0;
        K4[1] = eigenvalue[3];
        K4[2] = roe_averaged_tvelocity;
        K4[3] = (roe_averaged_total_enthalpy + (roe_averaged_nvelocity*roe_averaged_sound_speed));

        double delta_u5_bar = (convol[neighbourcellid].norm_cons_var[3] - convol[parentcellid].norm_cons_var[3]) - (((convol[neighbourcellid].norm_cons_var[2] - convol[parentcellid].norm_cons_var[2]) - (roe_averaged_tvelocity*(convol[neighbourcellid].norm_cons_var[0] - convol[parentcellid].norm_cons_var[0])))*roe_averaged_tvelocity);
        wave_strength[1] = (gamma - 1.0)*(((convol[neighbourcellid].norm_cons_var[0] - convol[parentcellid].norm_cons_var[0])*(roe_averaged_total_enthalpy - powf(roe_averaged_nvelocity, 2))) + (roe_averaged_nvelocity*(convol[neighbourcellid].norm_cons_var[1] - convol[parentcellid].norm_cons_var[1])) - delta_u5_bar)/powf(roe_averaged_sound_speed, 2);
        wave_strength[0] = ((eigenvalue[3]*(convol[neighbourcellid].norm_cons_var[0] - convol[parentcellid].norm_cons_var[0])) - (convol[neighbourcellid].norm_cons_var[1] - convol[parentcellid].norm_cons_var[1]) - (roe_averaged_sound_speed*wave_strength[1]))/(2.0*roe_averaged_sound_speed);
        wave_strength[2] = (convol[neighbourcellid].norm_cons_var[2] - convol[parentcellid].norm_cons_var[2]) - (roe_averaged_tvelocity*(convol[neighbourcellid].norm_cons_var[0] - convol[parentcellid].norm_cons_var[0]));
        wave_strength[3] = (convol[neighbourcellid].norm_cons_var[0] - convol[parentcellid].norm_cons_var[0]) - (wave_strength[0] + wave_strength[1]);

        for (j=0;j<4;j++)
        {
            double avg_flux = 0.5*(convol[parentcellid].norm_analytic_flux[j] + convol[neighbourcellid].norm_analytic_flux[j]);
            double diff_flux = 0.5*((wave_strength[0]*fabs(eigenvalue[0])*K1[j]) + (wave_strength[1]*fabs(eigenvalue[1])*K2[j]) + (wave_strength[2]*fabs(eigenvalue[2])*K3[j]) + (wave_strength[3]*fabs(eigenvalue[3])*K4[j]));
            E[i].normal_interface_flux[j] = (avg_flux - diff_flux)*E[i].L.length;
        }

        E[i].cart_interface_flux();
    }
}
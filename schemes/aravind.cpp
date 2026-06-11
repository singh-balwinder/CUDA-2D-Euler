#include "definitions.h"

// Aravind scheme
void aravind(vector <edge> &E, vector <cell> &convol)
{
    int i, j;
    int numEdge=E.size();
    double const gamma = 1.4;
    vector <double> positive_flux(4);
    vector <double> negative_flux(4);

    double soundspeedL, soundspeedR;
    double lambdaplusL, lambdaminusL;
    double lambdaplusR, lambdaminusR;
    double del_lambdaL, del_lambdaR;

    for (i=0;i<numEdge;i++)
    {
        double parentcellid = E[i].parentid;
        double neighbourcellid = E[i].neighbourid;

        convol[parentcellid].ncv(E[i].nx, E[i].ny);
        convol[neighbourcellid].ncv(E[i].nx, E[i].ny);

        soundspeedL = sqrt((gamma*convol[parentcellid].norm_prim_var[3])/convol[parentcellid].norm_prim_var[0]);
        soundspeedR = sqrt((gamma*convol[neighbourcellid].norm_prim_var[3])/convol[neighbourcellid].norm_prim_var[0]);

        lambdaplusL = powf(((convol[parentcellid].norm_prim_var[1] + soundspeedL) + fabs(convol[parentcellid].norm_prim_var[1] + soundspeedL)),2)/(4.0*soundspeedL);
        lambdaminusL = -powf(((convol[parentcellid].norm_prim_var[1] - soundspeedL) - fabs(convol[parentcellid].norm_prim_var[1] - soundspeedL)),2)/(4.0*soundspeedL);
        lambdaplusR = powf(((convol[neighbourcellid].norm_prim_var[1] + soundspeedR) + fabs(convol[neighbourcellid].norm_prim_var[1] + soundspeedR)),2)/(4.0*soundspeedR);
        lambdaminusR = -powf(((convol[neighbourcellid].norm_prim_var[1] - soundspeedR) - fabs(convol[neighbourcellid].norm_prim_var[1] - soundspeedR)),2)/(4.0*soundspeedR);

        /*
        lambdaplusL = ((E[i].norm_re_prim_var_par[1] + soundspeedL) + fabs(E[i].norm_re_prim_var_par[1] + soundspeedL))/2.0;
        lambdaminusL = ((E[i].norm_re_prim_var_par[1] - soundspeedL) - fabs(E[i].norm_re_prim_var_par[1] - soundspeedL))/2.0;
        lambdaplusR = ((E[i].norm_re_prim_var_neigh[1] + soundspeedR) + fabs(E[i].norm_re_prim_var_neigh[1] + soundspeedR))/2.0;
        lambdaminusR = ((E[i].norm_re_prim_var_neigh[1] - soundspeedR) - fabs(E[i].norm_re_prim_var_neigh[1] - soundspeedR))/2.0;
        */

        del_lambdaL = lambdaplusL - lambdaminusL;
        del_lambdaR = lambdaplusR - lambdaminusR;

        for (j=0;j<4;j++)
        {
            double positive_flux = (lambdaplusL/del_lambdaL)*(convol[parentcellid].norm_analytic_flux[j] - (lambdaminusL*convol[parentcellid].norm_cons_var[j]));
            double negative_flux = (lambdaminusR/del_lambdaR)*((lambdaplusR*convol[neighbourcellid].norm_cons_var[j]) - convol[neighbourcellid].norm_analytic_flux[j]);
            E[i].normal_interface_flux[j]=(positive_flux + negative_flux)*(E[i].L.length);
        }

        E[i].cart_interface_flux();
    }
}
#include "definitions.h"

void tempupdate(vector <cell> &convol, vector <edge> &E, double &del_t, bool &converged)
{
    int i, j;
    int numCell=convol.size();
    int numEdge=E.size();
    int numGhostcell = 0;
    int numInnercell;
    int edgeid;
    double massResidue, xmomResidue;
    double ymomResidue, energyResidue;
    double updated_cons_var[4];
    double L2norm[4];

    // number of ghost cells and inner cells
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
    numInnercell = numCell - numGhostcell;

    for (i=0;i<numCell;i++)
    {
        massResidue=0.0;
        xmomResidue=0.0;
        ymomResidue=0.0;
        energyResidue=0.0;
        if (convol[i].cellType == "GHOST")
        {
            continue;
        }
        else
        {
            // l1
            edgeid=convol[i].quad.l1.lineid;
            if (E[edgeid].parentid == convol[i].cellid)
            {
                massResidue=massResidue+E[edgeid].cartesian_interface_flux[0];
                xmomResidue=xmomResidue+E[edgeid].cartesian_interface_flux[1];
                ymomResidue=ymomResidue+E[edgeid].cartesian_interface_flux[2];
                energyResidue=energyResidue+E[edgeid].cartesian_interface_flux[3];
            }
            else
            {
                massResidue=massResidue-E[edgeid].cartesian_interface_flux[0];
                xmomResidue=xmomResidue-E[edgeid].cartesian_interface_flux[1];
                ymomResidue=ymomResidue-E[edgeid].cartesian_interface_flux[2];
                energyResidue=energyResidue-E[edgeid].cartesian_interface_flux[3];
            }

            // l2
            edgeid=convol[i].quad.l2.lineid;
            if (E[edgeid].parentid == convol[i].cellid)
            {
                massResidue=massResidue+E[edgeid].cartesian_interface_flux[0];
                xmomResidue=xmomResidue+E[edgeid].cartesian_interface_flux[1];
                ymomResidue=ymomResidue+E[edgeid].cartesian_interface_flux[2];
                energyResidue=energyResidue+E[edgeid].cartesian_interface_flux[3];
            }
            else
            {
                massResidue=massResidue-E[edgeid].cartesian_interface_flux[0];
                xmomResidue=xmomResidue-E[edgeid].cartesian_interface_flux[1];
                ymomResidue=ymomResidue-E[edgeid].cartesian_interface_flux[2];
                energyResidue=energyResidue-E[edgeid].cartesian_interface_flux[3];
            }

            // l3
            edgeid=convol[i].quad.l3.lineid;
            if (E[edgeid].parentid == convol[i].cellid)
            {
                massResidue=massResidue+E[edgeid].cartesian_interface_flux[0];
                xmomResidue=xmomResidue+E[edgeid].cartesian_interface_flux[1];
                ymomResidue=ymomResidue+E[edgeid].cartesian_interface_flux[2];
                energyResidue=energyResidue+E[edgeid].cartesian_interface_flux[3];
            }
            else
            {
                massResidue=massResidue-E[edgeid].cartesian_interface_flux[0];
                xmomResidue=xmomResidue-E[edgeid].cartesian_interface_flux[1];
                ymomResidue=ymomResidue-E[edgeid].cartesian_interface_flux[2];
                energyResidue=energyResidue-E[edgeid].cartesian_interface_flux[3];
            }

            // l4
            edgeid=convol[i].quad.l4.lineid;
            if (E[edgeid].parentid == convol[i].cellid)
            {
                massResidue=massResidue+E[edgeid].cartesian_interface_flux[0];
                xmomResidue=xmomResidue+E[edgeid].cartesian_interface_flux[1];
                ymomResidue=ymomResidue+E[edgeid].cartesian_interface_flux[2];
                energyResidue=energyResidue+E[edgeid].cartesian_interface_flux[3];
            }
            else
            {
                massResidue=massResidue-E[edgeid].cartesian_interface_flux[0];
                xmomResidue=xmomResidue-E[edgeid].cartesian_interface_flux[1];
                ymomResidue=ymomResidue-E[edgeid].cartesian_interface_flux[2];
                energyResidue=energyResidue-E[edgeid].cartesian_interface_flux[3];
            }
        }

        updated_cons_var[0]=convol[i].cons_var[0]-((del_t/convol[i].getArea())*(massResidue));
        updated_cons_var[1]=convol[i].cons_var[1]-((del_t/convol[i].getArea())*(xmomResidue));
        updated_cons_var[2]=convol[i].cons_var[2]-((del_t/convol[i].getArea())*(ymomResidue));
        updated_cons_var[3]=convol[i].cons_var[3]-((del_t/convol[i].getArea()))*(energyResidue);

        for (j=0;j<4;j++)
        {
            convol[i].error[j]=fabs(updated_cons_var[j])-fabs(convol[i].cons_var[j]);
        }

        for (j=0;j<4;j++)
        {
            convol[i].cons_var[j] = updated_cons_var[j];
        }
        convol[i].convert_c2p();
    }

    for (j=0;j<4;j++)
    {
        L2norm[j]=0.0;
        for (i=0;i<numCell;i++)
        {
            if (convol[i].cellType == "GHOST")
            {
                continue;
            }
            else
            {
                L2norm[j] = L2norm[j] + powf(convol[i].error[j],2);
            }
        }
        L2norm[j]=sqrt(L2norm[j])/numInnercell;
    }

    if ((L2norm[0]+L2norm[1]+L2norm[2]+L2norm[3])<=powf(10,-8))
    {
        converged = true;
    }

    
    // cout << convol[0].error[0] << endl;
    // cout << converged << endl;

    // cout << numInnercell << endl;
    
}
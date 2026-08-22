#include "definitions.h"

// kinked mach stem test case
void testcase2(vector <cell> &convol)
{
    int i;
    int numCell=convol.size();

    for (i=0;i<numCell;i++)
    {
        if (convol[i].cellType == "GHOST")
        {
            continue;
        }
        else
        {
            if (convol[i].quad.centroidy-(sqrt(3.0)*convol[i].quad.centroidx)+(1.0/(2.0*sqrt(3.0))) < 0.0)
            {
                convol[i].set(1.4, 0.0, 0.0, 1.0);
            }
            else
            {
                double den, ux, uy, pres;
                den = 7.999;
                ux = 7.14461;
                uy = -4.12494;
                pres = 116.5;
                convol[i].set(den, ux, uy, pres);
            }
        }
    }

    // print statements
    /*
    for (i=0;i<numCell;i++)
    {
        if (convol[i].cellType == "GHOST")
        {
            continue;
        }
        else
        {
            cout << convol[i].cons_var[1] << endl;
        }
    }
    */

    /*
    cout << convol[1].quad.l1.nx << "  " << convol[1].quad.l1.ny << endl;
    cout << convol[1].quad.l2.nx << "  " << convol[1].quad.l2.ny << endl;
    cout << convol[1].quad.l3.nx << "  " << convol[1].quad.l3.ny << endl;
    cout << convol[1].quad.l4.nx << "  " << convol[1].quad.l4.ny << endl;
    */
}
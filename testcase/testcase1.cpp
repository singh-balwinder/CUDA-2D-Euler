#include "definitions.h"

/*
// oblique shock test case
void testcase1(vector <cell> &convol)
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
            convol[i].set(1.0, 2.9, 0.0, (1/1.4));
        }
    }
}
*/

/*
// An oblique shock in the shock wave-boundary layer interaction
void testcase1(vector <cell> &convol)
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
            convol[i].set(1.0, 1.0, 0.0, 0.154);
        }
    }
}
*/

/*
// horizontal flow test case
void testcase1(vector <cell> &convol)
{
    int i;
    int numCell = convol.size();

    for (i=0;i<numCell;i++)
    {
        if (convol[i].cellType == "GHOST")
        {
            continue;
        }
        else
        {
            if (convol[i].quad.centroidy < 0.5)
            {
                convol[i].set(1.4, 2.0, 0.0, 1.0);
            }
            else
            {
                convol[i].set(1.4, 3.0, 0.0, 1.0);
            }
        }
    }
}
*/

/*
// supersonic flow over a ramp test case : M = 2.0
// supersonic flow over forward-facing step : M = 3.0
void testcase1(vector <cell> &convol)
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
            convol[i].set(1.4, 3.0, 0.0, 1.0);
        }
    }
}
*/

/*
void testcase1(vector <cell> &convol)
// odd-even decoupling test case
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
            convol[i].set(1.4, 0.0, 0.0, 1.0);
        }
    }
}
*/

/*
void testcase1(vector <cell> &convol)
// carbuncle
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
            convol[i].set(1.4, 20.0, 0.0, 1.0);
        }
    }
}
*/

/*
// 2-D enclosure test case (taken from Bernard parent 2012 JCP paper)
void testcase1(vector <cell> &convol)
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
            if (convol[i].quad.centroidx <= 0.0)
            {
                convol[i].set(0.999, 1183.80, -355.14, 10000.0);
            }
            else
            {
                convol[i].set(0.999, 1183.80, 236.76, 10000.0);
            }
        }
    }
}
*/

// shock diffraction test case
void testcase1(vector <cell> &convol)
{
    int i;
    int numCell=convol.size();

    for (i=0; i<numCell; i++)
    {
        if (convol[i].cellType == "GHOST")
        {
            continue;
        }
        else
        {
            if (convol[i].quad.centroidx < 0.05)
            {
                convol[i].set(7.041132, 4.07795, 0.0, 30.05945);
            }
            else
            {
                convol[i].set(1.4, 0.0, 0.0, 1.0);
            }
        }
    }
}
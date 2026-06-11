#include "definitions.h"
/*
    GMSH physical curve numbering
    1. WALL: 5
    2. OUTLET: 6
    3. INLET(S): 7/8
    4. DYNAMIC: 4
*/

void read_vtk(vector <point> &p, vector <edge> &E, vector <cell> &convol, int &numEdges, int &numIntcells)
{
    vector <line> l;
    string str1, str2, str3, str4;
    double x, y, z;
    int numPoints=1e7;
    int i, j;
    // to read polygons
    int a, b, c, d, e;
    // to read edges
    int f, g, h;
    // to read boundaries
    int m;
    int pointIter = 0;
    int cellIter = 0;
    int lineIter = 0;
    int edgeIter = 0;
    int loopIter1 = 0;
    int loopIter2 = 0;

    /* open the .vtk file */
    ifstream vtkFile("bck_step_200_200.vtk");

    /* check to see .vtk file is open or not */
    if (!vtkFile.is_open())
    {
        cerr << "Failed to open the .vtk file" << endl;
        // return 1;
    }

    bool readCoordinates = false;

    /* POINTS data */
    while (getline(vtkFile, str3))
    {
        loopIter1=loopIter1+1;
        if (str3.find("POINTS") != string::npos)
        {
            readCoordinates = true;
            istringstream iss(str3);
            iss >> str1 >> numPoints >> str2;
            continue;
        }
        if (readCoordinates)
        {
            istringstream iss(str3);
            iss >> x >> y >> z;
            /*
            // odd-even decoupling test case grid tweak
            if (fabs(y-10.0) <= powf(10,-8))
            {
                int int_ver = x;
                if (int_ver%2 == 0)
                {
                    y = y + powf(10,-3);
                }
                else
                {
                    y = y - powf(10,-3);
                }
            }
            */
            p.push_back(point(x, y, z, pointIter));
            pointIter=pointIter+1;
        }

        if (pointIter>=numPoints)
        {
            break;
        }

        vtkFile >> ws;
    }

    /* CELLS data */
    while (getline(vtkFile, str3))
    {
        loopIter2=loopIter2+1;
        istringstream iss1(str3);
        iss1 >> str4;
        if (str4 != "CELL_TYPES")
        {
            istringstream iss2(str3);
            iss2 >> f;
            if (f == 2)
            {
                iss1 >> g >> h;
                l.push_back(line(p[g], p[h], lineIter));
                lineIter = lineIter + 1;
            }

            if (f == 4)
            {
                a=f;
                iss1 >> b >> c >> d >> e;
                // cout << b << "  " << c << "  " << d << "  " << e << endl;
                convol.push_back(cell(p[b], p[c], p[d], p[e], cellIter, a));
                cellIter = cellIter+1;
            }
        }
        else
        {
            break;
        }

        vtkFile >> ws;
    }

    /* Boundary Type */
    bool readBoundary = false;
    while (getline(vtkFile, str3))
    {
        if (str3.find("LOOKUP_TABLE") != string::npos)
        {
            readBoundary = true;
            continue;
        }
        if (readBoundary)
        {
            istringstream iss3(str3);
            iss3 >> m;
            if (m == -1)
            {
                continue;
            }
            if (m == 5)
            {
                E.push_back(edge(l[edgeIter], edgeIter, "BOUNDARY"));
                E[edgeIter].boundaryType = "WALL";
                edgeIter=edgeIter+1;
            }
            if (m == 6)
            {
                E.push_back(edge(l[edgeIter], edgeIter, "BOUNDARY"));
                E[edgeIter].boundaryType = "OUTLET";
                edgeIter=edgeIter+1;
            }
            if (m == 7)
            {
                E.push_back(edge(l[edgeIter], edgeIter, "BOUNDARY"));
                E[edgeIter].boundaryType = "INLET";
                edgeIter=edgeIter+1;
            }
            /*
            if (m == 8)
            {
                E.push_back(edge(l[edgeIter], edgeIter, "BOUNDARY"));
                E[edgeIter].boundaryType = "INLET2";
                edgeIter=edgeIter+1;
            }
            if (m == 4)
            {
                E.push_back(edge(l[edgeIter], edgeIter, "BOUNDARY"));
                E[edgeIter].boundaryType = "DYNAMIC";
                edgeIter=edgeIter+1;
            }
            */
        }

        vtkFile >> ws;
    }

    /* close the .vtk file */
    vtkFile.close();

    numEdges = E.size();
    numIntcells = convol.size();

    /* Unique edge identification */
    for (i=0; i<numIntcells; i++)
    {
        convol[i].area = convol[i].getArea();

        /* procedure for l1 */
        bool flag1 = false;
        for (j=0; j<numEdges; j++)
        {
            if (convol[i].quad.l1.islinesame(E[j].L))
            {
                flag1 = true;
                if (E[j].edgeType == "BOUNDARY")
                {
                    convol[i].quad.l1.lineid=E[j].edgeid;
                    E[j].parentid=i;
                }
                if (E[j].edgeType == "INNER")
                {
                    convol[i].quad.l1.lineid=E[j].edgeid;
                    E[j].neighbourid=i;
                }
            }
        }
        if (flag1 == false)
        {
            E.push_back(edge(convol[i].quad.l1, edgeIter, "INNER"));
            E[edgeIter].parentid=i;
            convol[i].quad.l1.lineid=E[edgeIter].edgeid;
            edgeIter=edgeIter+1;
            numEdges = E.size();
        }

        /* procedure for l2 */
        bool flag2 = false;
        for (j=0; j<numEdges; j++)
        {
            if (convol[i].quad.l2.islinesame(E[j].L))
            {
                flag2 = true;
                if (E[j].edgeType == "BOUNDARY")
                {
                    convol[i].quad.l2.lineid=E[j].edgeid;
                    E[j].parentid=i;
                }
                if (E[j].edgeType == "INNER")
                {
                    convol[i].quad.l2.lineid=E[j].edgeid;
                    E[j].neighbourid=i;
                }
            }
        }
        if (flag2 == false)
        {
            E.push_back(edge(convol[i].quad.l2, edgeIter, "INNER"));
            E[edgeIter].parentid=i;
            convol[i].quad.l2.lineid=E[edgeIter].edgeid;
            edgeIter=edgeIter+1;
            numEdges = E.size();
        }

        /* procedure for l3 */
        bool flag3 = false;
        for (j=0; j<numEdges; j++)
        {
            if (convol[i].quad.l3.islinesame(E[j].L))
            {
                flag3 = true;
                if (E[j].edgeType == "BOUNDARY")
                {
                    convol[i].quad.l3.lineid=E[j].edgeid;
                    E[j].parentid=i;
                }
                if (E[j].edgeType == "INNER")
                {
                    convol[i].quad.l3.lineid=E[j].edgeid;
                    E[j].neighbourid=i;
                }
            }
        }
        if (flag3 == false)
        {
            E.push_back(edge(convol[i].quad.l3, edgeIter, "INNER"));
            E[edgeIter].parentid=i;
            convol[i].quad.l3.lineid=E[edgeIter].edgeid;
            edgeIter=edgeIter+1;
            numEdges = E.size();
        }

        /* procedure for l4 */
        bool flag4 = false;
        for (j=0; j<numEdges; j++)
        {
            if (convol[i].quad.l4.islinesame(E[j].L))
            {
                flag4 = true;
                if (E[j].edgeType == "BOUNDARY")
                {
                    convol[i].quad.l4.lineid=E[j].edgeid;
                    E[j].parentid=i;
                }
                if (E[j].edgeType == "INNER")
                {
                    convol[i].quad.l4.lineid=E[j].edgeid;
                    E[j].neighbourid=i;
                }
            }
        }
        if (flag4 == false)
        {
            E.push_back(edge(convol[i].quad.l4, edgeIter, "INNER"));
            E[edgeIter].parentid=i;
            convol[i].quad.l4.lineid=E[edgeIter].edgeid;
            edgeIter=edgeIter+1;
            numEdges = E.size();
        }
    }

    /* print statements
    cout << loopIter1 << endl;
    cout << numPoints << endl;

    for (i=0;i<numPoints;i++)
    {
        cout << p[i].pointid << "  " << p[i].x << "  " << p[i].y << endl;
    }

    cout << loopIter2 << endl;
    
    for (i=0;i<numCell;i++)
    {
        cout << i << "  " << convol[i].quad.p1.x << "  " << convol[i].quad.p1.y << endl;
    }

    cout << "Line data" << endl;

    for (i=0;i<edgeIter;i++)
    {
        cout << i << "  " << E[i].L.p1.x << "  " << E[i].L.p1.y << "  " << E[i].L.p2.x << "  " << E[i].L.p2.y << endl;
    }

    for (i=0;i<edgeIter;i++)
    {
        cout << i << "  " << E[i].parentid << "  " << E[i].neighbourid << endl;
    }

    cout << numCell << endl;
    cout << numEdge << endl;

    cout << "Polygon data" << endl;

    cout << convol[0].quad.l1.lineid << "  " << convol[0].quad.l2.lineid << "  " << convol[0].quad.l3.lineid << "  " << convol[0].quad.l4.lineid << endl;
    cout << convol[1].quad.l1.lineid << "  " << convol[1].quad.l2.lineid << "  " << convol[1].quad.l3.lineid << "  " << convol[1].quad.l4.lineid << endl;

    cout << "boundary type" << endl;
    for (i=0;i<10;i++)
    {
        cout << i << "  " << E[i].boundaryType << endl;
    }

    for (i=0;i<numEdge;i++)
    {
        cout << i << "  " << "nx = " << E[i].nx << "  " << "ny = " << E[i].ny << endl;
    }
    */

    /*
    cout << convol[1].quad.l1.nx << "  " << convol[1].quad.l1.ny << endl;
    cout << convol[1].quad.l2.nx << "  " << convol[1].quad.l2.ny << endl;
    cout << convol[1].quad.l3.nx << "  " << convol[1].quad.l3.ny << endl;
    cout << convol[1].quad.l4.nx << "  " << convol[1].quad.l4.ny << endl;
    */
    
    /*
    for (j=0;j<convol.size();j++)
    {
        cout << j << endl;
        cout << convol[j].quad.l4.nx << " " << convol[j].quad.l4.ny << endl;
    }
    */

    /*
    cout << convol[10000].quad.p1.x << " " << convol[10000].quad.p2.x << " " << convol[10000].quad.p3.x << " " << convol[10000].quad.p4.x << endl;
    */
}
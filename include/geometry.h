using namespace std;

class point
{
    public:
    int pointid;
    double x, y, z;

    // default constructor
    point()
    {
        
    }

    // parametrized constructor
    point(double a, double b, double c, int id)
    {
        pointid=id;
        x=a;
        y=b;
        z=c;
    }

    // methods
    void set(double a, double b, double c, int id)
    {
        pointid=id;
        x=a;
        y=b;
        z=c;
    }

    bool ispointsame(point pt)
    {
        double tol = 1e-9;
        if (fabs(x - pt.x) <= tol && fabs(y - pt.y) <= tol)
        {
            return true;
        }
        else
        {
            return false;
        }
    }
};

class line

{
    public:
    int lineid;
    point p1, p2;
    point midpoint;
    double length;
    double nx, ny;

    // default constructor
    line()
    {

    }

    // parametrized constructor 1
    line (point pt1, point pt2, int id)
    {
        lineid=id;
        p1=pt1;
        p2=pt2;
        midpoint.x = (pt1.x+pt2.x)/2.0;
        midpoint.y = (pt1.y+pt2.y)/2.0;
        midpoint.z = (pt1.z+pt2.z)/2.0;
        length=getLength();
        nx=(p2.y-p1.y)/length;
        ny=-(p2.x-p1.x)/length;
        lineTolerance();
    }

    // parametrized constructor 2
    line (point pt1, point pt2)
    {
        p1=pt1;
        p2=pt2;
        length=getLength();
        nx=(p2.y-p1.y)/length;
        ny=-(p2.x-p1.x)/length;
        lineTolerance();
    }

    // length calculation
    double getLength()
    {
        length=sqrt(powf((p2.x-p1.x),2)+powf((p2.y-p1.y),2));
        return length;
    }

    bool islinesame(line l1)
    {
        if (p1.ispointsame(l1.p1) && p2.ispointsame(l1.p2))
        {
            return true;
        }
        else if (p1.ispointsame(l1.p2) && p2.ispointsame(l1.p1))
        {
            return true;
        }
        else 
        {
            return false;
        }
    }

    void lineTolerance()
    {
        double tolerance = powf(10,-9);
        if (fabs(nx) <= tolerance)
        {
            nx=0;
        }
        if (fabs(ny) <= tolerance)
        {
            ny=0;
        }
    }

};

class edge
{
    public:
    int edgeid;
    int parentid, neighbourid;
    int numVertices;
    line L;
    string edgeType;
    string boundaryType;
    double nx, ny;
    double tx, ty;
    double normal_interface_flux[4];
    double cartesian_interface_flux[4];
    double alpha;

    // deafult constructor
    edge()
    {

    }
    
    // parametrized constructor
    edge(line l1, int id, string type)
    {
        edgeid=id;
        L=l1;
        edgeType=type;
        nx=(L.p2.y-L.p1.y)/L.length;
        ny=-(L.p2.x-L.p1.x)/L.length;
        edgeTolerance();
        tx=-ny;
        ty=nx;
    }

    void edgeTolerance()
    {
        double tolerance=powf(10,-10);
        if (fabs(nx) <= tolerance)
        {
            nx=0;
        }
        if (fabs(ny) <= tolerance)
        {
            ny=0;
        }
    }

    void cart_interface_flux()
    {
        cartesian_interface_flux[0]=normal_interface_flux[0];
        cartesian_interface_flux[1]=(normal_interface_flux[1]*nx)-(normal_interface_flux[2]*ny);
        cartesian_interface_flux[2]=(normal_interface_flux[1]*ny)+(normal_interface_flux[2]*nx);
        cartesian_interface_flux[3]=normal_interface_flux[3];
        fluxTolerance();
    }

    void fluxTolerance()
    {
        double fluxtolerance = powf(10,-10);
        int j;
        for (j=0;j<4;j++)
        {
            if (fabs(cartesian_interface_flux[j]) <= fluxtolerance)
            {
                cartesian_interface_flux[j] = 0.0;
            }
        }
    }
};

class polygon
{
    public:
    int polyid;
    double centroidx, centroidy;
    int numVertices;
    point p1, p2, p3, p4;
    line l1, l2, l3, l4;
    double v1, v2, v3, v4, v5, v6, v7, v8;

    // centroid calculation
    double getCentroidx()
    {
        centroidx=(p1.x+p2.x+p3.x+p4.x)/4.0;
        return centroidx;
    }

    double getCentroidy()
    {
        centroidy=(p1.y+p2.y+p3.y+p4.y)/4.0;
        return centroidy;
    }

    // default constructor
    polygon()
    {

    }

    // parametrized constructor
    polygon(point pt1, point pt2, point pt3, point pt4, int id, int num)
    {
        p1=pt1;
        p2=pt2;
        p3=pt3;
        p4=pt4;
        polyid=id;
        numVertices=num;


        // line details
        l1=line(p1, p2, 1);
        l2=line(p2, p3, 2);
        l3=line(p3, p4, 3);
        l4=line(p4, p1, 4);

        centroidx=getCentroidx();
        centroidy=getCentroidy();

        v1 = l1.midpoint.x - centroidx;
        v2 = l1.midpoint.y - centroidy;
        v3 = l2.midpoint.x - centroidx;
        v4 = l2.midpoint.y - centroidy;
        v5 = l3.midpoint.x - centroidx;
        v6 = l3.midpoint.y - centroidy;
        v7 = l4.midpoint.x - centroidx;
        v8 = l4.midpoint.y - centroidy;

        vTolerance();
    }

    void vTolerance()
    {
        double tolerance=powf(10,-12);
        if (fabs(v1) <= tolerance)
        {
            v1=0;
        }
        if (fabs(v2) <= tolerance)
        {
            v2=0;
        }
        if (fabs(v3) <= tolerance)
        {
            v3=0;
        }
        if (fabs(v4) <= tolerance)
        {
            v4=0;
        }
        if (fabs(v5) <= tolerance)
        {
            v5=0;
        }
        if (fabs(v6) <= tolerance)
        {
            v6=0;
        }
        if (fabs(v7) <= tolerance)
        {
            v7=0;
        }
        if (fabs(v8) <= tolerance)
        {
            v8=0;
        }
    }
};

struct host_edgeSoA
{
    vector <int> pCell, nCell;
    vector <double> nx, ny;
    vector <double> length;

    vector <double> cart_flux0, cart_flux1, cart_flux2, cart_flux3;

    host_edgeSoA(int n)
    {
        pCell.resize(n);
        nCell.resize(n);
        nx.resize(n);
        ny.resize(n);
        length.resize(n);

        cart_flux0.resize(n);
        cart_flux1.resize(n);
        cart_flux2.resize(n);
        cart_flux3.resize(n);
    }
};

struct host_cellSoA
{
    vector <double> area;
    vector <double> den, xMom, yMom, energy;

    host_cellSoA(int n)
    {
        area.resize(n);

        den.resize(n);
        xMom.resize(n);
        yMom.resize(n);
        energy.resize(n);
    }
};

struct device_edgeSoA
{
    int* pCell;
    int* nCell;

    double* nx;
    double* ny;

    double* length;

    double* denFlux;
    double* xMomFlux;
    double* yMomFlux;
    double* energyFlux;
};

struct device_cellSoA
{
    double* area;

    double* den;
    double* xMom;
    double* yMom;
    double* energy;
};
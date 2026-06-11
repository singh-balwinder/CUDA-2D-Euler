using namespace std;

class cell
{
    public:
    // general variables
    int cellid;
    int numPolyVertices;
    string cellType;

    // geometric variables
    polygon quad;
    double area;

    // physical system : inviscid gas dynamics
    double cons_var[4], prim_var[4];
    double norm_analytic_flux[4];
    double norm_cons_var[4], norm_prim_var[4];
    double norm_eigenvalue[4];

    // for convergence criteria
    double error[4];

    // default constructor
    cell()
    {

    }

    // parametrized constructor 1
    cell(point pt1, point pt2, point pt3, point pt4, int id, int num)
    {
        quad=polygon(pt1, pt2, pt3, pt4, id, num);
        cellid=id;
        numPolyVertices=num;
    }

    // parametrized constructor 3
    cell(int id, string type)
    {
        cellid=id;
        cellType=type;
    }

    // methods
    void set(double density, double u, double v, double pressure)
    {
        prim_var[0] = density;
        prim_var[1] = u;
        prim_var[2] = v;
        prim_var[3] = pressure;

        convert_p2c();
    }

    void convert_p2c()
    {
        const double gamma=1.4;
        cons_var[0] = prim_var[0];
        cons_var[1] = (prim_var[0]*prim_var[1]);
        cons_var[2] = (prim_var[0]*prim_var[2]);
        cons_var[3] = (prim_var[3]/(gamma-1.0))+(prim_var[0]*powf(prim_var[1],2)/2.0)+(prim_var[0]*powf(prim_var[2],2)/2.0);
    }

    void convert_c2p()
    {
        const double gamma = 1.4;
        prim_var[0] = cons_var[0];
        prim_var[1] = cons_var[1]/cons_var[0];
        prim_var[2] = cons_var[2]/cons_var[0];
        prim_var[3] = (gamma-1)*(cons_var[3]-(powf(cons_var[1],2)/(2.0*cons_var[0]))-(powf(cons_var[2],2)/(2.0*cons_var[0])));
    }

    void ncv(double nx, double ny)
    {
        // primitive variables (normal to edge)
        const double gamma=1.4;
        norm_prim_var[0] = prim_var[0];
        // normal velocity : u_n
        norm_prim_var[1] = (nx*prim_var[1])+(ny*prim_var[2]);
        // tangential velocity : u_t
        norm_prim_var[2] = -(ny*prim_var[1])+(nx*prim_var[2]);
        norm_prim_var[3] = prim_var[3];

        // conserved variables (normal to edge)
        norm_cons_var[0] = norm_prim_var[0];
        norm_cons_var[1] = norm_prim_var[0]*norm_prim_var[1];
        norm_cons_var[2] = norm_prim_var[0]*norm_prim_var[2];
        norm_cons_var[3] = (norm_prim_var[3]/(gamma-1.0))+(0.5*norm_prim_var[0]*powf(norm_prim_var[1],2))+(0.5*norm_prim_var[0]*powf(norm_prim_var[2],2));

        // analytic flux (normal to edge)
        norm_analytic_flux[0] = norm_prim_var[0]*norm_prim_var[1];
        norm_analytic_flux[1] = (norm_prim_var[0]*powf(norm_prim_var[1],2))+(norm_prim_var[3]);
        norm_analytic_flux[2] = norm_prim_var[0]*norm_prim_var[1]*norm_prim_var[2];
        norm_analytic_flux[3] = (norm_prim_var[3]+norm_cons_var[3])*(norm_prim_var[1]);

        // eigenvalues (normal to edge)
        norm_eigenvalue[0] = norm_prim_var[1]+sqrt((gamma*norm_prim_var[3]/norm_prim_var[0]));
        norm_eigenvalue[1] = norm_prim_var[1];
        norm_eigenvalue[2] = norm_prim_var[1];
        norm_eigenvalue[3] = norm_prim_var[1]-sqrt((gamma*norm_prim_var[3]/norm_prim_var[0]));
    }

    double findMaxeigen()
    {
        double maxEigen=fabs(norm_eigenvalue[0]);
        for (int i = 1; i < 4; i++) 
        {
            if (fabs(norm_eigenvalue[i]) > maxEigen) 
            {
                maxEigen = fabs(norm_eigenvalue[i]);
            }
        }

        return maxEigen;
    }

    double getArea()
    {
        double cellArea=fabs(((quad.p3.x-quad.p1.x)*(quad.p4.y-quad.p2.y))-((quad.p4.x-quad.p2.x)*(quad.p3.y-quad.p1.y)));
        return cellArea/2.0;
    }
};
// Zeroth order Raviart-Thomas vector element on quads, nonconforming, 2D

// the nodal functionals on the edges E_i are
// N_i(v) = int_{E_i} v.n 
// we use a 2-point Gauss quadrature on each edge (which is exact for 
// polynomials up to order 3)
static double NF_N_Q_RT0_2D_Xi[8] =
{ -sqrt(1./3.), sqrt(1./3.), 1, 1,
  -sqrt(1./3.), sqrt(1./3.), -1, -1 };
static double NF_N_Q_RT0_2D_Eta[8] =
{-1, -1, -sqrt(1./3.), sqrt(1./3.),
  1, 1, -sqrt(1./3.), sqrt(1./3.) };
static double NF_N_Q_RT0_2D_T[2] = {-sqrt(1./3.), sqrt(1./3.)};

void NF_N_Q_RT0_2D_EvalAll(TCollection *Coll, TBaseCell *Cell,
                           const double *PointValues, double *Functionals)
{
  // on the reference cell [-1,1]^2
  if(Cell == nullptr)
  {
    Functionals[0] = -PointValues[8] - PointValues[9];
    Functionals[1] = PointValues[2] + PointValues[3];
    Functionals[2] = PointValues[12] + PointValues[13];
    Functionals[3] = -PointValues[6] - PointValues[7];
  }
  else // on a real cell
  {
    double x0, x1, x2, x3, y0, y1, y2, y3;
    #ifdef __2D__
    Cell->GetVertex(0)->GetCoords(x0, y0);
    Cell->GetVertex(1)->GetCoords(x1, y1);
    Cell->GetVertex(2)->GetCoords(x2, y2);
    Cell->GetVertex(3)->GetCoords(x3, y3);
    #else
    ErrThrow("NF_N_Q_RT0_2D_EvalAll not implemented in 3D");
    #endif
    // length of edge, and outer normal
    double nx, ny;
    
    // first edge:
    nx = y1 - y0;
    ny = x0 - x1;
    Functionals[0] = 0.5*(PointValues[0] + PointValues[1])*nx
                    +0.5*(PointValues[8] + PointValues[9])*ny;
    Functionals[0] *= Cell->GetNormalOrientation(0);
    
    // second edge:
    nx = y2 - y1;
    ny = x1 - x2;
    Functionals[1] = 0.5*(PointValues[2] + PointValues[3])*nx
                    +0.5*(PointValues[10]+ PointValues[11])*ny;
    Functionals[1] *= Cell->GetNormalOrientation(1);
    
    // third edge:
    nx = y3 - y2;
    ny = x2 - x3;
    Functionals[2] = 0.5*(PointValues[4] + PointValues[5])*nx
                    +0.5*(PointValues[12]+ PointValues[13])*ny;
    Functionals[2] *= Cell->GetNormalOrientation(2);
    
    nx = y0 - y3;
    ny = x3 - x0;
    Functionals[3] = 0.5*(PointValues[6] + PointValues[7])*nx
                    +0.5*(PointValues[14]+ PointValues[15])*ny;
    Functionals[3] *= Cell->GetNormalOrientation(3);
  }
}

void NF_N_Q_RT0_2D_EvalEdge(TCollection *Coll, TBaseCell *Cell, int Joint,
                            const double *PointValues, double *Functionals)
{
  // this is needed for setting boundary conditions
  #ifdef __2D__
  double x0,x1,y0,y1;
  Cell->GetVertex(Joint)->GetCoords(x0,y0);
  Cell->GetVertex((Joint+1)%4)->GetCoords(x1,y1);// 4=number of edges
  double l = sqrt((x0-x1)*(x0-x1) + (y0-y1)*(y0-y1)); // length of joint
  Functionals[0] = 0.5*(PointValues[0] + PointValues[1])*l;
  #endif
}

TNodalFunctional2D *NF_N_Q_RT0_2D_Obj = new TNodalFunctional2D
        (NF_N_Q_RT0_2D, 4, 1, 8, 2, NF_N_Q_RT0_2D_Xi, NF_N_Q_RT0_2D_Eta,
         NF_N_Q_RT0_2D_T, NF_N_Q_RT0_2D_EvalAll, NF_N_Q_RT0_2D_EvalEdge);

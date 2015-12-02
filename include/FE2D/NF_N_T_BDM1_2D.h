/*
    TNodalFunctional2D(NodalFunctional2D id,
                       int n_allfunctionals, int n_edgefunctionals,
                       int n_pointsall, int n_pointsedge,
                       double *xi, double *eta, double *t,
                       DoubleFunctVect *evalall,
                       DoubleFunctVect *evaledge);
  used function names in this file must be unique in the entire program. 
  -> use the identifier as prefix
*/

// Tschebyschow-points
static double NF_N_T_BDM1_2D_a = -1./(2.*sqrt(2.)) + 0.5;

static double NF_N_T_BDM1_2D_Xi[] = 
{ NF_N_T_BDM1_2D_a, 1-NF_N_T_BDM1_2D_a,
  1-NF_N_T_BDM1_2D_a, NF_N_T_BDM1_2D_a,
  0, 0 };
static double NF_N_T_BDM1_2D_Eta[] = 
{ 0,   0,
  NF_N_T_BDM1_2D_a, 1-NF_N_T_BDM1_2D_a,
  1-NF_N_T_BDM1_2D_a, NF_N_T_BDM1_2D_a };
// NOTE: If you want to use other evaluation points for degress of freedom on
// the edges of a cell, you also have to change basis functions in 
// BF_N_T_BDM1_2D.h
//static double NF_N_T_BDM1_2D_T[] = {-0.333333333333,0.3333333333333};// equidistant points
//static double NF_N_T_BDM1_2D_T[] = {-0.577350269189626,0.577350269189626};//Gauss-points
static double NF_N_T_BDM1_2D_T[] = {-0.707106781186547,0.707106781186547};//Tschebyscheff-points

void NF_N_T_BDM1_2D_EvalAll(TCollection *Coll, TBaseCell *Cell, double *PointValues,
                          double *Functionals)
{
  // on the reference triangle with points (0,0), (1,0), (0,1) 
  if(Cell == nullptr)
  {
    Functionals[0] = -PointValues[6];
    Functionals[1] = -PointValues[7];
    
    Functionals[2] = PointValues[2] + PointValues[8];
    Functionals[3] = PointValues[3] + PointValues[9];
    
    Functionals[4] = -PointValues[4];
    Functionals[5] = -PointValues[5];
  }
  else // on a real cell
  {
    double x0, x1, x2, y0, y1, y2;
    #ifdef __2D__
    Cell->GetVertex(0)->GetCoords(x0, y0);
    Cell->GetVertex(1)->GetCoords(x1, y1);
    Cell->GetVertex(2)->GetCoords(x2, y2);
    #else
    ErrThrow("NF_N_T_BDM1_2D_EvalAll not implemented in 3D");
    #endif
    // length of edge, and outer normal
    double nx, ny;
    
    // first edge:
    nx = y1 - y0;
    ny = x0 - x1;
    Functionals[0] = PointValues[0]*nx + PointValues[6]*ny;
    Functionals[1] = PointValues[1]*nx + PointValues[7]*ny;
    
    // second edge:
    nx = y2 - y1;
    ny = x1 - x2;
    Functionals[2] = PointValues[2]*nx + PointValues[8]*ny;
    Functionals[3] = PointValues[3]*nx + PointValues[9]*ny;
    
    // third edge:
    nx = y0 - y2;
    ny = x2 - x0;
    Functionals[4] = PointValues[4]*nx + PointValues[10]*ny;
    Functionals[5] = PointValues[5]*nx + PointValues[11]*ny;
  }
}

void NF_N_T_BDM1_2D_EvalEdge(TCollection *Coll, TBaseCell *Cell, int Joint, double *PointValues,
                           double *Functionals)
{
// this is needed for setting boundary conditions
  /* the functional
   * int_Joint v.n
   * will be multiplied by the length of the Joint (edge). Otherwise one would
   * ensure int_Joint v.n=PointValues[0]. 
   * Example: If you would like to have u.n=1, then without multiplying by 
   *          the edge length l would result in having int_Joint u.n=1 on each
   *          boundary edge. This would mean one gets u.n=1/l on that 
   *          boundary. To avoid this, we introduce the factor l here. 
   * However I am not sure if this causes trouble elsewhere later. 
   * Be carefull!
   *                                            Ulrich Wilbrandt, 11.05.2012
  */
  double l; // length of joint
  double x0,x1,y0,y1;
  #ifdef __2D__
  Cell->GetVertex(Joint)->GetCoords(x0,y0);
  Cell->GetVertex((Joint+1)%3)->GetCoords(x1,y1);// 3=number of edges
  #endif
  l = sqrt((x0-x1)*(x0-x1) + (y0-y1)*(y0-y1));
  
  Functionals[0] = PointValues[0]*l; 
  Functionals[1] = PointValues[1]*l; 
}

TNodalFunctional2D *NF_N_T_BDM1_2D_Obj = new TNodalFunctional2D
        (NF_N_T_BDM1_2D, 6, 2, 6, 2, NF_N_T_BDM1_2D_Xi, NF_N_T_BDM1_2D_Eta,
         NF_N_T_BDM1_2D_T, NF_N_T_BDM1_2D_EvalAll, NF_N_T_BDM1_2D_EvalEdge);

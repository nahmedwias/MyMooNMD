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
static double NF_N_Q_BDM2_2D_a = sqrt(3./4.); 

static double NF_N_Q_BDM2_2D_Xi[] = 
{ -NF_N_Q_BDM2_2D_a, 0, NF_N_Q_BDM2_2D_a,
  1,1,1,
  NF_N_Q_BDM2_2D_a, 0,-NF_N_Q_BDM2_2D_a,
  -1,-1,-1,
  0,0 };
static double NF_N_Q_BDM2_2D_Eta[] = 
{ -1,-1,-1,
  -NF_N_Q_BDM2_2D_a,0,NF_N_Q_BDM2_2D_a,
  1,1, 1,
  NF_N_Q_BDM2_2D_a, 0,-NF_N_Q_BDM2_2D_a,
  0,0 };
// NOTE: If you want to use other evaluation points for degress of freedom on
// the edges of a cell, you also have to change basis functions in 
// BF_N_Q_BDM2_2D.h
//static double NF_N_Q_BDM2_2D_T[] = {-0.5,0,0.5};// equidistant points
//static double NF_N_Q_BDM2_2D_T[] = {-0.774596669241483,0,0.774596669241483};//Gauss-points
static double NF_N_Q_BDM2_2D_T[] = {-0.866025403784439,0,0.866025403784439};//Tschebyscheff-points

void NF_N_Q_BDM2_2D_EvalAll(TCollection *Coll, TBaseCell *Cell, double *PointValues,
                          double *Functionals)
{
  // on the reference cell [-1,1]^2
  if(Cell == nullptr)
  {
    Functionals[0] = -2*PointValues[14];
    Functionals[1] = -2*PointValues[15];
    Functionals[2] = -2*PointValues[16];
    
    Functionals[3] = 2*PointValues[3];
    Functionals[4] = 2*PointValues[4];
    Functionals[5] = 2*PointValues[5];
    
    Functionals[6] = 2*PointValues[20];
    Functionals[7] = 2*PointValues[21];
    Functionals[8] = 2*PointValues[22];
    
    Functionals[9] = -2*PointValues[9];
    Functionals[10]= -2*PointValues[10];
    Functionals[11]= -2*PointValues[11];
    
    Functionals[12]= PointValues[12];
    Functionals[13]= PointValues[27];
  }
  else
  {
    if(Cell->GetShapeDesc()->GetType() == Quadrangle) 
    {
      // not affine reference transform
      ErrThrow("NF_N_Q_BDM2_2D_EvalAll not tested for non affine ",
               "reference transformations");
    }
    double x0, x1, x2, x3, y0, y1, y2, y3;
    #ifdef __2D__
    Cell->GetVertex(0)->GetCoords(x0, y0);
    Cell->GetVertex(1)->GetCoords(x1, y1);
    Cell->GetVertex(2)->GetCoords(x2, y2);
    Cell->GetVertex(3)->GetCoords(x3, y3);
    #else
    ErrThrow("NF_N_Q_BDM2_2D_EvalAll not implemented in 3D");
    #endif
    
    // outer normal
    double nx, ny;
    
    // first edge:
    nx = y1 - y0;
    ny = x0 - x1;
    Functionals[0] = PointValues[0]*nx + PointValues[14]*ny;
    Functionals[1] = PointValues[1]*nx + PointValues[15]*ny;
    Functionals[2] = PointValues[2]*nx + PointValues[16]*ny;
    
    // second edge:
    nx = y2 - y1;
    ny = x1 - x2;
    Functionals[3] = PointValues[3]*nx + PointValues[17]*ny;
    Functionals[4] = PointValues[4]*nx + PointValues[18]*ny;
    Functionals[5] = PointValues[5]*nx + PointValues[19]*ny;
    
    // third edge:
    nx = y3 - y2;
    ny = x2 - x3;
    Functionals[6] = PointValues[6]*nx + PointValues[20]*ny;
    Functionals[7] = PointValues[7]*nx + PointValues[21]*ny;
    Functionals[8] = PointValues[8]*nx + PointValues[22]*ny;
    
    nx = y0 - y3;
    ny = x3 - x0;
    Functionals[9] = PointValues[9]*nx + PointValues[23]*ny;
    Functionals[10]= PointValues[10]*nx+ PointValues[24]*ny;
    Functionals[11]= PointValues[11]*nx+ PointValues[25]*ny;
    
    // the measure of the cell multiplied by the inverse measure of the 
    // refernce cell
    double measure = 0.25*Cell->GetMeasure();
    
    // we could use TQuadAffin if the Cell is a parallelogram, but this also
    // works for general quads
    TQuadBilinear referenceTransform;
    referenceTransform.SetCell(Cell);
    
    // transform the gradient of the (scalar) function phi(xi,eta) = xi
    // its gradient is (1,0) which is the vector with which we multiply to get
    // the correct dof
    
    double uref = 0., uxiref = 1., uetaref = 0., uorig, uxorig, uyorig;
    referenceTransform.GetOrigValues(NF_N_Q_BDM2_2D_Xi[12], 
                                     NF_N_Q_BDM2_2D_Eta[26], 1, &uref,  &uxiref,
                                     &uetaref, &uorig, &uxorig, &uyorig);
    Functionals[12] = (PointValues[12]*uxorig + PointValues[26]*uyorig) 
                         * measure;
    
    uxiref = 0.;
    uetaref = 1.;
    referenceTransform.GetOrigValues(NF_N_Q_BDM2_2D_Xi[13], 
                                     NF_N_Q_BDM2_2D_Eta[27], 1, &uref,  &uxiref,
                                     &uetaref, &uorig, &uxorig, &uyorig);
    Functionals[13] = (PointValues[13]*uxorig + PointValues[27]*uyorig) 
                         * measure;
  }
}

void NF_N_Q_BDM2_2D_EvalEdge(TCollection *Coll, TBaseCell *Cell, int Joint, double *PointValues,
                           double *Functionals)
{
// this is needed for setting boundary conditions.
  /* the functionals
   * int_Joint v.n q_1  and  int_Joint v.n q_2  and  int_Joint v.n q_3
   * (q_1, q2 and q_3 are two linearly independent polynomials of degree 2)
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
  Cell->GetVertex((Joint+1)%4)->GetCoords(x1,y1);// 4=number of edges
  #endif
  l = sqrt((x0-x1)*(x0-x1) + (y0-y1)*(y0-y1));
  Functionals[0] = PointValues[0]*l;
  Functionals[1] = PointValues[1]*l;
  Functionals[2] = PointValues[2]*l;
}

TNodalFunctional2D *NF_N_Q_BDM2_2D_Obj = new TNodalFunctional2D
        (NF_N_Q_BDM2_2D, 14, 3, 14, 3, NF_N_Q_BDM2_2D_Xi, NF_N_Q_BDM2_2D_Eta,
         NF_N_Q_BDM2_2D_T, NF_N_Q_BDM2_2D_EvalAll, NF_N_Q_BDM2_2D_EvalEdge);

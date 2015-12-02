/*
    TNodalFunctional2D(NodalFunctional2D id,
                       int n_allfunctionals, int n_edgefunctionals,
                       int n_pointsall, int n_pointsedge,
                       double *xi, double *eta, double *t,
                       DoubleFunctVect *evalall,
                       DoubleFunctVect *evaledge);
*/
// equidistant points on edges
//static double NF_N_Q_RT1_2D_a = 1./3.;
// Gauss points on edges
//static double NF_N_Q_RT1_2D_a = 1./sqrt(3.);
// Tschebyscheff points on edges
static double NF_N_Q_RT1_2D_a = 1./sqrt(2.);


static double NF_N_Q_RT1_2D_Xi[] 
 = {-NF_N_Q_RT1_2D_a, NF_N_Q_RT1_2D_a,
    1   ,1,
    NF_N_Q_RT1_2D_a,-NF_N_Q_RT1_2D_a,
    -1  ,-1,
    0, NF_N_Q_RT1_2D_a, 0, -NF_N_Q_RT1_2D_a };
static double NF_N_Q_RT1_2D_Eta[] 
 = {-1  ,-1,
   -NF_N_Q_RT1_2D_a,NF_N_Q_RT1_2D_a,
   1  , 1,
   NF_N_Q_RT1_2D_a,-NF_N_Q_RT1_2D_a,
   -NF_N_Q_RT1_2D_a, 0, NF_N_Q_RT1_2D_a, 0 };
// NOTE: If you want to use other evaluation points for degress of freedom on
// the edges of a cell, you also have to change basis functions in 
// BF_N_Q_RT1_2D.h
//static double NF_N_Q_RT1_2D_T[] = {-0.333333333333,0.3333333333333};// equidistant points
//static double NF_N_Q_RT1_2D_T[] = {-0.577350269189626,0.577350269189626};//Gauss-points
static double NF_N_Q_RT1_2D_T[] = {-0.707106781186547,0.707106781186547};//Tschebyscheff-points

void NF_N_Q_RT1_2D_EvalAll(TCollection *Coll, TBaseCell *Cell, double *PointValues,
                          double *Functionals)
{
  // on the reference cell [-1,1]^2
  if(Cell == nullptr)
  {
    Functionals[0] = -2*PointValues[12];
    Functionals[1] = -2*PointValues[13];
    Functionals[2] = 2*PointValues[2];
    Functionals[3] = 2*PointValues[3];
    Functionals[4] = 2*PointValues[16];
    Functionals[5] = 2*PointValues[17];
    Functionals[6] = -2*PointValues[6];
    Functionals[7] = -2*PointValues[7];
    
    Functionals[8] = PointValues[8];
    Functionals[9] = PointValues[21];
    Functionals[10]= PointValues[10];
    Functionals[11]= PointValues[23];
  }
  else // on a real cell
  {
    if(Cell->GetShapeDesc()->GetType() == Quadrangle) 
    {
      // not affine reference transform
      ErrThrow("NF_N_Q_RT1_2D_EvalAll not tested for non affine ",
               "reference transformations");
    }
    double x0, x1, x2, x3, y0, y1, y2, y3;
    #ifdef __2D__
    Cell->GetVertex(0)->GetCoords(x0, y0);
    Cell->GetVertex(1)->GetCoords(x1, y1);
    Cell->GetVertex(2)->GetCoords(x2, y2);
    Cell->GetVertex(3)->GetCoords(x3, y3);
    #else
    ErrThrow("NF_N_Q_RT1_2D_EvalAll not implemented in 3D");
    #endif
    // outer normal
    double nx, ny;
    
    // first edge:
    nx = y1 - y0;
    ny = x0 - x1;
    Functionals[0] = PointValues[0]*nx + PointValues[12]*ny;
    Functionals[1] = PointValues[1]*nx + PointValues[13]*ny;
    
    // second edge:
    nx = y2 - y1;
    ny = x1 - x2;
    Functionals[2] = PointValues[2]*nx + PointValues[14]*ny;
    Functionals[3] = PointValues[3]*nx + PointValues[15]*ny;
    
    // third edge:
    nx = y3 - y2;
    ny = x2 - x3;
    Functionals[4] = PointValues[4]*nx + PointValues[16]*ny;
    Functionals[5] = PointValues[5]*nx + PointValues[17]*ny;
    
    nx = y0 - y3;
    ny = x3 - x0;
    Functionals[6] = PointValues[6]*nx + PointValues[18]*ny;
    Functionals[7] = PointValues[7]*nx + PointValues[19]*ny;
    
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
    
    // first inner point
    double uref = 0., uxiref = 1., uetaref = 0., uorig, uxorig, uyorig;
    referenceTransform.GetOrigValues(NF_N_Q_RT1_2D_Xi[8], 
                                      NF_N_Q_RT1_2D_Eta[8], 1, &uref, &uxiref,
                                      &uetaref, &uorig, &uxorig, &uyorig);
    
    Functionals[8] = (PointValues[8]*uxorig + PointValues[20]*uyorig) * measure;
    
    // third inner point
    referenceTransform.GetOrigValues(NF_N_Q_RT1_2D_Xi[10], 
                                      NF_N_Q_RT1_2D_Eta[10], 1, &uref, &uxiref,
                                      &uetaref, &uorig, &uxorig, &uyorig);
    Functionals[10] = (PointValues[10]*uxorig + PointValues[22]*uyorig)*measure;
    
    // second inner point
    uxiref = 0.;
    uetaref = 1.;
    referenceTransform.GetOrigValues(NF_N_Q_RT1_2D_Xi[9], 
                                      NF_N_Q_RT1_2D_Eta[9], 1, &uref, &uxiref,
                                      &uetaref, &uorig, &uxorig, &uyorig);
    Functionals[9] = (PointValues[9]*uxorig + PointValues[21]*uyorig) * measure;
    
    // fourth inner point
    referenceTransform.GetOrigValues(NF_N_Q_RT1_2D_Xi[11], 
                                      NF_N_Q_RT1_2D_Eta[11], 1, &uref, &uxiref,
                                      &uetaref, &uorig, &uxorig, &uyorig);
    Functionals[11] = (PointValues[11]*uxorig + PointValues[23]*uyorig)*measure;
  }
}

void NF_N_Q_RT1_2D_EvalEdge(TCollection *Coll, TBaseCell *Cell, int Joint, double *PointValues,double *Functionals)
{
  // this is needed for setting boundary conditions.
  /* the functionals
   * int_Joint v.n q_1     and       int_Joint v.n q_2
   * (q_1 and q_2 are two linearly independent polynomials of degree 1)
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
}

TNodalFunctional2D *NF_N_Q_RT1_2D_Obj = new TNodalFunctional2D
        (NF_N_Q_RT1_2D, 12, 2, 12, 2, NF_N_Q_RT1_2D_Xi, NF_N_Q_RT1_2D_Eta,
         NF_N_Q_RT1_2D_T, NF_N_Q_RT1_2D_EvalAll, NF_N_Q_RT1_2D_EvalEdge);

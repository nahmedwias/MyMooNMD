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

static double NF_N_Q_RT3_2D_a = sqrt(2. + sqrt(2.)) / 2.;
static double NF_N_Q_RT3_2D_b = sqrt(2. - sqrt(2.)) / 2.;
static double NF_N_Q_RT3_2D_c = sqrt(3./4.);

static double NF_N_Q_RT3_2D_Xi[] = 
{ -NF_N_Q_RT3_2D_a,-NF_N_Q_RT3_2D_b, NF_N_Q_RT3_2D_b, NF_N_Q_RT3_2D_a,
  1, 1,1,1,
  NF_N_Q_RT3_2D_a,NF_N_Q_RT3_2D_b,-NF_N_Q_RT3_2D_b,-NF_N_Q_RT3_2D_a,
  -1,-1,-1,-1,
  -NF_N_Q_RT3_2D_a,-NF_N_Q_RT3_2D_b, NF_N_Q_RT3_2D_b, NF_N_Q_RT3_2D_a,
  -NF_N_Q_RT3_2D_a,-NF_N_Q_RT3_2D_b,NF_N_Q_RT3_2D_b,NF_N_Q_RT3_2D_a, 
  -NF_N_Q_RT3_2D_a,-NF_N_Q_RT3_2D_b,NF_N_Q_RT3_2D_b,NF_N_Q_RT3_2D_a,
  -NF_N_Q_RT3_2D_c,-NF_N_Q_RT3_2D_c,-NF_N_Q_RT3_2D_c,-NF_N_Q_RT3_2D_c,
  0, 0,0,0,
  NF_N_Q_RT3_2D_c, NF_N_Q_RT3_2D_c,NF_N_Q_RT3_2D_c,NF_N_Q_RT3_2D_c
};
static double NF_N_Q_RT3_2D_Eta[] = 
{ -1,-1,-1,-1,
  -NF_N_Q_RT3_2D_a,-NF_N_Q_RT3_2D_b,NF_N_Q_RT3_2D_b,NF_N_Q_RT3_2D_a,
  1,1, 1, 1,
  NF_N_Q_RT3_2D_a, NF_N_Q_RT3_2D_b,-NF_N_Q_RT3_2D_b,-NF_N_Q_RT3_2D_a,
  -NF_N_Q_RT3_2D_c,-NF_N_Q_RT3_2D_c,-NF_N_Q_RT3_2D_c,-NF_N_Q_RT3_2D_c,
  0, 0,0,0,
  NF_N_Q_RT3_2D_c, NF_N_Q_RT3_2D_c,NF_N_Q_RT3_2D_c,NF_N_Q_RT3_2D_c,
  -NF_N_Q_RT3_2D_a,-NF_N_Q_RT3_2D_b, NF_N_Q_RT3_2D_b, NF_N_Q_RT3_2D_a,
  -NF_N_Q_RT3_2D_a,-NF_N_Q_RT3_2D_b,NF_N_Q_RT3_2D_b,NF_N_Q_RT3_2D_a,
  -NF_N_Q_RT3_2D_a,-NF_N_Q_RT3_2D_b,NF_N_Q_RT3_2D_b,NF_N_Q_RT3_2D_a
};
// NOTE: If you want to use other evaluation points for degress of freedom on
// the edges of a cell, you also have to change basis functions in 
// BF_N_Q_RT3_2D.h
//static double NF_N_Q_RT3_2D_T[] = {-0.5,0,0.5};// equidistant points
//static double NF_N_Q_RT3_2D_T[] = {-0.774596669241483,0.774596669241483}//Gauss-points
static double NF_N_Q_RT3_2D_T[] = { -0.923879532511287,  -0.382683432365090,   0.382683432365090,   0.923879532511287};//Tschebyscheff-points

void NF_N_Q_RT3_2D_EvalAll(TCollection *Coll, TBaseCell *Cell, double *PointValues,
                          double *Functionals)
{
  // on the reference cell [-1,1]^2
  if(Cell == nullptr)
  {
    Functionals[0] = -2*PointValues[40];
    Functionals[1] = -2*PointValues[41];
    Functionals[2] = -2*PointValues[42];
    Functionals[3] = -2*PointValues[43];
    
    Functionals[4] = 2*PointValues[4];
    Functionals[5] = 2*PointValues[5];
    Functionals[6] = 2*PointValues[6];
    Functionals[7] = 2*PointValues[7];
    
    Functionals[8] = 2*PointValues[48];
    Functionals[9] = 2*PointValues[49];
    Functionals[10]= 2*PointValues[50];
    Functionals[11]= 2*PointValues[51];
    
    Functionals[12]= -2*PointValues[12];
    Functionals[13]= -2*PointValues[13];
    Functionals[14]= -2*PointValues[14];
    Functionals[15]= -2*PointValues[15];
    
    Functionals[16]= PointValues[56];
    Functionals[17]= PointValues[57];
    Functionals[18]= PointValues[58];
    Functionals[19]= PointValues[59];
    
    Functionals[20]= PointValues[60];
    Functionals[21]= PointValues[61];
    Functionals[22]= PointValues[62];
    Functionals[23]= PointValues[63];
    
    Functionals[24]= PointValues[64];
    Functionals[25]= PointValues[65];
    Functionals[26]= PointValues[66];
    Functionals[27]= PointValues[67];
    
    
    Functionals[28]= PointValues[28];
    Functionals[29]= PointValues[29];
    Functionals[30]= PointValues[30];
    Functionals[31]= PointValues[31];
    
    Functionals[32]= PointValues[32];
    Functionals[33]= PointValues[33];
    Functionals[34]= PointValues[34];
    Functionals[35]= PointValues[35];
    
    Functionals[36]= PointValues[36];
    Functionals[37]= PointValues[37];
    Functionals[38]= PointValues[38];
    Functionals[39]= PointValues[39];
  }
  else
  {
    if(Cell->GetShapeDesc()->GetType() == Quadrangle) 
    {
      // not affine reference transform
      ErrThrow("NF_N_Q_RT3_2D_EvalAll not tested for non affine ",
               "reference transformations");
    }
    double x0, x1, x2, x3, y0, y1, y2, y3;
    #ifdef __2D__
    Cell->GetVertex(0)->GetCoords(x0, y0);
    Cell->GetVertex(1)->GetCoords(x1, y1);
    Cell->GetVertex(2)->GetCoords(x2, y2);
    Cell->GetVertex(3)->GetCoords(x3, y3);
    #else
    ErrThrow("NF_N_Q_RT3_2D_EvalAll not implemented in 3D");
    #endif
    
    // outer normal
    double nx, ny;
    
    // first edge:
    nx = y1 - y0;
    ny = x0 - x1;
    Functionals[0] = PointValues[0]*nx + PointValues[40]*ny;
    Functionals[1] = PointValues[1]*nx + PointValues[41]*ny;
    Functionals[2] = PointValues[2]*nx + PointValues[42]*ny;
    Functionals[3] = PointValues[3]*nx + PointValues[43]*ny;
    
    // second edge:
    nx = y2 - y1;
    ny = x1 - x2;
    Functionals[4] = PointValues[4]*nx + PointValues[44]*ny;
    Functionals[5] = PointValues[5]*nx + PointValues[45]*ny;
    Functionals[6] = PointValues[6]*nx + PointValues[46]*ny;
    Functionals[7] = PointValues[7]*nx + PointValues[47]*ny;
    
    // third edge:
    nx = y3 - y2;
    ny = x2 - x3;
    Functionals[8] = PointValues[8]*nx + PointValues[48]*ny;
    Functionals[9] = PointValues[9]*nx + PointValues[49]*ny;
    Functionals[10]= PointValues[10]*nx+ PointValues[50]*ny;
    Functionals[11]= PointValues[11]*nx+ PointValues[51]*ny;
    
    nx = y0 - y3;
    ny = x3 - x0;
    Functionals[12]= PointValues[12]*nx+ PointValues[52]*ny;
    Functionals[13]= PointValues[13]*nx+ PointValues[53]*ny;
    Functionals[14]= PointValues[14]*nx+ PointValues[54]*ny;
    Functionals[15]= PointValues[15]*nx+ PointValues[55]*ny;
    
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
    
    // dofs in y-direction (y meaning in reference cell)
    for(unsigned int i = 0; i < 12; ++i)
    {
      double uref = 0., uxiref = 0., uetaref = 1., uorig, uxorig, uyorig;
      referenceTransform.GetOrigValues(NF_N_Q_RT3_2D_Xi[16+i], 
                                       NF_N_Q_RT3_2D_Eta[16+i], 1, &uref, 
                                       &uxiref, &uetaref, &uorig, &uxorig, 
                                       &uyorig);
      Functionals[16+i] = (PointValues[16+i]*uxorig + PointValues[56+i]*uyorig) 
                           * measure;
    }
    // dofs in x-direction (x meaning in reference cell)
    for(unsigned int i = 0; i < 12; ++i)
    {
      double uref = 0., uxiref = 1., uetaref = 0., uorig, uxorig, uyorig;
      referenceTransform.GetOrigValues(NF_N_Q_RT3_2D_Xi[28+i], 
                                       NF_N_Q_RT3_2D_Eta[28+i], 1, &uref, 
                                       &uxiref, &uetaref, &uorig, &uxorig, 
                                       &uyorig);
      Functionals[28+i] = (PointValues[28+i]*uxorig + PointValues[68+i]*uyorig) 
                           * measure;
    }
  }
}

void NF_N_Q_RT3_2D_EvalEdge(TCollection *Coll, TBaseCell *Cell, int Joint, double *PointValues,
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
  Functionals[3] = PointValues[3]*l;
}

TNodalFunctional2D *NF_N_Q_RT3_2D_Obj = new TNodalFunctional2D
        (NF_N_Q_RT3_2D, 40, 4, 40, 4, NF_N_Q_RT3_2D_Xi, NF_N_Q_RT3_2D_Eta,
         NF_N_Q_RT3_2D_T, NF_N_Q_RT3_2D_EvalAll, NF_N_Q_RT3_2D_EvalEdge);

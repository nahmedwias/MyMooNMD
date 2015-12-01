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

static double NF_N_T_RT3_2D_a = -sqrt(2. + sqrt(2.))/4. + 0.5;
static double NF_N_T_RT3_2D_b = -sqrt(2. - sqrt(2.))/4. + 0.5; 

static double NF_N_T_RT3_2D_Xi[] = 
{ NF_N_T_RT3_2D_a, NF_N_T_RT3_2D_b, 1-NF_N_T_RT3_2D_b, 1-NF_N_T_RT3_2D_a,
  1-NF_N_T_RT3_2D_a, 1-NF_N_T_RT3_2D_b, NF_N_T_RT3_2D_b, NF_N_T_RT3_2D_a, 
  0, 0, 0, 0,
  NF_N_T_RT3_2D_a, NF_N_T_RT3_2D_b,
  1-NF_N_T_RT3_2D_b, NF_N_T_RT3_2D_a,
  NF_N_T_RT3_2D_b, NF_N_T_RT3_2D_a };
static double NF_N_T_RT3_2D_Eta[] = 
{ 0,0,  0,  0,
  NF_N_T_RT3_2D_a,NF_N_T_RT3_2D_b,1-NF_N_T_RT3_2D_b,1-NF_N_T_RT3_2D_a,
  1-NF_N_T_RT3_2D_a,1-NF_N_T_RT3_2D_b,NF_N_T_RT3_2D_b,NF_N_T_RT3_2D_a,
  NF_N_T_RT3_2D_a, NF_N_T_RT3_2D_a,
  NF_N_T_RT3_2D_a, NF_N_T_RT3_2D_b,
  NF_N_T_RT3_2D_b, 1-NF_N_T_RT3_2D_b };
// NOTE: If you want to use other evaluation points for degress of freedom on
// the edges of a cell, you also have to change basis functions in 
// BF_N_T_RT3_2D.h
static double NF_N_T_RT3_2D_T[] = {-0.923879532511287,  -0.382683432365090,   0.382683432365090,   0.923879532511287};//Tschebyscheff-points

void NF_N_T_RT3_2D_EvalAll(TCollection *Coll, TBaseCell *Cell, double *PointValues,
                          double *Functionals)
{
  // on the reference triangle with points (0,0), (1,0), (0,1) 
  if(Cell == nullptr)
  {
    Functionals[0] = -PointValues[18];
    Functionals[1] = -PointValues[19];
    Functionals[2] = -PointValues[20];
    Functionals[3] = -PointValues[21];
    
    Functionals[4] = PointValues[4] + PointValues[22];
    Functionals[5] = PointValues[5] + PointValues[23];
    Functionals[6] = PointValues[6] + PointValues[24];
    Functionals[7] = PointValues[7] + PointValues[25];
    
    Functionals[8] = -PointValues[8];
    Functionals[9] = -PointValues[9];
    Functionals[10]= -PointValues[10];
    Functionals[11]= -PointValues[11];
    
    Functionals[12]= PointValues[12];
    Functionals[13]= PointValues[30];
    
    Functionals[14]= PointValues[13];
    Functionals[15]= PointValues[31];
    
    Functionals[16]= PointValues[14];
    Functionals[17]= PointValues[32];
    
    Functionals[18]= PointValues[15];
    Functionals[19]= PointValues[33];
    
    Functionals[20]= PointValues[16];
    Functionals[21]= PointValues[34];
    
    Functionals[22]= PointValues[17];
    Functionals[23]= PointValues[35];
  }
  else // on a real cell
  {
    double x0, x1, x2, y0, y1, y2;
    #ifdef __2D__
    Cell->GetVertex(0)->GetCoords(x0, y0);
    Cell->GetVertex(1)->GetCoords(x1, y1);
    Cell->GetVertex(2)->GetCoords(x2, y2);
    #else
    ErrThrow("NF_N_T_RT1_2D_EvalAll not implemented in 3D");
    #endif
    // length of edge, and outer normal
    double nx, ny;
    
    // first edge:
    nx = y1 - y0;
    ny = x0 - x1;
    Functionals[0] = PointValues[0]*nx + PointValues[18]*ny;
    Functionals[1] = PointValues[1]*nx + PointValues[19]*ny;
    Functionals[2] = PointValues[2]*nx + PointValues[20]*ny;
    Functionals[3] = PointValues[3]*nx + PointValues[21]*ny;
    
    // second edge:
    nx = y2 - y1;
    ny = x1 - x2;
    Functionals[4] = PointValues[4]*nx + PointValues[22]*ny;
    Functionals[5] = PointValues[5]*nx + PointValues[23]*ny;
    Functionals[6] = PointValues[6]*nx + PointValues[24]*ny;
    Functionals[7] = PointValues[7]*nx + PointValues[25]*ny;
    
    // third edge:
    nx = y0 - y2;
    ny = x2 - x0;
    Functionals[8] = PointValues[8]*nx + PointValues[26]*ny;
    Functionals[9] = PointValues[9]*nx + PointValues[27]*ny;
    Functionals[10]= PointValues[10]*nx+ PointValues[28]*ny;
    Functionals[11]= PointValues[11]*nx+ PointValues[29]*ny;
    
    // the measure of the cell multiplied by the inverse measure of the 
    // refernce cell
    double measure = 2*Cell->GetMeasure();
    
    TTriaAffin referenceTransform;
    referenceTransform.SetCell(Cell);
    // transform the gradient of the (scalar) function phi(xi,eta) = xi
    // its gradient is (1,0) which is the vector with which we multiply to get
    // the correct dof
    
    // dofs in x-direction (x meaning in reference cell)
    for(unsigned int i = 0; i < 6; ++i)
    {
      double uref = 0., uxiref = 1., uetaref = 0., uorig, uxorig, uyorig;
      referenceTransform.GetOrigValues(NF_N_Q_RT2_2D_Xi[12+i], 
                                       NF_N_Q_RT2_2D_Eta[12+i], 1, &uref, 
                                       &uxiref, &uetaref, &uorig, &uxorig, 
                                       &uyorig);
      Functionals[12+2*i] = ( PointValues[12+i]*uxorig 
                             +PointValues[30+i]*uyorig) * measure;
    }
    
    // dofs in y-direction (y meaning in reference cell)
    for(unsigned int i = 0; i < 6; ++i)
    {
      double uref = 0., uxiref = 0., uetaref = 1., uorig, uxorig, uyorig;
      referenceTransform.GetOrigValues(NF_N_Q_RT2_2D_Xi[12+i], 
                                       NF_N_Q_RT2_2D_Eta[12+i], 1, &uref, 
                                       &uxiref, &uetaref, &uorig, &uxorig, 
                                       &uyorig);
      Functionals[12+2*i+1] = ( PointValues[12+i]*uxorig 
                               +PointValues[30+i]*uyorig) * measure;
    }
  }
}

void NF_N_T_RT3_2D_EvalEdge(TCollection *Coll, TBaseCell *Cell, int Joint,
                            double *PointValues, double *Functionals)
{
// this is needed for setting boundary conditions
  /* the functionals
   * int_Joint v.n q_1  and  int_Joint v.n q_2  and  int_Joint v.n q_3
   * (q_1, q2 and q_3 are three linearly independent polynomials of degree 2)
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
  Functionals[2] = PointValues[2]*l;
  Functionals[3] = PointValues[3]*l;
}

TNodalFunctional2D *NF_N_T_RT3_2D_Obj = new TNodalFunctional2D
        (NF_N_T_RT3_2D, 24, 4, 18, 4, NF_N_T_RT3_2D_Xi, NF_N_T_RT3_2D_Eta,
         NF_N_T_RT3_2D_T, NF_N_T_RT3_2D_EvalAll, NF_N_T_RT3_2D_EvalEdge);

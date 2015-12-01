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

// equidistant points on edges
//static double NF_N_T_RT2_2D_a = 0.25;
// Gauss points on edges
//static double NF_N_T_RT2_2D_a = -sqrt(3./5.)/2+0.5;
// Tschebyscheff points on edges
static double NF_N_T_RT2_2D_a = -sqrt(3./4.)/2.+0.5;

static double NF_N_T_RT2_2D_Xi[] 
 = {NF_N_T_RT2_2D_a,0.5,1-NF_N_T_RT2_2D_a,
    1-NF_N_T_RT2_2D_a, 0.5, NF_N_T_RT2_2D_a,
    0,  0, 0,
    1./6., 2./3., 1./6. };
static double NF_N_T_RT2_2D_Eta[] 
 = {0,  0,  0,
    NF_N_T_RT2_2D_a, 0.5,1-NF_N_T_RT2_2D_a,
    1-NF_N_T_RT2_2D_a,0.5,NF_N_T_RT2_2D_a,
    1./6., 1./6., 2./3. };

// NOTE: If you want to use other evaluation points for degress of freedom on
// the edges of a cell, you also have to change basis functions in 
// BF_N_T_RT2_2D.h
//static double NF_N_T_RT2_2D_T[] = {-0.5,0,0.5};// equidistant points
//static double NF_N_T_RT2_2D_T[] = {-0.774596669241483,0.774596669241483};//Gauss-points
//Tschebyscheff-points
static double NF_N_T_RT2_2D_T[] = {-0.866025403784439,0,0.866025403784439};

void NF_N_T_RT2_2D_EvalAll(TCollection *Coll, TBaseCell *Cell, double *PointValues,
                          double *Functionals)
{
  // on the reference triangle with points (0,0), (1,0), (0,1) 
  if(Cell == nullptr)
  {
    Functionals[0] = -PointValues[12];
    Functionals[1] = -PointValues[13];
    Functionals[2] = -PointValues[14];
    
    Functionals[3] = PointValues[3] + PointValues[15];
    Functionals[4] = PointValues[4] + PointValues[16];
    Functionals[5] = PointValues[5] + PointValues[17];
    
    Functionals[6] = -PointValues[6];
    Functionals[7] = -PointValues[7];
    Functionals[8] = -PointValues[8];
    
    Functionals[9] = PointValues[9];
    Functionals[10]= PointValues[21];
    
    Functionals[11]= PointValues[10];
    Functionals[12]= PointValues[22];
    
    Functionals[13]= PointValues[11];
    Functionals[14]= PointValues[23];
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
    Functionals[0] = PointValues[0]*nx + PointValues[12]*ny;
    Functionals[1] = PointValues[1]*nx + PointValues[13]*ny;
    Functionals[2] = PointValues[2]*nx + PointValues[14]*ny;
    
    // second edge:
    nx = y2 - y1;
    ny = x1 - x2;
    Functionals[3] = PointValues[3]*nx + PointValues[15]*ny;
    Functionals[4] = PointValues[4]*nx + PointValues[16]*ny;
    Functionals[5] = PointValues[5]*nx + PointValues[17]*ny;
    
    // third edge:
    nx = y0 - y2;
    ny = x2 - x0;
    Functionals[6] = PointValues[6]*nx + PointValues[18]*ny;
    Functionals[7] = PointValues[7]*nx + PointValues[19]*ny;
    Functionals[8] = PointValues[8]*nx + PointValues[20]*ny;
    
    // the measure of the cell multiplied by the inverse measure of the 
    // refernce cell
    double measure = 2*Cell->GetMeasure();
    
    TTriaAffin referenceTransform;
    referenceTransform.SetCell(Cell);
    // transform the gradient of the (scalar) function phi(xi,eta) = xi
    // its gradient is (1,0) which is the vector with which we multiply to get
    // the correct dof
    
    // dofs in x-direction (x meaning in reference cell)
    for(unsigned int i = 0; i < 3; ++i)
    {
      double uref = 0., uxiref = 1., uetaref = 0., uorig, uxorig, uyorig;
      referenceTransform.GetOrigValues(NF_N_Q_RT2_2D_Xi[9+i], 
                                       NF_N_Q_RT2_2D_Eta[9+i], 1, &uref, 
                                       &uxiref, &uetaref, &uorig, &uxorig, 
                                       &uyorig);
      Functionals[9+2*i] = (PointValues[9+i]*uxorig + PointValues[21+i]*uyorig) 
                           * measure;
    }
    
    // dofs in y-direction (y meaning in reference cell)
    for(unsigned int i = 0; i < 3; ++i)
    {
      double uref = 0., uxiref = 0., uetaref = 1., uorig, uxorig, uyorig;
      referenceTransform.GetOrigValues(NF_N_Q_RT2_2D_Xi[9+i], 
                                       NF_N_Q_RT2_2D_Eta[9+i], 1, &uref, 
                                       &uxiref, &uetaref, &uorig, &uxorig, 
                                       &uyorig);
      Functionals[9+2*i+1] = ( PointValues[9+i]*uxorig 
                              +PointValues[21+i]*uyorig) * measure;
    }
  }
}

void NF_N_T_RT2_2D_EvalEdge(TCollection *Coll, TBaseCell *Cell, int Joint,
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
}

TNodalFunctional2D *NF_N_T_RT2_2D_Obj = new TNodalFunctional2D
        (NF_N_T_RT2_2D, 15, 3, 12, 3, NF_N_T_RT2_2D_Xi, NF_N_T_RT2_2D_Eta,
         NF_N_T_RT2_2D_T, NF_N_T_RT2_2D_EvalAll, NF_N_T_RT2_2D_EvalEdge);

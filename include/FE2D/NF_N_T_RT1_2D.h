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
//static double NF_N_T_RT1_2D_Xi[] =  { 1./3.,2./3.,2./3.,1./3.,0    , 0};
//static double NF_N_T_RT1_2D_Eta[] = { 0    ,0    ,1./3.,2./3.,2./3.,1./3.};
// Gauss points on edges
//static double NF_N_T_RT1_2D_Xi[] 
// =  { 0.211324865405187117745, 0.788675134594812882254,
//      0.788675134594812882254, 0.211324865405187117745,
//      0    , 0,
//      0.33333333333333333333};
//static double NF_N_T_RT1_2D_Eta[]
// = { 0    ,0,
//     0.211324865405187117745, 0.788675134594812882254,
//     0.788675134594812882254, 0.211324865405187117745,
//     0.33333333333333333333};
// Tschebyscheff points on edges
static double NF_N_T_RT1_2D_Xi[]
 =  { 0.14644660940672623779, 0.85355339059327376220,
      0.85355339059327376220, 0.14644660940672623779,
      0    , 0,
      0.33333333333333333333};
static double NF_N_T_RT1_2D_Eta[]
 = { 0    ,0,
     0.14644660940672623779, 0.85355339059327376220,
     0.85355339059327376220, 0.14644660940672623779,
     0.33333333333333333333};

// NOTE: If you want to use other evaluation points for degress of freedom on
// the edges of a cell, you also have to change basis functions in 
// BF_N_T_RT1_2D.h
//static double NF_N_T_RT1_2D_T[] = {-0.333333333333,0.3333333333333};// equidistant points
//static double NF_N_T_RT1_2D_T[] = {-0.577350269189626,0.577350269189626};//Gauss-points
static double NF_N_T_RT1_2D_T[] = {-0.707106781186547,0.707106781186547};//Tschebyscheff-points

void NF_N_T_RT1_2D_EvalAll(TCollection *Coll, TBaseCell *Cell, double *PointValues,
                          double *Functionals)
{
  // on the reference triangle with points (0,0), (1,0), (0,1) 
  if(Cell == nullptr)
  {
    Functionals[0] = -PointValues[7];
    Functionals[1] = -PointValues[8];
    Functionals[2] = PointValues[2] + PointValues[9];
    Functionals[3] = PointValues[3] + PointValues[10];
    Functionals[4] = -PointValues[4];
    Functionals[5] = -PointValues[5];
    Functionals[6] = PointValues[6];
    Functionals[7] = PointValues[13];
  }
  else // on a real cell
  {
    ErrThrow("NF_N_Q_RT0_2D_EvalAll not implemented on a real cell yet");
  }
}

void NF_N_T_RT1_2D_EvalEdge(TCollection *Coll, TBaseCell *Cell, int Joint, double *PointValues,
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

TNodalFunctional2D *NF_N_T_RT1_2D_Obj = new TNodalFunctional2D
        (NF_N_T_RT1_2D, 8, 2, 7, 2, NF_N_T_RT1_2D_Xi, NF_N_T_RT1_2D_Eta,
         NF_N_T_RT1_2D_T, NF_N_T_RT1_2D_EvalAll, NF_N_T_RT1_2D_EvalEdge);

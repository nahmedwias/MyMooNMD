/*
    TNodalFunctional2D(NodalFunctional2D id,
                       int n_allfunctionals, int n_edgefunctionals,
                       int n_pointsall, int n_pointsedge,
                       double *xi, double *eta, double *t,
                       DoubleFunctVect *evalall,
                       DoubleFunctVect *evaledge);
*/

static double NF_N_Q_RT0_2D_Xi[4] =  { 0, 1, 0, -1 };
static double NF_N_Q_RT0_2D_Eta[4] = {-1, 0, 1, 0  };
static double NF_N_Q_RT0_2D_T[2] = {0};

void NF_N_Q_RT0_2D_EvalAll(TCollection *Coll, TBaseCell *Cell, double *PointValues,
                          double *Functionals)
{
  // on the reference cell [-1,1]^2
  if(Cell == nullptr)
  {
    Functionals[0] = -2*PointValues[4];
    Functionals[1] = 2*PointValues[1];
    Functionals[2] = 2*PointValues[6];
    Functionals[3] = -2*PointValues[3];
  }
  else // on a real cell
  {
    ErrThrow("NF_N_Q_RT0_2D_EvalAll not implemented on a real cell yet");
  }
}

void NF_N_Q_RT0_2D_EvalEdge(TCollection *Coll, TBaseCell *Cell, int Joint, double *PointValues,
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
  Cell->GetVertex((Joint+1)%4)->GetCoords(x1,y1);// 4=number of edges
  #endif
  l = sqrt((x0-x1)*(x0-x1) + (y0-y1)*(y0-y1));
  Functionals[0] = PointValues[0]*l;
}

TNodalFunctional2D *NF_N_Q_RT0_2D_Obj = new TNodalFunctional2D
        (NF_N_Q_RT0_2D, 4, 1, 4, 1, NF_N_Q_RT0_2D_Xi, NF_N_Q_RT0_2D_Eta,
         NF_N_Q_RT0_2D_T, NF_N_Q_RT0_2D_EvalAll, NF_N_Q_RT0_2D_EvalEdge);

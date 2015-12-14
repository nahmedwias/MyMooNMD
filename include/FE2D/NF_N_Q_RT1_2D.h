/*
    TNodalFunctional2D(NodalFunctional2D id,
                       int n_allfunctionals, int n_edgefunctionals,
                       int n_pointsall, int n_pointsedge,
                       double *xi, double *eta, double *t,
                       DoubleFunctVect *evalall,
                       DoubleFunctVect *evaledge);
*/

static double NF_N_Q_RT1_2D_a = sqrt(3./5.);
static double NF_N_Q_RT1_2D_Xi[] =
{-NF_N_Q_RT1_2D_a, 0, NF_N_Q_RT1_2D_a,
 1, 1, 1,
 NF_N_Q_RT1_2D_a, 0, -NF_N_Q_RT1_2D_a,
 -1, -1, -1,
 -NF_N_Q_RT1_2D_a, 0, NF_N_Q_RT1_2D_a,
 -NF_N_Q_RT1_2D_a, 0, NF_N_Q_RT1_2D_a,
 -NF_N_Q_RT1_2D_a, 0, NF_N_Q_RT1_2D_a
};
static double NF_N_Q_RT1_2D_Eta[]  =
{-1, -1, -1,
 -NF_N_Q_RT1_2D_a, 0, NF_N_Q_RT1_2D_a,
  1, 1, 1,
  NF_N_Q_RT1_2D_a, 0, -NF_N_Q_RT1_2D_a,
  -NF_N_Q_RT1_2D_a, -NF_N_Q_RT1_2D_a, -NF_N_Q_RT1_2D_a,
  0, 0, 0, 
  NF_N_Q_RT1_2D_a, NF_N_Q_RT1_2D_a, NF_N_Q_RT1_2D_a
};

static double NF_N_Q_RT1_2D_T[] = { -NF_N_Q_RT1_2D_a, 0, NF_N_Q_RT1_2D_a };

void NF_N_Q_RT1_2D_EvalAll(TCollection *Coll, TBaseCell *Cell,
                           double *PointValues, double *Functionals)
{
  // on the reference cell [-1,1]^2
  if(Cell == nullptr)
  {
    Functionals[0] = -( 5*PointValues[21] + 8*PointValues[22]
                       +5*PointValues[23] )/9.;
    Functionals[1] = -NF_N_Q_RT1_2D_a*5*(-PointValues[21] + PointValues[23])/9.;
    Functionals[2] = ( 5*PointValues[3] + 8*PointValues[4] 
                     + 5*PointValues[5] )/9.;
    Functionals[3] = NF_N_Q_RT1_2D_a*5*(-PointValues[3] + PointValues[5] )/9.;
    Functionals[4] = ( 5*PointValues[27] + 8*PointValues[28]
                      +5*PointValues[29] )/9.;
    Functionals[5] = NF_N_Q_RT1_2D_a*5*(-PointValues[27] + PointValues[29] )/9.;
    Functionals[6] = -( 5*PointValues[9] + 8*PointValues[10]
                       +5*PointValues[11] )/9.;
    Functionals[7] = -NF_N_Q_RT1_2D_a*5*(-PointValues[9] + PointValues[11] )/9.;
    
    Functionals[8] = ( 25*PointValues[12]+40*PointValues[13]+25*PointValues[14]
                      +40*PointValues[15]+64*PointValues[16]+40*PointValues[17]
                      +25*PointValues[18]+40*PointValues[19]+25*PointValues[20]
                      )/81.;
    Functionals[9] = ( 25*PointValues[33]+40*PointValues[34]+25*PointValues[35]
                      +40*PointValues[36]+64*PointValues[37]+40*PointValues[38]
                      +25*PointValues[39]+40*PointValues[40]+25*PointValues[41]
                      )/81.;
    Functionals[10]= (-25*PointValues[12]-40*PointValues[13]-25*PointValues[14]
                      +25*PointValues[18]+40*PointValues[19]+25*PointValues[20]
                      )*NF_N_Q_RT1_2D_a/81.;
    Functionals[11]= (-25*PointValues[33]+25*PointValues[35]
                      -40*PointValues[36]+40*PointValues[38]
                      -25*PointValues[39]+25*PointValues[41]
                      )*NF_N_Q_RT1_2D_a/81.;
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
    Functionals[0] = ( ( 5*PointValues[0] + 8*PointValues[1] 
                        +5*PointValues[2] )*nx 
                      +( 5*PointValues[21] + 8*PointValues[22]
                        +5*PointValues[23])*ny )/18.;
    Functionals[1] = ( ( -5*PointValues[0] + 5*PointValues[2] )*nx 
                      +( -5*PointValues[21] + 5*PointValues[23] )*ny 
                      ) * NF_N_Q_RT1_2D_a/18.;
    Functionals[0] *= Cell->GetNormalOrientation(0);
    //Functionals[1] *= Cell->GetNormalOrientation(0);
    
    // second edge:
    nx = y2 - y1;
    ny = x1 - x2;
    Functionals[2] = ( ( 5*PointValues[3] + 8*PointValues[4] 
                        +5*PointValues[5] )*nx 
                      +( 5*PointValues[24] + 8*PointValues[25]
                        +5*PointValues[26])*ny )/18.;
    Functionals[3] = ( ( -5*PointValues[3] + 5*PointValues[5] )*nx 
                      +( -5*PointValues[24] + 5*PointValues[26] )*ny 
                      )*NF_N_Q_RT1_2D_a/18.;
    Functionals[2] *= Cell->GetNormalOrientation(1);
    //Functionals[3] *= Cell->GetNormalOrientation(1);
    
    // third edge:
    nx = y3 - y2;
    ny = x2 - x3;
    Functionals[4] = ( ( 5*PointValues[6] + 8*PointValues[7] 
                        +5*PointValues[8] )*nx 
                      +( 5*PointValues[27] + 8*PointValues[28]
                        +5*PointValues[29])*ny )/18.;
    Functionals[5] = ( ( -5*PointValues[6] + 5*PointValues[8] )*nx 
                      +( -5*PointValues[27] + 5*PointValues[29] )*ny 
                      )*NF_N_Q_RT1_2D_a/18.;
    Functionals[4] *= Cell->GetNormalOrientation(2);
    //Functionals[5] *= Cell->GetNormalOrientation(2);
    
    nx = y0 - y3;
    ny = x3 - x0;
    Functionals[6] = ( ( 5*PointValues[9] + 8*PointValues[10] 
                        +5*PointValues[11] )*nx 
                      +( 5*PointValues[30] + 8*PointValues[31]
                        +5*PointValues[32])*ny )/18.;
    Functionals[7] = ( ( -5*PointValues[9] + 5*PointValues[11] )*nx 
                      +( -5*PointValues[30] + 5*PointValues[32] )*ny 
                      )*NF_N_Q_RT1_2D_a/18.;
    Functionals[6] *= Cell->GetNormalOrientation(3);
    //Functionals[7] *= Cell->GetNormalOrientation(3);
    
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
    referenceTransform.GetOrigValues(NF_N_Q_RT1_2D_Xi[16], 
                                     NF_N_Q_RT1_2D_Eta[16], 1, &uref, &uxiref,
                                     &uetaref, &uorig, &uxorig, &uyorig);
    
    Functionals[8] = (
         uxorig * ( 25*PointValues[12]+40*PointValues[13]+25*PointValues[14]
                   +40*PointValues[15]+64*PointValues[16]+40*PointValues[17]
                   +25*PointValues[18]+40*PointValues[19]+25*PointValues[20]
                  )/81.
       + uyorig * ( 25*PointValues[33]+40*PointValues[34]+25*PointValues[35]
                   +40*PointValues[36]+64*PointValues[37]+40*PointValues[38]
                   +25*PointValues[39]+40*PointValues[40]+25*PointValues[41]
                  )/81.) * measure;
    
    // third inner point
    referenceTransform.GetOrigValues(NF_N_Q_RT1_2D_Xi[16], 
                                     NF_N_Q_RT1_2D_Eta[16], 1, &uref, &uxiref,
                                     &uetaref, &uorig, &uxorig, &uyorig);
    Functionals[10] = (
         uxorig * (-25*PointValues[12]-40*PointValues[13]-25*PointValues[14]
                   +25*PointValues[18]+40*PointValues[19]+25*PointValues[20]
                  )/81.
       + uyorig * (-25*PointValues[33]-40*PointValues[34]-25*PointValues[35]
                   +25*PointValues[39]+40*PointValues[40]+25*PointValues[41]
                  )/81.) * NF_N_Q_RT1_2D_a * measure;
    
    // second inner point
    uxiref = 0.;
    uetaref = 1.;
    referenceTransform.GetOrigValues(NF_N_Q_RT1_2D_Xi[16], 
                                     NF_N_Q_RT1_2D_Eta[16], 1, &uref, &uxiref,
                                     &uetaref, &uorig, &uxorig, &uyorig);
    Functionals[9] = (
         uxorig * ( 25*PointValues[12]+40*PointValues[13]+25*PointValues[14]
                   +40*PointValues[15]+64*PointValues[16]+40*PointValues[17]
                   +25*PointValues[18]+40*PointValues[19]+25*PointValues[20]
                  )/81.
       + uyorig * ( 25*PointValues[33]+40*PointValues[34]+25*PointValues[35]
                   +40*PointValues[36]+64*PointValues[37]+40*PointValues[38]
                   +25*PointValues[39]+40*PointValues[40]+25*PointValues[41]
                  )/81.) * measure;
    
    // fourth inner point
    referenceTransform.GetOrigValues(NF_N_Q_RT1_2D_Xi[16], 
                                     NF_N_Q_RT1_2D_Eta[16], 1, &uref, &uxiref,
                                     &uetaref, &uorig, &uxorig, &uyorig);
    Functionals[11] = (
         uxorig * (-25*PointValues[12]+25*PointValues[14]
                   -40*PointValues[15]+40*PointValues[17]
                   -25*PointValues[18]+25*PointValues[20]
                  )/81.
       + uyorig * (-25*PointValues[33]+25*PointValues[35]
                   -40*PointValues[36]+40*PointValues[38]
                   -25*PointValues[39]+25*PointValues[41]
                  )/81.) * NF_N_Q_RT1_2D_a * measure;
  }
}

void NF_N_Q_RT1_2D_EvalEdge(TCollection *Coll, TBaseCell *Cell, int Joint, double *PointValues,double *Functionals)
{
  #ifdef __2D__
  double x0, x1, y0, y1;
  Cell->GetVertex(Joint)->GetCoords(x0, y0);
  Cell->GetVertex((Joint+1)%4)->GetCoords(x1, y1); // 4=number of edges
  double l = sqrt((x0-x1)*(x0-x1) + (y0-y1)*(y0-y1)); // length of joint
  Functionals[0] = (5*PointValues[0]+8*PointValues[1]+5*PointValues[2])*l/18.;
  Functionals[1] = (-PointValues[0] + PointValues[2])*NF_N_Q_RT1_2D_a*l*5/18.;
  #endif
}

TNodalFunctional2D *NF_N_Q_RT1_2D_Obj = new TNodalFunctional2D
        (NF_N_Q_RT1_2D, 12, 2, 21, 3, NF_N_Q_RT1_2D_Xi, NF_N_Q_RT1_2D_Eta,
         NF_N_Q_RT1_2D_T, NF_N_Q_RT1_2D_EvalAll, NF_N_Q_RT1_2D_EvalEdge);

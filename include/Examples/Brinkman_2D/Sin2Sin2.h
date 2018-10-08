// Brinkman problem with sine2-sine2-velocity from Jarle Sogn resp. Mardal, Thai, Winther

double viscosity = -1;
double effective_viscosity = -1;
double permeability = -1;

void ExampleFile()
{
  Output::print<1>("Example: Sin2Sin2.h");
}

// ========================================================================
// exact solution
// ========================================================================
/* exact solution:
   u(x,y) = grad x (sin^2(Pi*x)*sin^2(Pi*y)) = (u1,u2)
          = (-d/dy (sin^2(Pi*x)*sin^2(Pi*y)) ,  d/dx (sin^2(Pi*x)*sin^2(Pi*y)))
   p(x,y) = -sin(2*Pi*x)
*/

void ExactU1(double x, double y, double *values)
{
    values[0] = - 2*Pi*cos(Pi*y)*sin(Pi*y)*sin(Pi*x)*sin(Pi*x);                   //u1
    values[1] = - 4*Pi*Pi*cos(Pi*y)*sin(Pi*y)*cos(Pi*x)*sin(Pi*x);                //u1_x
    values[2] = - 2*Pi*Pi*cos(Pi*y)*cos(Pi*y)*sin(Pi*x)*sin(Pi*x)
                + 2*Pi*Pi*sin(Pi*y)*sin(Pi*y)*sin(Pi*x)*sin(Pi*x);                //u1_y
    values[3] = 8*Pi*Pi*Pi*cos(Pi*y)*sin(Pi*y)*sin(Pi*x)*sin(Pi*x)
                + ( 4*Pi*Pi*Pi*sin(Pi*x)*sin(Pi*x)*cos(Pi*y)*sin(Pi*y)
                - 4*Pi*Pi*Pi*cos(Pi*x)*cos(Pi*x)*cos(Pi*y)*sin(Pi*y) ) ;          //Delta u1=u1_xx + u1_yy
}

void ExactU2(double x, double y, double *values)
{
    values[0] = 2*Pi*cos(Pi*x)*sin(Pi*x)*sin(Pi*y)*sin(Pi*y);                //u2
    values[1] = 2*Pi*Pi*cos(Pi*x)*cos(Pi*x)*sin(Pi*y)*sin(Pi*y)
                - 2*Pi*Pi*sin(Pi*x)*sin(Pi*x)*sin(Pi*y)*sin(Pi*y);           //u2_x
    values[2] = 4*Pi*Pi*cos(Pi*y)*sin(Pi*y)*cos(Pi*x)*sin(Pi*x);             //u2_y
    values[3] = - 8*Pi*Pi*Pi*cos(Pi*x)*sin(Pi*x)*sin(Pi*y)*sin(Pi*y)
                + 4*Pi*Pi*Pi*cos(Pi*y)*cos(Pi*y)*cos(Pi*x)*sin(Pi*x)
                - 4*Pi*Pi*Pi*sin(Pi*y)*sin(Pi*y)*cos(Pi*x)*sin(Pi*x);         //Delta u2=u2_xx + u2_yy
}

void ExactP(double x, double y, double *values)
{
    values[0] = -sin(2*Pi*x);                     //p
    values[1] = -2*Pi*cos(2*Pi*x);                //p_x
    values[2] = 0;                                //p_y
    values[3] = 4*Pi*Pi*sin(2*Pi*x)   +    0;     //Delta p=p_xx+p_yy
}


// ========================================================================
// boundary conditions (Parametrisierung des Randes); Param \in [0,1]
// ========================================================================

void BoundCondition(int i, double Param, BoundCond &cond)
{
    cond = DIRICHLET; // default

    for (int j = 0; j < TDatabase::ParamDB->n_neumann_boundary; j++)
    {
        if (i == TDatabase::ParamDB->neumann_boundary_id[j])
        {
            cond = NEUMANN;
            return;
        }
    }
    
    for (int j = 0; j < TDatabase::ParamDB->n_nitsche_boundary; j++)
    {
        if (i == TDatabase::ParamDB->nitsche_boundary_id[j])
        {
            cond = DIRICHLET_WEAK;
            return;
        }
    }
}

void U1BoundValue(int BdComp, double Param, double &value)
{
    
    // loop to impose Neumann boundary conditions
    for (int j = 0; j < TDatabase::ParamDB->n_neumann_boundary; j++)
    {
        if ( BdComp == TDatabase::ParamDB->neumann_boundary_id[j])
        {
            switch(BdComp)
            {
                    // case 0: value=0.;//TDatabase::ParamDB->neumann_boundary_value[j];
                    // break;
                case 1: value=0.;//TDatabase::ParamDB->neumann_boundary_value[j];
                    break;
                    // case 2: value=0.;//-TDatabase::ParamDB->neumann_boundary_value[j];
                    // break;
                case 3: value=0.;//-TDatabase::ParamDB->neumann_boundary_value[j];
                    break;
            }
            return;
        }
    }
    
    // loop to impose (strong or weak) Dirichlet
    switch(BdComp)
    {
        case 0: value= 0;
            break;
        case 1: value= 0;
            break;
        case 2: value= 0;
            break;
        case 3: value= 0;
            break;
        default: cout << "No boundary component with this number." << endl;
            break;
    }
}

void U2BoundValue(int BdComp, double Param, double &value)
{
    
    // loop to impose Neumann boundary conditions
    for (int j=0; j<TDatabase::ParamDB->n_neumann_boundary; j++)
    {
        if ( BdComp==TDatabase::ParamDB->neumann_boundary_id[j])
        {
            switch(BdComp)
            {
//                case 0:
//                    value=0.;//TDatabase::ParamDB->neumann_boundary_value[j];
//                    break;
                case 1:
                    value=0.;//TDatabase::ParamDB->neumann_boundary_value[j];
                    break;
//                case 2:
//                    value=0.;//-TDatabase::ParamDB->neumann_boundary_value[j];
//                    break;
                case 3:
                    value=0.;//-TDatabase::ParamDB->neumann_boundary_value[j];
                    break;
                default:
                    Output::print("I cannot impose Neumann boundary condition on component ", BdComp);
                    exit(1);
                    break;
            }
            return;
        }
    }
    
    // loop to impose (strong or weak) Dirichlet
    switch(BdComp)
    {
        case 0: value= 0;
            break;
        case 1: value= 0;
            break;
        case 2: value= 0;
            break;
        case 3: value= 0;
            break;
        default: cout << "No boundary component with this number." << endl;
            break;
    }

}



// FORMER VERSION OF ULRICH
//void BoundCondition(int i, double t, BoundCond &cond)
//{
//  cond = DIRICHLET;
//}
//
//void U1BoundValue(int BdComp, double Param, double &value)
//{
//  // Neumann boundary means setting T.n where n is outer normal and T is either
//  // 2nu D(u) - p  (for LAPLACETYPE==1) or nu grad(u) - p   (for LAPLACETYPE==0)
//  const int lt = TDatabase::ParamDB->LAPLACETYPE;
//  const double nu = 1./TDatabase::ParamDB->RE_NR;
//  BoundCond cond;
//  BoundCondition(BdComp, Param, cond);
//  switch(BdComp)
//  {
//    case 0:
//      if(cond == DIRICHLET)
//        value = sin(Pi*Param); // Dirichlet
//      else
//        value = lt == 0 ? 0. : -nu*Pi*Pi*sin(Pi*Param); // Neumann
//      break;
//    case 1:
//      if(cond == DIRICHLET)
//        value = 0.; // Dirichlet
//      else
//        value = lt == 0 ? -nu*Pi : -2*nu*Pi; // Neumann
//      break;
//    case 2:
//      if(cond == DIRICHLET)
//        value = sin(Pi*(1-Param)); // Dirichlet
//      else
//        value = lt == 0 ? 0. : nu*Pi*Pi*sin(Pi*(1-Param)); // Neumann
//      break;
//    case 3:
//      if(cond == DIRICHLET)
//        value = 0.; // Dirichlet
//      else
//        value = lt == 0 ? -nu*Pi : -2*nu*Pi; // Neumann
//      break;
//    default:
//      ErrThrow("wrong boundary part number", BdComp);
//      break;
//  }
//  return;
//}
//
//void U2BoundValue(int BdComp, double Param, double &value)
//{
//  // Neumann boundary means setting T.n where n is outer normal and T is either
//  // 2nu D(u) - p  (for LAPLACETYPE==1) or nu grad(u) - p   (for LAPLACETYPE==0)
//  const int lt = TDatabase::ParamDB->LAPLACETYPE;
//  const double nu = 1./TDatabase::ParamDB->RE_NR;
//  BoundCond cond;
//  BoundCondition(BdComp, Param, cond);
//  switch(BdComp)
//  {
//    case 0:
//      if(cond == DIRICHLET)
//        value = 0; // Dirichlet
//      else
//      {
//        value = nu*Pi*cos(Pi*Param) * (lt == 0 ? 1. : 2.); // Neumann
//        value += sin(Pi*Param);// Neumann (from pressure)
//      }
//      break;
//    case 1:
//      if(cond == DIRICHLET)
//        value = Pi*Param; // Dirichlet
//      else
//        value = 0.; // Neumann
//      break;
//    case 2:
//      if(cond == DIRICHLET)
//        value = -Pi*cos(Pi*(1-Param)); // Dirichlet
//      else
//      {
//        value = -nu*Pi*cos(Pi*(1-Param)) * (lt == 0 ? 1. : 2.); // Neumann
//        value += sin(Pi*(1-Param));// Neumann (from pressure)
//      }
//      break;
//    case 3:
//      if(cond == DIRICHLET)
//        value = -Pi*(1-Param); // Dirichlet
//      else
//        value = 0.; // Neumann
//      break;
//    default:
//      ErrThrow("wrong boundary part number", BdComp);
//      break;
//  }
//  return;
//}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y,
               double **parameters, double **coeffs)
{
    double val_u1[4];
    double val_u2[4];
    double val_p[4];
    
    for(int i = 0; i < n_points; i++)
    {
        ExactU1(X[i], Y[i], val_u1);
        ExactU2(X[i], Y[i], val_u2);
        ExactP(X[i], Y[i], val_p);
        
        coeffs[i][4] = viscosity;
        coeffs[i][5] = effective_viscosity;
        coeffs[i][6] = permeability;
         
        coeffs[i][0] = (coeffs[i][5] / coeffs[i][4]) * coeffs[i][6];
        coeffs[i][1] = -coeffs[i][5]*val_u1[3] + val_p[1] + (coeffs[i][4]/coeffs[i][6])*val_u1[0];//-coeffs[i][5]*val1[3] + val3[1] + (coeffs[i][4]/coeffs[i][6])*val1[0];  //0;// f1
        coeffs[i][2] = -coeffs[i][5]*val_u2[3] + val_p[2] + (coeffs[i][4]/coeffs[i][6])*val_u2[0];//-coeffs[i][5]*val2[3] + val3[2] + (coeffs[i][4]/coeffs[i][6])*val2[0];//0;//  // f2
        coeffs[i][3] = 0;// val_u1[1] + val_u2[2];  //                        // g (divergence)
        coeffs[i][7] = TDatabase::ParamDB->equal_order_stab_weight_PkPk;
        coeffs[i][8] = TDatabase::ParamDB->grad_div_stab_weight;
       }
}

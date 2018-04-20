// Brinkman problem, solution in ansatz space
// Poiseuille (exact solution in P2/P1)
// 
// u(x,y,z) = ( 0, 0 , 1-(x^2+y^2))^T=(u1,u2,u3)^T
// p(x,y,z) = -DP/Length*(z-zInlet)+DP/2;

// This example is adapted to the geometry cylinder_short_18Ktetra.mesh
// The pressure solution is set according to the Hagen-Poiseuille law (https://en.wikipedia.org/wiki/Hagen–Poiseuille_equation): DP= (8*mu*L)/(r^2) * v

double viscosity = -1;
double effective_viscosity = -1;
double permeability = -1;

void ExampleFile()
{
  Output::info<1>("EXAMPLE","Poiseuille.h");
}

// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y,  double z, double *values)
{
    values[0] = 0;           // u1
    values[1] = 0;           // u1_x
    values[2] = 0;           // u1_y
    values[3] = 0;           // u1_z
    values[4] = 0;           // Delta u1=u1_xx+u1_yy+u1_zz
}

void ExactU2(double x, double y,  double z, double *values)
{
    values[0] = 0;            // u2
    values[1] = 0;            // u2_x
    values[2] = 0;            // u2_y
    values[3] = 0;            // u2_z
    values[4] = 0;            // Delta u2=u2_xx+u2_yy+u2_zz
}

void ExactU3(double x, double y,  double z, double *values)
{
    values[0] = 1-(x*x+y*y);        // u3
    values[1] = -2*x;               // u3_x
    values[2] = -2*y;               // u3_y
    values[3] = 0;                  // u3_z
    values[4] = -4;                 // Delta u3=u3_xx+u3_yy+u3_zz

}

void ExactP(double x, double y,  double z, double *values)
{
    double Length = 1.;
    double Radius = 1.;
    double Umax = 1.;
    double DP = 4 * effective_viscosity * Length/(Radius * Radius) * Umax;
    //  double DP = 8*effective_viscosity*Length/(Radius*Radius)*Umax;
    //double Pinlet = DP/2;
    //double Poutlet = -DP/2;
    double zInlet = 1.;
    values[0] = -DP/Length*(z-zInlet)+DP/2;
    //  values[0] = -DP/3*z+5*DP/6;
    values[1] = 0;
    values[2] = 0;
    values[3] = -DP/Length;
    values[4] = 0;
}

// kind of boundary condition (for FE space needed);
// To all nodes (x,y,z) a bound condition has to be assigned to (Usually DIRICHLET or NEUMANN),
// The condition "DIRICHLET" assures the uniqueness of the pressure in case of Dirichlet b.c.
void BoundCondition(double x, double y, double z, BoundCond &cond)
{
    double zBottom = 0.;
    double zTop = 1.;
  
    cond = DIRICHLET;
    if (TDatabase::ParamDB->n_neumann_boundary==0)
    {
        return;
        // If n_neumann_boundary==1 this means in this example that we set Neumann bc. on top of the cylinder
    } else if (TDatabase::ParamDB->n_neumann_boundary==1){
        if (z==zTop) {
            cond = NEUMANN;
            return;
        }
        // If n_neumann_boundary==2 this means in this example that we set Neumann bc. on top and at the bottom of the cylinder
    } else if (TDatabase::ParamDB->n_neumann_boundary==2){
        if (z==zBottom || z==zTop) {
            cond = NEUMANN;
            return;
        }
    }
    
    
    
//  // first case: 2 Neumann boundaries  (label 1 and 2), one Nitsche (label 3)
//  if ((TDatabase::ParamDB->n_neumann_boundary + TDatabase::ParamDB->n_nitsche_boundary)==3 ) {
//    cond = NEUMANN;
//  } else {
//    // second case: all Dirichlet
//    cond = DIRICHLET; 
//    if (TDatabase::ParamDB->n_neumann_boundary==0) {
//      return;
//    } else {
//
//      // third case:  2 Neumann boundaries (label 1 and 2), Dirichlet on label 3
//      
//        
//        if (z==zTop) {
//            cond = NEUMANN;
//            return;
//        }
//  //    if (z==zBottom || z==zTop) {
//  //        cond = NEUMANN;
//   //       return;
//    //  }
//    }
//  }
}


// value of boundary condition
void U1BoundValue(double x, double y, double z, double &value)
{
    value = 0;
}

// value of boundary condition
void U2BoundValue(double x, double y, double z, double &value)
{
    value = 0;
}

// value of boundary condition
void U3BoundValue(double x, double y, double z, double &value)
{
    double zBottom = 0.;
    
    if (TDatabase::ParamDB->n_neumann_boundary==0 && TDatabase::ParamDB->n_nitsche_boundary==0)
        value = 1-(x*x+y*y);
    
    else if (TDatabase::ParamDB->n_neumann_boundary==1)
    {
        if (z==zBottom) {
            value = 1-(x*x+y*y);
        }
        else
            value = 0.;
    }
    else
        value = 0.;
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y, double *Z,
               double **parameters, double **coeffs)
{
    double *coeff;
    //double Length = 1;
    //double Radius = 1.;
    //double Umax = 1.;
    //double DP = 4* effective_viscosity * Length/(Radius*Radius)*Umax;
    //  double DP = 8 * effective_viscosity * Length/(Radius*Radius)*Umax;

    
    for(int i=0;i<n_points;i++)
    {
        coeff = coeffs[i];
        // Note: Setting f3=0, and PERMEABILITY=100000 gives the Stokes case. Here Delta u and nabla p balance each other for the functions chosen in the example;
        //coeff[0] = eps;
        coeff[5] = viscosity;//0.;
        coeff[6] = effective_viscosity;
        coeff[7] = permeability;
        coeff[1] = 0;//coeff[6]*(-12)+(-1)+(coeff[5]/coeff[7])*(3*(1-(Y[i]*Y[i]+Z[i]*Z[i]))); // f1
        coeff[2] = 0; // f2
        //coeff[3] = -coeff[6]*(-4)- DP/Length +(coeff[5]/coeff[7])*(1-(Y[i]*Y[i]+X[i]*X[i])); // 0; // f3
        coeff[3] = 0;//-coeff[6]*(-4)- DP/3 +(coeff[5]/coeff[7])*(1-(Y[i]*Y[i]+X[i]*X[i])); // 0; // f3
        coeff[4] = 0; // g
        coeff[8]=TDatabase::ParamDB->equal_order_stab_weight_PkPk;
        coeff[9]=TDatabase::ParamDB->equal_order_stab_weight_PkPk;
    }
 
}

// Brinkman problem, Poiseuille-Problem from Hannukainen


void ExampleFile()
{
    Output::print<1>("Example: Poiseuille_Hannukainen.h");
}
/*
 Square [0,1]x[0,1]
 t^2=K*(nu_eff/nu)
 u(x,y) = [(1+exp(1/t)-exp((1-y)/t) - exp(y/t)) / (1+exp(1/t)  ,  0], for t>0
 u(x,y) = [1  ,  0], for t=1
 p(x,y) = (nu/K)*(-x+1/2)
 (Note that f and g vanish)
 */

// ========================================================================
// exact solution
// ========================================================================


void ExactU1(double x, double y, double *values)
{
    double K = TDatabase::ParamDB->PERMEABILITY;
    double nu = TDatabase::ParamDB->VISCOSITY;
    double nu_eff = TDatabase::ParamDB->EFFECTIVE_VISCOSITY;
    double t = fabs(sqrt((nu_eff/nu)*K));
    
    if (K == 0)
    {
        ErrThrow("The permeability K is zero. Division by zero not allowed. ");
    }
    else if (nu == 0)
    {
        ErrThrow("The viscosity nu is zero. Division by zero not allowed. ");
    }
    
    if (t == 0)
    {
        values[0] = 1;                                                                  //u1
        values[1] = 0;                                                                  //u1_x
        values[2] = 0;                                                                  //u1_y
        values[3] = 0;                                                                  //Delta u1
    }
    else
    {
        values[0] = (1+exp(1/t)-exp((1-y)/t) - exp(y/t)) / (1+exp(1/t));             //u1
        values[1] = 0;                                                               //u1_x
        values[2] = (exp((1-y)/t)-exp(y/t))/(t*(1+exp(1/t)));                        //u1_y
        values[3] = (-exp((1-y)/t)-exp(y/t))/(t*t*(1+exp(1/t)));                     //Delta u1
    }
}

void ExactU2(double x, double y, double *values)
{
    values[0] = 0;            //u2
    values[1] = 0;            //u2_x
    values[2] = 0;            //u2_y
    values[3] = 0;            //Delta u2
}

void ExactP(double x, double y, double *values)
{
    values[0] = (0.5-x);                    //p
    values[1] = -1;                         //p_x
    values[2] = 0;                          //p_y
    values[3] = 0;                          //Delta p=p_xx+p_yy
}

// ========================================================================
// boundary conditions
// ========================================================================

void BoundCondition(int i, double Param, BoundCond &cond)
{
    cond = DIRICHLET; // default

    for (int j=0; j<TDatabase::ParamDB->n_neumann_boundary; j++)
    {
        if (i==TDatabase::ParamDB->neumann_boundary_id[j])
        {
            cond = NEUMANN;
            
            return;
        }
    }
    for (int j=0; j<TDatabase::ParamDB->n_nitsche_boundary; j++)
    {
        if (i==TDatabase::ParamDB->nitsche_boundary_id[j])
        {
            // Todo
            //cond = DIRICHLET_WEAK;
            return;
        }
    }
    
    if (TDatabase::ParamDB->n_neumann_boundary==0)
    {
        TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 1;
    }
    else
    {
        TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 0;
    }
}

void U1BoundValue(int BdComp, double Param, double &value)
{
    double K = TDatabase::ParamDB->PERMEABILITY;
    double nu = TDatabase::ParamDB->VISCOSITY;
    double nu_eff = TDatabase::ParamDB->EFFECTIVE_VISCOSITY;
    double t = fabs(sqrt((nu_eff/nu)*K));

    // loop to impose Neumann boundary conditions
    for (int j=0; j<TDatabase::ParamDB->n_neumann_boundary; j++)
    {
        if ( BdComp==TDatabase::ParamDB->neumann_boundary_id[j])
        {
            switch(BdComp)
            {
                case 1:
                    value=0.;//TDatabase::ParamDB->neumann_boundary_value[j];
                    break;
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
        case 0: value=0;
            break;
        case 1: value=(1+exp(1/t)-exp((1-Param)/t) - exp(Param/t)) / (1+exp(1/t));
            break;
        case 2: value=0;
            break;
        case 3: value= (1+exp(1/t)-exp((Param)/t) - exp((1-Param)/t)) / (1+exp(1/t)); //value=(1+exp(1/t)-exp((1-Param)/t) - exp(Param/t)) / ((nu/K)*(1+exp(1/t)));
            break;
            default: cout << "No boundary component with this number." << endl;
            break;
    }
}

void U2BoundValue(int BdComp, double Param, double &value)
{
    value=0;
}


// ========================================================================
// coefficients for Brinkman problem: viscosity, effective viscosity, permeability, f1, f2, g
// (the lhs of the Brinkman problem computed at quadrature points - for the error norms)
// ========================================================================
void LinCoeffs(int n_points, double *x, double *y,
               double **parameters, double **coeffs)
{
    static double eps = 1./TDatabase::ParamDB->RE_NR;
    static double nu     = TDatabase::ParamDB->VISCOSITY;
    static double nu_eff = TDatabase::ParamDB->EFFECTIVE_VISCOSITY;
    static double K      = TDatabase::ParamDB->PERMEABILITY;
    
    double *coeff;
    
    for(int i=0;i<n_points;i++)
    {
        coeff = coeffs[i];
        
        coeff[0] = eps;
        coeff[1] = 0;////(nu/K)-1;                                       // f1 (rhs of Brinkman problem for u1)
        coeff[2] = 0;                                       // f2 (rhs of Brinkman problem for u2)
        coeff[3] = 0;                                       // g (divergence term=u1_x+u2_y)
        coeff[4]=TDatabase::ParamDB->VISCOSITY;
        //    coeff[5]=TDatabase::ParamDB->EFFECTIVE_VISCOSITY;
        //    coeff[6]=TDatabase::ParamDB->PERMEABILITY;
        
        // effective viscosity unsteady
        if((x[i]-0.5)*(x[i]-0.5)+(y[i]-0.5)*(y[i]-0.5)> 0.09)
        {
            coeffs[i][5]= TDatabase::ParamDB->EFFECTIVE_VISCOSITY;
            coeffs[i][6]=TDatabase::ParamDB->PERMEABILITY;
        }
        else {
            coeffs[i][5]= TDatabase::ParamDB->EFFECTIVE_VISCOSITY/10000;
            //coeffs[i][5]= TDatabase::ParamDB->EFFECTIVE_VISCOSITY ;
            coeffs[i][6]=TDatabase::ParamDB->PERMEABILITY/10000;
        }
        coeff[7]=TDatabase::ParamDB->equal_order_stab_weight_P1P1;
        coeff[8]=TDatabase::ParamDB->equal_order_stab_weight_P2P2;
    }
}



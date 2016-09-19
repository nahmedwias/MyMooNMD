// ======================================================================
// @(#)Brinkman2D_Mixed.C        1.3 06/27/00
//
// common declaration for all Brinkman problems
//
// -nu \Delta (u1,u2) + \nabla p + K (u1,u2) = (f1,f2)
// \nabla \dcot (u1,u2) = g
//resp.
// -nu_eff \Delta (u1,u2) + \nabla p + nu/K (u1,u2) = (f1,f2)
// \nabla \dcot (u1,u2) = g
// ======================================================================
//#include <Convolution.h>
#include <Database.h>
#include <Brinkman2D_Mixed.h>

// ======================================================================
// Type 4, Standard Galerkin for Brinkman in [p div u] formulation
// Matrix Type:
// A11 0   B1^T
// 0   A22 B2^T
// B1  B2  0
// ======================================================================

// ======================================================================
// Type 1, Standard Galerkin for Brinkman in [p div v] formulation
// ======================================================================
void BrinkmanType1Galerkin(double Mult, double *coeff,
                           double *param, double hK,
                           double **OrigValues, int *N_BaseFuncts,
                           double ***LocMatrices, double **LocRhs)
{
    double ansatz00, ansatz10, ansatz01;    // ansatz functions
    double test00, test10, test01;          // test functions
    
    double ** MatrixA11 = LocMatrices[0];
    //double ** MatrixA12 = LocMatrices[1];
    //double ** MatrixA21 = LocMatrices[2];
    double ** MatrixA22 = LocMatrices[3];
    //double ** MatrixC = LocMatrices[4];
    double ** MatrixB1 = LocMatrices[5];
    double ** MatrixB2 = LocMatrices[6];
    double ** MatrixB1T = LocMatrices[7];
    double ** MatrixB2T = LocMatrices[8];
    
    double * Rhs_u1 = LocRhs[0];            // f_u1
    double * Rhs_u2 = LocRhs[1];            // f_u2
    // double * Rhs_p = LocRhs[2];          // g_q
    
    int N_U = N_BaseFuncts[0];              // number of basis functions for the velocity space
    int N_P = N_BaseFuncts[1];              // number of basis functions for the pressure space
    
    // values of fe functions at quadrature points: Origvalues = f(uk,vk) (Gauß-quadrature: \int f(u,v) = sum wk * f(uk,vk) at Gauß points)
    // Mult is the quadrature weight (wk)
    double *Orig0 = OrigValues[0];          // u_x
    double *Orig1 = OrigValues[1];          // u_y
    double *Orig2 = OrigValues[2];          // u
    double *Orig3 = OrigValues[3];          // p
    
    // values defined in the example
    //double c0 = coeff[0];                   // nu bzw eps (viscosity)
    double c1 = coeff[1];                   // f1
    double c2 = coeff[2];                   // f2
    //double c3 = coeff[3];                   // f3 (the rhs of incompressibility constraint)
    double nu = coeff[4];                   // viscosity
    double nu_eff = coeff[5];               // effective viscosity
    double K = coeff[6];                    // permeability
    
    double val;
    for(int i=0;i<N_U;i++)
    {
        test10 = Orig0[i];
        test01 = Orig1[i];
        test00 = Orig2[i];
        
        Rhs_u1[i] += Mult*test00*c1;
        Rhs_u2[i] += Mult*test00*c2;
        
        for(int j=0;j<N_U;j++)
        {
            double ansatz10 = Orig0[j]; // ansatz functions
            double ansatz01 = Orig1[j];
            double ansatz00 = Orig2[j];
            
            val  = nu_eff * (test10 * ansatz10 + test01 * ansatz01);        // nu_eff*(grad u,grad v)=nu_eff*(v_x*u_x + v_y*u_y)
            val += (nu/K) * (ansatz00 * test00);                            // nu/K*(u,v)
            MatrixA11[i][j] += Mult * val;
            
            val  = nu_eff* (test10 * ansatz10 + test01 * ansatz01);         // nu_eff*(grad u,grad v)=nu_eff*(v_x*u_x + v_y*u_y)
            val += (nu/K) * (ansatz00 * test00);                            // nu/K*(u,v)
            MatrixA22[i][j] += Mult * val;
        }
        
        for(int j=0;j<N_P;j++)
        {
            ansatz00 = Orig3[j];
            
            val = -Mult * ansatz00 * test10;                      // -(p,div v)=-p*v_x
            MatrixB1T[i][j] += val;
            
            val = -Mult * ansatz00 * test01;                      // -(p,div v)=-p*v_y
            MatrixB2T[i][j] += val;
        }
    }                              // endfor i
    
    for(int i=0;i<N_P;i++)
    {
        test00 = Orig3[i];
        
        for(int j=0;j<N_U;j++)
        {
            ansatz10 = Orig0[j];
            ansatz01 = Orig1[j];
            
            val = -Mult * test00 * ansatz10;                       // -(q,div u)=-q*u_x
            MatrixB1[i][j] += val;
            
            val = -Mult * test00 * ansatz01;                       // -(q,div u)=-q*u_y
            MatrixB2[i][j] += val;
        }                            // endfor j
        
        //for MatrixC: loop over U_P here
    }                                // endfor i
}

// ======================================================================
// Type 1b, Standard Galerkin for Scaled Brinkman in [p div v] formulation
// ======================================================================
void BrinkmanType1bGalerkin(double Mult, double *coeff,
                            double *param, double hK,
                            double **OrigValues, int *N_BaseFuncts,
                            double ***LocMatrices, double **LocRhs)
{
    double ansatz00, ansatz10, ansatz01;    // ansatz functions
    double test00, test10, test01;          // test functions
    
    double ** MatrixA11 = LocMatrices[0];
    //double ** MatrixA12 = LocMatrices[1];
    //double ** MatrixA21 = LocMatrices[2];
    double ** MatrixA22 = LocMatrices[3];
    //double ** MatrixC = LocMatrices[4];
    double ** MatrixB1 = LocMatrices[5];
    double ** MatrixB2 = LocMatrices[6];
    double ** MatrixB1T = LocMatrices[7];
    double ** MatrixB2T = LocMatrices[8];
    
    double * Rhs_u1 = LocRhs[0];            // f_u1
    double * Rhs_u2 = LocRhs[1];            // f_u2
    // double * Rhs_p = LocRhs[2];          // g_q
    
    int N_U = N_BaseFuncts[0];              // number of basis functions for the velocity space
    int N_P = N_BaseFuncts[1];              // number of basis functions for the pressure space
    
    // values of fe functions at quadrature points: Origvalues = f(uk,vk) (Gauß-quadrature: \int f(u,v) = sum wk * f(uk,vk) at Gauß points)
    // Mult is the quadrature weight (wk)
    double *Orig0 = OrigValues[0];          // u_x
    double *Orig1 = OrigValues[1];          // u_y
    double *Orig2 = OrigValues[2];          // u
    double *Orig3 = OrigValues[3];          // p
    
    // values defined in the example
    //double c0 = coeff[0];                   // nu bzw eps (viscosity)
    double c1 = coeff[1];                   // f1
    double c2 = coeff[2];                   // f2
    //double c3 = coeff[3];                   // f3 (the rhs of incompressibility constraint)
    double nu = coeff[4];                   // viscosity
    double nu_eff = coeff[5];               // effective viscosity
    double K = coeff[6];                    // permeability
    
    double val;
    for(int i=0;i<N_U;i++)
    {
        test10 = Orig0[i];
        test01 = Orig1[i];
        test00 = Orig2[i];
        
        Rhs_u1[i] += Mult*test00*c1;
        Rhs_u2[i] += Mult*test00*c2;
        
        for(int j=0;j<N_U;j++)
        {
            double ansatz10 = Orig0[j]; // ansatz functions
            double ansatz01 = Orig1[j];
            double ansatz00 = Orig2[j];
            
            val  = nu_eff*(test10*ansatz10+test01*ansatz01);      // nu*(v_x*u_x + v_y*u_y)
            val += (nu/K) *(ansatz00*test00);                       // K*(u * v)
            MatrixA11[i][j] += Mult * val;
            
            val  = nu_eff*(test10*ansatz10+test01*ansatz01);      // nu*(v_x*u_x + v_y*u_y)
            val += (nu/K)*(ansatz00*test00);                       // K*(u * v)
            MatrixA22[i][j] += Mult * val;
        }                            // endfor j
        
        for(int j=0;j<N_P;j++)
        {
            ansatz00 = Orig3[j];
            
            val = -Mult*ansatz00*test10;                      // -Mult* p*v_x
            MatrixB1T[i][j] += val;
            
            val = -Mult*ansatz00*test01;                      // -Mult*p*v_y
            MatrixB2T[i][j] += val;
        }
    }                              // endfor i
    
    for(int i=0;i<N_P;i++)
    {
        test00 = Orig3[i];
        
        for(int j=0;j<N_U;j++)
        {
            ansatz10 = Orig0[j];
            ansatz01 = Orig1[j];
            
            val = Mult*test00*ansatz10;                       // Mult*q*u_x
            MatrixB1[i][j] += val;
            
            val = Mult*test00*ansatz01;                       // Mult*q*u_y
            MatrixB2[i][j] += val;
        }                            // endfor j
        
        //for MatrixC: loop over U_P here
    }                                // endfor i
}

// ======================================================================
// Type 2, Standard Galerkin for Brinkman in [u gradq] formulation (macht nicht wirklich Sinn!!!!!!)
// ======================================================================
void BrinkmanType2Galerkin(double Mult, double *coeff,
                           double *param, double hK,
                           double **OrigValues, int *N_BaseFuncts,
                           double ***LocMatrices, double **LocRhs)
{
    double **MatrixA11 = LocMatrices[0];
    //double **MatrixA12 = LocMatrices[1];
    //double **MatrixA21 = LocMatrices[2];
    double **MatrixA22 = LocMatrices[3];
    ////int offset = TDatabase::ParamDB->NSTYPE == 14 ? 1 : 0;
    //double ** MatrixC = LocMatrices[4];
    int offset = 1;
    double **MatrixB1 = LocMatrices[4+offset];
    double **MatrixB2 = LocMatrices[5+offset];
    double **MatrixB1T = LocMatrices[6+offset];
    double **MatrixB2T = LocMatrices[7+offset];
    
    double *Rhs1 = LocRhs[0];                       // f1
    double *Rhs2 = LocRhs[1];                       // f2
    // double *Rhs_p = LocRhs[2];                   // g_q
    
    int N_U = N_BaseFuncts[0];
    int N_P = N_BaseFuncts[1];
    
    double *Orig0 = OrigValues[0];         // u_x
    double *Orig1 = OrigValues[1];         // u_y
    double *Orig2 = OrigValues[2];         // u
    //double *Orig3 = OrigValues[3];         // p
    double *Orig4 = OrigValues[4];         // p_x   MUSS NOCH DEFINIERT WERDEN !!!...eventuell
    double *Orig5 = OrigValues[5];         // p_y   MUSS NOCH DEFINIERT WERDEN!!!...eventuell
    
    //double c0 = coeff[0];                   // nu
    double c1 = coeff[1];                   // f1
    double c2 = coeff[2];                   // f2
    //double c3 = coeff[3];                   // f3 (the rhs of incompressibility constraint)
    double nu = coeff[4];                   // viscosity
    double nu_eff = coeff[5];               // effective viscosity
    double K = coeff[6];                    // permeability
    
    double val;
    double ansatz00, ansatz10, ansatz01;
    double test00, test10, test01;
    
    for(int i=0;i<N_U;i++)
    {
        double *Matrix11Row = MatrixA11[i];
        // double *Matrix12Row = MatrixA12[i];
        // double *Matrix21Row = MatrixA21[i];
        double *Matrix22Row = MatrixA22[i];
        
        test10 = Orig0[i];
        test01 = Orig1[i];
        test00 = Orig2[i];
        
        Rhs1[i] += Mult*test00*c1;
        Rhs2[i] += Mult*test00*c2;
        
        for(int j=0;j<N_U;j++)
        {
            ansatz10 = Orig0[j];
            ansatz01 = Orig1[j];
            ansatz00 = Orig2[j];
            
            val  = nu_eff*(test10*ansatz10+test01*ansatz01);
            val += (nu/K)* (ansatz00*test00);
            Matrix11Row[j] += Mult * val;
            
            // val  = 0;
            // Matrix12Row[j] += Mult * val;
            
            // val  = 0;
            // Matrix21Row[j] += Mult * val;
            
            val  = nu_eff*(test10*ansatz10+test01*ansatz01);
            val += (nu/K) *(ansatz00*test00);
            Matrix22Row[j] += Mult * val;
        }
        
        double *MatrixRow1 = MatrixB1T[i];
        double *MatrixRow2 = MatrixB2T[i];
        
        for(int j=0;j<N_P;j++)
        {
            ansatz10 = Orig4[j];
            ansatz01 = Orig5[j];
            
            val = Mult*ansatz10*test00;        // p_x*v
            MatrixRow1[j] += val;
            
            val = Mult*ansatz01*test00;        // p_y*v
            MatrixRow2[j] += val;
        }
    }
    
    for(int i=0;i<N_P;i++)
    {
        double *MatrixRow1 = MatrixB1[i];
        double *MatrixRow2 = MatrixB2[i];
        
        test10 = Orig4[i];
        test01 = Orig5[i];
        
        for(int j=0;j<N_U;j++)
        {
            ansatz00 = Orig2[j];
            
            val = -Mult*ansatz00*test10;
            MatrixRow1[j] += val;
            
            val = -Mult*ansatz00*test01;
            MatrixRow2[j] += val;
        }
    }
}

// ======================================================================
// Brinkman Problem with PSPG Stabilization (only P1/P1)
// assume Delta u = Delta v = 0
// according to Hannukainen: Computations with Finite element methods for the Brinkman Problem (2011)
// ======================================================================

void BrinkmanType1GalerkinResidualStab(double Mult, double *coeff,
                                       double *param, double hK,
                                       double **OrigValues, int *N_BaseFuncts,
                                       double ***LocMatrices, double **LocRhs)
{
    double ansatz00, ansatz10, ansatz01;    // ansatz functions
    double test00, test10, test01;          // test functions
    
    double ** MatrixA11 = LocMatrices[0];
    // double ** MatrixA12 = LocMatrices[1];
    // double ** MatrixA21 = LocMatrices[2];
    double ** MatrixA22 = LocMatrices[3];
    double ** MatrixC = LocMatrices[4];     //stabilization
    double ** MatrixB1 = LocMatrices[5];
    double ** MatrixB2 = LocMatrices[6];
    double ** MatrixB1T = LocMatrices[7];
    double ** MatrixB2T = LocMatrices[8];
    
    double * Rhs_u1 = LocRhs[0];            // f_u1
    double * Rhs_u2 = LocRhs[1];            // f_u2
    double * Rhs_p = LocRhs[2];             // g_q bzw f_q
    
    int N_U = N_BaseFuncts[0];              // number of basis functions for the velocity space
    int N_P = N_BaseFuncts[1];              // number of basis functions for the pressure space
    
    // values of fe functions at quadrature points: Origvalues = f(uk,vk) (Gauß-quadrature: \int f(u,v) = sum wk * f(uk,vk) at Gauß points)
    // Mult is the quadrature weight (wk)
    double *Orig0 = OrigValues[0];          // u_x
    double *Orig1 = OrigValues[1];          // u_y
    double *Orig2 = OrigValues[2];          // u
    double *Orig3 = OrigValues[3];          // p
    double *Orig4 = OrigValues[4];          // p_x
    double *Orig5 = OrigValues[5];          // p_y
    
    // values defined in the example
    // double c0 = coeff[0];                // t^2
    double c1 = coeff[1];                   // f1
    double c2 = coeff[2];                   // f2
    //double c3 = coeff[3];                   // f3 (the rhs of incompressibility constraint)
    double nu = coeff[4];                   // viscosity
    double nu_eff = coeff[5];               // effective viscosity
    double K = coeff[6];                    // permeability
    double alpha = coeff[7];
    double t = fabs(sqrt((nu_eff/nu) * K));
    double PSPGStab = -alpha * (hK*hK)/(t*t+hK*hK); // stabilization = (hK*hK)/(c0*c0+hK*hK) ///warum negativ, wenn positiv geht es schief-NOCHMAL TESTEN
    
    
    double val;
    for(int i=0;i<N_U;i++)
    {
        test10 = Orig0[i];
        test01 = Orig1[i];
        test00 = Orig2[i];
        
        Rhs_u1[i] += Mult * test00*c1;                                    // (f,v)
        Rhs_u2[i] += Mult * test00*c2;                                    // (f,v)
        Rhs_u1[i] += Mult * PSPGStab * (nu/K) * test00 * c1;               // stabilization: (nu/K)*(f,v)
        Rhs_u2[i] += Mult * PSPGStab * (nu/K) * test00 * c2;               // stabilization: (nu/K)*(f,v)
        
        for(int j=0;j<N_U;j++)
        {
            ansatz10 = Orig0[j];
            ansatz01 = Orig1[j];
            ansatz00 = Orig2[j];
            
            val  = nu_eff * (test10 * ansatz10 + test01 * ansatz01);                // nu_eff*(grad u,grad v)=nu_eff*(v_x*u_x + v_y*u_x)
            val += (nu/K) * (ansatz00 * test00);                                // nu/K *(u,v)
            val += PSPGStab * (nu/K) * (nu/K) * (ansatz00 * test00);                // stabilization: (nu/K)^2*(u,v)
            MatrixA11[i][j] += Mult * val;
            
            // val  = 0;
            // Matrix12Row[j] += Mult * val;
            
            // val  = 0;
            // Matrix21Row[j] += Mult * val;
            
            val  = nu_eff * (test10 * ansatz10 + test01 * ansatz01);            // nu_eff*(grad u,grad v)=nu_eff*(v_x*u_x + v_y*u_x)
            val += (nu/K) * (ansatz00 * test00);                            // nu/K *(u,v)
            val += PSPGStab * (nu/K) * (nu/K) * (ansatz00*test00);            // stabilization: (nu/K)^2*(u,v)
            MatrixA22[i][j] += Mult * val;
        }
        
        
        for(int j=0;j<N_P;j++)
        {
            ansatz00 = Orig3[j];
            ansatz10 = Orig4[j];
            ansatz01 = Orig5[j];
            
            val  = -ansatz00 * test10;                                     // -(p,div v)=-p*v_x
            val += PSPGStab * (nu/K) * (ansatz10 * test00);                    // stabilization: nu/K * (grad p,v)
            MatrixB1T[i][j] += Mult * val;
            
            val = -ansatz00*test01;                                     // -(p,div v)=-p*v_y
            val+= PSPGStab*(nu/K)*(ansatz01*test00);                    // stabilization: nu/K * (grad p,v)
            MatrixB2T[i][j] += Mult * val;
        }
    }
    
    for(int i=0;i<N_P;i++)
    {
        test00 = Orig3[i];
        test10 = Orig4[i];
        test01 = Orig5[i];
        
        Rhs_p[i] += Mult * PSPGStab * test10*c1;                    // stabilization: (f,grad q)
        Rhs_p[i] += Mult * PSPGStab * test01*c2;                    // stabilization: (f,grad q)
        
        for(int j=0;j<N_U;j++)
        {
            ansatz00 = Orig2[j];
            ansatz10 = Orig0[j];
            ansatz01 = Orig1[j];
            
            val  = -test00 * ansatz10;                                          // -(q,div u)=-q*u_x       WARUM MINUS??????
            val += PSPGStab * (nu/K) * (test10 * ansatz00);                     // stabilization: nu/K * (u,grad q)
            MatrixB1[i][j] += Mult * val;
            
            val  = -test00 * ansatz01;                                          // -(q,div u)=-q*u_y             WARUM MINUS????
            val += PSPGStab * (nu/K) * (test01 * ansatz00);                     // stabilization: nu/K * (u,grad q)
            MatrixB2[i][j] += Mult * val;
        }
        
        //for MatrixC: loop over N_P here
        
        for(int j=0;j<N_P;j++)
        {
            ansatz00 = Orig3[j];
            ansatz10 = Orig4[j];
            ansatz01 = Orig5[j];
            
            val = PSPGStab * (ansatz10 * test10 + ansatz01 * test01);               // stabilization: (grad p,grad q)
            MatrixC[i][j] += Mult * val;
        }
    }
}

// ======================================================================
// Brinkman Problem with PSPG Stabilization (for P2/P2)
// takes into account Delta u and Delta v
// according to Hannukainen:Computations with finite element methods for the Brinkman Problem (2011)
// ======================================================================

void BrinkmanType1GalerkinResidualStab2(double Mult, double *coeff,
                                        double *param, double hK,
                                        double **OrigValues, int *N_BaseFuncts,
                                        double ***LocMatrices, double **LocRhs)
{
    double ansatz00, ansatz10, ansatz01, ansatz20, ansatz02;    // ansatz functions
    double test00, test10, test01, test20, test02;              // test functions
    
    double ** MatrixA11 = LocMatrices[0];
    double ** MatrixA12 = LocMatrices[1];                       // stabilization P2/P2
    double ** MatrixA21 = LocMatrices[2];                       // stabilization P2/P2
    double ** MatrixA22 = LocMatrices[3];
    double ** MatrixC   = LocMatrices[4];                       // stabilization P1/P1
    double ** MatrixB1  = LocMatrices[5];
    double ** MatrixB2  = LocMatrices[6];
    double ** MatrixB1T = LocMatrices[7];
    double ** MatrixB2T = LocMatrices[8];
    
    double * Rhs_u1 = LocRhs[0];            // f_u1         (rhs tested with v1)
    double * Rhs_u2 = LocRhs[1];            // f_u2         (rhs tested with v2)
    double * Rhs_p  = LocRhs[2];            // g_q resp. f_q  (rhs tested with q)
    
    int N_U = N_BaseFuncts[0];              // number of basis functions for the velocity space (both for u1 and u2)
    int N_P = N_BaseFuncts[1];              // number of basis functions for the pressure space
    
    // values of fe functions at quadrature points: Origvalues = f(uk,vk) (Gauß-quadrature: \int f(u,v) = sum wk * f(uk,vk) at Gauß points)
    // Mult is the quadrature weight (wk)
    double *Orig0 = OrigValues[0];          // u_x
    double *Orig1 = OrigValues[1];          // u_y
    double *Orig2 = OrigValues[2];          // u
    double *Orig3 = OrigValues[3];          // p
    double *Orig4 = OrigValues[4];          // p_x
    double *Orig5 = OrigValues[5];          // p_y
    double *Orig6 = OrigValues[6];          // u_xx
    double *Orig7 = OrigValues[7];          // u_yy
    
    // values defined in the example
    //double c0 = coeff[0];                               // t^2= eps(:=1/Re) ; Re can be set in brinkman2d.dat
    double c1 = coeff[1];                               // f1
    double c2 = coeff[2];                               // f2
    //double c3 = coeff[3];                               // f3 (the rhs of incompressibility constraint)
    double nu = coeff[4];                               // viscosity
    double nu_eff = coeff[5];                           // effective viscosity
    double K = coeff[6];                                // permeability
    double t = fabs(sqrt((nu_eff/nu)*K));
    double alpha =coeff[8];                             // PSPG Stabilization Parameter
    double PSPGStab = -alpha * (hK*hK)/(t*t+hK*hK);     //stabilization = (hK*hK)/(c0*c0+hK*hK)
    
    double val;
    for(int i=0;i<N_U;i++)
    {
        test10 = Orig0[i];          // v_x
        test01 = Orig1[i];          // v_y
        test00 = Orig2[i];          // v
        test20 = Orig6[i];          // v_xx
        test02 = Orig7[i];          // v_yy
        
        Rhs_u1[i] += Mult * test00 * c1;                                            // (f,v)
        Rhs_u2[i] += Mult * test00 * c2;                                            // (f,v)
        Rhs_u1[i] += Mult * PSPGStab * nu/K * test00 * c1;                          // stabilization P1/P1: nu/K * (f,v)
        Rhs_u2[i] += Mult * PSPGStab * nu/K * test00 * c2;                          // stabilization P1/P1: nu/K * (f,v)
        Rhs_u1[i] += -Mult* PSPGStab * nu_eff * c1 * (test20 + test02);             // stabilization P2/P2: -nu_eff * (f, Delta v)
        Rhs_u2[i] += -Mult* PSPGStab * nu_eff * c2 * (test20 + test02);             // stabilization P2/P2: -nu_eff * (f, Delta v)
        
        for(int j=0;j<N_U;j++)
        {
            ansatz10 = Orig0[j];    // u_x
            ansatz01 = Orig1[j];    // u_y
            ansatz00 = Orig2[j];    // u
            ansatz20 = Orig6[j];    // u_xx
            ansatz02 = Orig7[j];    // u_yy
            
            val  = nu_eff * (test10 * ansatz10 + test01 * ansatz01);                                // nu_eff*(grad u,grad v)=nu_eff*(v_x*u_x + v_y*u_x)
            val += (nu/K) * (ansatz00 * test00);                                                    // nu/K * (u,v)
            val += PSPGStab * (nu/K) * (nu/K) * (ansatz00 * test00);                                // stabilization P1/P1: (nu/K)^2*(u,v)
            val += PSPGStab * nu_eff * nu_eff * ((ansatz20 + ansatz02) * (test20 + test02));        // stabilization P2/P2: (nu_eff)^2* (Delta u,Delta v)
            MatrixA11[i][j] += Mult * val;
            
            // val  = 0;
            // Matrix12Row[j] += Mult * val;
            
            val  = -PSPGStab * (nu_eff * nu)/K * (ansatz20 + ansatz02) * test00;                // stabilization P2/P2: -(nu_eff*nu)/K * (Delta u,v)
            val += -PSPGStab * (nu_eff * nu)/K * (test20 + test02) * ansatz00;                  // stabilization P2/P2: -(nu_eff*nu)/K * (Delta v,u)
            MatrixA12[i][j] += Mult * val;
            
            // val  = 0;
            // Matrix21Row[j] += Mult * val;
            
            val  = -PSPGStab * (nu_eff * nu)/K * (ansatz20 + ansatz02) * test00;                   // stabilization P2/P2: -(nu_eff*nu)/K * (Delta u,v)
            val += -PSPGStab * (nu_eff * nu)/K  * (test20 + test02) * ansatz00;                    // stabilization P2/P2: -(nu_eff*nu)/K * (Delta v,u)
            MatrixA21[i][j] += Mult * val;
            
            val  = nu_eff * (test10 * ansatz10 + test01 * ansatz01);                                // nu_eff*(grad u,grad v)=nu_eff*(v_x*u_x + v_y*u_x)
            val += (nu/K) * ansatz00 * test00;                                                      // nu/K * (u , v)
            val += PSPGStab * (nu/K) * (nu/K) * (ansatz00 * test00);                                // stabilization P1/P1: stabilization P1/P1: (nu/K)^2*(u,v)
            val += PSPGStab * nu_eff * nu_eff * ((ansatz20 + ansatz02) * (test20 + test02));        // stabilization P2/P2: (nu_eff)^2* (Delta u,Delta v)
            MatrixA22[i][j] += Mult * val;
        }
        
        for(int j=0;j<N_P;j++)
        {
            ansatz00 = Orig3[j];        // p
            ansatz10 = Orig4[j];        // p_x
            ansatz01 = Orig5[j];        // p_y
            
            val  = -ansatz00 * test10;                                                          // -(p,div v)=-p*v_x
            val += PSPGStab * nu/K * (ansatz10 * test00);                                       // stabilization P1/P1: nu/K * (v, grad p)
            val += -PSPGStab * nu_eff * (test20 + test02) * ansatz10;                           // stabilization P2/P2: -nu_eff * (Delta v, grad p)
            MatrixB1T[i][j] += Mult * val;
            
            val  = -ansatz00 * test01;                                                          // -(p,div v)=-p*v_y
            val += PSPGStab * nu/K * (ansatz01 * test00);                                       // stabilization P1/P1: nu/K * (v, grad p)
            val += -PSPGStab * nu_eff * (test20 + test02) * ansatz01;                           // stabilization P2/P2: -nu_eff * (Delta v, grad p)
            MatrixB2T[i][j] += Mult * val;
        }
    }
    
    for(int i=0;i<N_P;i++)
    {
        test00 = Orig3[i];          // q
        test10 = Orig4[i];          // q_x
        test01 = Orig5[i];          // q_y
        test20 = Orig6[i];
        test20 = Orig7[i];
        
        Rhs_p[i] += Mult * PSPGStab * test10 * c1;                                           // stabilization P1/P1: (f,grad q)
        Rhs_p[i] += Mult * PSPGStab * test01 * c2;                                           // stabilization P1/P1: (f,grad q)
        
        
        
        for(int j=0;j<N_U;j++)
        {
            ansatz00 = Orig2[j];        // u
            ansatz10 = Orig0[j];        // u_x
            ansatz01 = Orig1[j];        // u_y
            ansatz20 = Orig6[j];        // u_xx
            ansatz02 = Orig7[j];        // u_yy
            
            val  = -test00 * ansatz10;                                                      // -(q,div u)=-q*u_x     WARUM MINUS? (MASS CONSERVATION)
            val += PSPGStab * nu/K * (ansatz00 * test10);                                   // stabilization P1/P1: nu/K * (u,grad q)
            val += -PSPGStab * nu_eff *  (ansatz20 + ansatz02) * test10;                    // stabilization P2/P2: -nu_eff * (Delta u,grad q)
            MatrixB1[i][j] += Mult * val;
            
            val  = -test00 * ansatz01;                                                      // -(q,div u)=-q*u_y     WARUM MINUS? (MASS CONSERVATION)
            val += PSPGStab * nu/K * (ansatz00 * test01);                                   // stabilization P1/P1: nu/K * (u,grad q)
            val += -PSPGStab * nu_eff * (ansatz20 + ansatz02) * test01;                     // stabilization P2/P2: -nu_eff * (Delta u,grad q)
            MatrixB2[i][j] += Mult * val; 
        }
        
        // for MatrixC: loop over N_P here
        // (grad p, grad q)
        
        for(int j=0;j<N_P;j++)
        {
            ansatz00 = Orig3[j];            // p
            ansatz10 = Orig4[j];            // p_x
            ansatz01 = Orig5[j];            // p_y
            
            val = PSPGStab * (ansatz10 * test10+ansatz01*test01);                           // stabilization P1/P1: (grad p,grad q)
            MatrixC[i][j] += Mult * val;
        }
    }
}


//formulations with t^2=c0

// ======================================================================
// Brinkman Problem with PSPG Stabilization (only P1/P1)
// assume Delta u = Delta v = 0
// according to Hannukainen: Computations with Finite element methods for the Brinkman Problem (2011)
// ======================================================================
//
//void BrinkmanType1GalerkinResidualStab(double Mult, double *coeff,
//                                       double *param, double hK,
//                                       double **OrigValues, int *N_BaseFuncts,
//                                       double ***LocMatrices, double **LocRhs)
//{
//    double ansatz00, ansatz10, ansatz01;    // ansatz functions
//    double test00, test10, test01;          // test functions
//    
//    double ** MatrixA11 = LocMatrices[0];
//    //double ** MatrixA12 = LocMatrices[1];
//    //double ** MatrixA21 = LocMatrices[2];
//    double ** MatrixA22 = LocMatrices[3];
//    double ** MatrixC = LocMatrices[4];     //stabilization
//    double ** MatrixB1 = LocMatrices[5];
//    double ** MatrixB2 = LocMatrices[6];
//    double ** MatrixB1T = LocMatrices[7];
//    double ** MatrixB2T = LocMatrices[8];
//    
//    double * Rhs_u1 = LocRhs[0];            // f_u1
//    double * Rhs_u2 = LocRhs[1];            // f_u2
//    double * Rhs_p = LocRhs[2];             // g_q bzw f_q
//    
//    int N_U = N_BaseFuncts[0];              // number of basis functions for the velocity space
//    int N_P = N_BaseFuncts[1];              // number of basis functions for the pressure space
//    
//    // values of fe functions at quadrature points: Origvalues = f(uk,vk) (Gauß-quadrature: \int f(u,v) = sum wk * f(uk,vk) at Gauß points)
//    // Mult is the quadrature weight (wk)
//    double *Orig0 = OrigValues[0];          // u_x
//    double *Orig1 = OrigValues[1];          // u_y
//    double *Orig2 = OrigValues[2];          // u
//    double *Orig3 = OrigValues[3];          // p
//    double *Orig4 = OrigValues[4];          // p_x
//    double *Orig5 = OrigValues[5];          // p_y
//    
//    // values defined in the example
//    double c0 = coeff[0];                   // t^2
//    double c1 = coeff[1];                   // f1
//    double c2 = coeff[2];                   // f2
//    double c3 = coeff[3];                   // f3 (the rhs of incompressibility constraint)
//    double nu = coeff[4];                   //viscosity
//    double nu_eff = coeff[5];               // effective viscosity
//    double K = coeff[6];                    // permeability
//    double alpha = coeff[7];
//    double t = fabs(sqrt((nu_eff/nu)*K));
//    double PSPGStab = alpha*(hK*hK)/(t+hK*hK); //stabilization = (hK*hK)/(c0*c0+hK*hK) ///warum negativ, wenn positiv geht es schief-NOCHMAL TESTEN
//    
//    
//    double val;
//    for(int i=0;i<N_U;i++)
//    {
//        test10 = Orig0[i];
//        test01 = Orig1[i];
//        test00 = Orig2[i];
//        
//        Rhs_u1[i] += Mult*test00*c1;
//        Rhs_u2[i] += Mult*test00*c2;
//        Rhs_u1[i] += Mult* PSPGStab * test00 * c1;
//        Rhs_u2[i] += Mult* PSPGStab * test00 * c2;
//        
//        for(int j=0;j<N_U;j++)
//        {
//            ansatz10 = Orig0[j];
//            ansatz01 = Orig1[j];
//            ansatz00 = Orig2[j];
//            
//            val  = c0*(test10*ansatz10+test01*ansatz01);        // t^2*(v_x*u_x + v_y*u_x)
//            val += (ansatz00*test00);                           // (u * v)
//            val += PSPGStab*(ansatz00*test00);                  // stabilization
//            MatrixA11[i][j] += Mult * val;
//            
//            // val  = 0;
//            // Matrix12Row[j] += Mult * val;
//            
//            // val  = 0;
//            // Matrix21Row[j] += Mult * val;
//            
//            val  = c0*(test10*ansatz10+test01*ansatz01);        // t^2*(v_x*u_x + v_y*u_x)
//            val += (ansatz00*test00);                           // (u * v)
//            val += PSPGStab*(ansatz00*test00);                  // stabilization
//            MatrixA22[i][j] += Mult * val;
//        }                            // endfor j
//        
//        for(int j=0;j<N_P;j++)
//        {
//            ansatz00 = Orig3[j];
//            ansatz10 = Orig4[j];
//            ansatz01 = Orig5[j];
//            
//            val = -ansatz00*test10;                             // -Mult* p*v_x
//            val-= PSPGStab*(ansatz00*test10);                   // stabilization
//            
//            MatrixB1T[i][j] += Mult * val;
//            
//            val = -ansatz00*test01;                             // -Mult*p*v_y
//            val-= PSPGStab*(ansatz00*test01);                   // stabilization
//            
//            MatrixB2T[i][j] += Mult * val;
//        }
//    }                              // endfor i
//    
//    for(int i=0;i<N_P;i++)
//    {
//        test00 = Orig3[i];
//        test10 = Orig4[i];
//        test01 = Orig5[i];
//        
//        Rhs_p[i] += Mult*PSPGStab*test10*c1;                    // stabilization
//        Rhs_p[i] += Mult*PSPGStab*test01*c2;                    // stabilization
//        
//        for(int j=0;j<N_U;j++)
//        {
//            ansatz00 = Orig2[j];
//            ansatz10 = Orig0[j];
//            ansatz01 = Orig1[j];
//            
//            val = test00*ansatz10;                              // -Mult*q*u_x
//            val-= PSPGStab * (test00*ansatz10);                 // stabilization
//            MatrixB1[i][j] += Mult * val;
//            
//            val = test00*ansatz01;                              // -Mult*q*u_y
//            val-= PSPGStab * (test00*ansatz01);                 // stabilization
//            MatrixB2[i][j] += Mult * val;
//        }                            // endfor j
//        
//        //for MatrixC: loop over N_P here
//        
//        for(int j=0;j<N_P;j++)
//        {
//            ansatz00 = Orig3[j];
//            ansatz10 = Orig4[j];
//            ansatz01 = Orig5[j];
//            
//            val = PSPGStab*(ansatz10*test10);               // stabilization
//            MatrixC[i][j] += Mult * val;
//            
//            val= PSPGStab *(ansatz01*test01);               // stabilization
//            MatrixC[i][j] += Mult * val;
//        }
//    }
//}
//
//// ======================================================================
//// Brinkman Problem with PSPG Stabilization (for P2/P2)
//// takes into account Delta u and Delta v
//// according to Hannukainen:Computations with finite element methods for the Brinkman Problem (2011)
//// ======================================================================
//
//void BrinkmanType1GalerkinResidualStab2(double Mult, double *coeff,
//                                        double *param, double hK,
//                                        double **OrigValues, int *N_BaseFuncts,
//                                        double ***LocMatrices, double **LocRhs)
//{
//    double ansatz00, ansatz10, ansatz01, ansatz20, ansatz02;    // ansatz functions
//    double test00, test10, test01, test20, test02;              // test functions
//    
//    double ** MatrixA11 = LocMatrices[0];
//    double ** MatrixA12 = LocMatrices[1];                       // stabilization P2/P2
//    double ** MatrixA21 = LocMatrices[2];                       // stabilization P2/P2
//    double ** MatrixA22 = LocMatrices[3];
//    double ** MatrixC   = LocMatrices[4];                       //stabilization P1/P1
//    double ** MatrixB1  = LocMatrices[5];
//    double ** MatrixB2  = LocMatrices[6];
//    double ** MatrixB1T = LocMatrices[7];
//    double ** MatrixB2T = LocMatrices[8];
//    
//    double * Rhs_u1 = LocRhs[0];            // f_u1         (rhs tested with v1)
//    double * Rhs_u2 = LocRhs[1];            // f_u2         (rhs tested with v2)
//    double * Rhs_p  = LocRhs[2];            // g_q resp. f_q  (rhs tested with q)
//    
//    int N_U = N_BaseFuncts[0];              // number of basis functions for the velocity space (both for u1 and u2)
//    int N_P = N_BaseFuncts[1];              // number of basis functions for the pressure space
//    
//    // values of fe functions at quadrature points: Origvalues = f(uk,vk) (Gauß-quadrature: \int f(u,v) = sum wk * f(uk,vk) at Gauß points)
//    // Mult is the quadrature weight (wk)
//    double *Orig0 = OrigValues[0];          // u_x
//    double *Orig1 = OrigValues[1];          // u_y
//    double *Orig2 = OrigValues[2];          // u
//    double *Orig3 = OrigValues[3];          // p
//    double *Orig4 = OrigValues[4];          // p_x
//    double *Orig5 = OrigValues[5];          // p_y
//    double *Orig6 = OrigValues[6];          // u_xx
//    double *Orig7 = OrigValues[7];          // u_yy
//    
//    // values defined in the example
//    double c0 = coeff[0];                           // t^2= eps(:=1/Re) ; Re can be set in brinkman2d.dat
//    double c1 = coeff[1];                           // f1
//    double c2 = coeff[2];                           // f2
//    double c3 = coeff[3];                           // f3 (the rhs of incompressibility constraint)
//    double nu = coeff[4];
//    double nu_eff = coeff[5];
//    double K = coeff[6];                          // viscosity/permeability
//    double alpha =coeff[7];                       // PSPG Stabilization Parameter
//    double PSPGStab = -alpha*(hK*hK)/(c0*c0+hK*hK); //stabilization = (hK*hK)/(c0*c0+hK*hK)
//    
//    double val;
//    for(int i=0;i<N_U;i++)
//    {
//        test10 = Orig0[i];
//        test01 = Orig1[i];
//        test00 = Orig2[i];
//        test20 = Orig6[i];
//        test02 = Orig7[i];
//        
//        Rhs_u1[i] += Mult*test00*c1;
//        Rhs_u2[i] += Mult*test00*c2;
//        Rhs_u1[i] += Mult* PSPGStab * test00 * c1;              // stabilization P1/P1
//        Rhs_u2[i] += Mult* PSPGStab * test00 * c2;              // stabilization P1/P1
//        Rhs_u1[i] += -Mult*PSPGStab*c0*c1*(test20+test02);      // stabilization P2/P2
//        Rhs_u2[i] += -Mult*PSPGStab*c0*c2*(test20+test02);      // stabilization P2/P2
//        
//        for(int j=0;j<N_U;j++)
//        {
//            ansatz10 = Orig0[j];
//            ansatz01 = Orig1[j];
//            ansatz00 = Orig2[j];
//            ansatz20 = Orig6[j];
//            ansatz02 = Orig7[j];
//            
//            val  = c0*(test10*ansatz10+test01*ansatz01);                            // t^2*(v_x*u_x + v_y*u_x)
//            val += (ansatz00*test00);                                               // (u * v)
//            val += PSPGStab * (ansatz00*test00);                                    // stabilization P1/P1
//            val += PSPGStab * c0 * c0 * ((ansatz20+ansatz02)*(test20+test02));      // stabilization P2/P2
//            MatrixA11[i][j] += Mult * val;
//            
//            // val  = 0;
//            // Matrix12Row[j] += Mult * val;
//            
//            val  = -PSPGStab * c0 * (ansatz20+ansatz02) * test00;                  // stabilization P2/P2
//            val += -PSPGStab * c0 * (test20+test02) * ansatz00;                    // stabilization P2/P2
//            MatrixA12[i][j] += Mult * val;                                         // stabilization P2/P2
//            
//            // val  = 0;
//            // Matrix21Row[j] += Mult * val;
//            
//            val  = -PSPGStab * c0 * (ansatz20+ansatz02) * test00;                   // stabilization P2/P2
//            val += -PSPGStab * c0 * (test20+test02) * ansatz00;                     // stabilization P2/P2
//            MatrixA21[i][j] += Mult * val;                                          // stabilization P2/P2
//            
//            val  = c0*(test10*ansatz10+test01*ansatz01);                            // t^2*(v_x*u_x + v_y*u_x)
//            val += ansatz00*test00;                                                 // (u * v)
//            val += PSPGStab*(ansatz00*test00);                                      // stabilization P1/P1
//            val += PSPGStab * c0*c0* ((ansatz20+ansatz02)*(test20+test02));         // stabilization P2/P2
//            MatrixA22[i][j] += Mult * val;
//        }                            // endfor j
//        
//        for(int j=0;j<N_P;j++)
//        {
//            ansatz00 = Orig3[j];
//            ansatz10 = Orig4[j];
//            ansatz01 = Orig5[j];
//            
//            val = -ansatz00*test10;                                                 // -Mult* p*v_x
//            val+= PSPGStab*(ansatz10*test00);                                       // stabilization P1/P1
//            val+= -PSPGStab * c0 * (test20+test02) * ansatz10;                      // stabilization P2/P2
//            MatrixB1T[i][j] += Mult * val;
//            
//            val = -ansatz00*test01;                                                 // -Mult*p*v_y
//            val+= PSPGStab*(ansatz01*test00);                                       // stabilization P1/P1
//            val+= -PSPGStab * c0 * (test20+test02) * ansatz01;                      // stabilization P2/P2
//            MatrixB2T[i][j] += Mult * val;
//        }
//    }                              // endfor i
//    
//    for(int i=0;i<N_P;i++)
//    {
//        test00 = Orig3[i];
//        test10 = Orig4[i];
//        test01 = Orig5[i];
//        test20 = Orig6[i];
//        test20 = Orig7[i];
//        
//        Rhs_p[i] += Mult*PSPGStab*test10*c1;                                        // stabilization P1/P1
//        Rhs_p[i] += Mult*PSPGStab*test01*c2;                                        // stabilization P1/P1
//        
//        for(int j=0;j<N_U;j++)
//        {
//            ansatz00 = Orig2[j];
//            ansatz10 = Orig0[j];
//            ansatz01 = Orig1[j];
//            ansatz20 = Orig6[j];
//            ansatz02 = Orig7[j];
//            
//            val = -test00*ansatz10;                                                 // -Mult*q*u_x
//            val+= PSPGStab * (ansatz00*test10);                                     // stabilization P1/P1
//            val+= -PSPGStab * c0 * (ansatz20+ansatz02)*test10;                      // stabilization P2/P2
//            MatrixB1[i][j] += Mult * val;
//            
//            val = -test00*ansatz01;                                                 // -Mult*q*u_y
//            val+= PSPGStab * (ansatz00*test01);                                     // stabilization P1/P1
//            val+= -PSPGStab * c0 * (ansatz20+ansatz02)*test01;                      // stabilization P2/P2
//            MatrixB2[i][j] += Mult * val;
//        }                            // endfor j
//        
//        // for MatrixC: loop over N_P here
//        // (grad p, grad q)
//        
//        for(int j=0;j<N_P;j++)
//        {
//            ansatz00 = Orig3[j];
//            ansatz10 = Orig4[j];
//            ansatz01 = Orig5[j];
//            
//            val = PSPGStab*(ansatz10*test10);                                       // stabilization P1/P1
//            MatrixC[i][j] += Mult * val;
//            
//            val= PSPGStab *(ansatz01*test01);                                       // stabilization P1/P1
//            MatrixC[i][j] += Mult * val;
//        }
//    }
//}

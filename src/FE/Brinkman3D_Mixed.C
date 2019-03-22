// ======================================================================
// @(#)Brinkman3D_Mixed.C        1.3 11.10.2016
//
// common declaration for all Brinkman problems
//
// -nu \Delta (u1,u2,u3) + \nabla p + K (u1,u2,u3) = (f1,f2,f3)
// \nabla \dcot (u1,u2,u3) = g
//resp.
// -nu_eff \Delta (u1,u2,u3) + \nabla p + nu/K (u1,u2,u3) = (f1,f2,f3)
// \nabla \dcot (u1,u2,u3) = g = 0
// ======================================================================
//#include <Convolution.h>
#include <Database.h>
#include <Brinkman3D_Mixed.h>

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
// For NSTYPE: 14
void Brinkman3DType1Galerkin(double Mult, double *coeff,
                             double *, double,
                             double **OrigValues, int *N_BaseFuncts,
                             double ***LocMatrices, double **LocRhs)
{
    double ansatz000, ansatz100, ansatz010, ansatz001;      // ansatz functions
    double test000, test100, test010, test001;              // test functions
    
    double **MatrixA11 = LocMatrices[0];
    //double **MatrixA12 = LocMatrices[1];
    //double **MatrixA13 = LocMatrices[2];
    //double **MatrixA21 = LocMatrices[3];
    double **MatrixA22 = LocMatrices[4];
    //double **MatrixA23 = LocMatrices[5];
    //double **MatrixA31 = LocMatrices[6];
    //double **MatrixA32 = LocMatrices[7];
    double **MatrixA33 = LocMatrices[8];
    //double **MatrixC = LocMatrices[9];
    double **MatrixB1  = LocMatrices[10];
    double **MatrixB2  = LocMatrices[11];
    double **MatrixB3  = LocMatrices[12];
    double **MatrixB1T  = LocMatrices[13];
    double **MatrixB2T  = LocMatrices[14];
    double **MatrixB3T  = LocMatrices[15];
    
    double * Rhs_u1 = LocRhs[0];            // f_u1
    double * Rhs_u2 = LocRhs[1];            // f_u2
    double * Rhs_u3 = LocRhs[2];            // f_u3
    //double * Rhs_p = LocRhs[3];           // g_q
    
    int N_U = N_BaseFuncts[0];              // number of basis functions for the velocity space
    int N_P = N_BaseFuncts[1];              // number of basis functions for the pressure space
    
    // values of fe functions at quadrature points: Origvalues = f(uk,vk) (Gauß-quadrature: \int f(u,v) = sum wk * f(uk,vk) at Gauß points)
    // Mult is the quadrature weight (wk)
    double *Orig0 = OrigValues[0];          // u_x
    double *Orig1 = OrigValues[1];          // u_y
    double *Orig2 = OrigValues[2];          // u_z
    double *Orig3 = OrigValues[3];          // u
    double *Orig4 = OrigValues[4];          // p
    
    // values defined in the example
    //double c0 = coeff[0];                 // nu bzw eps (viscosity)
    double c1 = coeff[1];                   // f1
    double c2 = coeff[2];                   // f2
    double c3 = coeff[3];                   // f3 (the rhs of incompressibility constraint)
    //double c4 = coeff[4];                 // f4 (the rhs of incompressibility constraint)
    double nu = coeff[5];                   // viscosity
    double nu_eff = coeff[6];               // effective viscosity
    double K = coeff[7];                    // permeability
    
    double val;
    
    for(int i=0;i<N_U;i++)
    {
        // test functions
        test100 = Orig0[i]; //v_x
        test010 = Orig1[i]; //v_y
        test001 = Orig2[i]; //v_z
        test000 = Orig3[i]; //q
        
        Rhs_u1[i] += Mult*test000*c1;
        Rhs_u2[i] += Mult*test000*c2;
        Rhs_u3[i] += Mult*test000*c3;
        
        for(int j=0;j<N_U;j++)
        {
            // ansatz functions
            double ansatz100 = Orig0[j]; //u_x
            double ansatz010 = Orig1[j]; //u_y
            double ansatz001 = Orig2[j]; //u_z
            double ansatz000 = Orig3[j]; //p
            
            val  = nu_eff * (test100 * ansatz100 + test010 * ansatz010+ test001 * ansatz001);        // nu_eff*(grad u,grad v)=nu_eff*(v_x*u_x + v_y*u_y+ v_z*u_z)
            val += (nu/K) * (ansatz000 * test000);                            // nu/K*(u,v)
            MatrixA11[i][j] += Mult * val;
            
            val  = nu_eff * (test100 * ansatz100 + test010 * ansatz010+ test001 * ansatz001);        // nu_eff*(grad u,grad v)=nu_eff*(v_x*u_x + v_y*u_y+ v_z*u_z)
            val += (nu/K) * (ansatz000 * test000);                          // nu/K*(u,v)
            MatrixA22[i][j] += Mult * val;
            
            val  = nu_eff * (test100 * ansatz100 + test010 * ansatz010+ test001 * ansatz001);        // nu_eff*(grad u,grad v)=nu_eff*(v_x*u_x + v_y*u_y+ v_z*u_z)
            val += (nu/K) * (ansatz000 * test000);                            // nu/K*(u,v)
            MatrixA33[i][j] += Mult * val;
        }
        
        for(int j=0;j<N_P;j++)
        {
            ansatz000 = Orig4[j]; // p
            
            val = -Mult * ansatz000 * test100;                      // -(p,div v)=-p*v_x
            MatrixB1T[i][j] += val;
            
            val = -Mult * ansatz000 * test010;                      // -(p,div v)=-p*v_y
            MatrixB2T[i][j] += val;
            
            val = -Mult * ansatz000 * test001;                      // -(p,div v)=-p*v_z
            MatrixB3T[i][j] += val;
        }
    }// endfor i
    
    
    for(int i=0;i<N_P;i++)
    {
        test000 = Orig4[i]; // q
        
        for(int j=0;j<N_U;j++)
        {
            ansatz100 = Orig0[j]; // u_x
            ansatz010 = Orig1[j]; // u_y
            ansatz001 = Orig2[j]; // u_z
            
            val = -Mult * test000 * ansatz100;                       // -(q,div u)=-q*u_x
            MatrixB1[i][j] += val;
            
            val = -Mult * test000 * ansatz010;                       // -(q,div u)=-q*u_y
            MatrixB2[i][j] += val;
            
            val = -Mult * test000 * ansatz001;                       // -(q,div u)=-q*u_z
            MatrixB3[i][j] += val;
        }                            // endfor j
        //for MatrixC: loop over U_P here
    }                                // endfor i
}

// ======================================================================
// Type 2, Standard Galerkin for Brinkman in [p div v] formulation
// ======================================================================
// For NSTYPE: 4
void Brinkman3DType2Galerkin(double Mult, double *coeff,
                             double *, double,
                             double **OrigValues, int *N_BaseFuncts,
                             double ***LocMatrices, double **LocRhs)
{
    double ansatz000, ansatz100, ansatz010, ansatz001;    // ansatz functions
    double test000, test100, test010, test001;          // test functions
    
    double **MatrixA11 = LocMatrices[0];
    //double **MatrixA12 = LocMatrices[1];
    //double **MatrixA13 = LocMatrices[2];
    //double **MatrixA21 = LocMatrices[3];
    double **MatrixA22 = LocMatrices[4];
    //double **MatrixA23 = LocMatrices[5];
    //double **MatrixA31 = LocMatrices[6];
    //double **MatrixA32 = LocMatrices[7];
    double **MatrixA33 = LocMatrices[8];
    double **MatrixB1  = LocMatrices[9];
    double **MatrixB2  = LocMatrices[10];
    double **MatrixB3  = LocMatrices[11];
    double **MatrixB1T  = LocMatrices[12];
    double **MatrixB2T  = LocMatrices[13];
    double **MatrixB3T  = LocMatrices[14];
    
    double * Rhs_u1 = LocRhs[0];            // f_u1
    double * Rhs_u2 = LocRhs[1];            // f_u2
    double * Rhs_u3 = LocRhs[2];            // f_u3
    //double * Rhs_p = LocRhs[3];           // g_q
    
    int N_U = N_BaseFuncts[0];              // number of basis functions for the velocity space
    int N_P = N_BaseFuncts[1];              // number of basis functions for the pressure space
    
    // values of fe functions at quadrature points: Origvalues = f(uk,vk) (Gauß-quadrature: \int f(u,v) = sum wk * f(uk,vk) at Gauß points)
    // Mult is the quadrature weight (wk)
    double *Orig0 = OrigValues[0];          // u_x
    double *Orig1 = OrigValues[1];          // u_y
    double *Orig2 = OrigValues[2];          // u_z
    double *Orig3 = OrigValues[3];          // u
    double *Orig4 = OrigValues[4];          // p
    
    // values defined in the example
    //double c0 = coeff[0];                 // nu bzw eps (viscosity)
    double c1 = coeff[1];                   // f1
    double c2 = coeff[2];                   // f2
    double c3 = coeff[3];                   // f3 (the rhs of incompressibility constraint)
    //   double c4 = coeff[4];              // f4 (the rhs of incompressibility constraint)
    double nu = coeff[5];                   // viscosity
    double nu_eff = coeff[6];               // effective viscosity
    double K = coeff[7];                    // permeability
    
    double val;
    
    for(int i=0;i<N_U;i++)
    {
        // test functions
        test100 = Orig0[i]; //v_x
        test010 = Orig1[i]; //v_y
        test001 = Orig2[i]; //v_z
        test000 = Orig3[i]; //v
        
        Rhs_u1[i] += Mult*test000*c1;
        Rhs_u2[i] += Mult*test000*c2;
        Rhs_u3[i] += Mult*test000*c3;
        
        for(int j=0;j<N_U;j++)
        {
            // ansatz functions
            double ansatz100 = Orig0[j]; //u_x
            double ansatz010 = Orig1[j]; //u_y
            double ansatz001 = Orig2[j]; //u_z
            double ansatz000 = Orig3[j]; //u
            
            val  = nu_eff * (test100 * ansatz100 + test010 * ansatz010+ test001 * ansatz001);        // nu_eff*(grad u,grad v)=nu_eff*(v_x*u_x + v_y*u_y+ v_z*u_z)
            val += (nu/K) * (ansatz000 * test000);                            // nu/K*(u,v)
            MatrixA11[i][j] += Mult * val;
            
            val  = nu_eff * (test100 * ansatz100 + test010 * ansatz010+ test001 * ansatz001);        // nu_eff*(grad u,grad v)=nu_eff*(v_x*u_x + v_y*u_y+ v_z*u_z)
            val += (nu/K) * (ansatz000 * test000);                          // nu/K*(u,v)
            MatrixA22[i][j] += Mult * val;
            
            val  = nu_eff * (test100 * ansatz100 + test010 * ansatz010+ test001 * ansatz001);        // nu_eff*(grad u,grad v)=nu_eff*(v_x*u_x + v_y*u_y+ v_z*u_z)
            val += (nu/K) * (ansatz000 * test000);                            // nu/K*(u,v)
            MatrixA33[i][j] += Mult * val;
        }
        
        for(int j=0;j<N_P;j++)
        {
            ansatz000 = Orig4[j]; // p
            
            val = -Mult * ansatz000 * test100;                      // -(p,div v)=-p*v_x
            MatrixB1T[i][j] += val;
            
            val = -Mult * ansatz000 * test010;                      // -(p,div v)=-p*v_y
            MatrixB2T[i][j] += val;
            
            val = -Mult * ansatz000 * test001;                      // -(p,div v)=-p*v_z
            MatrixB3T[i][j] += val;
        }
    }
    // endfor i
    
    for(int i=0;i<N_P;i++)
    {
        test000 = Orig4[i]; // q
        
        for(int j=0;j<N_U;j++)
        {
            ansatz100 = Orig0[j]; // u_x
            ansatz010 = Orig1[j]; // u_y
            ansatz001 = Orig2[j]; // u_z
            
            val = -Mult * test000 * ansatz100;                       // -(q,div u)=-q*u_x
            MatrixB1[i][j] += val;
            
            val = -Mult * test000 * ansatz010;                       // -(q,div u)=-q*u_y
            MatrixB2[i][j] += val;
            
            val = -Mult * test000 * ansatz001;                       // -(q,div u)=-q*u_z
            MatrixB3[i][j] += val;
        }                            // endfor j
        //for MatrixC: loop over U_P here
    }                                // endfor i
}


// ======================================================================
// PSPG Stabilization (P2/P2 and P1/P1) for Brinkman Problem
// according to Hannukainen: Computations with Finite element methods for the Brinkman Problem (2011)
// ======================================================================
// For NSTYPE: 14
// Bem.: A_ij contains integrands build by u_j and v_i, i=1,2,3
//       B_i^T contains integrands build by p and v_i, i=1,2,3
//       B_i contains integrands build by u_i and q, i=1,2,3
//       C contains integrands build by p and q!

void ResidualStabPkPk_for_Brinkman3DType1Galerkin(double Mult, double *coeff,
                                                double *, double hK,
                                                double **OrigValues, int *N_BaseFuncts,
                                                double ***LocMatrices, double **LocRhs)
{
    double ansatz000, ansatz100, ansatz010, ansatz001, ansatz200, ansatz020, ansatz002;     // ansatz functions
    double test000, test100, test010, test001, test200, test020, test002;                   // test functions
    
    double **MatrixA11 = LocMatrices[0];
    // double **MatrixA12 = LocMatrices[1];
    // double **MatrixA13 = LocMatrices[2];
    // double **MatrixA21 = LocMatrices[3];
    double **MatrixA22 = LocMatrices[4];
    // double **MatrixA23 = LocMatrices[5];
    // double **MatrixA31 = LocMatrices[6];
    // double **MatrixA32 = LocMatrices[7];
    double **MatrixA33 = LocMatrices[8];
    double **MatrixC = LocMatrices[9];
    double **MatrixB1  = LocMatrices[10];
    double **MatrixB2  = LocMatrices[11];
    double **MatrixB3  = LocMatrices[12];
    double **MatrixB1T  = LocMatrices[13];
    double **MatrixB2T  = LocMatrices[14];
    double **MatrixB3T  = LocMatrices[15];
    
    double * Rhs_u1 = LocRhs[0];            // f_u1
    double * Rhs_u2 = LocRhs[1];            // f_u2
    double * Rhs_u3 = LocRhs[2];            // f_u3
    double * Rhs_p = LocRhs[3];             // g_q
    
    int N_U = N_BaseFuncts[0];              // number of basis functions for the velocity space
    int N_P = N_BaseFuncts[1];              // number of basis functions for the pressure space
    
    // values of fe functions at quadrature points: Origvalues = f(uk,vk) (Gauß-quadrature: \int f(u,v) = sum wk * f(uk,vk) at Gauß points)
    // Mult is the quadrature weight (wk)
    double *Orig0 = OrigValues[0];          // u_x
    double *Orig1 = OrigValues[1];          // u_y
    double *Orig2 = OrigValues[2];          // u_z
    double *Orig3 = OrigValues[3];          // u
    double *Orig4 = OrigValues[4];          // p
    double *Orig5 = OrigValues[5];          // p_x
    double *Orig6 = OrigValues[6];          // p_y
    double *Orig7 = OrigValues[7];          // p_z
    double *Orig8 = OrigValues[8];          // u_xx
    double *Orig9 = OrigValues[9];          // u_yy
    double *Orig10 = OrigValues[10];        // u_zz
    
    
    // values defined in the example
    // double c0 = coeff[0];                 // nu bzw eps (viscosity)
    double c1 = coeff[1];                   // f1
    double c2 = coeff[2];                   // f2
    double c3 = coeff[3];                   // f3
    // double c4 = coeff[4];                 // g (the rhs of incompressibility constraint)
    double nu = coeff[5];                   // viscosity
    double nu_eff = coeff[6];               // effective viscosity
    double K = coeff[7];                    // permeability
    double alpha = coeff[9];
    double t = fabs(sqrt((nu_eff/nu) * K));
    /////double PSPGStab = -alpha * (K/nu) * (hK*hK)/(t*t+hK*hK); // stabilization = (hK*hK)/(c0*c0+hK*hK) ///warum negativ, wenn positiv geht es schief-NOCHMAL TESTEN
    double PSPGStab = -alpha * (hK*hK)/(t*t+hK*hK); // stabilization = (hK*hK)/(c0*c0+hK*hK) ///warum negativ, wenn positiv geht es schief-NOCHMAL TESTEN
    double val;
    
    for(int i=0;i<N_U;i++)
    {
        test100 = Orig0[i];     // v_x
        test010 = Orig1[i];     // v_y
        test001 = Orig2[i];     // v_z
        test000 = Orig3[i];     // v
        test200 = Orig8[i];     // v_xx
        test020 = Orig9[i];     // v_yy
        test002 = Orig10[i];    // v_zz
        
        
        Rhs_u1[i] = Mult * PSPGStab  * (nu/K) * test000 * c1;                          // stabilization P1/P1: (nu/K)*(f1,v)
        Rhs_u2[i] = Mult * PSPGStab  * (nu/K) * test000 * c2;                          // stabilization P1/P1: (nu/K)*(f2,v)
        Rhs_u3[i] = Mult * PSPGStab  * (nu/K) * test000 * c3;                          // stabilization P1/P1: (nu/K)*(f3,v)
        Rhs_u1[i] -= Mult * PSPGStab  * nu_eff * (test200+test020+test002) * c1;               // stabilization P2/P2: nueff * (f1,Delta v)
        Rhs_u2[i] -= Mult * PSPGStab  * nu_eff * (test200+test020+test002) * c2;               // stabilization P2/P2: nueff *(f2,Delta v)
        Rhs_u3[i] -= Mult * PSPGStab  * nu_eff * (test200+test020+test002) * c3;               // stabilization P2/P2: nueff *(f3,Delta v)
        
        
        for(int j=0;j<N_U;j++)
        {
            ansatz100 = Orig0[j];   // u_x
            ansatz010 = Orig1[j];   // u_y
            ansatz001 = Orig2[j];   // u_z
            ansatz000 = Orig3[j];   // u
            ansatz200 = Orig8[j];   // u_xx
            ansatz020 = Orig9[j];   // u_yy
            ansatz002 = Orig10[j];  // u_zz
            
            
            val  = PSPGStab * (nu/K) * (nu/K) * (ansatz000 * test000);                              // stabilization P1/P1: (nu/K)^2 * (u,v)
            val += PSPGStab * nu_eff * nu_eff * (ansatz200+ansatz020+ansatz002)*(test200+test020+test002);   // stabilization P2/P2: nu_eff^2* (u_xx+u_yy+u_zz)*(v_xx+v_yy+v_zz)
            val -= PSPGStab * nu_eff * (nu/K) * (ansatz200+ansatz020+ansatz002)*test000;            // stabilization P2/P2: nu_eff* (nu/K) * (u_xx+u_yy+u_zz)*(v)
            val -= PSPGStab * nu_eff * (nu/K) * ansatz000 * (test200+test020+test002);                       // stabilization P2/P2: nu_eff* (nu/K) * (u)*(v_xx+v_yy+v_zz)
            MatrixA11[i][j] += Mult * val;
            
            
            
            val  = PSPGStab * (nu/K) * (nu/K) * (ansatz000 * test000);                              // stabilization P1/P1: (nu/K)^2 * (u,v)
            val += PSPGStab * nu_eff * nu_eff * (ansatz200+ansatz020+ansatz002)*(test200+test020+test002);    // stabilization P2/P2: nu_eff^2 *(u_xx+u_yy+u_zz)*(v_xx+v_yy+v_zz)
            val -= PSPGStab * nu_eff * (nu/K) * (ansatz200+ansatz020+ansatz002)*test000;            // stabilization P2/P2: nu_eff* (nu/K) * (u_xx+u_yy+u_zz)*(v)
            val -= PSPGStab * nu_eff * (nu/K) * ansatz000 * (test200+test020+test002);                       // stabilization P2/P2: nu_eff* (nu/K) * (u)*(v_xx+v_yy+v_zz)
            MatrixA22[i][j] += Mult * val;
            
            
            val  = PSPGStab * (nu/K) * (nu/K) * (ansatz000 * test000);                              // stabilization P1/P1: (nu/K)^2*(u,v)
            val += PSPGStab * nu_eff * nu_eff * (ansatz200+ansatz020+ansatz002)*(test200+test020+test002);    // stabilization P2/P2: nu_eff^2* (u_xx+u_yy+u_zz)*(v_xx+v_yy+v_zz)
            val -= PSPGStab * nu_eff * (nu/K) * (ansatz200+ansatz020+ansatz002)*test000;            // stabilization P2/P2: nu_eff* (nu/K) * (u_xx+u_yy+u_zz)*(v)
            val -= PSPGStab * nu_eff * (nu/K) * ansatz000 * (test200+test020+test002);                       // stabilization P2/P2: nu_eff* (nu/K) * (u)*(v_xx+v_yy+v_zz)
            MatrixA33[i][j] += Mult * val;
        }
        
        
        for(int j=0;j<N_P;j++)
        {
            ansatz000 = Orig4[j];
            ansatz100 = Orig5[j];
            ansatz010 = Orig6[j];
            ansatz001 = Orig7[j];
            
            val  = PSPGStab * (nu/K) * (ansatz100 * test000);                           // stabilization P1/P1: nu/K * (grad p,v) --> nu/K * p_x*v
            //val += PSPGStab * (nu/K) * (ansatz000 * test100);                         // stabilization P1/P1: nu/K * (p,div v) --> nu/K * p*v_x
            val -= PSPGStab * nu_eff * ansatz100 * (test200+test020+test002);           // stabilization P2/P2: nu_eff* p_x* Delta v
            MatrixB1T[i][j] += Mult * val;
            
            
            val  = PSPGStab * (nu/K) * (ansatz010 * test000);                            // stabilization P1/P1: nu/K * (grad p,v) --> nu/K * p_y*v
            //val+= PSPGStab * (nu/K) * (ansatz000 * test010);                          // stabilization P1/P1: nu/K * (p,div v) --> nu/K * p*v_y
            val -= PSPGStab * nu_eff * ansatz010 * (test200+test020+test002);           // stabilization P2/P2: nu_eff* p_x* Delta v
            MatrixB2T[i][j] += Mult * val;
            
            
            val  = PSPGStab * (nu/K) * (ansatz001 * test000);                                     // stabilization P1/P1: nu/K * (grad p,v)  --> nu/K * p_z*v
            //val+= PSPGStab * (nu/K) * (ansatz000 * test001);                                   // stabilization P1/P1: nu/K * (p,div v)  --> nu/K * p*v_z
            val -= PSPGStab * nu_eff * ansatz001 * (test200+test020+test002);           // stabilization P2/P2: nu_eff* p_x* Delta v
            MatrixB3T[i][j] += Mult * val;
        }
    }
    
    for(int i=0;i<N_P;i++)
    {
        test000 = Orig4[i];
        test100 = Orig5[i];
        test010 = Orig6[i];
        test001 = Orig7[i];
        
        Rhs_p[i]  = Mult * PSPGStab * test100*c1;                     // stabilization: (f,grad q) --> (f1,q_x)
        Rhs_p[i] += Mult * PSPGStab * test010*c2;                    // stabilization: (f,grad q) --> (f2,q_y)
        Rhs_p[i] += Mult * PSPGStab * test001*c3;                    // stabilization: (f,grad q) --> (f3,q_z)
        
        
        for(int j=0;j<N_U;j++)
        {
            ansatz100 = Orig0[j];
            ansatz010 = Orig1[j];
            ansatz001 = Orig2[j];
            ansatz000 = Orig3[j];
            ansatz200 = Orig8[j];
            ansatz020 = Orig9[j];
            ansatz002 = Orig10[j];
            
            
            val  = PSPGStab * (nu/K) * (test100 * ansatz000);                               // stabilization P1/P1: nu/K * (u,grad q) --> nu/K * q_x*u_x
            val -= PSPGStab * nu_eff * test100 * (ansatz200+ansatz020+ansatz002);           // stabilization P2/P2: nu_eff* q_x* Delta u
            MatrixB1[i][j] += Mult * val;
            
            
            val  = PSPGStab * (nu/K) * (test010 * ansatz000);                               // stabilization P1/P1: nu/K * (u,grad q) --> nu/K * q_y*u_y
            val -= PSPGStab * nu_eff * test010 * (ansatz200+ansatz020+ansatz002);           // stabilization P2/P2: nu_eff* q_x* Delta u
            MatrixB2[i][j] += Mult * val;
            
            val  = PSPGStab * (nu/K) * (test001 * ansatz000);                               // stabilization P1/P1: nu/K * (u,grad q) --> nu/K * q_z*u_z
            val -= PSPGStab * nu_eff * test001 * (ansatz200+ansatz020+ansatz002);           // stabilization P2/P2: nu_eff* q_x* Delta u
            MatrixB3[i][j] += Mult * val;
        }
        
        //for Matrix C: loop over N_P here
        for(int j=0;j<N_P;j++)
        {
            ansatz000 = Orig4[j];
            ansatz100 = Orig5[j];
            ansatz010 = Orig6[j];
            ansatz001 = Orig7[j];
            
            val = PSPGStab * (ansatz100 * test100 + ansatz010 * test010+ ansatz001 * test001);               // stabilization: (grad p,grad q)
            MatrixC[i][j] += Mult * val;
        }
    }
}


// ======================================================================
// Brinkman Problem with GradDiv Stabilization
// ======================================================================
// For NSTYPE: 14
void GradDivStab_for_Brinkman3DType1Galerkin(double, double *, double *, double,
                                             double **, int *, double ***,
                                             double **)
{
}

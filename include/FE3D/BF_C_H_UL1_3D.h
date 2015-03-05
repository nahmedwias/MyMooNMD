// ***********************************************************************
// Q1 + bubble element for LPS, conforming, 3D
// Author : Sashi
// History: 25.06.2011 
// ***********************************************************************

static void C_H_UL1_3D_Funct(double xi, double eta, double zeta,
                          double *values)
{

 double t1, t2, t3, t4, t6, t8, t11, t13, t15, t16, t19, t20,t23, t26;
 
  t1 = 1.0-xi;
  t2 = 1.0-eta;
  t3 = t1*t2;
  t4 = 1.0-zeta;
  t6 = xi*xi;
  t8 = eta*eta;
  t11 = zeta*zeta;
  t13 = (1.0-t6)*(1.0-t8)*(1.0-t11);
  t15 = 1.0+xi;
  t16 = t15*t2;
  t19 = 1.0+eta;
  t20 = t1*t19;
  t23 = t15*t19;
  t26 = 1.0+zeta;
      
  values[0] = t3*t4/8.0-t13/8.0;
  values[1] = t16*t4/8.0-t13/8.0;
  values[2] = t20*t4/8.0-t13/8.0;
  values[3] = t23*t4/8.0-t13/8.0;
  values[4] = t3*t26/8.0-t13/8.0;
  values[5] = t16*t26/8.0-t13/8.0;
  values[6] = t20*t26/8.0-t13/8.0;
  values[7] = t23*t26/8.0-t13/8.0;
  values[8] = t13;
 
}

static void C_H_UL1_3D_DeriveXi(double xi, double eta, double zeta,
                             double *values)
{
  double t1, t2, t4, t5, t8, t10, t11, t14,t16, t19, t21,  t25;
  
  t1 = 1.0-eta;
  t2 = 1.0-zeta;
  t4 = t1*t2/8.0;
  t5 = eta*eta;
  t8 = zeta*zeta;
  t10 = xi*(1.0-t5)*(1.0-t8);
  t11 = t10/4.0;
  t14 = 1.0+eta;
  t16 = t14*t2/8.0;
  t19 = 1.0+zeta;
  t21 = t1*t19/8.0;
  t25 = t14*t19/8.0;
      
  values[0] = -t4+t11;
  values[1] = t4+t11;
  values[2] = -t16+t11;
  values[3] = t16+t11;
  values[4] = -t21+t11;
  values[5] = t21+t11;
  values[6] = -t25+t11;
  values[7] = t25+t11;
  values[8] = -2.0*t10;   
}

static void C_H_UL1_3D_DeriveEta(double xi, double eta, double zeta,
                             double *values)
{
 double t1,t2,  t4, t5, t8,t10,  t11, t13, t15, t19, t21,t24; 
 
  t1 = 1.0-xi;
  t2 = 1.0-zeta;
  t4 = t1*t2/8.0;
  t5 = xi*xi;
  t8 = zeta*zeta;
  t10 = (1.0-t5)*eta*(1.0-t8);
  t11 = t10/4.0;
  t13 = 1.0+xi;
  t15 = t13*t2/8.0;
  t19 = 1.0+zeta;
  t21 = t1*t19/8.0;
  t24 = t13*t19/8.0;
      
  values[0] = -t4+t11;
  values[1] = -t15+t11;
  values[2] = t4+t11;
  values[3] = t15+t11;
  values[4] = -t21+t11;
  values[5] = -t24+t11;
  values[6] = t21+t11;
  values[7] = t24+t11;
  values[8] = -2.0*t10; 
 
}

static void C_H_UL1_3D_DeriveZeta(double xi, double eta, double zeta,
                             double *values)
{
  double t1, t2,  t4, t5, t7, t10, t11,t13, t15, t17, t19,  t22;
 
  t1 = 1.0-xi;
  t2 = 1.0-eta;
  t4 = t1*t2/8.0;
  t5 = xi*xi;
  t7 = eta*eta;
  t10 = (1.0-t5)*(1.0-t7)*zeta;
  t11 = t10/4.0;
  t13 = 1.0+xi;
  t15 = t13*t2/8.0;
  t17 = 1.0+eta;
  t19 = t1*t17/8.0;
  t22 = t13*t17/8.0;
      
  values[0] = -t4+t11;
  values[1] = -t15+t11;
  values[2] = -t19+t11;
  values[3] = -t22+t11;
  values[4] = t4+t11;
  values[5] = t15+t11;
  values[6] = t19+t11;
  values[7] = t22+t11;
  values[8] = -2.0*t10;   
}

static void C_H_UL1_3D_DeriveXiXi(double xi, double eta, double zeta,
                             double *values)
{
  double t1, t3, t5, t6;
  
  t1 = eta*eta;
  t3 = zeta*zeta;
  t5 = (1.0-t1)*(1.0-t3);
  t6 = t5/4.0;
      
  values[0] = t6;
  values[1] = t6;
  values[2] = t6;
  values[3] = t6;
  values[4] = t6;
  values[5] = t6;
  values[6] = t6;
  values[7] = t6;
  values[8] = -2.0*t5; 
}

static void C_H_UL1_3D_DeriveXiEta(double xi, double eta, double zeta,
                             double *values)
{
  double t1, t3, t5, t6, t7, t8, t9, t10;
  
  t1 = zeta/8.0;
  t3 = zeta*zeta;
  t5 = xi*eta*(1.0-t3);
  t6 = t5/2.0;
  t7 = 1.0/8.0-t1-t6;
  t8 = -1.0/8.0+t1-t6;
  t9 = 1.0/8.0+t1-t6;
  t10 = -1.0/8.0-t1-t6;
      
  values[0] = t7;
  values[1] = t8;
  values[2] = t8;
  values[3] = t7;
  values[4] = t9;
  values[5] = t10;
  values[6] = t10;
  values[7] = t9;
  values[8] = 4.0*t5; 
}

static void C_H_UL1_3D_DeriveXiZeta(double xi, double eta, double zeta,
                             double *values)
{
  double t1, t2, t5, t6, t7, t8, t9,t10 ;

  t1 = eta/8.0;
  t2 = eta*eta;
  t5 = xi*(1.0-t2)*zeta;
  t6 = t5/2.0;
  t7 = 1.0/8.0-t1-t6;
  t8 = -1.0/8.0+t1-t6;
  t9 = 1.0/8.0+t1-t6;
  t10 = -1.0/8.0-t1-t6;
      
  values[0] = t7;
  values[1] = t8;
  values[2] = t9;
  values[3] = t10;
  values[4] = t8;
  values[5] = t7;
  values[6] = t10;
  values[7] = t9;
  values[8] = 4.0*t5; 
}

static void C_H_UL1_3D_DeriveEtaEta(double xi, double eta, double zeta,
                             double *values)
{
  double t1, t3, t5, t6;
  
  t1 = xi*xi;
  t3 = zeta*zeta;
  t5 = (1.0-t1)*(1.0-t3);
  t6 = t5/4.0;
      
  values[0] = t6;
  values[1] = t6;
  values[2] = t6;
  values[3] = t6;
  values[4] = t6;
  values[5] = t6;
  values[6] = t6;
  values[7] = t6;
  values[8] = -2.0*t5;  
}

static void C_H_UL1_3D_DeriveEtaZeta(double xi, double eta, double zeta,
                             double *values)
{
  double t1, t2, t5, t6,t7, t8,t9, t10;

  t1 = xi/8.0;
  t2 = xi*xi;
  t5 = (1.0-t2)*eta*zeta;
  t6 = t5/2.0;
  t7 = 1.0/8.0-t1-t6;
  t8 = 1.0/8.0+t1-t6;
  t9 = -1.0/8.0+t1-t6;
  t10 = -1.0/8.0-t1-t6;
      
  values[0] = t7;
  values[1] = t8;
  values[2] = t9;
  values[3] = t10;
  values[4] = t9;
  values[5] = t10;
  values[6] = t7;
  values[7] = t8;
  values[8] = 4.0*t5; 
}

static void C_H_UL1_3D_DeriveZetaZeta(double xi, double eta, double zeta,
                             double *values)
{
  double t1, t3, t5, t6;
  
  t1 = xi*xi;
  t3 = eta*eta;
  t5 = (1.0-t1)*(1.0-t3);
  t6 = t5/4.0;
      
  values[0] = t6;
  values[1] = t6;
  values[2] = t6;
  values[3] = t6;
  values[4] = t6;
  values[5] = t6;
  values[6] = t6;
  values[7] = t6;
  values[8] = -2.0*t5;   
}

TBaseFunct3D *BF_C_H_UL1_3D_Obj = 
new TBaseFunct3D(9, BF_C_H_UL1_3D, BFUnitHexahedron, 
                 C_H_UL1_3D_Funct, C_H_UL1_3D_DeriveXi,
                 C_H_UL1_3D_DeriveEta, C_H_UL1_3D_DeriveZeta,
                 C_H_UL1_3D_DeriveXiXi, C_H_UL1_3D_DeriveXiEta,
                 C_H_UL1_3D_DeriveXiZeta, C_H_UL1_3D_DeriveEtaEta,
                 C_H_UL1_3D_DeriveEtaZeta, C_H_UL1_3D_DeriveZetaZeta,
                 2, 1,
                 0, NULL);

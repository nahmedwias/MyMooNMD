// =======================================================================
// @(#)HexaAffin.C      1.6 02/22/00
//
// Class:      THexaAffin
//
// Purpose:    affin reference transformations for hexahedron
//
// Author:     Daniel Quoos  
//
// History:    08.07.97 start implementation
// 
// =======================================================================

#include <MooNMD_Io.h>
#include <HexaAffin.h>
#include <FEDatabase3D.h>
#include <LinAlg.h>
#include <string.h>
#include <stdlib.h>

/** constuctor */
THexaAffin::THexaAffin()
{
}

/** transfer from reference element to original element */
void THexaAffin::GetOrigFromRef(double xi, double eta, double zeta, double &X, double &Y, double &Z)
{
  X = xc0 + xc1*xi + xc2*eta + xc3*zeta;
  Y = yc0 + yc1*xi + yc2*eta + yc3*zeta;
  Z = zc0 + zc1*xi + zc2*eta + zc3*zeta;
}

/** transfer a set of point from reference to original element */
void THexaAffin::GetOrigFromRef(int N_Points, double *xi, double *eta, double *zeta,
                                double *X, double *Y, double *Z, double *absdetjk)
{
  int i;
  double Xi, Eta, Zeta;
  double absdet=fabs(detjk);
  
  for(i=0;i<N_Points;i++)
  {
    Xi = xi[i];
    Eta = eta[i];
    Zeta = zeta[i];
    X[i] = xc0 + xc1*Xi + xc2*Eta + xc3*Zeta;
    Y[i] = yc0 + yc1*Xi + yc2*Eta + yc3*Zeta;
    Z[i] = zc0 + zc1*Xi + zc2*Eta + zc3*Zeta;
    absdetjk[i] = absdet;
  }
}

/** transfer from reference element to original element */
void THexaAffin::GetOrigFromRef(double *ref, double *orig)
{
  orig[0]=xc0 + xc1*ref[0] + xc2*ref[1] + xc3*ref[2];
  orig[1]=yc0 + yc1*ref[0] + yc2*ref[1] + yc3*ref[2];
  orig[2]=zc0 + zc1*ref[0] + zc2*ref[1] + zc3*ref[2];
}

/** transfer from original element to reference element */
void THexaAffin::GetRefFromOrig(double X, double Y, double Z, double &xi, double &eta, double &zeta)
{
  double xt=(X - xc0)/detjk;
  double yt=(Y - yc0)/detjk;
  double zt=(Z - zc0)/detjk;

  xi  = (yc2*zc3 - yc3*zc2)*xt - (xc2*zc3 - xc3*zc2)*yt + (xc2*yc3 - xc3*yc2)*zt;
  eta = -(yc1*zc3 - yc3*zc1)*xt + (xc1*zc3 - xc3*zc1)*yt - (xc1*yc3 - xc3*yc1)*zt;
  zeta = (yc1*zc2 - yc2*zc1)*xt - (xc1*zc2 - xc2*zc1)*yt + (xc1*yc2 - xc2*yc1)*zt;

  //  xi   = (yc2*zc3 - yc3*zc2)*xt + (yc3*zc1 - yc1*zc3)*yt + (yc1*zc2 - yc2*zc1)*zt;
  // eta  = (xc3*zc2 - xc2*zc3)*xt + (xc1*zc3 - xc3*zc1)*yt + (xc2*zc1 - xc1*zc2)*zt;
  // zeta = (xc2*yc3 - xc3*yc2)*xt + (xc3*yc1 - xc1*yc3)*yt + (xc1*yc2 - xc2*yc1)*zt;
}

/** transfer from original element to reference element */
void THexaAffin::GetRefFromOrig(double *orig, double *ref)
{
  double xt=(orig[0] - xc0)/detjk;
  double yt=(orig[1] - yc0)/detjk;
  double zt=(orig[2] - zc0)/detjk;

  ref[0]  = (yc2*zc3 - yc3*zc2)*xt - (xc2*zc3 - xc3*zc2)*yt + (xc2*yc3 - xc3*yc2)*zt;
  ref[1] = -(yc1*zc3 - yc3*zc1)*xt + (xc1*zc3 - xc3*zc1)*yt - (xc1*yc3 - xc3*yc1)*zt;
  ref[2] = (yc1*zc2 - yc2*zc1)*xt - (xc1*zc2 - xc2*zc1)*yt + (xc1*yc2 - xc2*yc1)*zt;
 
  //  ref[0] = (yc2*zc3 - yc3*zc2)*xt + (yc3*zc1 - yc1*zc3)*yt + (yc1*zc2 - yc2*zc1)*zt;
  // ref[1] = (xc3*zc2 - xc2*zc3)*xt + (xc1*zc3 - xc3*zc1)*yt + (xc2*zc1 - xc1*zc2)*zt;
  // ref[2] = (xc2*yc3 - xc3*yc2)*xt + (xc3*yc1 - xc1*yc3)*yt + (xc1*yc2 - xc2*yc1)*zt;
}

/** calculate functions and derivatives from reference element
    to original element */
void THexaAffin::GetOrigValues(BaseFunct3D BaseFunct,
                               int N_Points, double *xi, double *eta, double *zeta,
                               int N_Functs, QuadFormula3D QuadFormula)
{
  int i,j,k;
  double **refvaluesD000, **origvaluesD000;
  double **refvaluesD100, **origvaluesD100;
  double **refvaluesD010, **origvaluesD010;
  double **refvaluesD001, **origvaluesD001;
  double **refvaluesD200, **origvaluesD200;
  double **refvaluesD110, **origvaluesD110;
  double **refvaluesD020, **origvaluesD020;
  double **refvaluesD101, **origvaluesD101;
  double **refvaluesD011, **origvaluesD011;
  double **refvaluesD002, **origvaluesD002;
  double *refD000, *origD000;
  double *refD100, *origD100;
  double *refD010, *origD010;
  double *refD001, *origD001;
  double *refD200, *origD200;
  double *refD110, *origD110;
  double *refD020, *origD020;
  double *refD101, *origD101;
  double *refD011, *origD011;
  double *refD002, *origD002;
  double *aux;
  double AllData[MaxN_BaseFunctions3D][5];
  double GeoData[5][5];
  
  refvaluesD000=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula, D000);
  if(refvaluesD000==NULL)
    {
      TFEDatabase3D::GetBaseFunct3D(BaseFunct)->MakeRefElementData(QuadFormula);
      refvaluesD000=TFEDatabase3D::GetRefElementValues(BaseFunct, 
                                                       QuadFormula, D000);
    }
  
  origvaluesD000=TFEDatabase3D::GetOrigElementValues(BaseFunct, D000);
  if(origvaluesD000==NULL)
    {
      origvaluesD000 = new double* [MaxN_QuadPoints_3D];
      aux = new double [MaxN_QuadPoints_3D*MaxN_BaseFunctions3D];
      for(i=0;i<MaxN_QuadPoints_3D;i++)
        origvaluesD000[i] = aux+i*MaxN_BaseFunctions3D;
      TFEDatabase3D::RegisterOrigElementValues(BaseFunct, D000, origvaluesD000);
    }
  
  refvaluesD100=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula, D100);
  origvaluesD100=TFEDatabase3D::GetOrigElementValues(BaseFunct, D100);
  if(origvaluesD100==NULL)
    {
      origvaluesD100 = new double* [MaxN_QuadPoints_3D];
      aux = new double [MaxN_QuadPoints_3D*MaxN_BaseFunctions3D];
      for(i=0;i<MaxN_QuadPoints_3D;i++)
        origvaluesD100[i] = aux+i*MaxN_BaseFunctions3D;
      TFEDatabase3D::RegisterOrigElementValues(BaseFunct, D100, origvaluesD100);
    }
  
  refvaluesD010=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula, D010);
  origvaluesD010=TFEDatabase3D::GetOrigElementValues(BaseFunct, D010);
  if(origvaluesD010==NULL)
    {
      origvaluesD010 = new double* [MaxN_QuadPoints_3D];
      aux = new double [MaxN_QuadPoints_3D*MaxN_BaseFunctions3D];
      for(i=0;i<MaxN_QuadPoints_3D;i++)
        origvaluesD010[i] = aux+i*MaxN_BaseFunctions3D;
      TFEDatabase3D::RegisterOrigElementValues(BaseFunct, D010, origvaluesD010);
    } 
  
  refvaluesD001=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula, D001);
  origvaluesD001=TFEDatabase3D::GetOrigElementValues(BaseFunct, D001);
  if(origvaluesD001==NULL)
    { 
      origvaluesD001 = new double* [MaxN_QuadPoints_3D];
      aux = new double [MaxN_QuadPoints_3D*MaxN_BaseFunctions3D];
      for(i=0;i<MaxN_QuadPoints_3D;i++)
        origvaluesD001[i] = aux+i*MaxN_BaseFunctions3D;
      TFEDatabase3D::RegisterOrigElementValues(BaseFunct, D001, origvaluesD001);
    } 
  
  refvaluesD200=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula, D200);
  origvaluesD200=TFEDatabase3D::GetOrigElementValues(BaseFunct, D200);
  if(origvaluesD200==NULL)
    { 
      origvaluesD200 = new double* [MaxN_QuadPoints_3D];
      aux = new double [MaxN_QuadPoints_3D*MaxN_BaseFunctions3D];
      for(i=0;i<MaxN_QuadPoints_3D;i++)
        origvaluesD200[i] = aux+i*MaxN_BaseFunctions3D;
      TFEDatabase3D::RegisterOrigElementValues(BaseFunct, D200, origvaluesD200);
    } 
  
  refvaluesD110=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula, D110);
  origvaluesD110=TFEDatabase3D::GetOrigElementValues(BaseFunct, D110);
  if(origvaluesD110==NULL)
    { 
      origvaluesD110 = new double* [MaxN_QuadPoints_3D];
      aux = new double [MaxN_QuadPoints_3D*MaxN_BaseFunctions3D];
      for(i=0;i<MaxN_QuadPoints_3D;i++)
        origvaluesD110[i] = aux+i*MaxN_BaseFunctions3D;
      TFEDatabase3D::RegisterOrigElementValues(BaseFunct, D110, origvaluesD110);
    } 
  
  refvaluesD101=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula, D101);
  origvaluesD101=TFEDatabase3D::GetOrigElementValues(BaseFunct, D101);
  if(origvaluesD101==NULL)
    { 
      origvaluesD101 = new double* [MaxN_QuadPoints_3D];
      aux = new double [MaxN_QuadPoints_3D*MaxN_BaseFunctions3D];
      for(i=0;i<MaxN_QuadPoints_3D;i++)
        origvaluesD101[i] = aux+i*MaxN_BaseFunctions3D;
      TFEDatabase3D::RegisterOrigElementValues(BaseFunct, D101, origvaluesD101);
    } 
  
  refvaluesD011=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula, D011);
  origvaluesD011=TFEDatabase3D::GetOrigElementValues(BaseFunct, D011);
  if(origvaluesD011==NULL)
    { 
      origvaluesD011 = new double* [MaxN_QuadPoints_3D];
      aux = new double [MaxN_QuadPoints_3D*MaxN_BaseFunctions3D];
      for(i=0;i<MaxN_QuadPoints_3D;i++)
        origvaluesD011[i] = aux+i*MaxN_BaseFunctions3D;
      TFEDatabase3D::RegisterOrigElementValues(BaseFunct, D011, origvaluesD011);
    } 
  
  refvaluesD020=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula, D020);
  origvaluesD020=TFEDatabase3D::GetOrigElementValues(BaseFunct, D020);
  if(origvaluesD020==NULL)
    {
      origvaluesD020 = new double* [MaxN_QuadPoints_3D];
      aux = new double [MaxN_QuadPoints_3D*MaxN_BaseFunctions3D];
      for(i=0;i<MaxN_QuadPoints_3D;i++)
        origvaluesD020[i] = aux+i*MaxN_BaseFunctions3D;
      TFEDatabase3D::RegisterOrigElementValues(BaseFunct, D020, origvaluesD020);
    } 
  
  refvaluesD002=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula, D002);
  origvaluesD002=TFEDatabase3D::GetOrigElementValues(BaseFunct, D002);
  if(origvaluesD002==NULL)
    {
      origvaluesD002 = new double* [MaxN_QuadPoints_3D];
      aux = new double [MaxN_QuadPoints_3D*MaxN_BaseFunctions3D];
      for(i=0;i<MaxN_QuadPoints_3D;i++)
        origvaluesD002[i] = aux+i*MaxN_BaseFunctions3D;
      TFEDatabase3D::RegisterOrigElementValues(BaseFunct, D002, origvaluesD002);
    } 
  // D000
  for(i=0;i<N_Points;i++)
    {
      refD000 = refvaluesD000[i];
      origD000 = origvaluesD000[i];
      
      for(j=0;j<N_Functs;j++)
        {
          origD000[j] = refD000[j];
        } // endfor j
    } // endfor i
  
  // /*
  // D100, D010 and D001
  for(i=0;i<N_Points;i++)
    { 
      refD100 = refvaluesD100[i];
      origD100 = origvaluesD100[i];
      
      refD010 = refvaluesD010[i];
      origD010 = origvaluesD010[i];
      
      refD001 = refvaluesD001[i];
      origD001 = origvaluesD001[i];
      
      for(j=0;j<N_Functs;j++)
        {
          origD100[j] = ((yc2*zc3 - yc3*zc2)*refD100[j] + (yc3*zc1 - yc1*zc3)*refD010[j] + (yc1*zc2 - yc2*zc1)*refD001[j]) *rec_detjk;
          origD010[j] = ((xc3*zc2 - xc2*zc3)*refD100[j] + (xc1*zc3 - xc3*zc1)*refD010[j] + (xc2*zc1 - xc1*zc2)*refD001[j]) *rec_detjk;
          origD001[j] = ((xc2*yc3 - yc2*xc3)*refD100[j] + (xc3*yc1 - xc1*yc3)*refD010[j] + (xc1*yc2 - xc2*yc1)*refD001[j]) *rec_detjk;
          
          
          // cout << "10: " << origD10[j] << endl;
          // cout << "01: " << origD01[j] << endl;
          // cout << endl;
        } // endfor j
      // cout << "----------" << endl;
    } // endfor i
  
  /*
    for(i=0;i<N_Points;i++)
    {
    // reset matrix
    memset(GeoData, 0, 25*SizeOfDouble);
    
    refD10 = refvaluesD10[i];
    refD01 = refvaluesD01[i];
    refD20 = refvaluesD20[i];
    refD11 = refvaluesD11[i];
    refD02 = refvaluesD02[i];
    
    origD10 = origvaluesD10[i];
    origD01 = origvaluesD01[i];
    origD20 = origvaluesD20[i];
    origD11 = origvaluesD11[i];
    origD02 = origvaluesD02[i];
    
    GeoData[0][0] = xc1;
    GeoData[0][1] = yc1;

    GeoData[1][0] = xc2;
    GeoData[1][1] = yc2;

    GeoData[2][2] = xc1*xc1;
    GeoData[2][3] = 2*xc1*yc1;
    GeoData[2][4] = yc1*yc1;
    
    GeoData[3][2] = xc1*xc2;
    GeoData[3][3] = yc1*xc2+xc1*yc2; 
    GeoData[3][4] = yc1*yc2;
    
    GeoData[4][2] = xc2*xc2;
    GeoData[4][3] = 2*xc2*yc2;
    GeoData[4][4] = yc2*yc2;

    for(j=0;j<N_Functs;j++)
    {
    AllData[j][0] = refD10[j];
    AllData[j][1] = refD01[j];
    AllData[j][2] = refD20[j];
    AllData[j][3] = refD11[j];
    AllData[j][4] = refD02[j];
    } // endfor j
    
    // subroutine for solving a multiple systems of linear equations
    // void SolveMultipleSystems(double *a, double *b, int N_Eqn,
    //                        int LDA, int LDB, int N_Rhs);
    SolveMultipleSystems((double *)GeoData, (double *)AllData, 5,
    5, 5, N_Functs);
    
    for(j=0;j<N_Functs;j++)
    {
    origD10[j] = AllData[j][0];
    origD01[j] = AllData[j][1];
    origD20[j] = AllData[j][2];
    origD11[j] = AllData[j][3];
    origD02[j] = AllData[j][4];
    
    } // endfor j
    
    } // endfor i
     */
      }

/** calculate functions and derivatives from reference element
    to original element, for all given elements */
void THexaAffin::GetOrigValues(int N_Sets, BaseFunct3D *BaseFuncts,
                               int N_Points, double *xi, double *eta, double *zeta, 
                               QuadFormula3D QuadFormula,
                               bool *Needs2ndDer)
{
  int i,j,k,N_, start, end;
  double **refvaluesD000, **origvaluesD000;
  double **refvaluesD100, **origvaluesD100;
  double **refvaluesD010, **origvaluesD010;
  double **refvaluesD001, **origvaluesD001;
  double **refvaluesD200, **origvaluesD200;
  double **refvaluesD110, **origvaluesD110;
  double **refvaluesD101, **origvaluesD101;
  double **refvaluesD011, **origvaluesD011;
  double **refvaluesD020, **origvaluesD020;
  double **refvaluesD002, **origvaluesD002;
  double *refD000, *origD000;
  double *refD100, *origD100;
  double *refD010, *origD010;
  double *refD001, *origD001;
  double *refD200, *origD200;
  double *refD110, *origD110;
  double *refD101, *origD101;
  double *refD011, *origD011;
  double *refD020, *origD020;
  double *refD002, *origD002;
  double r20, r11, r02, o20, o11, o02;
  double *aux;
  double GeoData[3][3];
  double Eye[3][3];
  BaseFunct3D BaseFunct;
  int N_Functs;
  bool SecondDer;
  int ii,ij,ik;
  double tmp,Eye1[3][3];
  
  SecondDer = FALSE;
  for(i=0;i<N_Sets;i++)
  {
    BaseFunct=BaseFuncts[i];
    N_Functs = TFEDatabase3D::GetBaseFunct3D(BaseFunct)->GetDimension();
      
    refvaluesD000=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula, D000);
    if(refvaluesD000==NULL)
    {
      TFEDatabase3D::GetBaseFunct3D(BaseFunct)->MakeRefElementData(QuadFormula);
      refvaluesD000=TFEDatabase3D::GetRefElementValues(BaseFunct, 
                                                       QuadFormula, D000);
    }
      
    origvaluesD000=TFEDatabase3D::GetOrigElementValues(BaseFunct, D000);
    if(origvaluesD000==NULL)
    {
      origvaluesD000 = new double* [MaxN_QuadPoints_3D];
      aux = new double [MaxN_QuadPoints_3D*MaxN_BaseFunctions3D];
      for(j=0;j<MaxN_QuadPoints_3D;j++)
        origvaluesD000[j] = aux+j*MaxN_BaseFunctions3D;
      TFEDatabase3D::RegisterOrigElementValues(BaseFunct, D000, origvaluesD000);
    }
      
    for(j=0;j<N_Points;j++)
    {
      refD000 = refvaluesD000[j];
      origD000 = origvaluesD000[j];
          
      memcpy(origD000, refD000, N_Functs*SizeOfDouble);
    } // endfor j
      
    refvaluesD100=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula, D100);
    origvaluesD100=TFEDatabase3D::GetOrigElementValues(BaseFunct, D100);
    if(origvaluesD100==NULL)
    {
      origvaluesD100 = new double* [MaxN_QuadPoints_3D];
      aux = new double [MaxN_QuadPoints_3D*MaxN_BaseFunctions3D];
      for(j=0;j<MaxN_QuadPoints_3D;j++)
        origvaluesD100[j] = aux+j*MaxN_BaseFunctions3D;
      TFEDatabase3D::RegisterOrigElementValues(BaseFunct, D100, origvaluesD100);
    }
      
    refvaluesD010=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula, D010);
    origvaluesD010=TFEDatabase3D::GetOrigElementValues(BaseFunct, D010);
    if(origvaluesD010==NULL)
    {
       origvaluesD010 = new double* [MaxN_QuadPoints_3D];
       aux = new double [MaxN_QuadPoints_3D*MaxN_BaseFunctions3D];
       for(j=0;j<MaxN_QuadPoints_3D;j++)
         origvaluesD010[j] = aux+j*MaxN_BaseFunctions3D;
       TFEDatabase3D::RegisterOrigElementValues(BaseFunct, D010, origvaluesD010);
    }
      
    refvaluesD001=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula, D001);
    origvaluesD001=TFEDatabase3D::GetOrigElementValues(BaseFunct, D001);
    if(origvaluesD001==NULL)
    {
      origvaluesD001 = new double* [MaxN_QuadPoints_3D];
      aux = new double [MaxN_QuadPoints_3D*MaxN_BaseFunctions3D];
      for(j=0;j<MaxN_QuadPoints_3D;j++)
        origvaluesD001[j] = aux+j*MaxN_BaseFunctions3D;
      TFEDatabase3D::RegisterOrigElementValues(BaseFunct, D001, origvaluesD001);
    }
  }
      /*
      if(Needs2ndDer[i])
        {
          SecondDer = TRUE;
          
          refvaluesD200=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula, D200);
          origvaluesD200=TFEDatabase3D::GetOrigElementValues(BaseFunct, D200);
          if(origvaluesD200==NULL)
            {
              origvaluesD200 = new double* [MaxN_QuadPoints_3D];
              aux = new double [MaxN_QuadPoints_3D*MaxN_BaseFunctions3D];
              for(j=0;j<MaxN_QuadPoints_3D;j++)
                origvaluesD200[j] = aux+j*MaxN_BaseFunctions3D;
              TFEDatabase3D::RegisterOrigElementValues(BaseFunct, D200, origvaluesD200);
            }
          
          refvaluesD110=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula, D110);
          origvaluesD110=TFEDatabase3D::GetOrigElementValues(BaseFunct, D110);
          if(origvaluesD110==NULL)
            {
              origvaluesD110 = new double* [MaxN_QuadPoints_3D];
              aux = new double [MaxN_QuadPoints_3D*MaxN_BaseFunctions3D];
              for(j=0;j<MaxN_QuadPoints_3D;j++)
                origvaluesD110[j] = aux+j*MaxN_BaseFunctions3D;
              TFEDatabase3D::RegisterOrigElementValues(BaseFunct, D110, origvaluesD110);
            }
          
          refvaluesD101=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula, D101);
          origvaluesD101=TFEDatabase3D::GetOrigElementValues(BaseFunct, D101);
          if(origvaluesD101==NULL)
            {
              origvaluesD101 = new double* [MaxN_QuadPoints_3D];
              aux = new double [MaxN_QuadPoints_3D*MaxN_BaseFunctions3D];
              for(j=0;j<MaxN_QuadPoints_3D;j++)
                origvaluesD101[j] = aux+j*MaxN_BaseFunctions3D;
              TFEDatabase3D::RegisterOrigElementValues(BaseFunct, D101, origvaluesD101);
            }
          
          refvaluesD011=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula, D011);
          origvaluesD011=TFEDatabase3D::GetOrigElementValues(BaseFunct, D011);
          if(origvaluesD011==NULL)
            {
              origvaluesD011 = new double* [MaxN_QuadPoints_3D];
              aux = new double [MaxN_QuadPoints_3D*MaxN_BaseFunctions3D];
              for(j=0;j<MaxN_QuadPoints_3D;j++)
                origvaluesD011[j] = aux+j*MaxN_BaseFunctions3D;
              TFEDatabase3D::RegisterOrigElementValues(BaseFunct, D011, origvaluesD011);
            }
          
          refvaluesD020=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula, D020);
          origvaluesD020=TFEDatabase3D::GetOrigElementValues(BaseFunct, D020);
          if(origvaluesD020==NULL)
            {
              origvaluesD020 = new double* [MaxN_QuadPoints_3D];
              aux = new double [MaxN_QuadPoints_3D*MaxN_BaseFunctions3D];
              for(j=0;j<MaxN_QuadPoints_3D;j++)
                origvaluesD020[j] = aux+j*MaxN_BaseFunctions3D;
              TFEDatabase3D::RegisterOrigElementValues(BaseFunct, D020, origvaluesD020);
            }
          
          refvaluesD002=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula, D002);
          origvaluesD002=TFEDatabase3D::GetOrigElementValues(BaseFunct, D002);
          if(origvaluesD002==NULL)
            {
              origvaluesD002 = new double* [MaxN_QuadPoints_3D];
              aux = new double [MaxN_QuadPoints_3D*MaxN_BaseFunctions3D];
              for(j=0;j<MaxN_QuadPoints_3D;j++)
                origvaluesD002[j] = aux+j*MaxN_BaseFunctions3D;
              TFEDatabase3D::RegisterOrigElementValues(BaseFunct, D002, origvaluesD002);
            } // endfor Needs2ndDer[i]
        } // endfor i
      */
    // D100, D010 and D001
    for(i=0;i<N_Sets;i++)
    {
      BaseFunct=BaseFuncts[i];
      N_Functs = TFEDatabase3D::GetBaseFunct3D(BaseFunct)->GetDimension();
          
      refvaluesD100=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula, D100);
      origvaluesD100=TFEDatabase3D::GetOrigElementValues(BaseFunct, D100);
      
      refvaluesD010=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula, D010);
      origvaluesD010=TFEDatabase3D::GetOrigElementValues(BaseFunct, D010);
          
      refvaluesD001=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula, D001);
      origvaluesD001=TFEDatabase3D::GetOrigElementValues(BaseFunct, D001);
          
      for(j=0;j<N_Points;j++)
      {
        refD100 = refvaluesD100[j];
        origD100 = origvaluesD100[j];
              
        refD010 = refvaluesD010[j];
        origD010 = origvaluesD010[j];
              
        refD001 = refvaluesD001[j];
        origD001 = origvaluesD001[j];
              
        for(k=0;k<N_Functs;k++)
        {
          origD100[k] = ((yc2*zc3 - yc3*zc2)*refD100[k] + (yc3*zc1 - yc1*zc3)*refD010[k] + (yc1*zc2 - yc2*zc1)*refD001[k]) *rec_detjk;
          origD010[k] = ((xc3*zc2 - xc2*zc3)*refD100[k] + (xc1*zc3 - xc3*zc1)*refD010[k] + (xc2*zc1 - xc1*zc2)*refD001[k]) *rec_detjk;
          origD001[k] = ((xc2*yc3 - xc3*yc2)*refD100[k] + (xc3*yc1 - xc1*yc3)*refD010[k] + (xc1*yc2 - xc2*yc1)*refD001[k]) *rec_detjk;                   
        } // endfor k
      } // endfor j
    } // endfor i
   
    // leave if no second derivatives are needed
    if(!SecondDer) return;
      /*
        // find transformation matrix for second derivatives
        // do this only once since matrix is constant for an affine mapping
        
        // reset matrices
        memset(GeoData, 0, 9*SizeOfDouble);
        memset(Eye, 0, 9*SizeOfDouble);
        Eye[0][0] = 1;
        Eye[1][1] = 1;
        Eye[2][2] = 1;
        
        GeoData[0][0] = xc1*xc1;
        GeoData[0][1] = 2*xc1*yc1;
        GeoData[0][2] = yc1*yc1;
        
        GeoData[1][0] = xc1*xc2;
        GeoData[1][1] = yc1*xc2+xc1*yc2; 
        GeoData[1][2] = yc1*yc2;
        
        GeoData[2][0] = xc2*xc2;
        GeoData[2][1] = 2*xc2*yc2;
        GeoData[2][2] = yc2*yc2;
        
        // subroutine for solving a multiple systems of linear equations
        // void SolveMultipleSystems(double *a, double *b, int N_Eqn,
        //                        int LDA, int LDB, int N_Rhs);
        //SolveMultipleSystems((double *)GeoData, (double *)Eye, 3,
        //                 3, 3, 3);
        //SolveLinearSystem((double *)GeoData, (double *)Eye, 3, 3);
                    
        SolveMultipleSystemsNew((double *)GeoData, (double *)Eye, 3,
        3, 3, 3);
        
        for(i=0;i<N_Sets;i++)
        {
  if(Needs2ndDer[i])
  {
  BaseFunct=BaseFuncts[i];
  N_Functs = TFEDatabase3D::GetBaseFunct3D(BaseFunct)->GetDimension();
  
  refvaluesD20=TFEDatabase3D::GetRefElementValues
  (BaseFunct, QuadFormula, D20);
  origvaluesD20=TFEDatabase3D::GetOrigElementValues
  (BaseFunct, D20);
  
  refvaluesD11=TFEDatabase3D::GetRefElementValues
  (BaseFunct, QuadFormula, D11);
  origvaluesD11=TFEDatabase3D::GetOrigElementValues
  (BaseFunct, D11);
  
  refvaluesD02=TFEDatabase3D::GetRefElementValues
  (BaseFunct, QuadFormula, D02);
  origvaluesD02=TFEDatabase3D::GetOrigElementValues
  (BaseFunct, D02);
  
  for(j=0;j<N_Points;j++)
  {
  refD20 = refvaluesD20[j];
  refD11 = refvaluesD11[j];
  refD02 = refvaluesD02[j];
  
  origD20 = origvaluesD20[j];
  origD11 = origvaluesD11[j];
  origD02 = origvaluesD02[j];
  
  for(k=0;k<N_Functs;k++)
  {
    r20 = refD20[k];
    r11 = refD11[k];
    r02 = refD02[k];
  
    o20 = Eye[0][0]*r20+Eye[0][1]*r11+Eye[0][2]*r02;
    o11 = Eye[1][0]*r20+Eye[1][1]*r11+Eye[1][2]*r02;
    o02 = Eye[2][0]*r20+Eye[2][1]*r11+Eye[2][2]*r02;
  
    origD20[k] = o20;
    origD11[k] = o11;
    origD02[k] = o02;
  } // endfor k
  } // endif
  } // endfor j
  } // endfor i
*/
}
/** calculate functions and derivatives from reference element
    to original element */
void THexaAffin::GetOrigValues(double xi, double eta, double zeta,
                               int N_BaseFunct,
                               double *uref, double *uxiref, double *uetaref, double *uzetaref,
                               double *uorig, double *uxorig, double *uyorig, double *uzorig)
{
  int i;
  
  // D000
  for(i=0;i<N_BaseFunct;i++)
    uorig[i] = uref[i];

  // D100, D010 and D001
  for(i=0;i<N_BaseFunct;i++)
  {
    uxorig[i] = ((yc2*zc3 - yc3*zc2)*uxiref[i] + (yc3*zc1 - yc1*zc3)*uetaref[i] + (yc1*zc2 - yc2*zc1)*uzetaref[i]) *rec_detjk;
    uyorig[i] = ((xc3*zc2 - xc2*zc3)*uxiref[i] + (xc1*zc3 - xc3*zc1)*uetaref[i] + (xc2*zc1 - xc1*zc2)*uzetaref[i]) *rec_detjk;
    uzorig[i] = ((xc2*yc3 - xc3*yc2)*uxiref[i] + (xc3*yc1 - xc1*yc3)*uetaref[i] + (xc1*yc2 - xc2*yc1)*uzetaref[i]) *rec_detjk;
  } // endfor i
}

void THexaAffin::SetCell(TBaseCell *cell)
{
  int i;

  Cell = cell;

  Cell->GetVertex(0)->GetCoords(x0, y0, z0);
  Cell->GetVertex(1)->GetCoords(x1, y1, z1);
  Cell->GetVertex(2)->GetCoords(x2, y2, z2);
  Cell->GetVertex(3)->GetCoords(x3, y3, z3);
  Cell->GetVertex(4)->GetCoords(x4, y4, z4);
  Cell->GetVertex(5)->GetCoords(x5, y5, z5);
  Cell->GetVertex(6)->GetCoords(x6, y6, z6);
  Cell->GetVertex(7)->GetCoords(x7, y7, z7);

  xc0 = (x1 + x3 + x4 - x0) * 0.5;
  xc1 = (x1 - x0) * 0.5;
  xc2 = (x3 - x0) * 0.5;
  xc3 = (x4 - x0) * 0.5;

  yc0 = (y1 + y3 + y4 - y0) * 0.5;
  yc1 = (y1 - y0) * 0.5;
  yc2 = (y3 - y0) * 0.5;
  yc3 = (y4 - y0) * 0.5;

  zc0 = (z1 + z3 + z4 - z0) * 0.5;
  zc1 = (z1 - z0) * 0.5;
  zc2 = (z3 - z0) * 0.5;
  zc3 = (z4 - z0) * 0.5;

  detjk=xc1*yc2*zc3 + xc2*yc3*zc1 + xc3*yc1*zc2 - xc3*yc2*zc1 - xc2*yc1*zc3 - xc1*yc3*zc2;
  rec_detjk=1/detjk;
}

/** return outer normal unit vector */
void THexaAffin::GetOuterNormal(int j, double s, double t,
                                double &n1, double &n2, double &n3)
{
//   double len;
 
  switch(j)
  {
    case 0:
      n1 = yc2*zc1 - zc2*yc1;
      n2 = zc2*xc1 - xc2*zc1;
      n3 = xc2*yc1 - yc2*xc1;
    break;
 
    case 1:
      n1 = yc1*zc3 - zc1*yc3;
      n2 = zc1*xc3 - xc1*zc3;
      n3 = xc1*yc3 - yc1*xc3;
    break;

    case 2:
      n1 = yc2*zc3 - zc2*yc3;  //OLD  yc2*zc3 - zc2*xc3;   sashi 20.1.12
      n2 = zc2*xc3 - xc2+zc3;
      n3 = xc2*yc3 - yc2*xc3;
    break;
 
    case 3:
      n1 = zc1*yc3 - yc1*zc3;
      n2 = xc1*zc3 - zc1*xc3;
      n3 = yc1*xc3 - xc1*yc3;
    break;

    case 4:
      n1 = zc2*yc3 - yc2*zc3;
      n2 = xc2*zc3 - zc2*xc3;
      n3 = yc2*xc3 - xc2*yc3;
    break;

    case 5:
      n1 = zc2*yc1 - yc2*zc1;
      n2 = xc2*zc1 - zc2*xc1;
      n3 = yc2*xc1 - xc2*yc1;
    break;

    default:
      Error("Wrong local joint number" << endl);
  }

  
  // without scaling needed for finding area, Sashi, 20. Jan 2012
//   len = sqrt(n1*n1 + n2*n2 + n3*n3);
// 
//   n1 /= len;
//   n2 /= len;
//   n3 /= len;
}

/** return two tangent vectors */
void THexaAffin::GetTangentVectors(int j, double p1, double p2,
        double &t11, double &t12, double &t13,
        double &t21, double &t22, double &t23)
{
  switch(j)
  {
    case 0:
      t11 = xc2; t12 = yc2; t13 = zc2;
      t21 = xc1; t22 = yc1; t23 = zc1;
    break;

    case 1:
      t11 = xc1; t12 = yc1; t13 = zc1;
      t21 = xc3; t22 = yc3; t23 = zc3;
    break;

    case 2:
      t11 = xc2; t12 = yc2; t13 = zc2;
      t21 = xc3; t22 = yc3; t23 = zc3;
    break;

    case 3:
      t11 = xc3; t12 = yc3; t13 = zc3;
      t21 = xc1; t22 = yc1; t23 = zc1;
    break;

    case 4:
      t11 = xc3; t12 = yc3; t13 = zc3;
      t21 = xc2; t22 = yc2; t23 = zc2;
    break;

    case 5:
      t11 = xc1; t12 = yc1; t13 = zc1;
      t21 = xc2; t22 = yc2; t23 = zc2;
    break;

    default:
      Error("Wrong joint number for hexahedron! " << endl);
      return;
  }
} // end THexaAffin::GetTangentVectors


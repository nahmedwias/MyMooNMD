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
#include "BaseCell.h"
#include <string.h>
#include <stdlib.h>

/** constuctor */
THexaAffin::THexaAffin()
{
}

/** transfer from reference element to original element */
void THexaAffin::GetOrigFromRef(double xi, double eta, double zeta, double &X,
                                double &Y, double &Z)
{
  X = xc0 + xc1*xi + xc2*eta + xc3*zeta;
  Y = yc0 + yc1*xi + yc2*eta + yc3*zeta;
  Z = zc0 + zc1*xi + zc2*eta + zc3*zeta;
}

/** transfer a set of point from reference to original element */
void THexaAffin::GetOrigFromRef(int N_Points, const double *xi,
                                const double *eta, const double *zeta,
                                double *X, double *Y, double *Z,
                                double *absdetjk)
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
void THexaAffin::GetOrigFromRef(const double *ref, double *orig)
{
  orig[0]=xc0 + xc1*ref[0] + xc2*ref[1] + xc3*ref[2];
  orig[1]=yc0 + yc1*ref[0] + yc2*ref[1] + yc3*ref[2];
  orig[2]=zc0 + zc1*ref[0] + zc2*ref[1] + zc3*ref[2];
}

/** transfer from original element to reference element */
void THexaAffin::GetRefFromOrig(double X, double Y, double Z, double &xi,
                                double &eta, double &zeta)
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
void THexaAffin::GetRefFromOrig(const double *orig, double *ref)
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
void THexaAffin::GetOrigValues(BaseFunct3D BaseFunct, int N_Points,
                               const double *, const double *, const double *,
                               int N_Functs, QuadFormula3D QuadFormula)
{
  int i,j;
  double **refvaluesD000, **origvaluesD000;
  double **refvaluesD100, **origvaluesD100;
  double **refvaluesD010, **origvaluesD010;
  double **refvaluesD001, **origvaluesD001;
  double **origvaluesD200;
  double **origvaluesD110;
  double **origvaluesD020;
  double **origvaluesD101; //**refvaluesD101,
  double **origvaluesD011; //**refvaluesD011;
  double **origvaluesD002; //**refvaluesD002;
  double *refD000, *origD000;
  double *refD100, *origD100;
  double *refD010, *origD010;
  double *refD001, *origD001;
  double *aux;
//  double AllData[MaxN_BaseFunctions3D][5];
//  double GeoData[5][5];
  
  {
    TBaseFunct3D* bf = TFEDatabase3D::GetBaseFunct3D(BaseFunct);
    if(bf->GetBaseVectDim() != 1)
    {
      ErrMsg("Piola map for HexaAffin not yet implemented");
      exit(1);
    }
  }
  
  refvaluesD000=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula, D000);
  if(refvaluesD000==nullptr)
    {
      TFEDatabase3D::GetBaseFunct3D(BaseFunct)->MakeRefElementData(QuadFormula);
      refvaluesD000=TFEDatabase3D::GetRefElementValues(BaseFunct, 
                                                       QuadFormula, D000);
    }
  
  origvaluesD000=TFEDatabase3D::GetOrigElementValues(BaseFunct, D000);
  if(origvaluesD000==nullptr)
    {
      origvaluesD000 = new double* [MaxN_QuadPoints_3D];
      aux = new double [MaxN_QuadPoints_3D*MaxN_BaseFunctions3D];
      for(i=0;i<MaxN_QuadPoints_3D;i++)
        origvaluesD000[i] = aux+i*MaxN_BaseFunctions3D;
      TFEDatabase3D::RegisterOrigElementValues(BaseFunct, D000, origvaluesD000);
    }
  
  refvaluesD100=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula, D100);
  origvaluesD100=TFEDatabase3D::GetOrigElementValues(BaseFunct, D100);
  if(origvaluesD100==nullptr)
    {
      origvaluesD100 = new double* [MaxN_QuadPoints_3D];
      aux = new double [MaxN_QuadPoints_3D*MaxN_BaseFunctions3D];
      for(i=0;i<MaxN_QuadPoints_3D;i++)
        origvaluesD100[i] = aux+i*MaxN_BaseFunctions3D;
      TFEDatabase3D::RegisterOrigElementValues(BaseFunct, D100, origvaluesD100);
    }
  
  refvaluesD010=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula, D010);
  origvaluesD010=TFEDatabase3D::GetOrigElementValues(BaseFunct, D010);
  if(origvaluesD010==nullptr)
    {
      origvaluesD010 = new double* [MaxN_QuadPoints_3D];
      aux = new double [MaxN_QuadPoints_3D*MaxN_BaseFunctions3D];
      for(i=0;i<MaxN_QuadPoints_3D;i++)
        origvaluesD010[i] = aux+i*MaxN_BaseFunctions3D;
      TFEDatabase3D::RegisterOrigElementValues(BaseFunct, D010, origvaluesD010);
    } 
  
  refvaluesD001=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula, D001);
  origvaluesD001=TFEDatabase3D::GetOrigElementValues(BaseFunct, D001);
  if(origvaluesD001==nullptr)
    { 
      origvaluesD001 = new double* [MaxN_QuadPoints_3D];
      aux = new double [MaxN_QuadPoints_3D*MaxN_BaseFunctions3D];
      for(i=0;i<MaxN_QuadPoints_3D;i++)
        origvaluesD001[i] = aux+i*MaxN_BaseFunctions3D;
      TFEDatabase3D::RegisterOrigElementValues(BaseFunct, D001, origvaluesD001);
    } 
  
//  refvaluesD200=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula, D200);
  origvaluesD200=TFEDatabase3D::GetOrigElementValues(BaseFunct, D200);
  if(origvaluesD200==nullptr)
    { 
      origvaluesD200 = new double* [MaxN_QuadPoints_3D];
      aux = new double [MaxN_QuadPoints_3D*MaxN_BaseFunctions3D];
      for(i=0;i<MaxN_QuadPoints_3D;i++)
        origvaluesD200[i] = aux+i*MaxN_BaseFunctions3D;
      TFEDatabase3D::RegisterOrigElementValues(BaseFunct, D200, origvaluesD200);
    } 
  
//  refvaluesD110=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula, D110);
  origvaluesD110=TFEDatabase3D::GetOrigElementValues(BaseFunct, D110);
  if(origvaluesD110==nullptr)
    { 
      origvaluesD110 = new double* [MaxN_QuadPoints_3D];
      aux = new double [MaxN_QuadPoints_3D*MaxN_BaseFunctions3D];
      for(i=0;i<MaxN_QuadPoints_3D;i++)
        origvaluesD110[i] = aux+i*MaxN_BaseFunctions3D;
      TFEDatabase3D::RegisterOrigElementValues(BaseFunct, D110, origvaluesD110);
    } 
  
//  refvaluesD101=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula, D101);
  origvaluesD101=TFEDatabase3D::GetOrigElementValues(BaseFunct, D101);
  if(origvaluesD101==nullptr)
    { 
      origvaluesD101 = new double* [MaxN_QuadPoints_3D];
      aux = new double [MaxN_QuadPoints_3D*MaxN_BaseFunctions3D];
      for(i=0;i<MaxN_QuadPoints_3D;i++)
        origvaluesD101[i] = aux+i*MaxN_BaseFunctions3D;
      TFEDatabase3D::RegisterOrigElementValues(BaseFunct, D101, origvaluesD101);
    } 
  
//  refvaluesD011=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula, D011);
  origvaluesD011=TFEDatabase3D::GetOrigElementValues(BaseFunct, D011);
  if(origvaluesD011==nullptr)
    { 
      origvaluesD011 = new double* [MaxN_QuadPoints_3D];
      aux = new double [MaxN_QuadPoints_3D*MaxN_BaseFunctions3D];
      for(i=0;i<MaxN_QuadPoints_3D;i++)
        origvaluesD011[i] = aux+i*MaxN_BaseFunctions3D;
      TFEDatabase3D::RegisterOrigElementValues(BaseFunct, D011, origvaluesD011);
    } 
  
//  refvaluesD020=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula, D020);
  origvaluesD020=TFEDatabase3D::GetOrigElementValues(BaseFunct, D020);
  if(origvaluesD020==nullptr)
    {
      origvaluesD020 = new double* [MaxN_QuadPoints_3D];
      aux = new double [MaxN_QuadPoints_3D*MaxN_BaseFunctions3D];
      for(i=0;i<MaxN_QuadPoints_3D;i++)
        origvaluesD020[i] = aux+i*MaxN_BaseFunctions3D;
      TFEDatabase3D::RegisterOrigElementValues(BaseFunct, D020, origvaluesD020);
    } 
  
//  refvaluesD002=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula, D002);
  origvaluesD002=TFEDatabase3D::GetOrigElementValues(BaseFunct, D002);
  if(origvaluesD002==nullptr)
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
                               int N_Points, const double *xi,
                               const double *eta, const double *zeta, 
                               QuadFormula3D QuadFormula,
                               bool *Needs2ndDer)
{
  int i,j;
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
  double *aux;
  BaseFunct3D BaseFunct;
  int N_Functs;
  bool SecondDer;
  
  SecondDer = false;
  for(i=0;i<N_Sets;i++)
  {
    BaseFunct=BaseFuncts[i];
    N_Functs = TFEDatabase3D::GetBaseFunct3D(BaseFunct)->GetDimension();
    int BaseVectDim =
        TFEDatabase3D::GetBaseFunct3D(BaseFunct)->GetBaseVectDim();
      
    refvaluesD000=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula, D000);
    if(refvaluesD000==nullptr)
    {
      TFEDatabase3D::GetBaseFunct3D(BaseFunct)->MakeRefElementData(QuadFormula);
      refvaluesD000=TFEDatabase3D::GetRefElementValues(BaseFunct, 
                                                       QuadFormula, D000);
    }
      
    origvaluesD000=TFEDatabase3D::GetOrigElementValues(BaseFunct, D000);
    if(origvaluesD000==nullptr)
    {
      origvaluesD000 = new double* [MaxN_QuadPoints_3D];
      aux = new double [MaxN_QuadPoints_3D*N_Functs*BaseVectDim];
      for(j=0;j<MaxN_QuadPoints_3D;j++)
        origvaluesD000[j] = aux+j*N_Functs*BaseVectDim;
      TFEDatabase3D::RegisterOrigElementValues(BaseFunct, D000, origvaluesD000);
    }
      
    for(j=0;j<N_Points;j++)
    {
      refD000 = refvaluesD000[j];
      origD000 = origvaluesD000[j];
          
      if(BaseVectDim == 1)
      {
        // simply copy values
        memcpy(origD000, refD000, N_Functs*BaseVectDim*SizeOfDouble);
      }
      else
      {
        // do Piola transformation
        PiolaMapOrigFromRef(xi[j], eta[j], zeta[j], N_Functs, refD000,
                            origD000);
      }
    } // endfor j
      
    refvaluesD100=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula, D100);
    origvaluesD100=TFEDatabase3D::GetOrigElementValues(BaseFunct, D100);
    if(origvaluesD100==nullptr)
    {
      origvaluesD100 = new double* [MaxN_QuadPoints_3D];
      aux = new double [MaxN_QuadPoints_3D*N_Functs*BaseVectDim];
      for(j=0;j<MaxN_QuadPoints_3D;j++)
        origvaluesD100[j] = aux+j*N_Functs*BaseVectDim;
      TFEDatabase3D::RegisterOrigElementValues(BaseFunct, D100, origvaluesD100);
    }
      
    refvaluesD010=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula, D010);
    origvaluesD010=TFEDatabase3D::GetOrigElementValues(BaseFunct, D010);
    if(origvaluesD010==nullptr)
    {
       origvaluesD010 = new double* [MaxN_QuadPoints_3D];
       aux = new double [MaxN_QuadPoints_3D*N_Functs*BaseVectDim];
       for(j=0;j<MaxN_QuadPoints_3D;j++)
         origvaluesD010[j] = aux+j*N_Functs*BaseVectDim;
       TFEDatabase3D::RegisterOrigElementValues(BaseFunct, D010, origvaluesD010);
    }
      
    refvaluesD001=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula, D001);
    origvaluesD001=TFEDatabase3D::GetOrigElementValues(BaseFunct, D001);
    if(origvaluesD001==nullptr)
    {
      origvaluesD001 = new double* [MaxN_QuadPoints_3D];
      aux = new double [MaxN_QuadPoints_3D*N_Functs*BaseVectDim];
      for(j=0;j<MaxN_QuadPoints_3D;j++)
        origvaluesD001[j] = aux+j*N_Functs*BaseVectDim;
      TFEDatabase3D::RegisterOrigElementValues(BaseFunct, D001, origvaluesD001);
    }
    
    if(Needs2ndDer[i])
    {
      SecondDer = true;
      
      refvaluesD200=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula,
                                                       D200);
      origvaluesD200=TFEDatabase3D::GetOrigElementValues(BaseFunct, D200);
      if(origvaluesD200==nullptr)
      {
        origvaluesD200 = new double* [MaxN_QuadPoints_3D];
        aux = new double [MaxN_QuadPoints_3D*N_Functs*BaseVectDim];
        for(j=0;j<MaxN_QuadPoints_3D;j++)
          origvaluesD200[j] = aux+j*N_Functs*BaseVectDim;
        TFEDatabase3D::RegisterOrigElementValues(BaseFunct, D200,
                                                 origvaluesD200);
      }
          
      refvaluesD110=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula,
                                                       D110);
      origvaluesD110=TFEDatabase3D::GetOrigElementValues(BaseFunct, D110);
      if(origvaluesD110==nullptr)
      {
        origvaluesD110 = new double* [MaxN_QuadPoints_3D];
        aux = new double [MaxN_QuadPoints_3D*N_Functs*BaseVectDim];
        for(j=0;j<MaxN_QuadPoints_3D;j++)
          origvaluesD110[j] = aux+j*N_Functs*BaseVectDim;
        TFEDatabase3D::RegisterOrigElementValues(BaseFunct, D110,
                                                 origvaluesD110);
      }
          
      refvaluesD101=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula, 
                                                       D101);
      origvaluesD101=TFEDatabase3D::GetOrigElementValues(BaseFunct, D101);
      if(origvaluesD101==nullptr)
      {
        origvaluesD101 = new double* [MaxN_QuadPoints_3D];
        aux = new double [MaxN_QuadPoints_3D*N_Functs*BaseVectDim];
        for(j=0;j<MaxN_QuadPoints_3D;j++)
          origvaluesD101[j] = aux+j*N_Functs*BaseVectDim;
        TFEDatabase3D::RegisterOrigElementValues(BaseFunct, D101,
                                                 origvaluesD101);
      }
      
      refvaluesD011=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula,
                                                       D011);
      origvaluesD011=TFEDatabase3D::GetOrigElementValues(BaseFunct, D011);
      if(origvaluesD011==nullptr)
      {
        origvaluesD011 = new double* [MaxN_QuadPoints_3D];
        aux = new double [MaxN_QuadPoints_3D*N_Functs*BaseVectDim];
        for(j=0;j<MaxN_QuadPoints_3D;j++)
          origvaluesD011[j] = aux+j*N_Functs*BaseVectDim;
        TFEDatabase3D::RegisterOrigElementValues(BaseFunct, D011,
                                                 origvaluesD011);
      }
      
      refvaluesD020=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula, 
                                                       D020);
      origvaluesD020=TFEDatabase3D::GetOrigElementValues(BaseFunct, D020);
      if(origvaluesD020==nullptr)
      {
        origvaluesD020 = new double* [MaxN_QuadPoints_3D];
        aux = new double [MaxN_QuadPoints_3D*N_Functs*BaseVectDim];
        for(j=0;j<MaxN_QuadPoints_3D;j++)
          origvaluesD020[j] = aux+j*N_Functs*BaseVectDim;
        TFEDatabase3D::RegisterOrigElementValues(BaseFunct, D020, 
                                                 origvaluesD020);
      }
      
      refvaluesD002=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula,
                                                       D002);
      origvaluesD002=TFEDatabase3D::GetOrigElementValues(BaseFunct, D002);
      if(origvaluesD002==nullptr)
      {
        origvaluesD002 = new double* [MaxN_QuadPoints_3D];
        aux = new double [MaxN_QuadPoints_3D*N_Functs*BaseVectDim];
        for(j=0;j<MaxN_QuadPoints_3D;j++)
          origvaluesD002[j] = aux+j*N_Functs*BaseVectDim;
        TFEDatabase3D::RegisterOrigElementValues(BaseFunct, D002,
                                                 origvaluesD002);
      } // endfor Needs2ndDer[i]
    } // endfor i
  }

  // D100, D010 and D001
  for(i=0;i<N_Sets;i++)
  {
    BaseFunct=BaseFuncts[i];
    int BaseVectDim =
        TFEDatabase3D::GetBaseFunct3D(BaseFunct)->GetBaseVectDim();
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
            
      // inverse of transformation (DF^{-1}) multiplied by determinant^{-1}
      // for the inverse you can use maxima and type
      // A:matrix([x1, x2, x3],[y1, y2, y3],[z1, z2, z3]);
      // B:invert(A);
      // ratsimp(A.B); // check
      double i11 = (yc2*zc3 - yc3*zc2) * rec_detjk;
      double i12 = (yc3*zc1 - yc1*zc3) * rec_detjk;
      double i13 = (yc1*zc2 - yc2*zc1) * rec_detjk;
      double i21 = (xc3*zc2 - xc2*zc3) * rec_detjk;
      double i22 = (xc1*zc3 - xc3*zc1) * rec_detjk;
      double i23 = (xc2*zc1 - xc1*zc2) * rec_detjk;
      double i31 = (xc2*yc3 - xc3*yc2) * rec_detjk;
      double i32 = (xc3*yc1 - xc1*yc3) * rec_detjk;
      double i33 = (xc1*yc2 - xc2*yc1) * rec_detjk;
      
      if(BaseVectDim == 1)
      {
        for(int k = 0; k < N_Functs; k++)
        {
          origD100[k] = i11 * refD100[k] + i12 * refD010[k] + i13 * refD001[k];
          origD010[k] = i21 * refD100[k] + i22 * refD010[k] + i23 * refD001[k];
          origD001[k] = i31 * refD100[k] + i32 * refD010[k] + i33 * refD001[k];
        } // endfor k
      }
      else
      {
        // Piola transformation
        // this is DF divided by determinant
        double a11 = xc1 * rec_detjk;
        double a12 = xc2 * rec_detjk;
        double a13 = xc3 * rec_detjk;
        double a21 = yc1 * rec_detjk;
        double a22 = yc2 * rec_detjk;
        double a23 = yc3 * rec_detjk;
        double a31 = zc1 * rec_detjk;
        double a32 = zc2 * rec_detjk;
        double a33 = zc3 * rec_detjk;
        for(int k = 0; k < N_Functs; k++)
        {
          double p11 = a11 * refD100[k             ] 
                     + a12 * refD100[k + N_Functs  ] 
                     + a13 * refD100[k + 2*N_Functs];
          double p21 = a11 * refD010[k             ] 
                     + a12 * refD010[k + N_Functs  ] 
                     + a13 * refD010[k + 2*N_Functs];
          double p31 = a11 * refD001[k             ] 
                     + a12 * refD001[k + N_Functs  ] 
                     + a13 * refD001[k + 2*N_Functs];
          origD100[k] = i11*p11 + i12*p21 + i13*p31;
          origD010[k] = i21*p11 + i22*p21 + i23*p31;
          origD001[k] = i31*p11 + i32*p21 + i33*p31;
          
          double p12 = a21 * refD100[k             ] 
                     + a22 * refD100[k + N_Functs  ] 
                     + a23 * refD100[k + 2*N_Functs];
          double p22 = a21 * refD010[k             ] 
                     + a22 * refD010[k + N_Functs  ] 
                     + a23 * refD010[k + 2*N_Functs];
          double p32 = a21 * refD001[k             ] 
                     + a22 * refD001[k + N_Functs  ] 
                     + a23 * refD001[k + 2*N_Functs];
          origD100[k + N_Functs] = i11*p12 + i12*p22 + i13*p32;
          origD010[k + N_Functs] = i21*p12 + i22*p22 + i23*p32;
          origD001[k + N_Functs] = i31*p12 + i32*p22 + i33*p32;
          
          double p13 = a31 * refD100[k             ] 
                     + a32 * refD100[k + N_Functs  ] 
                     + a33 * refD100[k + 2*N_Functs];
          double p23 = a31 * refD010[k             ] 
                     + a32 * refD010[k + N_Functs  ] 
                     + a33 * refD010[k + 2*N_Functs];
          double p33 = a31 * refD001[k             ] 
                     + a32 * refD001[k + N_Functs  ] 
                     + a33 * refD001[k + 2*N_Functs];
          origD100[k + 2*N_Functs] = i11*p13 + i12*p23 + i13*p33;
          origD010[k + 2*N_Functs] = i21*p13 + i22*p23 + i23*p33;
          origD001[k + 2*N_Functs] = i31*p13 + i32*p23 + i33*p33;
        }
      }
    } // endfor j
  } // endfor i
   
  // leave if no second derivatives are needed
  if(!SecondDer) return;
  
  double GeoData[6][6];
  double Eye[6][6];
  // reset matrices
  memset(GeoData, 0, 36*SizeOfDouble);
  memset(Eye, 0, 36*SizeOfDouble);
  Eye[0][0] = 1;
  Eye[1][1] = 1;
  Eye[2][2] = 1;
  Eye[3][3] = 1;
  Eye[4][4] = 1;
  Eye[5][5] = 1;
  
  GeoData[0][0] = xc1*xc1;
  GeoData[0][1] = 2*xc1*yc1;
  GeoData[0][2] = 2*xc1*zc1;
  GeoData[0][3] = yc1*yc1;
  GeoData[0][4] = 2*yc1*zc1;
  GeoData[0][5] = zc1*zc1;
  
  GeoData[1][0] = xc1*xc2;
  GeoData[1][1] = yc1*xc2 + xc1*yc2; 
  GeoData[1][2] = xc1*zc2 + xc2*zc1;
  GeoData[1][3] = yc1*yc2;
  GeoData[1][4] = yc1*zc2 + yc2*zc1;
  GeoData[1][5] = zc1*zc2;
  
  GeoData[2][0] = xc1*xc3;
  GeoData[2][1] = xc1*yc3 + xc3*yc1;
  GeoData[2][2] = xc1*zc3 + xc3*zc1;
  GeoData[2][3] = yc1*yc3;
  GeoData[2][4] = yc1*zc3 + yc3*zc1;
  GeoData[2][5] = zc1*zc3;
  
  GeoData[3][0] = xc2*xc2;
  GeoData[3][1] = 2*xc2*yc2;
  GeoData[3][2] = 2*xc2*zc2;
  GeoData[3][3] = yc2*yc2;
  GeoData[3][4] = 2*yc2*zc2;
  GeoData[3][5] = zc2*zc2;
  
  GeoData[4][0] = xc2*xc3;
  GeoData[4][1] = xc2*yc3 + xc3*yc2;
  GeoData[4][2] = xc2*zc3 + xc3*zc2;
  GeoData[4][3] = yc2*yc3;
  GeoData[4][4] = yc2*zc3 + yc3*zc2;
  GeoData[4][5] = zc2*zc3;
  
  GeoData[5][0] = xc3*xc3;
  GeoData[5][1] = 2*xc3*yc3;
  GeoData[5][2] = 2*xc3*zc3;
  GeoData[5][3] = yc3*yc3;
  GeoData[5][4] = 2*yc3*zc3;
  GeoData[5][5] = zc3*zc3;
  
  // subroutine for solving a multiple systems of linear equations
  // void SolveMultipleSystems(double *a, double *b, int N_Eqn,
  //                        int LDA, int LDB, int N_Rhs);
  //SolveMultipleSystems((double *)GeoData, (double *)Eye, 3,
  //                 3, 3, 3);
  //SolveLinearSystem((double *)GeoData, (double *)Eye, 3, 3);
              
  SolveMultipleSystemsNew((double *)GeoData, (double *)Eye, 6, 6, 6, 6);
  
  for(i=0;i<N_Sets;i++)
  {
    if(Needs2ndDer[i])
    {
      BaseFunct=BaseFuncts[i];
      N_Functs = TFEDatabase3D::GetBaseFunct3D(BaseFunct)->GetDimension();
      int BaseVectDim =
        TFEDatabase3D::GetBaseFunct3D(BaseFunct)->GetBaseVectDim();
      if(BaseVectDim != 1)
      {
        ErrMsg("second derivatives of vector valued basis function not " <<
               "supported");
        exit(0);
      }
      
      refvaluesD200 = TFEDatabase3D::GetRefElementValues
      (BaseFunct, QuadFormula, D200);
      origvaluesD200 = TFEDatabase3D::GetOrigElementValues
      (BaseFunct, D200);
      
      refvaluesD110 = TFEDatabase3D::GetRefElementValues
      (BaseFunct, QuadFormula, D110);
      origvaluesD110 = TFEDatabase3D::GetOrigElementValues
      (BaseFunct, D110);
      
      refvaluesD101 = TFEDatabase3D::GetRefElementValues
      (BaseFunct, QuadFormula, D101);
      origvaluesD101 = TFEDatabase3D::GetOrigElementValues
      (BaseFunct, D101);
      
      refvaluesD020 = TFEDatabase3D::GetRefElementValues
      (BaseFunct, QuadFormula, D020);
      origvaluesD020 = TFEDatabase3D::GetOrigElementValues
      (BaseFunct, D020);
      
      refvaluesD011 = TFEDatabase3D::GetRefElementValues
      (BaseFunct, QuadFormula, D011);
      origvaluesD011 = TFEDatabase3D::GetOrigElementValues
      (BaseFunct, D011);
      
      refvaluesD002 = TFEDatabase3D::GetRefElementValues
      (BaseFunct, QuadFormula, D002);
      origvaluesD002 = TFEDatabase3D::GetOrigElementValues
      (BaseFunct, D002);
      
      for(j=0;j<N_Points;j++)
      {
        double *refD200 = refvaluesD200[j];
        double *refD110 = refvaluesD110[j];
        double *refD101 = refvaluesD101[j];
        double *refD020 = refvaluesD020[j];
        double *refD011 = refvaluesD011[j];
        double *refD002 = refvaluesD002[j];
        
        double *origD200 = origvaluesD200[j];
        double *origD110 = origvaluesD110[j];
        double *origD101 = origvaluesD101[j];
        double *origD020 = origvaluesD020[j];
        double *origD011 = origvaluesD011[j];
        double *origD002 = origvaluesD002[j];
        
        for(int k = 0; k < N_Functs; k++)
        {
          double r200 = refD200[k];
          double r110 = refD110[k];
          double r101 = refD101[k];
          double r020 = refD020[k];
          double r011 = refD011[k];
          double r002 = refD002[k];
          
          double o200, o110, o101, o020, o011, o002;
          o200 =  Eye[0][0]*r200 + Eye[0][1]*r110 + Eye[0][2]*r101 
                + Eye[0][3]*r020 + Eye[0][4]*r011 + Eye[0][5]*r002;
          o110 =  Eye[1][0]*r200 + Eye[1][1]*r110 + Eye[1][2]*r101 
                + Eye[1][3]*r020 + Eye[1][4]*r011 + Eye[1][5]*r002;
          o101 =  Eye[2][0]*r200 + Eye[2][1]*r110 + Eye[2][2]*r101 
                + Eye[2][3]*r020 + Eye[2][4]*r011 + Eye[2][5]*r002;
          o020 =  Eye[3][0]*r200 + Eye[3][1]*r110 + Eye[3][2]*r101 
                + Eye[3][3]*r020 + Eye[3][4]*r011 + Eye[3][5]*r002;
          o011 =  Eye[4][0]*r200 + Eye[4][1]*r110 + Eye[4][2]*r101 
                + Eye[4][3]*r020 + Eye[4][4]*r011 + Eye[4][5]*r002;
          o002 =  Eye[5][0]*r200 + Eye[5][1]*r110 + Eye[5][2]*r101 
                + Eye[5][3]*r020 + Eye[5][4]*r011 + Eye[5][5]*r002;
          
          origD200[k] = o200;
          origD110[k] = o110;
          origD101[k] = o101;
          origD020[k] = o020;
          origD011[k] = o011;
          origD002[k] = o002;
        } // endfor k
      } // endif
    } // endfor j
  } // endfor i
}
/** calculate functions and derivatives from reference element
    to original element */
void THexaAffin::GetOrigValues(double xi, double eta, double zeta,
                               int N_BaseFunct,
                               const double *uref, const double *uxiref,
                               const double *uetaref,const  double *uzetaref,
                               double *uorig, double *uxorig, double *uyorig, double *uzorig,
                               int _BaseVectDim)
{
  if(_BaseVectDim == 1)
  {
    int i;
    // D000
    for(i=0;i<N_BaseFunct;i++)
      uorig[i] = uref[i];
    
    // D100, D010 and D001
    for(i=0;i<N_BaseFunct;i++)
    {
      uxorig[i] = ( (yc2*zc3 - yc3*zc2)*uxiref[i] 
                  + (yc3*zc1 - yc1*zc3)*uetaref[i] 
                  + (yc1*zc2 - yc2*zc1)*uzetaref[i] ) *rec_detjk;
      uyorig[i] = ( (xc3*zc2 - xc2*zc3)*uxiref[i] 
                  + (xc1*zc3 - xc3*zc1)*uetaref[i] 
                  + (xc2*zc1 - xc1*zc2)*uzetaref[i] ) *rec_detjk;
      uzorig[i] = ( (xc2*yc3 - xc3*yc2)*uxiref[i] 
                  + (xc3*yc1 - xc1*yc3)*uetaref[i] 
                  + (xc1*yc2 - xc2*yc1)*uzetaref[i] ) *rec_detjk;
    } // endfor i
  }
  else if(_BaseVectDim == 3)
  {
    // D000
    PiolaMapOrigFromRef(xi, eta, zeta, N_BaseFunct, uref, uorig);
    
    // D100, D010, D001
    // not yet implemented
  }
  else
  {
    ErrMsg("unknown basis function dimension");
    exit(0);
  }
}

void THexaAffin::SetCell(const TBaseCell *cell)
{

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
void THexaAffin::GetOuterNormal(int j, double, double, double &n1, double &n2,
                                double &n3) const
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
void THexaAffin::GetTangentVectors(int j, double, double,
        double &t11, double &t12, double &t13,
        double &t21, double &t22, double &t23) const
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


/** Piola transformation for vectorial basis functions */
void THexaAffin::PiolaMapOrigFromRef(double, double, double, int N_Functs,
                                     const double *refD000, double *origD000)
{
  double a11 = xc1 * rec_detjk;
  double a12 = xc2 * rec_detjk;
  double a13 = xc3 * rec_detjk;
  double a21 = yc1 * rec_detjk;
  double a22 = yc2 * rec_detjk;
  double a23 = yc3 * rec_detjk;
  double a31 = zc1 * rec_detjk;
  double a32 = zc2 * rec_detjk;
  double a33 = zc3 * rec_detjk;
  for(int k = 0; k < N_Functs; k++)
  {
    // three components:
    origD000[k             ] = a11 * refD000[k             ] 
                             + a12 * refD000[k +   N_Functs] 
                             + a13 * refD000[k + 2*N_Functs]; 
    origD000[k + N_Functs  ] = a21 * refD000[k             ] 
                             + a22 * refD000[k +   N_Functs] 
                             + a23 * refD000[k + 2*N_Functs]; 
    origD000[k + 2*N_Functs] = a31 * refD000[k             ] 
                             + a32 * refD000[k +   N_Functs] 
                             + a33 * refD000[k + 2*N_Functs];
  }
}



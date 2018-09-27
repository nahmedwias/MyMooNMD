#include "DarcyMixed.h"
#include <Database.h>

/* ======================================================================
   depending on the global orientation of the normals at inner edges the sign 
   of the basis functions need to be changed. This is only important for 
   H-div conforming finite elements (such as Raviart-Thomas and 
   Brezzi-Douglas-Marini).
   
   This information is stored in the FE-Descriptor. However the FE-Descriptor 
   is not accessible in the assembling routine DarcyRaviartThomas. Therefore we
   use global parameters in TDatabase. These need to be set for every cell.
   
   Here it is assumed that Raviart-Thomas elements of order 0,1,2, or 3 on
   triangles or quadrilaterals are used
*/
template <int d> int GetSignOfThisDOF(int N_DOF, int DOF);

template <>
int GetSignOfThisDOF<2>(int N_DOF, int DOF)
{
  switch (N_DOF)
  {
  case 3:// Raviart-Thomas zeroth order, triangles
    return TDatabase::ParamDB->NORMAL_ORIENTATION_TRIA[DOF];
    break;
    // note that RT2 on quadrilaterals and RT3 on triangles have 24 basis functions
  case 4:// Raviart-Thomas zeroth order, quadrilaterals
    return TDatabase::ParamDB->NORMAL_ORIENTATION_QUAD[DOF];
    break;
  case 6: // BDM first order, triangles
    // degree of freedom on an edge, no inner degree of freedom
    return TDatabase::ParamDB->NORMAL_ORIENTATION_TRIA[(DOF-DOF%2)/2];
    break;
  case 8:
    // Raviart-Thomas first order, triangles
    if(TDatabase::ParamDB->VELOCITY_SPACE == 1001)
    {
      if(DOF<6) // degree of freedom on an edge
        return TDatabase::ParamDB->NORMAL_ORIENTATION_TRIA[(DOF-DOF%2)/2];
      else // inner degree of freedom
        return 1;
    }
    else if (TDatabase::ParamDB->VELOCITY_SPACE == 1011) //BDM first order, quadrilaterals
    {
      return TDatabase::ParamDB->NORMAL_ORIENTATION_QUAD[(DOF-DOF%2)/2];
    }
    break;
  case 12:
    // Raviart-Thomas first order, quadrilaterals
    if(TDatabase::ParamDB->VELOCITY_SPACE == 1001)
    {
    if(DOF<8) // degree of freedom on an edge
      return TDatabase::ParamDB->NORMAL_ORIENTATION_QUAD[(DOF-DOF%2)/2];
    else // inner degree of freedom
      return 1;
    }
    else if (TDatabase::ParamDB->VELOCITY_SPACE == 1012) //BDM second order, triangles
    {
      if(DOF<9) // degree of freedom on an edge
        return TDatabase::ParamDB->NORMAL_ORIENTATION_TRIA[(DOF-DOF%3)/3];
      else // inner degree of freedom
        return 1;
    }
    break;
  case 14://BDM second order, quadrilaterals
    if(DOF<12) // degree of freedom on an edge
      return TDatabase::ParamDB->NORMAL_ORIENTATION_QUAD[(DOF-DOF%3)/3];
    else // inner degree of freedom
      return 1;
    break;
  case 15:// Raviart-Thomas second order, triangles
    if(DOF<9) // degree of freedom on an edge
      return TDatabase::ParamDB->NORMAL_ORIENTATION_TRIA[(DOF-DOF%3)/3];
    else // inner degree of freedom
      return 1;
    break;
  case 20:// BDM third order. triangles
    if(DOF<12) // degree of freedom on an edge
      return TDatabase::ParamDB->NORMAL_ORIENTATION_TRIA[(DOF-DOF%4)/4];
    else // inner degree of freedom
      return 1;
    break;
  case 22://BDM third order, quadrilaterals
    if(DOF<16) // degree of freedom on an edge
      return TDatabase::ParamDB->NORMAL_ORIENTATION_QUAD[(DOF-DOF%4)/4];
    else // inner degree of freedom
      return 1;
    break;
  case 24:
    if(TDatabase::ParamDB->VELOCITY_SPACE == 1002)
    {
      // Raviart-Thomas second order, quadrilaterals
      if(DOF<12) // degree of freedom on an edge
        return TDatabase::ParamDB->NORMAL_ORIENTATION_QUAD[(DOF-DOF%3)/3];
      else // inner degree of freedom
        return 1;
    }
    else if(TDatabase::ParamDB->VELOCITY_SPACE == 1003)
    {
      // Raviart-Thomas third order, triangles
      if(DOF<12) // degree of freedom on an edge
      return TDatabase::ParamDB->NORMAL_ORIENTATION_TRIA[(DOF-DOF%4)/4];
      else // inner degree of freedom
        return 1;
    }
    else
    {
      ErrThrow("VELOCITY_SPACE has to be set to either 1002 or 1003\n");
    }
    break;
  case 40:// Raviart-Thomas third order, quadrilaterals
    if(DOF<16) // degree of freedom on an edge
      return TDatabase::ParamDB->NORMAL_ORIENTATION_QUAD[(DOF-DOF%4)/4];
    else
      return 1;
    break;
  default:
       ErrThrow("WARNING: Unknown Raviart-Thomas or BDM element !", N_DOF);
       break;
  }
  //dummy return, should not be reached
  return -4711;
}

// ======================================================================
// matrix and rhs for mixed formulation 
// ======================================================================
// this method returns the correct sign of the normal to the i-th face
// of the current cell. N is the total number of local degrees of freedom.
// This is a really dirty hack which is necessary because we don't have
// access to the cell (where such information is stored) during the local
// assembling routine. Therefore this information is written into the 
// database and used here.
template <>
int GetSignOfThisDOF<3>(int N, int i)
{
  switch(N)
  {
    case 4: // Raviart-Thomas zeroth order on tetrahedra
      return TDatabase::ParamDB->NORMAL_ORIENTATION_TETRA[i];
      break;
    case 6: // Raviart-Thomas zeroth order on hexahedra
      return TDatabase::ParamDB->NORMAL_ORIENTATION_HEXA[i];
      break;
    case 12: //Brezzi-Douglas-Duran-Fortin first order on tetrahedra
      return TDatabase::ParamDB->NORMAL_ORIENTATION_TETRA[(i-i%3)/3];
      break;
    case 15: // Raviart-Thomas first order on tetrahedra
      if(i<12)// degree of freedom on the faces
        return TDatabase::ParamDB->NORMAL_ORIENTATION_TETRA[(i-i%3)/3];
      else // inner degree of freedom
        return 1;
      break;
    case 18: //Brezzi-Douglas-Duran-Fortin first order on hexahedra
      return TDatabase::ParamDB->NORMAL_ORIENTATION_HEXA[(i-i%3)/3];
      break;
    case 30: //Brezzi-Douglas-Duran-Fortin second order on tetrahedra
      if(i<24) // degree of freedom on the faces
        return TDatabase::ParamDB->NORMAL_ORIENTATION_TETRA[(i-i%6)/6];
      else // inner degree of freedom
        return 1;
      break;
    case 36:
      if(TDatabase::ParamDB->VELOCITY_SPACE == 1001)
      {// Raviart-Thomas first order on hexahedra
      if(i<24) // degree of freedom on the faces
        return TDatabase::ParamDB->NORMAL_ORIENTATION_HEXA[(i-i%4)/4];
      else // inner degree of freedom
        return 1;
      }else if (TDatabase::ParamDB->VELOCITY_SPACE == 1002)
      {// Raviart-Thomas second order on tetrahedra
      if(i<24) // degree of freedom on the faces
        return TDatabase::ParamDB->NORMAL_ORIENTATION_TETRA[(i-i%6)/6];
      else // inner degree of freedom
        return 1;
      }
      break;
    case 39: //Brezzi-Douglas-Duran-Fortin second order on hexahedra
      if(i<36) //degree of freedom on the faces
        return TDatabase::ParamDB->NORMAL_ORIENTATION_HEXA[(i-i%6)/6];
      else //inner degree of freedom
        return 1;
      break;
    case 60: //Brezzi-Douglas-Duran-Fortin third order on tetrahedra
      if(i<40) // degree of freedom on the faces
        return TDatabase::ParamDB->NORMAL_ORIENTATION_TETRA[(i-i%10)/10];
      else // inner degree of freedom
        return 1;
      break;
    case 70: //Raviart-Thomas third order on tetrahedra
      if(i<40) // degree of freedom on the faces
        return TDatabase::ParamDB->NORMAL_ORIENTATION_TETRA[(i-i%10)/10];
      else // inner degree of freedom
        return 1;
      break;
    case 72: //Brezzi-Douglas-Duran-Fortin third order on hexahedra
      if(i<60) //degree of freedom on the faces
        return TDatabase::ParamDB->NORMAL_ORIENTATION_HEXA[(i-i%10)/10];
      else //inner degree of freedom
        return 1;
      break;
    case 108: // Raviart-Thomas second order on hexahedra
        if(i<54) // degree of freedom on the faces
                return TDatabase::ParamDB->NORMAL_ORIENTATION_HEXA[(i-i%9)/9];
        else // inner degree of freedom
                return 1;
        break;
    default:
      ErrThrow("unsupported number of degrees of freedom for mixed elements ",
               N);
      break;
  }
  return 1;
}

// ======================================================================
// (DarcyType 1)
// Standard Galerkin with Raviart-Thomas (RT) or Brezzi-Douglas-Marini (BDM)
// elements
// ======================================================================
template <int d>
void BilinearAssembleDarcyGalerkin(double Mult, double *coeff, double *param,
                                   double hK, double **OrigValues,
                                   int *N_BaseFuncts, double ***LocMatrices,
                                   double **LocRhs)
{
  double val;
  double ansatz, ansatz_x, ansatz_y, ansatz_z;
  double test, test_x, test_y, test_z;
  double test_x_100, test_y_010, test_z_001;
  double test_div;
  
  // ( A  B1 )   ( 0 2 )
  // ( B2 C  )   ( 3 1 )
  
  double **MatrixA = LocMatrices[0];
//  double **MatrixC = LocMatrices[1];
  double **MatrixB1 = LocMatrices[2];
  double **MatrixB2 = LocMatrices[3];
  
  double *Rhs0 = LocRhs[0];
  double *Rhs1 = LocRhs[1];
  
  int N_U = N_BaseFuncts[0];
  int N_P = N_BaseFuncts[1];

  double *Orig0 = OrigValues[0];   // u
  double *Orig1 = OrigValues[1];   // p
  double *Orig2 = OrigValues[2];   // u_x
  double *Orig3 = OrigValues[3];   // u_y
  double *Orig4 = OrigValues[4];   // u_z (unused in 2D)
  

  double c0 = coeff[0];  // sigma
  double f1 = coeff[1];  // f1
  double f2 = coeff[2];  // f2
  double f3 = coeff[3];  // unused 2D, f3 in 3D
  double g = d == 2 ? coeff[3] : coeff[4];  // g(x,y)

  int *newsign = new int[N_U];
  for(int i=0;i<N_U;i++)
  {
    // here check whether signs should be inverted
    newsign[i] = GetSignOfThisDOF<d>(N_U,i);
  }
  // A, B1, B2
  for(int i=0;i<N_U;i++)
  {
    // A:
    test_x = newsign[i]*Orig0[i];
    test_y = newsign[i]*Orig0[N_U+i];
    test_z = d == 2 ? 0. : newsign[i]*Orig0[2*N_U+i];
    
    Rhs0[i] += Mult*(f1*test_x + f2*test_y + f3*test_z);
    
    for(int j=0;j<N_U;j++)
    {
      ansatz_x = newsign[j]*Orig0[j];
      ansatz_y = newsign[j]*Orig0[N_U+j];
      ansatz_z = d == 2 ? 0. : newsign[j]*Orig0[2*N_U+j];

      // A: u_x v_x + u_y v_y
      val  = c0*(test_x*ansatz_x + test_y*ansatz_y + test_z*ansatz_z);
      MatrixA[i][j] += Mult * val;
    }
    // B1, B2:
    test_x_100 = newsign[i]*Orig2[i];
    test_y_010 = newsign[i]*Orig3[N_U+i];
    test_z_001 = d == 2 ? 0. : newsign[i]*Orig4[2*N_U+i];
    test_div = test_x_100 + test_y_010 + test_z_001;
    for(int j=0;j<N_P;j++)
    {
      ansatz = Orig1[j];
      val = Mult*test_div*ansatz;
      // (p div v)
      MatrixB1[i][j] -= val;
      // (q, div u)
      MatrixB2[j][i] -= val; // slow (consider moving this into another loop)
    }
  }
  
  for(int i=0;i<N_P;i++)
  {
    test = Orig1[i];
    // assemble rhs: div u = g
    // rhs: -(g,q)
    Rhs1[i] -= Mult*test*g;
  }
  delete [] newsign;
}

#ifdef __3D__
template void BilinearAssembleDarcyGalerkin<3>(double Mult, double *coeff, double *param,
                                   double hK, double **OrigValues,
                                   int *N_BaseFuncts, double ***LocMatrices,
                                   double **LocRhs);
#else
template void BilinearAssembleDarcyGalerkin<2>(double Mult, double *coeff, double *param,
                                   double hK, double **OrigValues,
                                   int *N_BaseFuncts, double ***LocMatrices,
                                   double **LocRhs);
#endif



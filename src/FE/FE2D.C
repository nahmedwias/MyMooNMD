// =======================================================================
// @(#)FE2D.C        1.1 10/30/98
//
// Class:       TFE2D
// Purpose:     store all information for one finite element class
//
// Author:      Gunar Matthies  09.07.98
//
// =======================================================================
#ifdef _MPI
#  include "mpi.h"
# endif

#include <FE2D.h>
#include <FEDatabase2D.h>
#include <NodalFunctional2D.h>
#include <Database.h>


/** constructor */
TFE2D::TFE2D()
{
}

/** constructor with data */
TFE2D::TFE2D(BaseFunct2D basefunct_id, NodalFunctional2D nodalfunctional_id,
         RefTrans2D reftransid, FEDesc2D fedesc_id, int n_info)
{
  BaseFunct_ID = basefunct_id;
  BaseFunct = TFEDatabase2D::GetBaseFunct2D(BaseFunct_ID);

  NodalFunctional_ID = nodalfunctional_id;
  NodalFunctional  = TFEDatabase2D::GetNodalFunctional2D(NodalFunctional_ID);

  RefTransID = reftransid;

  FEDesc_ID = fedesc_id;
  FEDesc = TFEDatabase2D::GetFEDesc2D(FEDesc_ID);

  N_Info = n_info;
  N_DOF = FEDesc->GetN_DOF();

  Size = N_Info + N_DOF;
}

/** check N[i](b[j]) = delta[ij] */
void TFE2D::CheckNFandBF() const
{
  #ifdef _MPI
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if(rank == TDatabase::ParamDB->Par_P0)
  #endif 
  {
   Output::print<3>("CheckNFandBF: BaseFunct_ID: ", this->BaseFunct_ID,
                 " NodalFunctional_ID: ", this->NodalFunctional_ID);
  }
  
  int N_Points;
  double *xi, *eta;
  this->NodalFunctional->GetPointsForAll(N_Points, xi, eta);
  
  // dimension of the basis function (usually 1, for H(div) elements it is 2)
  int baseVectDim = this->BaseFunct->GetBaseVectDim();
  // number of basis functions, this is the length of the array needed to 
  // evaluate the basis functions (therefore the factor baseVectDim)
  int nBaseFunc = this->GetN_DOF() * baseVectDim;
  
  double AllPointValues[N_Points][nBaseFunc];
  for(int k = 0; k < N_Points; k++)
  {
    BaseFunct->GetDerivatives(D00, xi[k], eta[k], AllPointValues[k]);
  }
  
  double PointValues[N_Points * baseVectDim];
  double FunctionalValues[this->GetN_DOF()];
  for(int k = 0; k < this->N_DOF; k++)
  {
    for(int l = 0; l < N_Points; l++)
    {
      for(int i = 0; i < baseVectDim; ++i)
      {
        PointValues[l + i * N_Points] 
          = AllPointValues[l][k + i*this->GetN_DOF()];
      }
    }
    
    NodalFunctional->GetAllFunctionals(nullptr, nullptr, PointValues,
                                       FunctionalValues);
    
    for(int i = 0; i < this->N_DOF; i++)
    {
      if( fabs(FunctionalValues[i]) < 1e-10 )
      {
        FunctionalValues[i] = 0;
      }
      //Output::print(k, " ", i, " ", FunctionalValues[i]);
      if( i == k && fabs(FunctionalValues[i]-1) > 1e-8 )
      {
        Output::print<3>("BF: ", k, " NF: ", i, " ", FunctionalValues[i]);
      }
      if( i != k && fabs(FunctionalValues[i]-0) > 1e-8 )
      {
        Output::print<3>("BF: ", k, " NF: ", i, " ", FunctionalValues[i]);
      }
    }
  }
}

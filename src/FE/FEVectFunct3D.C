// =======================================================================
// @(#)FEVectFunct3D.C        1.3 07/22/99
// 
// Class:       TFEVectFunct3D
// Purpose:     a function from a finite element space in 3D
//
// Author:      Gunar Matthies (13.07.2000)
//
// History:     start of implementation 13.07.2000 (Gunar Matthies)
//
//              WriteSol/ReadSol    13.12.10 (Sashikumaar Ganesan)
// =======================================================================
#ifdef _MPI
# include "mpi.h"
#endif

#include <FEVectFunct3D.h>
#include <FEDatabase3D.h>
#include <NodalFunctional3D.h>
#include <HexaAffin.h>
#include <HexaTrilinear.h>
#include <HexaIsoparametric.h>
#include <TetraAffin.h>
#include <TetraIsoparametric.h>
#include <Database.h>
#include <BoundFace.h>
#include "BaseCell.h"

#include <fstream>
#include <stdlib.h>
#include <sstream>
#include <MooNMD_Io.h>
// #include <malloc.h>
#include <dirent.h> 
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>

/// Default constructor. Constructs an empty object.
TFEVectFunct3D::TFEVectFunct3D()
{
  N_Components = 0;
}



/** constructor with vector initialization */
TFEVectFunct3D::TFEVectFunct3D(std::shared_ptr<const TFESpace3D> fespace3D,
                               const std::string& name,
                               const std::string& description, double *values,
                               int length, int n_components)
: TFEFunction3D(fespace3D, name, description, values, length)
{
  N_Components = n_components;
}

TFEVectFunct3D& TFEVectFunct3D::operator=( const TFEVectFunct3D & other)
{
  //call base class copy assignment
  TFEFunction3D::operator=(other);

  this->N_Components  = other.N_Components;

  return *this;
}

/** compute integral and measure */
void TFEVectFunct3D::compute_flux(int surface_id, double& flux) const
{
  flux = 0.;

  auto coll = FESpace3D->GetCollection();

  for(int i=0; i< coll->GetN_Cells(); i++) {
    auto cell = coll->GetCell(i); //boundaryCells[i];

    int *DOF = FESpace3D->GetGlobalDOF(cell->GetCellIndex());
    for(size_t joint_id=0; joint_id< (size_t) cell->GetN_Faces(); joint_id++) {
      TJoint* joint = cell->GetJoint(joint_id);
            
      if (joint->GetType() == BoundaryFace ||
	  joint->GetType() == IsoBoundFace) {
                
	// convert the joint to an object of BoundFace type
	TBoundFace *boundface = (TBoundFace *)joint;

	/// check if the face is on the desired component
	if (boundface->GetBoundComp()->get_physical_id()==surface_id) {


	  // ===================
	  // get quadrature data
	  // ===================
	  int nFaceVertices = cell->getNumberOfFaceVertices(joint_id);
	  // set quadrature formula and compute quadrature info
	  FE3D FEId = FESpace3D->GetFE3D(cell->GetCellIndex(),cell);
	  int fe_degree = TFEDatabase3D::GetPolynomialDegreeFromFE3D(FEId);
                        
	  QuadFormula2D FaceQuadFormula; //=BaryCenterTria;
	  switch(nFaceVertices) {
	  case 3:
	    // triangular face
	    FaceQuadFormula = TFEDatabase3D::GetQFTriaFromDegree(2*fe_degree);
	    FaceQuadFormula = Gauss3Tria;
	    break;
	  case 4:
	    // quadrilateral face
	    FaceQuadFormula = TFEDatabase3D::GetQFQuadFromDegree(2*fe_degree);
	    break;
    default:
      ErrThrow("wrong number of face vertices ", nFaceVertices);
	  }
	  int N_Points;
	  const double* faceWeights;
	  const double *t,*s;
	  // get a quadrature formula good enough for the velocity FE space
	  TQuadFormula2D *qf2 = TFEDatabase3D::GetQuadFormula2D(FaceQuadFormula);
	  qf2->GetFormulaData(N_Points, faceWeights, t, s);
  
	  // ====================================
	  // generate data on reference mesh cell for the 2d face of 3d cell
	  TFEDatabase3D::GetBaseFunct3DFromFE3D(FEId)
	    ->MakeRefElementData(FaceQuadFormula);
                        
	  BaseFunct3D *BaseFuncts = TFEDatabase3D::GetBaseFunct3D_IDFromFE3D();
	  int* N_BaseFunct = TFEDatabase3D::GetN_BaseFunctFromFE3D();
                        
                        
	  // values of base functions in all quadrature points on face
	  double **JointValues = TFEDatabase3D::GetJointValues3D
	    (BaseFuncts[FEId], FaceQuadFormula, joint_id);
                        
	  TFEDatabase3D::GetBaseFunct3D(BaseFuncts[FEId])->ChangeBF(FESpace3D->GetCollection(),
								    cell, N_Points, JointValues);
                        
  
	  // compute normal vector
	  std::vector<double> normal;
	  double transformationDeterminant;

	  cell->computeNormalAndTransformationData(joint_id,
						   normal,
						   transformationDeterminant);
	  
	  // note: the normal computed above is not always directed outward (for boundary cells)
	  normal.resize(3);
	  double x,y,z;
	  ///@attention we assume that the bound.face is planar
	  boundface->GetXYZofTS(t[0], s[0], x, y, z);
	  boundface->get_normal_vector(x,y,z,normal[0],normal[1],normal[2]);
	  Output::print<4>(" ** computed normal vector on (", x ,",", y , "," , z,
			   ") => n = (",normal[0],",",normal[1],",",normal[2],")");
	    
	  // compute \int_F u.n = sum_{gauss pt} w_k \sum_j u_j.n phi_j(x_k)
	  double value = 0;
	  for(int l=0; l < N_Points; l++) {
	    double u_n = 0;
	    
	    // compute u.n on l-th Gauss point
	    for(int k=0; k<N_BaseFunct[FEId]; k++) {
	      int global_dof_from_local = DOF[k];
	      for(size_t icoor=0; icoor<3; icoor++) {
		double *u_icoor_values = this->GetComponent(icoor)->GetValues();
		double u_icoor_on_x_k = u_icoor_values[ global_dof_from_local ];
		u_n += JointValues[l][k]* (u_icoor_on_x_k*normal[icoor]) ;
	      }
	    }
	      
	    value += faceWeights[l] * transformationDeterminant * u_n;
	  }
	  
	  flux += value;
	  
	}
	
      }
    }
  }

}

/** calculate errors to given vector function */
void TFEVectFunct3D::GetDeformationTensorErrors( 
  DoubleFunct3D *Exact, DoubleFunct3D *Exact1,
  DoubleFunct3D *Exact2,
  int N_Derivatives,
  MultiIndex3D *NeededDerivatives,
  int N_Errors, TFEFunction3D::ErrorMethod *ErrorMeth, 
  const CoeffFct3D& Coeff, 
  TAuxParam3D *Aux,
  int n_fespaces, TFESpace3D **fespaces,
  double *errors)
{
  int i,j,k,l, N_LocalUsedElements;
  int N_Cells, N_Points, N_Parameters, N_;
  int Used[N_FEs3D], *N_BaseFunct;
  TFESpace3D *fespace;
  FE3D LocalUsedElements[N_FEs3D], CurrentElement;
  BaseFunct3D BaseFunct, *BaseFuncts;
  const double *weights, *xi, *eta, *zeta;
  double X[MaxN_QuadPoints_3D], Y[MaxN_QuadPoints_3D], Z[MaxN_QuadPoints_3D];
  double AbsDetjk[MaxN_QuadPoints_3D];
  double *Param[MaxN_QuadPoints_3D], *aux, *aux1, *aux2, *aux3;
  double *Derivatives[3*MaxN_QuadPoints_3D];
  double *ExactVal[3*MaxN_QuadPoints_3D];
  double *AuxArray[MaxN_QuadPoints_3D];
  int *DOF;
  double **OrigFEValues, *Orig, value, value1, value2;
  double FEFunctValues[MaxN_BaseFunctions3D];
  double FEFunctValues1[MaxN_BaseFunctions3D];
  double FEFunctValues2[MaxN_BaseFunctions3D];
  int *GlobalNumbers, *BeginIndex;
  double LocError[4], *Values0,*Values1, *Values2;
  double hK;
  bool *SecondDer;

  BaseFuncts = TFEDatabase3D::GetBaseFunct3D_IDFromFE3D();
  N_BaseFunct = TFEDatabase3D::GetN_BaseFunctFromFE3D();

  SecondDer = new bool[n_fespaces];
  for(i=0;i<n_fespaces;i++)
    SecondDer[i] = false;

  N_Parameters = Aux->GetN_Parameters();
  aux1 = new double [MaxN_QuadPoints_3D*N_Parameters];
  for(j=0;j<MaxN_QuadPoints_3D;j++)
    Param[j] = aux1 + j*N_Parameters;

  aux2 = new double [3*MaxN_QuadPoints_3D*N_Derivatives];
  for(j=0;j<3*MaxN_QuadPoints_3D;j++)
    Derivatives[j] = aux2 + j*N_Derivatives;
  
  aux3 = new double [3*MaxN_QuadPoints_3D * 4];
  for(j=0;j<3*MaxN_QuadPoints_3D;j++)
    ExactVal[j] = aux3 + j*4;

  // 20 <= number of term
  aux = new double [MaxN_QuadPoints_3D*20]; 
  for(j=0;j<MaxN_QuadPoints_3D;j++)
    AuxArray[j] = aux + j*20;

  fespace = fespaces[0];
  GlobalNumbers = fespace->GetGlobalNumbers();
  BeginIndex = fespace->GetBeginIndex();
  Values0 = Values;
  Values1 = Values+Length;
  Values2 = Values1+Length;

  for(i=0;i<N_Errors;i++)
    errors[i] = 0.0;

// ########################################################################
// loop over all cells
// ########################################################################
  auto Coll = fespaces[0]->GetCollection(); // all spaces use same Coll
  N_Cells = Coll->GetN_Cells();
  for(i=0;i<N_Cells;i++)
  {
    auto cell = Coll->GetCell(i);
    hK = cell->GetDiameter();

    // ####################################################################
    // find local used elements on this cell
    // ####################################################################
    memset(Used, 0, N_FEs3D*SizeOfInt);
    for(j=0;j<n_fespaces;j++)
    {
      CurrentElement = fespaces[j]->GetFE3D(i, cell);
      Used[CurrentElement] = 1;
    }

    N_LocalUsedElements = 0;
    memset(LocalUsedElements, 0, SizeOfInt*N_FEs3D);
    j = 0;
    for(k=0;k<N_FEs3D;k++)
      if(Used[k])
      {
        LocalUsedElements[j] = (FE3D)k;
        j++;
      }
    N_LocalUsedElements = j;

    // ####################################################################
    // calculate values on original element
    // ####################################################################
    TFEDatabase3D::GetOrig(
      N_LocalUsedElements, LocalUsedElements, 
      Coll, cell, SecondDer,
      N_Points, xi, eta, zeta, weights, X, Y, Z, AbsDetjk);

    if(N_Parameters>0)
      Aux->GetParameters(N_Points, cell, i, xi, eta, zeta, X, Y, Z, Param); 

    // calculate all needed derivatives of this FE function
    CurrentElement = fespace->GetFE3D(i, cell);
    BaseFunct = BaseFuncts[CurrentElement];
    N_ = N_BaseFunct[CurrentElement];

    DOF = GlobalNumbers + BeginIndex[i];
    for(l=0;l<N_;l++)
    {
      FEFunctValues[l] = Values0[DOF[l]];
      FEFunctValues1[l] = Values1[DOF[l]];
      FEFunctValues2[l] = Values2[DOF[l]];
    }

    // for all needed derivatives
    for(k=0;k<N_Derivatives;k++)
    {
      OrigFEValues = TFEDatabase3D::GetOrigElementValues(BaseFunct,
                                                         NeededDerivatives[k]);
      // for all quadrature points
      for(j=0;j<N_Points;j++)
      {
        Orig = OrigFEValues[j];
        value = 0;
        value1 = 0;
        value2 = 0;
        for(l=0;l<N_;l++)
        {
          value += FEFunctValues[l] * Orig[l];
          value1 += FEFunctValues1[l] * Orig[l];
          value2 += FEFunctValues2[l] * Orig[l];
        } // endfor l
        Derivatives[j][k] = value;
        Derivatives[j+N_Points][k] = value1;
        Derivatives[j+2*N_Points][k] = value2;
      } // endfor j
    } // endfor k

    // exact value for first component
    for(j=0;j<N_Points;j++)
      Exact(X[j], Y[j], Z[j], ExactVal[j]);

    // exact value for second component
    for(j=0;j<N_Points;j++)
      Exact1(X[j], Y[j], Z[j], ExactVal[j+N_Points]);

    // exact value for third component
    for(j=0;j<N_Points;j++)
      Exact2(X[j], Y[j], Z[j], ExactVal[j+2*N_Points]);

    if(Coeff)
      Coeff(N_Points, X, Y, Z, Param, AuxArray);      

    ErrorMeth(N_Points, {{X, Y, Z}}, AbsDetjk, weights, hK, Derivatives, 
              ExactVal, AuxArray, LocError);

    for(j=0;j<N_Errors;j++)
      errors[j] += LocError[j];

  } // endfor i

  for(j=0;j<N_Errors;j++)
  {
    if (errors[j]>0)
      errors[j] = sqrt(errors[j]);
  }

  delete aux;
  delete aux1;
  delete aux2;
  delete aux3;
  delete SecondDer;
} // TFEFunction3D::GetDeformationTensorErrors


double TFEVectFunct3D::GetL2NormDivergenceError(DoubleFunct3D* Exact_u1,
                                                DoubleFunct3D* Exact_u2,
                                                DoubleFunct3D* Exact_u3)
{
  auto Coll = FESpace3D->GetCollection();
  int N_Cells = Coll->GetN_Cells();

  double FEFunctValues0[MaxN_BaseFunctions3D];
  double FEFunctValues1[MaxN_BaseFunctions3D];
  double FEFunctValues2[MaxN_BaseFunctions3D];
  double X[MaxN_QuadPoints_3D], Y[MaxN_QuadPoints_3D], Z[MaxN_QuadPoints_3D];
  double AbsDetJK[MaxN_QuadPoints_3D];
  bool SecondDer = false;

  double div_error = 0.;
  // loop over all cells
  for(int i = 0; i < N_Cells; i++)
  {
    TBaseCell *cell = Coll->GetCell(i);
    FE3D FEid = FESpace3D->GetFE3D(i, cell);
    TFE3D * fe = TFEDatabase3D::GetFE3D(FEid);
    fe->GetBaseFunct3D_ID();

    // compute transformation to reference cell
    const double *xi, *eta, *zeta, *weights;
    int N_Points;
    TFEDatabase3D::GetOrig(1, &FEid, Coll, cell, &SecondDer,
        N_Points, xi, eta, zeta, weights, X, Y, Z, AbsDetJK);

    // calculate all needed derivatives of this FE function
    int N_Bf = fe->GetN_DOF();
    int *DOF = FESpace3D->GetGlobalDOF(i);
    for(int jj = 0; jj < N_Bf; jj++)
    {
      int k = DOF[jj];
      FEFunctValues0[jj] = Values[k];
      FEFunctValues1[jj] = Values[k+Length];
      FEFunctValues2[jj] = Values[k+2*Length];
    }
    BaseFunct3D BaseFunct = fe->GetBaseFunct3D_ID();
    auto OrigFEValuesX = TFEDatabase3D::GetOrigElementValues(BaseFunct, D100);
    auto OrigFEValuesY = TFEDatabase3D::GetOrigElementValues(BaseFunct, D010);
    auto OrigFEValuesZ = TFEDatabase3D::GetOrigElementValues(BaseFunct, D001);

    // loop over all quadrature points
    for(int j = 0; j < N_Points; j++)
    {
      double local_divergence_fe = 0;
      for(int l = 0; l < N_Bf; l++)
      {
        local_divergence_fe += FEFunctValues0[l] * OrigFEValuesX[j][l];
        local_divergence_fe += FEFunctValues1[l] * OrigFEValuesY[j][l];
        local_divergence_fe += FEFunctValues2[l] * OrigFEValuesZ[j][l];
      }
      double local_divergence_exact = 0;
      double exact_val[5];
      Exact_u1(X[j], Y[j], Z[j], exact_val);
      local_divergence_exact += exact_val[1];
      Exact_u2(X[j], Y[j], Z[j], exact_val);
      local_divergence_exact += exact_val[2];
      Exact_u3(X[j], Y[j], Z[j], exact_val);
      local_divergence_exact += exact_val[3];
      
      auto local_div_error = local_divergence_fe - local_divergence_exact;
      div_error += AbsDetJK[j] * weights[j] * local_div_error * local_div_error;
    } // endfor j
  } // endfor i
  div_error = std::sqrt(div_error);
  return div_error;
} // TFEVectFunct3D::GetL2NormDivergenceError


void TFEVectFunct3D::FindValueLocal(const TBaseCell* cell, int cell_no, 
				    double x, double y, double z, 
				    double* values) const
{
 this->TFEFunction3D::FindValueLocal(cell, cell_no, x, y, z, values);
 auto u2 = this->GetComponent(1);
 u2->FindValueLocal(cell, cell_no, x, y, z, values+1);
 auto u3 = this->GetComponent(2);
 u3->FindValueLocal(cell, cell_no, x, y, z, values+2);
 delete u2;
 delete u3;
}


/** write the solution into a data file - written by Sashi **/
void TFEVectFunct3D::WriteSol(double t, const std::string& directory,
                              const std::string& basename)
{
  int i, N_Joints, N_Cells;
  char Dquot;

  #ifdef _MPI
  int rank;
  MPI_Comm_rank(TDatabase::ParamDB->Comm, &rank);
  #endif

  Dquot = 34; //  see ASCII Chart
  auto Coll = FESpace3D->GetCollection();
  N_Cells = Coll->GetN_Cells();

  i=0;
  auto cell =  Coll->GetCell(i);
  N_Joints = cell->GetN_Joints();

  std::ostringstream os;
  os << " ";

  #ifdef _MPI
  OutPut("Writing solution into "<< directory << "/" << basename << rank
         << ".Sol MooNMD file"<< endl);
  os.seekp(std::ios::beg);
  os << directory << "/" << basename<<rank<<".Sol" << ends;
  #else
  OutPut("Writing solution into "<< directory << "/" << basename << t
         << ".Sol MooNMD file"<< endl);
  os.seekp(std::ios::beg);
  os << directory << "/" << basename << t<<".Sol" << ends;
  #endif  
  
  std::ofstream dat(os.str().c_str());

  if (!dat)
   {
    cerr << "cannot open file for output" << endl;
    exit(0);
   }

    dat << "# Solution of the vector "<<Dquot<<Name<<Dquot<<", written by MooNMD"<< endl;
    dat << "# N_Cells, Cell_Type, N_Dim, N_Dof" <<  endl;
    dat <<N_Cells << " " << N_Joints << " " << N_Components << " " << Length<< endl;
    dat <<  endl;

    dat << "# Dof "<< " Nodal Values"<< endl;

    for(i=0;i<Length;i++)
     dat << i << " " << Values[i] << " " << Values[Length + i] <<" " << Values[2*Length + i] << endl;
  dat.close();  
  
  
}//WriteSol




/** Read the solution from a given data file - written by Sashi **/
void TFEVectFunct3D::ReadSol(const std::string& BaseName)
{
 int i, j, rank, N_Joints, N_Cells, N_cells, N_joints, N_components, length;
 char line[100];

#ifdef _MPI 
   MPI_Comm_rank(TDatabase::ParamDB->Comm, &rank);
#else
   rank = 0;
#endif

  auto Coll = FESpace3D->GetCollection();
  N_Cells = Coll->GetN_Cells();

  i=0;
  auto cell =  Coll->GetCell(i);
  N_Joints = cell->GetN_Joints();

  std::ifstream dat(BaseName);
  if (!dat)
   {
    cerr << "cannot open '" <<  BaseName << "' for input" << endl;
#ifdef _MPI
    MPI_Finalize();
#endif
    exit(0);
   }
   
  dat.getline (line, 99);
  dat.getline (line, 99);
  dat >> N_cells >> N_joints >> N_components >> length;
  dat.getline (line, 99);

  if(N_cells!=N_Cells || N_joints!=N_Joints || N_components!=N_Components || length!=Length )
   {
//     printf("Given data file does not match with this FE Vector function !\n" );
    printf("Rank %d, N_cells %d, N_joints %d,  N_components %d, length %d\n",rank, N_cells, N_joints, N_components,length);
   printf("Rank %d, Needed N_cells %d, N_joints %d,  N_components %d, length %d\n",rank, N_Cells, N_Joints, N_Components,Length);
//     OutPut(N_Cells <<", "<< N_Joints<< ", " << N_Components<< " , "<< Length<< " , " <<endl);
#ifdef _MPI
    MPI_Finalize();
#endif
    exit(0);
   }

  dat.getline (line, 99);
  
#ifdef _MPI
  printf("Reading nodal values of the FE Vector function:  \n");  
#else
  OutPut("Reading nodal values of the FE Vector function !"<<endl);
#endif

  for(i=0;i<Length;i++)
   {
    dat.getline (line, 99);
    dat >> j >> Values[i] >> Values[Length + i] >> Values[2*Length + i];
   }

  dat.close();   
  
} // ReadSol


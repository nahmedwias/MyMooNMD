#include <Domain.h>
#include <FEDatabase2D.h>
#include <iostream>
#include <Database.h>
#include <FESpace2D.h>
#include <MainUtilities.h>
#include <FEMatrix.h>
#include <BlockMatrix.h>
#include <PrRobustNSE2D.h>

#include <BlockFEMatrix.h>
#include <list>
#include <sys/stat.h>
#include <sys/types.h>

int main(int argc, char* argv[])
{
  TDatabase Database;
  TFEDatabase2D FEDatabase;
  TDomain domain;
  
  unsigned int nRefinements = 2;
  
  Output::set_outfile(TDatabase::ParamDB->OUTFILE);
  
  // possibly change parameters in the database, if they are not meaningful now
  Database.CheckParameterConsistencyNSE();
  // write all Parameters to the OUTFILE (not to console) for later reference
  Database.WriteParamDB(argv[0]);
  
   domain.Init((char*)"Default_UnitSquare", (char*)"UnitSquare");
  // domain.Init((char*)"Default_UnitSquare", (char*)"TwoTriangles");
  
  for(int i = 0; i < nRefinements; i++)
  {
    domain.RegRefineAll();
  }
  
  TDatabase::ParamDB->EXAMPLE = 41;
  TDatabase::ParamDB->VELOCITY_SPACE = 2;
  TDatabase::ParamDB->PRESSURE_SPACE = 1;
  TDatabase::ParamDB->PROJECTION_SPACE = 1012;
  TDatabase::ParamDB->NSTYPE = 4;
  TDatabase::ParamDB->SOLVER_TYPE = 2;
  TDatabase::ParamDB->PROBLEM_TYPE = 3;
  TDatabase::ParamDB->LAPLACETYPE = 0;
  TDatabase::ParamDB->NSE_NONLINEAR_FORM = 0;
  TDatabase::ParamDB->MEASURE_ERRORS= 1;
  TDatabase::ParamDB->WRITE_VTK = 1;
  TDatabase::ParamDB->USE_ISOPARAMETRIC = 0;
  TDatabase::ParamDB->RE_NR = 1e5;
  TDatabase::ParamDB->DISCTYPE = 300;
  
  //===================================================
  Example_NSE2D example;  
  PrRobustNSE2D prRoubustFEM(domain, example);  
  prRoubustFEM.assembleMatrixRhs();
  prRoubustFEM.solve();
  //===================================================
  double errors[5];
  TAuxParam2D aux;
  MultiIndex2D allDerivs[3] = {D00, D10, D01};
  const TFESpace2D *v_space = &prRoubustFEM.get_velocity_space();
  const TFESpace2D *p_space = &prRoubustFEM.get_pressure_space();
  TFEFunction2D *u1 = prRoubustFEM.get_velocity_component(0);
  TFEFunction2D *u2 = prRoubustFEM.get_velocity_component(1);
  TFEFunction2D &p = prRoubustFEM.get_pressure();
  
  u1->GetErrors(example.get_exact(0), 3, allDerivs, 2, L2H1Errors, 
                  nullptr, &aux, 1, &v_space, errors);
  
  u2->GetErrors(example.get_exact(1), 3, allDerivs, 2, L2H1Errors, 
                  nullptr, &aux, 1, &v_space, errors+2);

  double l2 = sqrt(errors[0]*errors[0] + errors[2]*errors[2]);
  double h1 = sqrt(errors[1]*errors[1] + errors[3]*errors[3]);
  
  p.GetErrors(example.get_exact(2), 3, allDerivs, 2, L2H1Errors, 
                  nullptr, &aux, 1, &p_space, errors);
  double l2p = errors[0];
  double h1p = errors[1];
  
  if(fabs(l2 - 0.00017658) > 1e-6)
  {
    ErrThrow("Program 1: L2 error of velocity is not correct.");
  }
  if(fabs(h1-0.00470131) > 1e-6)
  {
    ErrThrow("Program 1: H1 error of velocity is not correct.");
  }
  
  if(fabs(l2p - 0.00272494) > 1e-6)
  {
    ErrThrow("Program 1: L2 error of pressure is  not correct.");
  }
  
  if(fabs(h1p - 0.0705915) > 1e-6)
  {
    ErrThrow("Program 1: H1 error of pressure is  not correct.");
  }
  
  cout<<"test passed: " << endl;
  return 0;
}

// /*
// 
// #include <Domain.h>
// #include <Database.h>
// #include <FEDatabase2D.h>
// #include <FESpace2D.h>
// #include <MainUtilities.h>
// #include <FEMatrix.h>
// #include <BlockMatrix.h>
// #include <list>
// 
// void BC(int BdComp, double t, BoundCond &cond)
// { cond = NEUMANN; }
// 
// int main(int argc, char* argv[])
// {
//   TDatabase Database;
//   TFEDatabase2D FEDatabase;
//   TDomain domain;
//   
//   unsigned int nRefinements = 0;
//   
//   // domain.Init((char*)"Default_UnitSquare", (char*)"TwoTriangles");
//   domain.Init((char*)"Default_UnitSquare", (char*)"UnitSquare");
//   
//   for(int i = 0; i < nRefinements; i++)
//   {
//     domain.RegRefineAll();
//   }
//   TCollection *coll = domain.GetCollection(It_Finest, 0);
//   
//   int vel_ord, proj_ord;
//   TDatabase::ParamDB->VELOCITY_SPACE = 1;
//   switch(TDatabase::ParamDB->VELOCITY_SPACE)
//   {
//     case 1:
//       vel_ord = 2;
//       proj_ord = 1012;
//       break;
//     case 2:
//       vel_ord = 2;
//       proj_ord = 1001;
//       break;
//     case 3:
//       vel_ord = 3;
//       proj_ord = 1002;
//       break;
//     case 4:
//       vel_ord = 4;
//       proj_ord = 1003;
//       break;
//     default:
//       ErrMsg("velocity space " << TDatabase::ParamDB->VELOCITY_SPACE << "is not allowed:");
//       exit(0);
//   }
//   
//   /**
//    * \Pi_h: V_h \rightarrow X_h(RT or BDM) 
//    */
//   
//   // test space V_h
//   TFESpace2D *V_h = new TFESpace2D(coll, (char*) "v_h", (char*)"v_h", BC,
//           vel_ord, nullptr);
//   // anzats space X_h (vector valued)
//   TFESpace2D *X_h = new TFESpace2D(coll, (char*) "x_h", (char*)"x_h", 
//                                    BoundConditionNoBoundCondition, proj_ord, nullptr);
//   
//   std::shared_ptr<TMatrix> M1(new FEMatrix(V_h, X_h));
//   std::shared_ptr<TMatrix> M2(new FEMatrix(V_h, X_h));
//   BlockMatrix A(2,1, {M1, M2});
//   // cout<<V_h->GetN_DegreesOfFreedom()<<"  " << X_h->GetN_DegreesOfFreedom()<<endl;
//   
//   
//   unsigned int nCells = coll->GetN_Cells();
//   for(unsigned int i = 0; i<nCells; ++i)
//     coll->GetCell(i)->SetNormalOrientation();
//   
//   
//   TBaseCell *cell;
//   
//   // declartion for Vh
//   TFE2D *elementVh;
//   TBaseFunct2D *bfVh;
//   int nBFVh;
//   // declartion for Xh
//   TFE2D *elementXh;
//   TBaseFunct2D *bfXh;
//   TNodalFunctional2D *nfXh;
//   int nBFXh, baseVectDim;
//   double *xi, *eta;
//   int nPoints;
//   // cout<<nCells<<"  " << nPoints<<endl;  
//   for(unsigned int i = 0; i<nCells; ++i)
//   {
//     cell = coll->GetCell(i);
//     
//     elementVh = TFEDatabase2D::GetFE2D(V_h->GetFE2D(i,cell));
//     bfVh = elementVh->GetBaseFunct2D(); // basis function
//     int nDofVh = bfVh->GetDimension(); // number of basis 
//     
//     
//     elementXh = TFEDatabase2D::GetFE2D(X_h->GetFE2D(i, cell));
//     bfXh = elementXh->GetBaseFunct2D(); // basis function object
//     int nDofXh = bfXh->GetDimension(); // no of basis ftns in Xh
//     baseVectDim= bfXh->GetBaseVectDim();
//     
//     // nodal functional object
//     nfXh = elementXh->GetNodalFunctional2D();
//     // points on the reference element Xh
//     nfXh->GetPointsForAll(nPoints,xi, eta);
//     
//     // number of basis functions, this is the length of the array needed to 
//     // evaluate the basis functions (therefore the factor baseVectDim)    
//     int nBaseFunct = nDofVh*baseVectDim;    
//     
//      // the id of the reference transformation
//     RefTrans2D refTransID = elementXh->GetRefTransID();
//     TFEDatabase2D::SetCellForRefTrans(cell, refTransID);
//     
//     double uxorig[nPoints][nBaseFunct];
//     double AllPointValues[nDofVh];
//     
//     for(int k = 0; k < nPoints; k++)
//     {
//       bfVh->GetDerivatives(D00, xi[k], eta[k], AllPointValues);
//       
//       TFEDatabase2D::GetOrigValues(refTransID, xi[k], eta[k], bfVh, coll,
//                                    (TGridCell *) cell, AllPointValues,
//                                    nullptr, nullptr, uxorig[k], nullptr,
//                                    nullptr);
//     }
//     Output::print("cell ", i, "  verices:\n  ", cell->GetVertex(0), "\n  ", 
//                   cell->GetVertex(1), "\n  ", cell->GetVertex(2));
//     
//     for(int k = 0; k < nPoints; k++)
//     {
//       cout<<xi[k] << " " << eta[k] << " "<<endl;
//       
//       for(int l = 0; l < nDofVh; l++)
//       {
//         cout<< "a["<<k<<"]["<<l<<"] " << uxorig[k][l] << " " ;
//       }
//       cout<<endl;
//     }
//     
//     double PointValues[nPoints* baseVectDim];
//     double PointValues1[nPoints* baseVectDim];
//     memset(PointValues,0,nPoints* baseVectDim*sizeof(double));
//     memset(PointValues1,0,nPoints* baseVectDim*sizeof(double));
//     double FunctionalValues[nDofXh];
//     double FunctionalValues1[nDofXh];
//     
//     for(int j=0; j<nDofVh; ++j)
//     {
//       for(int k = 0; k < nDofXh; k++)
//       {
//         Output::print("local: Vh ", j, " Xh ", k, "\t global: Vh ", 
//                       V_h->GetGlobalDOF(i)[j], " Xh ",
//                       X_h->GetGlobalDOF(i)[k]);
//       }
//     }
//     
//     for(int j=0; j<nDofVh; ++j)
//     {
//       cout<<"dof: " << j<<endl;
//       for(int k=0; k<nPoints; ++k)
//       {
//         PointValues[k] = uxorig[k][j];
//         cout<<" point: " << k << " x( " << xi[k] << " " 
//               << eta[k] << ")  "<< PointValues[k] << "  " <<endl;;
//         PointValues1[k + nPoints] = uxorig[k][j];
//         //cout<< PointValues1[k + nPoints] << " \t" ;
//         //cout<<endl;
//       }      
//       nfXh->GetAllFunctionals(coll, cell, PointValues, FunctionalValues);
//       nfXh->GetAllFunctionals(coll, cell, PointValues1, FunctionalValues1);
//       
//       
//       for(int k = 0; k < nDofXh; k++)
//       {
//         int e = elementXh->GetFEDesc2D()->GetJointOfThisDOF(k);
//         int n = 1; // normal orientation
//         if(e != -1)
//         {
//           n = cell->GetNormalOrientation(e);
//         }
//         if( fabs(FunctionalValues[k]) < 1e-10 )
//         {
//           FunctionalValues[k] = 0;
//         }
//         Output::print<1>("xdof: " , k, "\t", FunctionalValues[k], " ", FunctionalValues1[k],
//                          "  -> ", V_h->GetGlobalDOF(i)[j], ", ", X_h->GetGlobalDOF(i)[k]);
//         M1->set(V_h->GetGlobalDOF(i)[j],X_h->GetGlobalDOF(i)[k], FunctionalValues[k]*n);
//         M2->set(V_h->GetGlobalDOF(i)[j],X_h->GetGlobalDOF(i)[k], FunctionalValues1[k]*n);
//       }
//       Output::print<1>("\n");
//     }    
//   }
//   //A.get_combined_matrix()->Print();
//   A.get_combined_matrix()->PrintFull();
// }*/

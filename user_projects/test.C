#include <Domain.h>
#include <FEDatabase2D.h>
#include <iostream>
#include <Database.h>
#include <FESpace2D.h>
#include <MainUtilities.h>
#include <FEMatrix.h>
#include <BlockMatrix.h>

#include <BlockFEMatrix.h>
#include <list>

// int main(int argc, char* argv[])
// {
//   TDatabase database;
//   TFEDatabase2D FEDatabase;
//   
//   TDomain domain;
//   domain.Init((char*) "Default_UnitSquare", (char*)"PeriodicSquares");
//   
//   TCollection *coll = domain.GetCollection(It_Finest,0);
//   
//   char *PsBaseName = TDatabase::ParamDB->BASENAME;
//   std::stringstream os;
//   os.seekp(std::ios::beg);
//   os << PsBaseName << 0 << ".ps" << ends;
//   domain.PS(os.str().c_str(),coll);
// }

void BC(int BdComp, double t, BoundCond &cond)
{ cond = DIRICHLET; }

int main(int argc, char* argv[])
{
  TDatabase Database;
  TFEDatabase2D FEDatabase;
  TDomain domain;
  
  unsigned int nRefinements = 1;
  
  domain.Init((char*)"Default_UnitSquare", (char*)"TwoTriangles");
  //domain.Init((char*)"Default_UnitSquare", (char*)"UnitSquare");
  
  for(int i = 0; i < nRefinements; i++)
  {
    domain.RegRefineAll();
  }
  TCollection *coll = domain.GetCollection(It_Finest, 0);
  
  int vel_ord, proj_ord;
  TDatabase::ParamDB->VELOCITY_SPACE = 2;
  switch(TDatabase::ParamDB->VELOCITY_SPACE)
  {
    case 1:
      vel_ord = 1;
      proj_ord = 1000;
      break;
    case 2:
      vel_ord = 2;
      proj_ord = 1001;
      break;
    case 3:
      vel_ord = 3;
      proj_ord = 1002;
      break;
    case 4:
      vel_ord = 4;
      proj_ord = 1003;
      break;
    default:
      ErrMsg("velocity space " << TDatabase::ParamDB->VELOCITY_SPACE << "is not allowed:");
      exit(0);
  }
//   
//   /**
//    * \Pi_h: V_h \rightarrow X_h(RT or BDM) 
//    */
//   
  // test space V_h
  TFESpace2D *V_h = new TFESpace2D(coll, (char*) "v_h", (char*)"v_h", BC,
          vel_ord, nullptr);
  // anzats space X_h (vector valued)
  TFESpace2D *X_h = new TFESpace2D(coll, (char*) "x_h", (char*)"x_h", 
                                   BoundConditionNoBoundCondition, proj_ord, nullptr);
  
  BlockFEMatrix bfem(V_h, X_h);
  
  return 0;
}
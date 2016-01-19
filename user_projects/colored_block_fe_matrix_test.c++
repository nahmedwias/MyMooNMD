/*
 * Unit testing of the colored block matrix. Since this requires
 * the presence of FEspaces, this starts as a usual ParMooN program,
 * which is needed to build up FESpaces which can than be used for
 * the ColoredBlockFEMatrix.
 *
 *
 *  Created on: Dec 17, 2015
 *      Author: bartsch
 */


#include <ColoredBlockFEMatrix.h>
#include <FEMatrix.h>

#include <MooNMD_Io.h>
#include <FEDatabase2D.h>
#include <Database.h>
#include <FESpace2D.h>
#include <MainUtilities.h>

#include <vector>
#include <tuple>


int main(int argc, char* argv[])
{

   //  declaration of databases
   TDatabase Database;
   TFEDatabase2D FEDatabase;

   // default construct a domain object
   TDomain domain;

   // Set Database values (this is what is usually done by the input-file)
   TDatabase::ParamDB->UNIFORM_STEPS = 1;
   TDatabase::ParamDB->LEVELS = 1;

   // the domain is initialised with default description and default
   // initial mesh
   domain.Init((char*)"Default_UnitSquare", (char*)"UnitSquare");

   // refine grid up to the coarsest level
   for(int i=0; i<TDatabase::ParamDB->UNIFORM_STEPS + TDatabase::ParamDB->LEVELS; i++)
   {
     domain.RegRefineAll();
   }

   //collection from finest cell level
   TCollection *coll = domain.GetCollection(It_EQ, 0);

   //create two fespaces2D to fiddle around with
   size_t first_ansatz_order = 2;
   size_t second_ansatz_order = 1;
   TFESpace2D first_fe_space(coll, (char*)"first_fe_space", (char*)"first_fe_space", //may act as velo space dummy
                             BoundConditionNSE, first_ansatz_order, nullptr);
   TFESpace2D second_fe_space(coll, (char*)"second_fe_space", (char*)"second_fe_space", //may act as pressuer space dummy
                              BoundCondition_FEM_FCT, second_ansatz_order, nullptr);

//TODO test all methods on a customized matrix!
//   //create four FE Matrices
//   FEMatrix fe_matrix_1(&first_fe_space);
//   FEMatrix fe_matrix_2(&first_fe_space, &second_fe_space);
//   FEMatrix fe_matrix_3(&second_fe_space, &first_fe_space);
//   FEMatrix fe_matrix_4(&second_fe_space);
//
//
//   size_t dim1 = first_fe_space.GetN_DegreesOfFreedom();
//   size_t dim2 = second_fe_space.GetN_DegreesOfFreedom();
//
//   TMatrix t_matrix_1(dim1,dim2);
//
//   //colored block fe matrix
//   ColoredBlockFEMatrix myMatrix({&first_fe_space,&second_fe_space});
//   myMatrix.check_pointer_types();
//   myMatrix.print_and_check();
//
//   myMatrix.replace_blocks(fe_matrix_1, {{0,0}}, { false } );
//
//   myMatrix.replace_blocks(fe_matrix_2, {{0,1}, {1,0}}, { false, true });
//
//   myMatrix.replace_blocks(fe_matrix_4, {{1,1}}, { false });
//
//   myMatrix.check_pointer_types();
//
//   myMatrix.print_coloring_pattern("full fe matrix",true);
//   myMatrix.print_and_check();
//
//   // try combination method of base
//   std::shared_ptr<TMatrix> combined_base = myMatrix.ColoredBlockMatrix::get_combined_matrix();
//   combined_base->PrintFull("combined base");
//
//   //try combination method of derived
//   std::shared_ptr<TMatrix> combined_derived = myMatrix.ColoredBlockFEMatrix::get_combined_matrix();
//   combined_derived->PrintFull("combined derived");
//
//   //check copying
//   ColoredBlockFEMatrix hisMatrix(myMatrix);
//   hisMatrix.print_coloring_pattern("copy constructed fe matrix",true);
//   hisMatrix.check_pointer_types();
//   hisMatrix.print_and_check();
//
//   ColoredBlockFEMatrix herMatrix = myMatrix;
//   herMatrix.print_coloring_pattern("copy assigned fe matrix",true);
//   herMatrix.check_pointer_types();
//
//   //check adding of actives
//   fe_matrix_4.setEntries(std::vector<double>(81,1.0));
//
//   herMatrix.add_scaled_matrix_to_blocks(fe_matrix_4, 1.0 ,{{1,1}}, {false});
//   herMatrix.get_combined_matrix()->PrintFull("added standard");
//   herMatrix.check_pointer_types();
//
//   herMatrix.add_scaled_actives(fe_matrix_4, - 2.0 ,{{1,1}}, {false});
//   herMatrix.get_combined_matrix()->PrintFull("subtracted active");
//   herMatrix.check_pointer_types();
//
//   //check scaling of entries
//
//   herMatrix.scale_blocks(10, {{1,1}});
//   herMatrix.get_combined_matrix()->PrintFull("scaled");
//   herMatrix.check_pointer_types();
//   herMatrix.print_and_check();
//
//   herMatrix.scale_blocks_actives(100, {{1,1}});
//   herMatrix.get_combined_matrix()->PrintFull("scaled actives");
//   herMatrix.check_pointer_types();
//   herMatrix.print_and_check();


   // try out some named constructors
   {//CD2D
   TFESpace2D third_fe_space(coll, (char*)"third_fe_space", (char*)"third_fe_space",
                                  BoundConditionNSE, 3, nullptr);

   ColoredBlockFEMatrix cd2d_blockmat = ColoredBlockFEMatrix::CD2D(third_fe_space);

   cd2d_blockmat.print_and_check("matrix for cd2d");
   cd2d_blockmat.check_pointer_types();
   }
   {//NSE2D Typ 1
     ColoredBlockFEMatrix nse2d_1_blockmat=
         ColoredBlockFEMatrix::NSE2D_Type1(first_fe_space, second_fe_space);
     nse2d_1_blockmat.print_and_check("matrix for nse2d, nstype 1");
     nse2d_1_blockmat.check_pointer_types();
     //   std::shared_ptr<TMatrix> combi = nse2d_1_blockmat.get_combined_matrix();
     //   combi->PrintFull("nstype 1 combi");
     //   std::shared_ptr<TMatrix> combi_base = nse2d_1_blockmat.ColoredBlockMatrix::get_combined_matrix();
     //   combi_base->PrintFull("nstype 1 combi_base");
   }
   {//NSE2D Typ 2
     ColoredBlockFEMatrix nse2d_2_blockmat=
         ColoredBlockFEMatrix::NSE2D_Type2(first_fe_space, second_fe_space);
     nse2d_2_blockmat.print_and_check("matrix for nse2d, nstype 2");
     nse2d_2_blockmat.check_pointer_types();
   }
   {//NSE2D Typ 3
     ColoredBlockFEMatrix nse2d_3_blockmat=
         ColoredBlockFEMatrix::NSE2D_Type3(first_fe_space, second_fe_space);
     nse2d_3_blockmat.print_and_check("matrix for nse2d, nstype 3");
     nse2d_3_blockmat.check_pointer_types();
   }
   {//NSE2D Typ 4
     ColoredBlockFEMatrix nse2d_4_blockmat=
         ColoredBlockFEMatrix::NSE2D_Type4(first_fe_space, second_fe_space);
     nse2d_4_blockmat.print_and_check("matrix for nse2d, nstype 4");
     nse2d_4_blockmat.check_pointer_types();
   }
   {//NSE2D Typ 14
     ColoredBlockFEMatrix nse2d_14_blockmat=
         ColoredBlockFEMatrix::NSE2D_Type14(first_fe_space, second_fe_space);
     nse2d_14_blockmat.print_and_check("matrix for nse2d, nstype 14");
     nse2d_14_blockmat.check_pointer_types();
   }
   {//Darcy2D
     ColoredBlockFEMatrix darcy_blockmat=
         ColoredBlockFEMatrix::Darcy2D(first_fe_space, second_fe_space);
     darcy_blockmat.print_and_check("matrix for darcy 2d");
     darcy_blockmat.check_pointer_types();
   }

}


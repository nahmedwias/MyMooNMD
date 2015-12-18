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
   //(BoundConditionNSE chosen as bc because its available in global scope and sets only diri bdry)
   size_t first_ansatz_order = 1;
   size_t second_ansatz_order = 2;
   TFESpace2D first_fe_space(coll, (char*)"first_fe_space", (char*)"first_fe_space",
                             BoundConditionNSE, first_ansatz_order, nullptr);
   TFESpace2D second_fe_space(coll, (char*)"second_fe_space", (char*)"second_fe_space",
                              BoundConditionNSE, second_ansatz_order, nullptr);

   //create four FE Matrices
   FEMatrix fe_matrix_1(&first_fe_space);
   FEMatrix fe_matrix_2(&first_fe_space, &second_fe_space);
   FEMatrix fe_matrix_3(&second_fe_space, &first_fe_space);
   FEMatrix fe_matrix_4(&second_fe_space);


   size_t dim1 = first_fe_space.GetN_DegreesOfFreedom();
   size_t dim2 = second_fe_space.GetN_DegreesOfFreedom();

   TMatrix t_matrix_1(dim1,dim2);

   //colored block fe matrix
   ColoredBlockFEMatrix myMatrix({dim1,dim2},{dim1,dim2});
   myMatrix.check_pointer_types();

   myMatrix.replace_blocks(fe_matrix_1, {{0,0}}, { false } );

   myMatrix.replace_blocks(fe_matrix_2, {{0,1}, {1,0}}, { false, true });

   myMatrix.replace_blocks(fe_matrix_4, {{1,1}}, { false });

   myMatrix.check_pointer_types();

   myMatrix.print_coloring_pattern("full fe matrix",true);

   std::shared_ptr<TMatrix> combined = myMatrix.get_combined_matrix();
   combined->PrintFull("combined");

   //check copying
   ColoredBlockFEMatrix hisMatrix(myMatrix);
   hisMatrix.print_coloring_pattern("copy constructed fe matrix",true);
   hisMatrix.check_pointer_types();

   ColoredBlockFEMatrix herMatrix = myMatrix;
   herMatrix.print_coloring_pattern("copy assigned fe matrix",true);
   herMatrix.check_pointer_types();

   //check adding of actives
   fe_matrix_4.setEntries(std::vector<double>(17,1.0));

   herMatrix.add_scaled_matrix_to_blocks(fe_matrix_4, 1.0 ,{{1,1}}, {false});
   herMatrix.get_combined_matrix()->PrintFull("added standard");
   herMatrix.check_pointer_types();

   herMatrix.add_scaled_actives(fe_matrix_4, - 2.0 ,{{1,1}}, {false});
   herMatrix.get_combined_matrix()->PrintFull("subtracted active");
   herMatrix.check_pointer_types();

   //check scaling of entries

   herMatrix.scale_blocks(10, {{1,1}});
   herMatrix.get_combined_matrix()->PrintFull("scaled");
   herMatrix.check_pointer_types();
   herMatrix.print_and_check();

   herMatrix.scale_blocks_actives(100, {{1,1}});
   herMatrix.get_combined_matrix()->PrintFull("scaled actives");
   herMatrix.check_pointer_types();
   herMatrix.print_and_check();



}


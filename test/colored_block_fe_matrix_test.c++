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
#include <BlockVector.h>

#include <vector>
#include <tuple>
#include <type_traits>


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

  Output::setVerbosity(1);

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

  // Create FeSpaces2D to fiddle around with

  size_t first_ansatz_order = 2;
  size_t second_ansatz_order = 1;
  size_t third_ansatz_order = 3;

  TFESpace2D first_fe_space(coll, (char*)"first_fe_space",
                            (char*)"first_fe_space", //may act as velo space dummy
                            BoundConditionNSE, first_ansatz_order, nullptr);

  TFESpace2D second_fe_space(coll, (char*)"second_fe_space",
                             (char*)"second_fe_space", //may act as pressure space dummy
                             BoundCondition_FEM_FCT, second_ansatz_order, nullptr);

  TFESpace2D third_fe_space(coll, (char*)"third_fe_space",
                            (char*)"third_fe_space", //yet another space
                            BoundConditionNSE, third_ansatz_order, nullptr);

  {
    //test default constructor
    ColoredBlockFEMatrix zero_matrix;
    zero_matrix.check_coloring();
    zero_matrix.check_pointer_types();
  }


  { // test standard methods with custom-made 2x2 FEMatrix, including
    // one transposed storage and one transposed-storage memory hack

    //create four FE Matrices to fiddle around with
    FEMatrix fe_matrix_1(&first_fe_space);
    FEMatrix fe_matrix_2(&first_fe_space, &second_fe_space);
    FEMatrix fe_matrix_3(&second_fe_space, &first_fe_space);
    FEMatrix fe_matrix_4(&second_fe_space);

    // one TMatrix, too (copied from one FEMatrix for the sake of a non-emtpy structure)
    TMatrix t_matrix_1(fe_matrix_3);
    t_matrix_1.setEntries(std::vector<double>(t_matrix_1.GetN_Entries(),1.0));

   //custom construct and check
   ColoredBlockFEMatrix myMatrix({&first_fe_space,&second_fe_space});
   myMatrix.check_pointer_types(); //casts to FEMatrix work?
   myMatrix.check_coloring(); //coloring is unbroken?

   //replace method and check
   myMatrix.replace_blocks(fe_matrix_1, {{0,0}}, { false } );
   myMatrix.replace_blocks(fe_matrix_4, {{1,1}}, { false });

   myMatrix.replace_blocks(fe_matrix_3, {{0,1}, {1,0}}, { true, false }); //replace which leads to color merge
   myMatrix.replace_blocks(fe_matrix_2, {{0,1}}, {false} ); //replace which leads to color split

   myMatrix.check_pointer_types(); //casts to FEMatrix work?
   myMatrix.check_coloring(); //coloring is unbroken?

   myMatrix.replace_blocks(fe_matrix_3, {{0,1}, {1,0}}, { true, false }); //replace which leads to color merge
   myMatrix.check_pointer_types(); //casts to FEMatrix work?
   myMatrix.check_coloring(); //coloring is unbroken?

   // try combination method of base
   {
     std::shared_ptr<TMatrix> combined_base
     = myMatrix.ColoredBlockMatrix::get_combined_matrix();
//     combined_base->PrintFull("combined base");
     if(combined_base->GetNorm(-2) != 0)
     {
       ErrThrow("Wrong norm of base-class combined matrix!");
     }
   }

   //try combination method of derived
   {
     std::shared_ptr<TMatrix> combined_derived
     = myMatrix.ColoredBlockFEMatrix::get_combined_matrix();
//     combined_derived->PrintFull("combined derived");
     if( fabs(combined_derived->GetNorm(-2) - 2.82843) > 1e-5)
     {
       ErrThrow("Wrong norm of derived-class combined matrix!");
     }
   }

   // test the space getter methods
   if(&myMatrix.get_ansatz_space(0,0) != &first_fe_space )
   {
     ErrThrow("get_ansatz_space not working correctly.")
   }
   if(&myMatrix.get_test_space(0,1) != &first_fe_space )
   {
     ErrThrow("get_ansatz_space not working correctly.")
   }
   if(&myMatrix.get_row_space(1) != &second_fe_space )
   {
     ErrThrow("get_row_space not working correctly.")
   }
   if(&myMatrix.get_column_space(1) != &second_fe_space )
   {
     ErrThrow("get_column_space not working correctly.")
   }

   // check standard adding (no color split intended), here with a simple TMatrix
   myMatrix.add_matrix(t_matrix_1, 1.0 ,{{1,0}, {0,1}}, {false, true});
   myMatrix.check_pointer_types(); //casts to FEMatrix work?
   myMatrix.check_coloring(); //coloring is unbroken?
   if( fabs( myMatrix.get_combined_matrix()->GetNorm(-2) - 6.9282 ) > 1e-5)
   {//check norm of dirichlet-row corrected matrix
     ErrThrow("Wrong matrix norm after adding of a TMatrix!")
   }
   if( fabs( myMatrix.ColoredBlockMatrix::get_combined_matrix()->GetNorm(-2) - 8.48528 ) > 1e-5 )
   {//check norm of rough, algebraic matrix
     ErrThrow("Wrong matrix norm after adding of a TMatrix!")
   }

   //check actives adding (color split intended) with an FEMatrix
   fe_matrix_1.setEntries(std::vector<double>(17, 2.0));
   myMatrix.add_matrix_actives(fe_matrix_1, -2.0 ,{{0,0}}, { false });
   myMatrix.check_pointer_types(); //casts to FEMatrix work?
   myMatrix.check_coloring(); //coloring is unbroken?

   if( fabs( myMatrix.get_combined_matrix()->GetNorm(-2) - 13.8564 ) > 1e-4 )
   {//check norm of dirichlet-row corrected matrix
     ErrThrow("Wrong matrix norm after actives adding!")
   }
   if( fabs( myMatrix.ColoredBlockMatrix::get_combined_matrix()->GetNorm(-2) - 14.6969 ) > 1e-4 )
   {//check norm of rough, algebraic matrix
     ErrThrow("Wrong matrix norm after actives adding!")
   }


   //check scaling of entries
   myMatrix.scale_blocks(0.5, {{0,0},{1,1}});
   myMatrix.check_pointer_types(); //casts to FEMatrix work?
   myMatrix.check_coloring(); //coloring is unbroken?

   if( fabs( myMatrix.get_combined_matrix()->GetNorm(-2) - 9.16515 ) > 1e-5 )
   {//check norm of dirichlet-row corrected matrix
     ErrThrow("Wrong matrix norm after block scaling!")
   }
   if( fabs( myMatrix.ColoredBlockMatrix::get_combined_matrix()->GetNorm(-2) - 10.3923 ) > 1e-5 )
   {//check norm of rough, algebraic matrix
     ErrThrow("Wrong matrix norm after block scaling!")
   }

   //check scaling of active entries
   myMatrix.replace_blocks(fe_matrix_1, {{0,0}}, { false });//give us a new block in {0,0}
   myMatrix.scale_blocks_actives(-2, {{0,0}});
   myMatrix.check_pointer_types(); //casts to FEMatrix work?
   myMatrix.check_coloring(); //coloring is unbroken?

   if( fabs( myMatrix.get_combined_matrix()->GetNorm(-2) - 13.8564 ) > 1e-4 )
   {//check norm of dirichlet-row corrected matrix
     ErrThrow("Wrong matrix norm after block scaling!")
   }
   if( fabs( myMatrix.ColoredBlockMatrix::get_combined_matrix()->GetNorm(-2) - 15.748 ) > 1e-4 )
   {//check norm of rough, algebraic matrix
     ErrThrow("Wrong matrix norm after block scaling!")
   }


  }

  { //make a default NSE Matrix and two vectors
    ColoredBlockFEMatrix blockmat=
            ColoredBlockFEMatrix::NSE2D_Type1(first_fe_space, second_fe_space);

    //check matrix-vector multiplication (incl. actives)
    BlockVector preimage_act(blockmat, false);
    BlockVector image_act(blockmat,true);

    for(size_t i = 0 ; i< preimage_act.length() ; ++i)
    { //fill one vector with with ones
      preimage_act.at(i) = 1;
    }

    blockmat.apply(preimage_act, image_act);

    if(image_act.norm() != 4)
    {
      ErrThrow("Norm of BlockVector from multiplication is not correct!");
    }

    //check usage in std::vector
    std::vector<ColoredBlockFEMatrix> myMatrices;
    myMatrices.push_back(blockmat);
    myMatrices.push_back(blockmat);
    myMatrices.push_back(blockmat); //three pushes
    myMatrices.pop_back(); //one pop
    myMatrices.at(0).check_pointer_types(); //random access and check element
    myMatrices.at(0).check_coloring();


    //check copying

    //copy construction
    ColoredBlockFEMatrix hisMatrix(blockmat);
    //hisMatrix.print_coloring_pattern("copy constructed fe matrix",true);
    hisMatrix.check_pointer_types();
    hisMatrix.check_coloring();

    //copy assignment
    ColoredBlockFEMatrix herMatrix({&first_fe_space});
    herMatrix = blockmat;
    //herMatrix.print_coloring_pattern("copy assigned fe matrix",true);
    herMatrix.check_pointer_types();
    herMatrix.check_coloring();

    //check moving

    Output::setVerbosity(5);
    Output::print("Test moving.");

    //move constructor
    //ColoredBlockFEMatrix moveConstructedMatrix(ColoredBlockFEMatrix::NSE2D_Type14(first_fe_space, second_fe_space));
    ColoredBlockFEMatrix moveConstructedMatrix(std::move(ColoredBlockFEMatrix({&first_fe_space, &second_fe_space})));
    moveConstructedMatrix.check_pointer_types();
    moveConstructedMatrix.check_coloring();

    ColoredBlockFEMatrix moveAssignedMatrix({&first_fe_space});
    moveAssignedMatrix = ColoredBlockFEMatrix::Darcy2D(first_fe_space, second_fe_space);
    moveAssignedMatrix.check_pointer_types();
    moveAssignedMatrix.check_coloring();

    Output::setVerbosity(0);
  }


  // try out some named constructors, plus the block getter
  // methods for solvers and assemblers
  {//CD2D
    ColoredBlockFEMatrix blockmat = ColoredBlockFEMatrix::CD2D(third_fe_space);
    //blockmat.print_and_check("matrix for cd2d");
    blockmat.check_pointer_types(); //casts to FEMatrix work?
    blockmat.check_coloring(); //coloring is unbroken?
    // assemble blocks getter
    std::vector<std::shared_ptr<FEMatrix>> blocks_for_assembler
    = blockmat.get_blocks_uniquely();
    if (blocks_for_assembler.size() != 1)
    {
      ErrThrow("Incorrect blocks_for_assembler.size() !");
    }
    //solver blocks getter
    std::vector<std::shared_ptr<const FEMatrix>> blocks_for_solver
    = blockmat.get_blocks();
    if (blocks_for_assembler.size() != 1)
    {
      ErrThrow("Incorrect blocks_for_solver.size() !");
    }

  }
   {//NSE2D Typ 1
     ColoredBlockFEMatrix blockmat=
         ColoredBlockFEMatrix::NSE2D_Type1(first_fe_space, second_fe_space);
     //blockmat.print_and_check("matrix for nse2d, nstype 1");
     blockmat.check_pointer_types(); //casts to FEMatrix work?
     blockmat.check_coloring(); //coloring is unbroken?
     // assemble blocks getter
     std::vector<std::shared_ptr<FEMatrix>> blocks_for_assembler
     = blockmat.get_blocks_uniquely();
     if (blocks_for_assembler.size() != 3)
     {
       ErrThrow("Incorrect blocks_for_assembler.size() !");
     }
     //solver blocks getter
     std::vector<std::shared_ptr<const FEMatrix>> blocks_for_solver
     = blockmat.get_blocks();
     if (blocks_for_solver.size() != 9)
     {
       ErrThrow("Incorrect blocks_for_solver.size() !");
     }
   }
   {//NSE2D Typ 2
     ColoredBlockFEMatrix blockmat=
         ColoredBlockFEMatrix::NSE2D_Type2(first_fe_space, second_fe_space);
     //blockmat.print_and_check("matrix for nse2d, nstype 2");
     blockmat.check_pointer_types(); //casts to FEMatrix work?
     blockmat.check_coloring(); //coloring is unbroken?
     // assemble blocks getter
     std::vector<std::shared_ptr<FEMatrix>> blocks_for_assembler
     = blockmat.get_blocks_uniquely();
     if (blocks_for_assembler.size() != 5)
     {
       ErrThrow("Incorrect blocks_for_assembler.size() !");
     }
     //solver blocks getter
     std::vector<std::shared_ptr<const FEMatrix>> blocks_for_solver
     = blockmat.get_blocks();
     if (blocks_for_solver.size() != 9)
     {
       ErrThrow("Incorrect blocks_for_solver.size() !");
     }
   }
   {//NSE2D Typ 3
     ColoredBlockFEMatrix blockmat=
         ColoredBlockFEMatrix::NSE2D_Type3(first_fe_space, second_fe_space);
     //blockmat.print_and_check("matrix for nse2d, nstype 3");
     blockmat.check_pointer_types(); //casts to FEMatrix work?
     blockmat.check_coloring(); //coloring is unbroken?
     // assemble blocks getter
     std::vector<std::shared_ptr<FEMatrix>> blocks_for_assembler
     = blockmat.get_blocks_uniquely();
     if (blocks_for_assembler.size() != 6)
     {
       ErrThrow("Incorrect blocks_for_assembler.size() !");
     }
     //solver blocks getter
     std::vector<std::shared_ptr<const FEMatrix>> blocks_for_solver
     = blockmat.get_blocks();
     if (blocks_for_solver.size() != 9)
     {
       ErrThrow("Incorrect blocks_for_solver.size() !");
     }
   }
   {//NSE2D Typ 4
     ColoredBlockFEMatrix blockmat=
         ColoredBlockFEMatrix::NSE2D_Type4(first_fe_space, second_fe_space);
     //blockmat.print_and_check("matrix for nse2d, nstype 4");
     blockmat.check_pointer_types(); //casts to FEMatrix work?
     blockmat.check_coloring(); //coloring is unbroken?
     // assemble blocks getter
     std::vector<std::shared_ptr<FEMatrix>> blocks_for_assembler
     = blockmat.get_blocks_uniquely();
     if (blocks_for_assembler.size() != 8)
     {
       ErrThrow("Incorrect blocks_for_assembler.size() !");
     }
     //solver blocks getter
     std::vector<std::shared_ptr<const FEMatrix>> blocks_for_solver
     = blockmat.get_blocks();
     if (blocks_for_solver.size() != 9)
     {
       ErrThrow("Incorrect blocks_for_solver.size() !");
     }
   }
   {//NSE2D Typ 14
     ColoredBlockFEMatrix blockmat=
         ColoredBlockFEMatrix::NSE2D_Type14(first_fe_space, second_fe_space);
     //blockmat.print_and_check("matrix for nse2d, nstype 14");
     blockmat.check_pointer_types(); //casts to FEMatrix work?
     blockmat.check_coloring(); //coloring is unbroken?
     // assemble blocks getter
     std::vector<std::shared_ptr<FEMatrix>> blocks_for_assembler
     = blockmat.get_blocks_uniquely();
     if (blocks_for_assembler.size() != 9)
     {
       ErrThrow("Incorrect blocks_for_assembler.size() !");
     }
     //solver blocks getter
     std::vector<std::shared_ptr<const FEMatrix>> blocks_for_solver
     = blockmat.get_blocks();
     if (blocks_for_solver.size() != 9)
     {
       ErrThrow("Incorrect blocks_for_solver.size() !");
     }
   }
   {//Darcy2D
     ColoredBlockFEMatrix blockmat=
         ColoredBlockFEMatrix::Darcy2D(first_fe_space, second_fe_space);
     //blockmat.print_and_check("matrix for darcy 2d");
     blockmat.check_pointer_types(); //casts to FEMatrix work?
     blockmat.check_coloring(); //coloring is unbroken?
     // assemble blocks getter
     std::vector<std::shared_ptr<FEMatrix>> blocks_for_assembler
     = blockmat.get_blocks_uniquely();
     if (blocks_for_assembler.size() != 4)
     {
       ErrThrow("Incorrect blocks_for_assembler.size() !");
     }
     //solver blocks getter
     std::vector<std::shared_ptr<const FEMatrix>> blocks_for_solver
     = blockmat.get_blocks();
     if (blocks_for_solver.size() != 4)
     {
       ErrThrow("Incorrect blocks_for_solver.size() !");
     }
   }

   {//more tests on block getter methods (assemble parts getter, "true" for zero blocks)
     ColoredBlockFEMatrix blockmat= //use NSE Type 2 for the checks
         ColoredBlockFEMatrix::NSE2D_Type2(first_fe_space, second_fe_space);

     std::vector<std::vector<size_t>> positions_A = {{0,0},{0,1},{1,0},{1,1}};
     std::vector<std::vector<size_t>> positions_B = {{0,2},{1,2},{2,0},{2,1}};
     // adapt and use that if you want to check the handling of not-so clever input
     //std::vector<std::vector<size_t>> stupid_input = {{1,1}};

     // positionwise assemble getter
     std::vector<std::shared_ptr<FEMatrix>> A_blocks_for_assembler
       = blockmat.get_blocks_uniquely(positions_A);
     if(A_blocks_for_assembler.size() != 1)
     {
       ErrThrow("Incorrect A_blocks_for_assembler.size() !");
     }
     std::vector<std::shared_ptr<FEMatrix>> A_blocks_for_assembler_with_zeroes
       = blockmat.get_blocks_uniquely(positions_A, true); //with "true" for zeroes
     if(A_blocks_for_assembler_with_zeroes.size() != 2)
     {
       ErrThrow("Incorrect A_blocks_for_assembler_with_zeroes.size() !");
     }
     std::vector<std::shared_ptr<FEMatrix>> B_blocks_for_assembler
       = blockmat.get_blocks_uniquely(positions_B);
     if(B_blocks_for_assembler.size() != 4)
     {
       ErrThrow("Incorrect B_blocks_for_assembler.size() !");
     }
   }


}


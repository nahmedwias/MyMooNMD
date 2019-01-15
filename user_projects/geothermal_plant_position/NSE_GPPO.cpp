
#include "NavierStokes.h"
#include "NSE_GPPO.hpp"
#include "Multigrid.h"
#include "LocalAssembling.h"
#ifdef __2D__
#include "Assemble2D.h"
#else
#include "Assemble3D.h"
#endif




#ifdef __2D__
///************************************************************************** */
template<int d>
NSE_GPPO<d>::NSE_GPPO(const TDomain &domain,
		      const ParameterDatabase& param_db,
		      const Example_NSE& example)
  :NavierStokes<d>(domain, param_db, example),
 coefficient_function_FEspace(new FESpace(this->get_pressure_space().GetCollection(),
         "coefficient_function_FEspace", "s",  BoundConditionNoBoundCondition, 1))
{
  if (false) //param_db["variable_sigma_fct_type"])
  {
    cout <<" ***** coefficient_function_type 2 detected **** "<< endl;
    Output::print("It is assumed that a mesh and a fitting TFEFunction ReadSol() file are provided.");

    /// TCollection *coll = brinkman2d.get_pressure_space().GetCollection();
    coefficient_function_vector = BlockVector(coefficient_function_FEspace->GetN_DegreesOfFreedom());

    coefficient_function = FEFunction(coefficient_function_FEspace.get(), "coefficient_function", "coefficient_function",
            &coefficient_function_vector.at(0), coefficient_function_FEspace->GetN_DegreesOfFreedom());

    coefficient_function.ReadSol(param_db["read_coefficient_function_directory"]);
    // coefficient_function.WriteSol("/Users/blank/ParMooN/Tests/Geothermal_Plants_Position_Optimization", "fe_function_permeability_mod_out.txt");
    /////// coefficient_function.Interpolate(&read_coefficient_function);
  }
}

/*
NSE2D_GPPO::~NSE2D_GPPO()
{
  //  delete coefficient_function_FEspace;
}
 */

// ********************// START READ IN PERMEABILITY // ********************* //

/** ************************************************************************ */
void mapping_local_parameters_NSE(const double *in, double *out)
{
  // coordinates:  x at in[0], y at in[1]
  out[0] = in[2];
}

/** ************************************************************************ */
template<int d>
void NSE_GPPO<d>::assemble_with_coefficient_fct(bool read_coeff_fct)
{
  if(this->NavierStokes<d>::systems.size() > 1)
    {
      ErrThrow("NavierStokes_Adjoint::assemble_additional_terms not yet implemented for "
               "multigrid");
    }

  TFEFunction2D* coefficient_function =nullptr;
  if (read_coeff_fct) {
    coefficient_function = &(this->get_coefficient_function());
  }
  
  cout << "coefficient_function: "<< coefficient_function << endl;
  const TFESpace2D *coefficient_function_space;

  for(auto& s : this->NavierStokes<d>::systems) //OLD: (System_per_grid& s : this->NavierStokes<d>::systems)
  {
    s.rhs.reset(); //right hand side reset (TODO: is that necessary?)
    s.matrix.reset(); // reset matrix (needed for mdml where this is called)

    const FESpace * v_space = s.velocity_space.get();
    const FESpace * p_space = s.pressure_space.get();

    // declare the variables which Assemble2D needs and each nstype has to fill
    size_t N_FESpaces = 2;

    const FESpace *fespmat[3] = {v_space, p_space, nullptr};

    if (coefficient_function)
    {
      fespmat[2] = coefficient_function_space;
    }

    size_t n_sq_mat;
    SquareMatrixD *sq_matrices[5]{nullptr};//it's five pointers maximum (Type14)

    size_t n_rect_mat;
    MatrixD *rect_matrices[4]{nullptr};//it's four pointers maximum (Types 2, 4, 14)

    size_t N_Rhs = 3;
    double *RHSs[3] = {s.rhs.block(0), s.rhs.block(1), s.rhs.block(2)};
    const FESpace *fesprhs[3] = {v_space, v_space, p_space};

    BoundaryConditionFunction* boundary_conditions[3] = {
            v_space->get_boundary_condition(), v_space->get_boundary_condition(),
            p_space->get_boundary_condition() };

    std::array<BoundaryValuesFunction*, 3> non_const_bound_values;
    non_const_bound_values[0] = this->example.get_bd()[0];
    non_const_bound_values[1] = this->example.get_bd()[1];
    non_const_bound_values[2] = this->example.get_bd()[2];


    //same for all: the local assembling object
    TFEFunction2D *fe_functions[4] =
    { s.u.GetComponent(0), s.u.GetComponent(1), &s.p, NULL };

    if (coefficient_function)
    {
      coefficient_function_space = coefficient_function->GetFESpace2D(); // ->GetFESpace2D();
      fe_functions[3] = coefficient_function;
    }

    LocalAssembling<d> la(this->db, LocalAssembling_type::NSE3D_Linear,
            fe_functions, this->example.get_coeffs());

    if (coefficient_function)
    {
      // modify la such that it includes the TFEFunction2D coefficient_function
      la.setBeginParameter({0});
      la.SetN_Parameters(1);
      la.setN_ParamFct(1);
      la.setParameterFct({mapping_local_parameters_NSE});
      la.setN_FeValues(1);
      la.setFeValueFctIndex({3});
      la.setFeValueMultiIndex({D00});
    }

    std::vector<std::shared_ptr<FEMatrix>> blocks =
            s.matrix.get_blocks_uniquely();

    n_sq_mat = 5;
    n_rect_mat = 4;
    switch(TDatabase::ParamDB->NSTYPE)
    {// switch over known Block Matrix types, treat each one individually,
    // using a priori knowledge about the structure and the way it fits
    // to the LocalAssembling2D object
    // TODO remove all reinterpret_casts as soon as Assembling process takes only FEMatrices
    // we have to use reinterpret_casts because dynamic downcasting won't work here
    // FIXME replace global switch by local checking of blockmatrix type!
    case 1:
      sq_matrices[0] =  reinterpret_cast<SquareMatrixD*>(blocks.at(0).get());
      rect_matrices[0] = reinterpret_cast<MatrixD*>(blocks.at(1).get());
      rect_matrices[1] = reinterpret_cast<MatrixD*>(blocks.at(2).get());
      break;
    case 2:
      sq_matrices[0] =  reinterpret_cast<SquareMatrixD*>(blocks.at(0).get());
      rect_matrices[0] = reinterpret_cast<MatrixD*>(blocks.at(3).get()); //first the lying B blocks
      rect_matrices[1] = reinterpret_cast<MatrixD*>(blocks.at(4).get());
      rect_matrices[2] = reinterpret_cast<MatrixD*>(blocks.at(1).get()); //than the standing B blocks
      rect_matrices[3] = reinterpret_cast<MatrixD*>(blocks.at(2).get());
      break;
    case 3:
      sq_matrices[0] = reinterpret_cast<SquareMatrixD*>(blocks.at(0).get());
      sq_matrices[1] = reinterpret_cast<SquareMatrixD*>(blocks.at(1).get());
      sq_matrices[2] = reinterpret_cast<SquareMatrixD*>(blocks.at(3).get());
      sq_matrices[3] = reinterpret_cast<SquareMatrixD*>(blocks.at(4).get());
      rect_matrices[0] = reinterpret_cast<MatrixD*>(blocks.at(2).get());
      rect_matrices[1] = reinterpret_cast<MatrixD*>(blocks.at(5).get());
      break;
    case 4:
      sq_matrices[0] = reinterpret_cast<SquareMatrixD*>(blocks.at(0).get());
      sq_matrices[1] = reinterpret_cast<SquareMatrixD*>(blocks.at(1).get());
      sq_matrices[2] = reinterpret_cast<SquareMatrixD*>(blocks.at(3).get());
      sq_matrices[3] = reinterpret_cast<SquareMatrixD*>(blocks.at(4).get());
      rect_matrices[0] = reinterpret_cast<MatrixD*>(blocks.at(6).get()); //first the lying B blocks
      rect_matrices[1] = reinterpret_cast<MatrixD*>(blocks.at(7).get());
      rect_matrices[2] = reinterpret_cast<MatrixD*>(blocks.at(2).get()); //than the standing B blocks
      rect_matrices[3] = reinterpret_cast<MatrixD*>(blocks.at(5).get());
      break;
    case 14:
      sq_matrices[0] = reinterpret_cast<SquareMatrixD*>(blocks.at(0).get());
      sq_matrices[1] = reinterpret_cast<SquareMatrixD*>(blocks.at(1).get());
      sq_matrices[2] = reinterpret_cast<SquareMatrixD*>(blocks.at(3).get());
      sq_matrices[3] = reinterpret_cast<SquareMatrixD*>(blocks.at(4).get());
      sq_matrices[4] = reinterpret_cast<SquareMatrixD*>(blocks.at(8).get());
      rect_matrices[0] = reinterpret_cast<MatrixD*>(blocks.at(6).get()); //first the lying B blocks
      rect_matrices[1] = reinterpret_cast<MatrixD*>(blocks.at(7).get());
      rect_matrices[2] = reinterpret_cast<MatrixD*>(blocks.at(2).get()); //than the standing B blocks
      rect_matrices[3] = reinterpret_cast<MatrixD*>(blocks.at(5).get());
      break;
    default:
      ErrThrow("Sorry, the structure of that BlockMatrix is unknown to class NSE2D. "
              "I don't know how to pass its blocks to Assemble2D.");
    }

    //HOTFIX: Check the documentation!
    assemble_nse = Hotfixglobal_AssembleNSE::WITHOUT_CONVECTION;


    // call the assemble method with the information that has been patched together
#ifdef __3D__
    Assemble3D(
#else
    Assemble2D(
#endif
            N_FESpaces, fespmat, n_sq_mat, sq_matrices,
            n_rect_mat, rect_matrices, N_Rhs, RHSs, fesprhs,
            boundary_conditions, non_const_bound_values.data(), la);

    // assemble on the boundary if needed
    this->assemble_boundary_terms();


    // copy Dirichlet values from right hand side into solution
    s.solution.copy_nonactive(s.rhs);

    // TODO Maybe we have to explicitely set non-actives in non-diagonal blocks
    // to zero here, that was done in former code, but maybe we can move it to the solver part

    //tidy up
    delete fe_functions[0];
    delete fe_functions[1];

    if(this->db["space_discretization_type"].is("local_projection"))
    {
      LPS_parameter_set lp{this->db["lps_coeff_type"], this->db["lps_delta0"],
        this->db["lps_delta1"]};
      auto C = LPS_for_pressure_Scott_Zhang(blocks.at(8), false,
              this->example.get_nu(), lp);
      s.matrix.replace_blocks(*C, {{2,2}}, {false}); // creates a copy
    }
  }
}



#else

/*
#ifdef _MPI
NSE3D_GPPO::NSE3D_GPPO( TDomain &domain, const ParameterDatabase& param_db, const Example_NSE3D& example)
:NSE3D(domain.refine_and_get_hierarchy_of_collections(param_db, 0), param_db, example)
{
}
#else*/
template<int d>
NSE_GPPO<d>::NSE_GPPO( const TDomain &domain, const ParameterDatabase& param_db, const Example_NSE& example)
:NavierStokes<d>(domain, param_db, example) //(domain.refine_and_get_hierarchy_of_collections(param_db), param_db, example)
{
}
//#endif

/*
void read_coefficient_function(std::string read_coefficient_function_directory, TCollection *coll, TFEFunction3D *coefficient_function_ptr)
{

  cout <<" ***** coefficient_function_type 2 detected **** "<< endl;
          Output::print("It is assumed that a mesh and a fitting TFEFunction2D ReadSol() file are provided.");

          // fe space of piecewise constant functions
         // TCollection *coll = brinkman2d.get_pressure_space().GetCollection();

          //Domain.GetCollection();
          //   cout << "*********** number of cells:" << coll->GetN_Cells()<<endl;
          /// TCollection tmp_collection = TCollection(brinkman2d.get_pressure_space().GetCollection()->GetN_Cells(),
          /// brinkman2d.get_pressure_space().GetCollection()->GetCells());
          /// read_coll = *tmp_collection;

          TFESpace3D coefficient_function_FEspace(coll, "coefficient_function_FEspace", "s",
              BoundConditionNoBoundCondition, 0);
          BlockVector coefficient_function_vector(coefficient_function_FEspace.GetN_DegreesOfFreedom());

          TFEFunction3D coefficient_function(&coefficient_function_FEspace, "coefficient_function", "coefficient_function",
              &coefficient_function_vector.at(0), coefficient_function_FEspace.GetN_DegreesOfFreedom());

          coefficient_function.ReadSol(parmoon_db["read_coefficient_function_directory"]);

          ////coefficient_function.Interpolate(&read_coefficient_function);

          //TFEFunction3D *coefficient_function_ptr;
          coefficient_function_ptr = &coefficient_function;

         /// brinkman2d.assemble(i, coefficient_function_ptr);

}












*/
/*
void assemble(TFEFunction3D* coefficient_function)
{
   // assemble all matrices and right hand side
 NSE3D::assemble_linear_terms();

}
*/


/*
void assemble(TFEFunction3D* coefficient_function)

{
  const TFESpace3D *coefficient_function_space;

  for(System_per_grid& s : this->systems)
  {
    s.rhs.reset(); //right hand side reset (TODO: is that necessary?)
    s.matrix.reset(); // reset matrix (needed for mdml where this is called)

    const TFESpace3D * v_space = s.velocity_space.get();
    const TFESpace3D * p_space = s.pressure_space.get();

    // declare the variables which Assemble3D needs and each nstype has to fill
    size_t N_FESpaces = 2;

    const TFESpace3D *fespmat[3] = {v_space, p_space, nullptr};

    if (coefficient_function)
    {
      fespmat[2] = coefficient_function_space;
    }

    size_t n_sq_mat;
    TSquareMatrix3D *sq_matrices[10]{nullptr};//it's 10 pointers maximum (Type14)

    size_t n_rect_mat;
    TMatrix3D *rect_matrices[6]{nullptr};//it's four pointers maximum (Types 2, 4, 14)

    size_t N_Rhs = 4;
    double *RHSs[4] = {s.rhs.block(0), s.rhs.block(1), s.rhs.block(2), s.rhs.block(3)};
    const TFESpace3D *fesprhs[4] = {v_space, v_space, v_space, p_space};

    BoundCondFunct3D * boundary_conditions[4] = {
      v_space->get_boundary_condition(), v_space->get_boundary_condition()
      , v_space->get_boundary_condition(), p_space->get_boundary_condition() };

    std::array<BoundValueFunct3D*, 4> non_const_bound_values;
    non_const_bound_values[0] = example.get_bd()[0];
    non_const_bound_values[1] = example.get_bd()[1];
    non_const_bound_values[2] = example.get_bd()[2];
    non_const_bound_values[3] = example.get_bd()[3];


    //same for all: the local assembling object
    TFEFunction3D *fe_functions[5] =
    { s.u.GetComponent(0), s.u.GetComponent(1), s.u.GetComponent(2), &s.p, NULL };

    if (coefficient_function)
    {
      coefficient_function_space = coefficient_function->GetFESpace3D();
      fe_functions[4] = coefficient_function;
    }

    LocalAssembling3D la(this->db, LocalAssembling_type::NSE3D_Linear,
                        feFunction.data(), example_.get_coeffs());

    if (coefficient_function)
 {
   // modify la such that it includes the TFEFunction2D coefficient_function
   la.setBeginParameter({0});
   la.setN_ParamFct(1);
   la.setParameterFct({mapping_local_parameters});
   la.setN_FeValues(1);
   la.setFeValueFctIndex({3});
   la.setFeValueMultiIndex({D00});
 }

    std::vector<std::shared_ptr<FEMatrix>> blocks =
        s.matrix.get_blocks_uniquely();

    n_sq_mat = 10;
    n_rect_mat = 6;
    switch(TDatabase::ParamDB->NSTYPE)
    {// switch over known Block Matrix types, treat each one individually,
      // using a priori knowledge about the structure and the way it fits
      // to the LocalAssembling2D object
      // TODO remove all reinterpret_casts as soon as Assembling process takes only FEMatrices
      // we have to use reinterpret_casts because dynamic downcasting won't work here
      // FIXME replace global switch by local checking of blockmatrix type!
    case 1:
      sqMatrices[0]=reinterpret_cast<TSquareMatrix3D*>(blocks[0].get());

      reMatrices[0]=reinterpret_cast<TMatrix3D*>(blocks[1].get());
      reMatrices[1]=reinterpret_cast<TMatrix3D*>(blocks[2].get());
      reMatrices[2]=reinterpret_cast<TMatrix3D*>(blocks[3].get());
      break;
    case 2:
      sqMatrices[0]=reinterpret_cast<TSquareMatrix3D*>(blocks[0].get());
      reMatrices[0]=reinterpret_cast<TMatrix3D*>(blocks[4].get());
      reMatrices[1]=reinterpret_cast<TMatrix3D*>(blocks[5].get());
      reMatrices[2]=reinterpret_cast<TMatrix3D*>(blocks[6].get());
      reMatrices[3]=reinterpret_cast<TMatrix3D*>(blocks[1].get());
      reMatrices[4]=reinterpret_cast<TMatrix3D*>(blocks[2].get());
      reMatrices[5]=reinterpret_cast<TMatrix3D*>(blocks[3].get());
      break;
    case 3:
      sqMatrices[0]=reinterpret_cast<TSquareMatrix3D*>(blocks[0].get());
      sqMatrices[1]=reinterpret_cast<TSquareMatrix3D*>(blocks[1].get());
      sqMatrices[2]=reinterpret_cast<TSquareMatrix3D*>(blocks[2].get());
      sqMatrices[3]=reinterpret_cast<TSquareMatrix3D*>(blocks[4].get());
      sqMatrices[4]=reinterpret_cast<TSquareMatrix3D*>(blocks[5].get());
      sqMatrices[5]=reinterpret_cast<TSquareMatrix3D*>(blocks[6].get());
      sqMatrices[6]=reinterpret_cast<TSquareMatrix3D*>(blocks[8].get());
      sqMatrices[7]=reinterpret_cast<TSquareMatrix3D*>(blocks[9].get());
      sqMatrices[8]=reinterpret_cast<TSquareMatrix3D*>(blocks[10].get());
      reMatrices[0]=reinterpret_cast<TMatrix3D*>(blocks[3].get());
      reMatrices[1]=reinterpret_cast<TMatrix3D*>(blocks[7].get());
      reMatrices[2]=reinterpret_cast<TMatrix3D*>(blocks[11].get());
      break;
    case 4:
      sqMatrices[0]=reinterpret_cast<TSquareMatrix3D*>(blocks[0].get());
      sqMatrices[1]=reinterpret_cast<TSquareMatrix3D*>(blocks[1].get());
      sqMatrices[2]=reinterpret_cast<TSquareMatrix3D*>(blocks[2].get());
      sqMatrices[3]=reinterpret_cast<TSquareMatrix3D*>(blocks[4].get());
      sqMatrices[4]=reinterpret_cast<TSquareMatrix3D*>(blocks[5].get());
      sqMatrices[5]=reinterpret_cast<TSquareMatrix3D*>(blocks[6].get());
      sqMatrices[6]=reinterpret_cast<TSquareMatrix3D*>(blocks[8].get());
      sqMatrices[7]=reinterpret_cast<TSquareMatrix3D*>(blocks[9].get());
      sqMatrices[8]=reinterpret_cast<TSquareMatrix3D*>(blocks[10].get());
      reMatrices[0]=reinterpret_cast<TMatrix3D*>(blocks[12].get()); //first the lying B blocks
      reMatrices[1]=reinterpret_cast<TMatrix3D*>(blocks[13].get());
      reMatrices[2]=reinterpret_cast<TMatrix3D*>(blocks[14].get());
      reMatrices[3]=reinterpret_cast<TMatrix3D*>(blocks[3].get()); //than the standing B blocks
      reMatrices[4]=reinterpret_cast<TMatrix3D*>(blocks[7].get());
      reMatrices[5]=reinterpret_cast<TMatrix3D*>(blocks[11].get());
      break;
    case 14:
      sqMatrices[0]=reinterpret_cast<TSquareMatrix3D*>(blocks[0].get());
      sqMatrices[1]=reinterpret_cast<TSquareMatrix3D*>(blocks[1].get());
      sqMatrices[2]=reinterpret_cast<TSquareMatrix3D*>(blocks[2].get());
      sqMatrices[3]=reinterpret_cast<TSquareMatrix3D*>(blocks[4].get());
      sqMatrices[4]=reinterpret_cast<TSquareMatrix3D*>(blocks[5].get());
      sqMatrices[5]=reinterpret_cast<TSquareMatrix3D*>(blocks[6].get());
      sqMatrices[6]=reinterpret_cast<TSquareMatrix3D*>(blocks[8].get());
      sqMatrices[7]=reinterpret_cast<TSquareMatrix3D*>(blocks[9].get());
      sqMatrices[8]=reinterpret_cast<TSquareMatrix3D*>(blocks[10].get());
      sqMatrices[9]=reinterpret_cast<TSquareMatrix3D*>(blocks[15].get());
      reMatrices[0]=reinterpret_cast<TMatrix3D*>(blocks[12].get());
      reMatrices[1]=reinterpret_cast<TMatrix3D*>(blocks[13].get());
      reMatrices[2]=reinterpret_cast<TMatrix3D*>(blocks[14].get());
      reMatrices[3]=reinterpret_cast<TMatrix3D*>(blocks[3].get());
      reMatrices[4]=reinterpret_cast<TMatrix3D*>(blocks[7].get());
      reMatrices[5]=reinterpret_cast<TMatrix3D*>(blocks[11].get());
      break;
      default:
        ErrThrow("Sorry, the structure of that BlockMatrix is unknown to class NSE2D. "
            "I don't know how to pass its blocks to Assemble2D.");
    }

    //HOTFIX: Check the documentation!
    assemble_nse = Hotfixglobal_AssembleNSE::WITHOUT_CONVECTION;


    // call the assemble method with the information that has been patched together
    Assemble3D(N_FESpaces, fespmat, n_sq_mat, sq_matrices,
               n_rect_mat, rect_matrices, N_Rhs, RHSs, fesprhs,
               boundary_conditions, non_const_bound_values.data(), la);

    // assemble on the boundary if needed
      assemble_boundary_terms();


    // copy Dirichlet values from right hand side into solution
    s.solution.copy_nonactive(s.rhs);

    // TODO Maybe we have to explicitely set non-actives in non-diagonal blocks
    // to zero here, that was done in former code, but maybe we can move it to the solver part

    //tidy up
    delete fe_functions[0];
    delete fe_functions[1];

    if(db["space_discretization_type"].is("local_projection"))
    {
      LPS_parameter_set lp{db["lps_coeff_type"], db["lps_delta0"],
                           db["lps_delta1"]};
      auto C = LPS_for_pressure_Scott_Zhang(blocks.at(11), false,
                                            this->example.get_nu(), lp);
      s.matrix.replace_blocks(*C, {{3,3}}, {false}); // creates a copy
    }
  }
}

*/


#endif



#ifdef __3D__
template class NSE_GPPO<3>;
#else
template class NSE_GPPO<2>;
#endif


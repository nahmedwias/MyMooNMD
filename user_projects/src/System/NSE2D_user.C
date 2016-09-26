#include <NSE2D_user.h>
#include <MainUtilities.h> // GetVelocityAndPressureSpace
#include <Database.h>
#include <LinAlg.h> // Ddot, IntoL20FEFunction
#include <Upwind.h>
#include <GridTransfer.h>
#include <Multigrid.h>
#include<Assemble2D.h>

/** ************************************************************************ */
void NSE2D::assemble(TFEFunction2D* rho_field,
                     TFEFunction2D* mu_field)
{
  for(System_per_grid& s : this->systems)
  {
    s.rhs.reset(); //right hand side reset (TODO: is that necessary?)
    s.matrix.reset(); // reset matrix (needed for mdml where this is called)

    const TFESpace2D * v_space = &s.velocity_space;
    const TFESpace2D * p_space = &s.pressure_space;

    // declare the variables which Assemble2D needs and each nstype has to fill
    size_t N_FESpaces = 2;

    const TFESpace2D *fespmat[2] = {v_space, p_space};
    size_t n_sq_mat;
    TSquareMatrix2D *sq_matrices[5]{nullptr};//it's five pointers maximum (Type14)

    size_t n_rect_mat;
    TMatrix2D *rect_matrices[4]{nullptr};//it's four pointers maximum (Types 2, 4, 14)

    size_t N_Rhs = 2; //is 3 if NSE type is 4 or 14
    double *RHSs[3] = {s.rhs.block(0), s.rhs.block(1), nullptr}; //third place gets only filled
    const TFESpace2D *fesprhs[3] = {v_space, v_space, nullptr};  // if NSE type is 4 or 14

    BoundCondFunct2D * boundary_conditions[3] = {
      v_space->GetBoundCondition(), v_space->GetBoundCondition(),
      p_space->GetBoundCondition() };
    std::array<BoundValueFunct2D*, 3> non_const_bound_values;
    non_const_bound_values[0] = example.get_bd()[0];
    non_const_bound_values[1] = example.get_bd()[1];
    non_const_bound_values[2] = example.get_bd()[2];

//    double resultats[3];
//    rho_field->FindGradient(0.187, 0.42987, resultats);
//    cout << resultats[0] << endl;
//    cout << resultats[1] << endl;
//    cout << resultats[2] << endl;

    //same for all: the local asembling object
    TFEFunction2D *fe_functions[5] =
      { s.u.GetComponent(0), s.u.GetComponent(1), &s.p, nullptr, nullptr };
    LocalAssembling2D la(NSE2D_Galerkin, fe_functions,
                     this->example.get_coeffs());
    if (rho_field != nullptr && mu_field != nullptr)
    {
      fe_functions[3] = rho_field;
      fe_functions[4] = mu_field;
      la.setBeginParameter({0});
      la.setFeFunctions2D(fe_functions); //reset - now velo comp included
      la.setFeValueFctIndex({0,1,3,4});
      la.setFeValueMultiIndex({D00,D00,D00,D00});
      la.setN_Parameters(4);
      la.setN_FeValues(4);
      la.setN_ParamFct(1);
      // la.setParameterFct({NSParamsVelo});
      //...this should do the trick
    }

    std::vector<std::shared_ptr<FEMatrix>> blocks = s.matrix.get_blocks_uniquely();

    switch(TDatabase::ParamDB->NSTYPE)
    {// switch over known Block Matrix types, treat each one individually,
      // using a priori knowledge about the structure and the way it fits
      // to the LocalAssembling2D object
      // TODO remove all reinterpret_casts as soon as Assembling process takes only FEMatrices
      // we have to use reinterpret_casts because dynamic downcasting won't work here
      // FIXME replace global switch by local checking of blockmatrix type!
      case 1:
        n_sq_mat = 1;
        sq_matrices[0] =  reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());

        n_rect_mat = 2;
        rect_matrices[0] = reinterpret_cast<TMatrix2D*>(blocks.at(1).get());
        rect_matrices[1] = reinterpret_cast<TMatrix2D*>(blocks.at(2).get());

        break;
      case 2:
        n_sq_mat = 1;
        sq_matrices[0] =  reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());

        n_rect_mat = 4;
        rect_matrices[0] = reinterpret_cast<TMatrix2D*>(blocks.at(3).get()); //first the lying B blocks
        rect_matrices[1] = reinterpret_cast<TMatrix2D*>(blocks.at(4).get());
        rect_matrices[2] = reinterpret_cast<TMatrix2D*>(blocks.at(1).get()); //than the standing B blocks
        rect_matrices[3] = reinterpret_cast<TMatrix2D*>(blocks.at(2).get());

        break;
      case 3:
        n_sq_mat = 4;
        sq_matrices[0] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());
        sq_matrices[1] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(1).get());
        sq_matrices[2] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(3).get());
        sq_matrices[3] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(4).get());

        n_rect_mat = 2;
        rect_matrices[0] = reinterpret_cast<TMatrix2D*>(blocks.at(2).get());
        rect_matrices[1] = reinterpret_cast<TMatrix2D*>(blocks.at(5).get());
        break;

      case 4:
        n_sq_mat = 4;
        sq_matrices[0] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());
        sq_matrices[1] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(1).get());
        sq_matrices[2] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(3).get());
        sq_matrices[3] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(4).get());

        n_rect_mat = 4;
        rect_matrices[0] = reinterpret_cast<TMatrix2D*>(blocks.at(6).get()); //first the lying B blocks
        rect_matrices[1] = reinterpret_cast<TMatrix2D*>(blocks.at(7).get());
        rect_matrices[2] = reinterpret_cast<TMatrix2D*>(blocks.at(2).get()); //than the standing B blocks
        rect_matrices[3] = reinterpret_cast<TMatrix2D*>(blocks.at(5).get());

        RHSs[2] = s.rhs.block(2); // NSE type 4 includes pressure rhs
        fesprhs[2]  = p_space;
        N_Rhs = 3;

        break;

      case 14:
        n_sq_mat = 5;
        sq_matrices[0] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());
        sq_matrices[1] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(1).get());
        sq_matrices[2] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(3).get());
        sq_matrices[3] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(4).get());
        sq_matrices[4] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(8).get());

        n_rect_mat = 4;
        rect_matrices[0] = reinterpret_cast<TMatrix2D*>(blocks.at(6).get()); //first the lying B blocks
        rect_matrices[1] = reinterpret_cast<TMatrix2D*>(blocks.at(7).get());
        rect_matrices[2] = reinterpret_cast<TMatrix2D*>(blocks.at(2).get()); //than the standing B blocks
        rect_matrices[3] = reinterpret_cast<TMatrix2D*>(blocks.at(5).get());

        RHSs[2] = s.rhs.block(2); // NSE type 14 includes pressure rhs
        fesprhs[2]  = p_space;
        N_Rhs = 3;

        break;
      default:
        ErrThrow("Sorry, the structure of that BlockMatrix is unknown to class NSE2D. "
            "I don't know how to pass its blocks to Assemble2D.");
    }

    // call the assemble method with the information that has been patched together
    Assemble2D(N_FESpaces, fespmat, n_sq_mat, sq_matrices,
               n_rect_mat, rect_matrices, N_Rhs, RHSs, fesprhs,
               boundary_conditions, non_const_bound_values.data(), la);

    // copy Dirichlet values from right hand side into solution
    s.solution.copy_nonactive(s.rhs);

    // TODO Maybe we have to explicitely set non-actives in non-diagonal blocks
    // to zero here, that was done in former code, but maybe we can move it to the solver part

    //tidy up
    delete fe_functions[0];
    delete fe_functions[1];
  }
}

/** ************************************************************************ */
void NSE2D::assemble_nonlinear_term(TFEFunction2D* rho_field,
                                    TFEFunction2D* mu_field)
{
  // the class LocalAssembling2D which we will need next, requires an array of
  // pointers to finite element functions, i.e. TFEFunction2D **.

  //Nonlinear assembling requires an approximate velocity solution on every grid!
  if(systems.size() > 1)
  {
    for( int block = 0; block < 2 ;++block)
    {
      std::vector<const TFESpace2D*> spaces;
      std::vector<double*> u_entries;
      std::vector<size_t> u_ns_dofs;
      for(auto &s : systems )
      {
        spaces.push_back(&s.velocity_space);
        u_entries.push_back(s.solution.block(block));
        u_ns_dofs.push_back(s.solution.length(block));
      }
      GridTransfer::RestrictFunctionRepeatedly(spaces, u_entries, u_ns_dofs);
    }
  }
  
  bool mdml =  this->solver.is_using_multigrid() 
            && this->solver.get_multigrid()->is_using_mdml();
  bool is_stokes = this->db["problem_type"].is(3); // otherwise Navier-Stokes
  
  if(mdml && !is_stokes)
  {
    // in case of upwinding we only assemble the linear terms. The nonlinear
    // term is not assembled but replaced by a call to the upwind method.
    // Note that we assemble the same terms over and over again here. Not 
    // nice, but otherwise we would have to store the linear parts in a 
    // separate BlockFEMatrix.
    this->assemble();
  }

  for(System_per_grid& s : this->systems)
  {
    //hold the velocity space, we'll need it...
    const TFESpace2D * v_space = &s.velocity_space;

    //the variables we will have to fill for the call to Assemble2D
    size_t n_fe_spaces = 1;
    const TFESpace2D* fe_spaces[1]{v_space};

    size_t n_sq_mat;
    TSquareMatrix2D* sq_mat[2]{nullptr};//two pointers maximum

    size_t n_rect_mat = 0;
    TMatrix2D** rect_mat = nullptr;

    size_t n_rhs = 0;
    double** rhs = nullptr;
    const TFESpace2D** fe_rhs = nullptr;

    BoundCondFunct2D * boundary_conditions[1] = { v_space->GetBoundCondition() };
    std::array<BoundValueFunct2D*, 3> non_const_bound_values;
    non_const_bound_values[0] = example.get_bd()[0];
    non_const_bound_values[1] = example.get_bd()[1];
    non_const_bound_values[2] = example.get_bd()[2];

    //same for all: the local asembling object
    TFEFunction2D *fe_functions[5] =
    { s.u.GetComponent(0), s.u.GetComponent(1), &s.p, nullptr, nullptr };
    LocalAssembling2D la_nonlinear(NSE2D_Galerkin_Nonlinear, fe_functions,
                         this->example.get_coeffs());
    if (rho_field != nullptr && mu_field != nullptr)
    {
      fe_functions[3] = rho_field;
      fe_functions[4] = mu_field;
      la_nonlinear.setBeginParameter({0});
      la_nonlinear.setFeFunctions2D(fe_functions); //reset - now velo comp included
      la_nonlinear.setFeValueFctIndex({0,1,3,4});
      la_nonlinear.setFeValueMultiIndex({D00,D00,D00,D00});
      la_nonlinear.setN_Parameters(4);
      la_nonlinear.setN_FeValues(4);
      la_nonlinear.setN_ParamFct(1);
      // la.setParameterFct({NSParamsVelo});
      //...this should do the trick
    }


    //fetch us (a) pointer(s) to the diagonal A block(s)
    std::vector<std::shared_ptr<FEMatrix>> blocks = s.matrix.get_blocks_uniquely({{0,0},{1,1}});

    switch(TDatabase::ParamDB->NSTYPE)
    {// switch over known Block Matrix types, treat each one individually,
      // using a priori knowledge about the structure and the way it fits
      // to the LocalAssembling2D object
      // TODO remove all reinterpret casts as soon as Assembling process takes only FEMatrices
      // FIXME replace global switch by local checking of blockmatrix type!
      case 1:
        n_sq_mat = 1;
        sq_mat[0] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());
        break;
      case 2:
        n_sq_mat = 1;
        sq_mat[0] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());
        break;
      case 3:
        n_sq_mat = 2;
        sq_mat[0] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());
        sq_mat[1] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(1).get());
        break;
      case 4:
        n_sq_mat = 2;
        sq_mat[0] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());
        sq_mat[1] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(1).get());
        break;
      case 14:
        n_sq_mat = 2;
        sq_mat[0] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());
        sq_mat[1] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(1).get());
        break;
      default:
        ErrThrow("Sorry, the structure of that BlockMatrix is unknown to class NSE2D. "
            "I don't know how to pass its blocks to Assemble2D.");
    }
    
    // do upwinding TODO remove dependency of global values
    bool on_finest_grid = &systems.front() == &s;
    bool do_upwinding = (TDatabase::ParamDB->DISCTYPE == UPWIND 
                         || (mdml && !on_finest_grid))
                        && !is_stokes;

    // assemble the nonlinear part of NSE
    if(!do_upwinding)
    {
      for(size_t i =0; i < n_sq_mat; ++i)
      {
        //reset the matrices, linear part is assembled anew
        sq_mat[i]->reset();
      }
      //do the actual assembling
      Assemble2D(n_fe_spaces, fe_spaces, n_sq_mat, sq_mat, n_rect_mat, rect_mat,
                 n_rhs, rhs, fe_rhs, boundary_conditions,
                 non_const_bound_values.data(), la_nonlinear);
    }
    else
    {
      switch(TDatabase::ParamDB->NSTYPE)
      {
        case 1:
        case 2:
          // do upwinding with one matrix
          UpwindForNavierStokes(la_nonlinear.GetCoeffFct(), sq_mat[0],
                                la_nonlinear.get_fe_function(0),
                                la_nonlinear.get_fe_function(1));
          Output::print<3>("UPWINDING DONE : level ");
          break;

        case 3:
        case 4:
        case 14:
          // do upwinding with two matrices
          Output::print<3>("UPWINDING DONE : level ");
          UpwindForNavierStokes(la_nonlinear.GetCoeffFct(), sq_mat[0],
                                la_nonlinear.get_fe_function(0),
                                la_nonlinear.get_fe_function(1));
          UpwindForNavierStokes(la_nonlinear.GetCoeffFct(), sq_mat[1],
                                la_nonlinear.get_fe_function(0),
                                la_nonlinear.get_fe_function(1));
          break;
      } // endswitch
    } // endif

    //end nonlinear assembling

    // copy Dirichlet values from right hand side into solution
    //(is this necessary here? solution has not been touched!)
    s.solution.copy_nonactive(s.rhs);

    //tidy up
    delete fe_functions[0];
    delete fe_functions[1];

  }


}


/** ************************************************************************ */
void NSE2D::update_solution(BlockVector weight_vector)
{
  const unsigned int length0 = this->get_solution().length(0);
  const unsigned int length1 = this->get_solution().length(1);
  const unsigned int length2 = weight_vector.length();

  if (length0 != length2 || length1 != length2)
  {
    Output::warn<1>("MULTIPHASE: ATTENTION, BlockVectors  "
        "don't have same length! ", length0, " ", length1, " ", length2);
    ErrThrow("See Warning");
  }
  else
  {
    Output::print<3>("NSE2D: ", "Ok, blocks are same length and can be multiplied");
  }

  for (unsigned int i = 0; i < length0; i++)
  {
    this->get_solution().at(i) *= weight_vector.at(i);
  }

  for (unsigned int i = 0; i < length0; i++)
  {
    this->get_solution().at(i + length0) *= weight_vector.at(i);
  }
}


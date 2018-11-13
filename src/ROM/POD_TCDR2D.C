/** ************************************************************************ 
*
* @name       POD_TCDR2D
* @brief      POD for time-dep. scalar convection-diffusion problems in 2D
*
* @author     Swetlana Giere
* @date       08.03.2017 (start of implementation)
*
****************************************************************************/

#include <POD_TCDR2D.h>
#include <Database.h>
#include <MainUtilities.h>
#include <LocalProjection.h>
#include <Assemble2D.h>


/**************************************************************************** */
POD_TCDR2D::POD_TCDR2D(TCollection& coll, const ParameterDatabase& param_db,
		const Example_TimeCD2D& ex)
 : POD(param_db),
   example(ex),
   fe_space(&coll, (char*)"space", (char*)"cd2d space", example.get_bc(0), TDatabase::ParamDB->ANSATZ_ORDER),
   // TODO CB: Building the matrix here and rebuilding later is due to the
   // highly non-functional class TFEVectFunction2D (and TFEFunction2D,
   // which do neither provide default constructors nor working copy assignments.)
   gramian_matrix({&fe_space}),
   pod_mode(this->gramian_matrix, false),
   fe_function(&this->fe_space, (char*)"c", (char*)"c",
   this->pod_mode.get_entries(), this->pod_mode.length()),
   outputWriter(param_db)
{
  this->gramian_matrix = BlockFEMatrix::CD2D(fe_space);
  this->set_parameters();
  int ndof = this->fe_space.GetN_DegreesOfFreedom();

  double hmin, hmax;
  coll.GetHminHmax(&hmin, &hmax);
  Output::print<1>("N_Cells    : ", setw(12), coll.GetN_Cells());
  Output::print<1>("h(min,max) : ", setw(12), hmin, " ", setw(12), hmax);
  Output::print<1>("dof        : ", setw(12), ndof);
  Output::print<1>("active dof : ", setw(12), this->fe_space.GetN_ActiveDegrees());

  // Initialize the output object, add the fe function to it.
  outputWriter.add_fe_function(&this->fe_function);
}

/**************************************************************************** */
POD_TCDR2D::~POD_TCDR2D()
{
}

/**************************************************************************** */
void POD_TCDR2D::set_parameters()
{
  //set problem_type to Time_CD if not yet set
  if(!rom_db["problem_type"].is(2))
  {
    if (rom_db["problem_type"].is(0))
    {
      rom_db["problem_type"] = 2;
    }
    else
    {
      Output::warn<2>("The parameter problem_type doesn't correspond to Time_CD."
          "It is now reset to the correct value for Time_CD (=2).");
      rom_db["problem_type"] = 2;
    }
  }

  // an error when using ansatz order 0
  if(TDatabase::ParamDB->ANSATZ_ORDER == 0)
  {
    throw std::runtime_error("Ansatz order 0 is no use in convection diffusion "
        "reaction problems! (Vanishing convection and diffusion term).");
  }

  // For assembling of the Gramian matrix for the computation of
  // the POD basis only the "Galerkin" part of the matrices is needed.
  // TDatabase::ParamDB->DISCTYPE = GALERKIN;
  if(!rom_db["space_discretization_type"].is("galerkin"))
    rom_db["space_discretization_type"] = "galerkin";

}

/**************************************************************************** */
void POD_TCDR2D::assemble_gramian()
{
  Output::print<1>("Assembling the gramian matrix for POD computation...");
  LocalAssembling_type type;
  if (rom_db["pod_inner_product"].get<std::string>() == "l2")
  {
    type = LocalAssembling_type::TCDMassOnly;
  }
  else if (rom_db["pod_inner_product"].get<std::string>() == "eucl")
  {
    Output::print<1>("For given parameter 'pod_inner_product' no assembling needed.");
    return;
  }

  TFEFunction2D * pointer_to_function = &this->fe_function;
  LocalAssembling2D la_mat(this->get_db(), type, &pointer_to_function, this->example.get_coeffs());
  const TFESpace2D * _fe_space = &this->fe_space;
  BoundCondFunct2D * boundary_conditions = _fe_space->get_boundary_condition();
  int N_Matrices = 1;
  BoundValueFunct2D * non_const_bound_value[1] {this->example.get_bd()[0]};

  //fetch gramian matrix as block
  std::vector<std::shared_ptr<FEMatrix>> mat_blocks = this->gramian_matrix.get_blocks_uniquely();
  TSquareMatrix2D * mat_block[1]{reinterpret_cast<TSquareMatrix2D*>(mat_blocks.at(0).get())};

  //do the assembling
  mat_block[0]->reset();
  Assemble2D(1, &_fe_space, N_Matrices, mat_block, 0, NULL, 0, NULL,
                 NULL, &boundary_conditions, non_const_bound_value, la_mat);
}

/**************************************************************************** */
void POD_TCDR2D::compute_pod_basis()
{
  if(this->rom_db["pod_inner_product"] != "eucl")
  {
    assemble_gramian();
    //cout<<"norm: "<< gramian_matrix.get_blocks().at(0)->GetNorm()<<endl;exit(0);
    POD::set_gramian(this->gramian_matrix.get_combined_matrix());
  }
  POD::compute_basis();

  if(rom_db["pod_inner_product"].get<std::string>() == "eucl")
  {
    Output::print<1>("POD w.r.t euclidean inner product.");
    POD::write_pod(rom_db["pod_basename"].get<std::string>()+"eucl.");
  }
  else if(rom_db["pod_inner_product"].get<std::string>() == "l2")
  {
    Output::print<1>("POD w.r.t L2 inner product.");
    POD::write_pod(rom_db["pod_basename"].get<std::string>()+"l2.");
  }
  else
  {
    ErrThrow("Unknown value of parameter 'pod_inner_product'. Choose eucl or l2.");
  }
}

/**************************************************************************** */
void POD_TCDR2D::read_basis()
{
  POD::read_basis();
}

/**************************************************************************** */
void POD_TCDR2D::output()
{
  if (this->fe_space.GetN_DegreesOfFreedom() != POD::get_basis().size1())
  {
	ErrThrow("Current FE space dimension does not coincide with the number "
			 "of dof of POD basis.\n"
			 "Dimension of FE space : ", this->fe_space.GetN_DegreesOfFreedom(),
			 "\nDOF of POD basis      : ", POD::get_basis().size1());
  }

  /* write averages of snapshots into a vtk-file */
  if(this->rom_db["pod_fluct"])
  {
    for(int j=0; j< this->fe_space.GetN_DegreesOfFreedom(); ++j)
    {
      this->pod_mode.get_entries()[j] = POD::get_snaps_avr()(j);
    }
    Output::print<1>("Writing vtk for snapshots average into ",
    		         this->rom_db["output_basename"], 0 ,".vtk");
    this->outputWriter.write(0);
  }

  for(int i=0; i < std::min(POD::get_rank(),10); ++i)
  {
	for(int j=0; j<this->fe_space.GetN_DegreesOfFreedom(); ++j)
	{
	  this->pod_mode.get_entries()[j] = POD::get_basis()(j,i);
	}
	Output::print<1>("Writing vtk for POD basis into ",
			         this->rom_db["output_basename"], i+1 ,".vtk");
	this->outputWriter.write(i+1);
  }
}





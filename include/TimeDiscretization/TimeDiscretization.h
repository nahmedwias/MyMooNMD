#ifndef TIMEDISCRETIZATION_H
#define TIMEDISCRETIZATION_H

#include <RungeKuttaTable.h>
#include <ParameterDatabase.h>
#include <vector>
#include <memory>
#include <BlockFEMatrix.h>
#include <BlockVector.h>

class TimeDiscretization
{
  protected:
    std::shared_ptr<RungeKuttaTable> rk;

    std::vector<double> bdf_coefficients;

    ParameterDatabase db;

    double current_time_step_length;


    static ParameterDatabase default_TimeDiscretization_database();

public:
    // 
    TimeDiscretization(const ParameterDatabase& param_db);

    //
    unsigned int current_step_;
    
    /// parameter to be set for scaling and descaling of the 
    /// pressure-velocity blocks
    unsigned int n_scale_block;

    void set_time_disc_parameters();

    // a function which takes the matrices, right-hand side and solution vector
    // and perform the ...
    void prepare_timestep(BlockFEMatrix& system_matrix, 
        const BlockFEMatrix& mass_matrix, std::vector<BlockVector> & rhs, 
        const std::vector<BlockVector> old_solutions);
    //
    void prepare_rhs_from_time_disc(BlockFEMatrix& system_matrix, 
        const BlockFEMatrix& mass_matrix, std::vector<BlockVector> & rhs, 
        const std::vector<BlockVector> old_solutions);
    
    /// @brief this function will prepare the BlockFEMatrix 
    /// which passes to the solver. In details, the system_matrix
    /// will be scaled depending on the time stepping scheme
    /// and the mass matrix will be added to it. The B-blocks
    /// will be only scaled with the corresponding factor:
    /// @param system_matrix stiffness matrix
    /// @param mass_matrix  the mass matrix
    /// @param it_counter nonlinear iteration counter 
    /// @param nl_b_block the BT and B-blocks are nonlinear
    ///        nl_b_block = 0, 2 means B and BT's blocks have 
    ///        to be scaled 
    ///        nl_b_block = 1 means only the BT's have to be scaled
    void prepare_system_matrix(BlockFEMatrix& system_matrix, 
        const BlockFEMatrix& mass_matrix, 
        unsigned int it_counter);
    
    // 
    void reset_linear_matrices(BlockFEMatrix& matrix, 
                               const BlockFEMatrix& mass);
    
    // pre-stage of the bdf-schemes
    bool pre_stage_bdf;
    
    // a method which returns how many old solutions I need
    unsigned int n_old_solutions() const
    {
      if(rk)
      {
         return 1;
      }
      else if(bdf_coefficients.size() > 0)
      {
        return bdf_coefficients.size() - 1;
      }
      else
      {
        ErrThrow("I am not implemented");
      }
    }
};

#endif // TIMEDISCRETIZATION_H

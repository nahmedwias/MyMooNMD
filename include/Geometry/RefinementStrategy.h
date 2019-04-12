// =============================================================================
//
// Purpose: This class provides refinement strategies. Provided eta_K's from
//          an error estimator, cells are marked in a boolean vector which are
//          to be refined.
//
// Author: Moritz Hoffmann, 06/07/15.
//
// =============================================================================

#ifndef PARMOON_REFINEMENTSTRATEGY_H
#define PARMOON_REFINEMENTSTRATEGY_H


#include "ParameterDatabase.h"
#include "ErrorEstimator.h"
#include <array>
#include <functional>
#include <vector>

enum class RefinementStrategyType
{
  // this refinement strategy marks cells if the error is above some tolerance
  // and relaxes the tolerance if not enough cells have been marked yet
  MaximalEstimatedLocalError = 0,
  // this refinement strategy marks cells if the error is above some tolerance
  // and relaxes the tolerance if not enough cells have been marked yet or the
  // accumulated error of the cells is not some prescribed portion of the global
  // error
  PortionOfEstimatedGlobalError,
  // this refinement strategy is computationally more heavy than the other
  // strategies and aims at refining as few cells as possible while sustaining
  // uniform convergence of the adaptive process
  EquilibrationStrategy
};

// enable output of refinement strategy type
std::ostream &operator<<(std::ostream &os, RefinementStrategyType type);

template <int d>
class RefinementStrategy 
{
  protected:
    RefinementStrategyType refine_strategy;
    double reftol;
    double coarsetol;
    double min_fraction_to_change;
    double decrease_reftol_factor;
    double increase_coarsetol_factor;
    double fraction_of_error;
    int max_cell_level;
    std::vector<bool> refinements;
    int current_estimator;
  public:
    
    RefinementStrategy(const ParameterDatabase& db);
    ~RefinementStrategy() = default;
    
    /// @brief default parameters used in this class
    static ParameterDatabase default_refinement_strategy_database();

    /// @brief With a given ErrorEstimator find out which cells should be 
    /// refined/coarsened.
    /// This fills the member vector 'refinements'
    void apply_estimator(ErrorEstimator<d> &estimator);
    
    /// @brief Find out which cells should be refined using an indicator.
    /// The indicator function maps coordinates (x,y,z) (z=0 in 2D) to double.
    /// All cells which have vertices with both positive and negative indicator
    /// are marked to be refined. This fills the member vector 'refinements'.
    void apply_indicator(const TCollection* coll,
                         std::function<double(std::array<double, 3>)> indicator);

    /// @brief find out if a given cell (by index) should be refined.
    /// @note You have to call 'apply_estimator' before this makes any sense.
    bool should_refine_cell(size_t cell_index) const;
};


#endif //PARMOON_REFINEMENTSTRATEGY_H

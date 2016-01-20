// =============================================================================
// 
// Purpose: 
// 
// Author: Moritz Hoffmann, 15/12/15.
//
// =============================================================================

#ifndef PARMOON_STODA2DREFINEMENTSTRATEGY_H
#define PARMOON_STODA2DREFINEMENTSTRATEGY_H


#include <RefinementStrategy.h>
#include <CDErrorEstimator2D.h>
#include <NSEErrorEstimator2D.h>
#include <memory>

class StoDa2DRefinementStrategy : public RefinementStrategy {

protected:
    std::unique_ptr<RefinementStrategy> d_strategy;
    std::unique_ptr<RefinementStrategy> nse_strategy;

    // values to be set in applyEstimator(...)
    int reference_id_darcy;
    int reference_id_stokes;
    
    TCollection *finest_collection; // stokes and darcy domain together
    TCollection *stokes_collection;
    TCollection *darcy_collection;

    void applyEstimator(ErrorEstimator2D &estimator2D);
public:
    StoDa2DRefinementStrategy();

    StoDa2DRefinementStrategy(int refine_strategy,                             // TDatabase::ParamDB->REFINE_STRATEGY
                              double reftol,                                   // TDatabase::ParamDB->REFTOL
                              double coarsetol,                                // TDatabase::ParamDB->COARSETOL
                              double min_fraction_to_change,                   // TDatabase::ParamDB->MIN_FRACTION_TO_CHANGE
                              double decrease_reftol_factor,                   // TDatabase::ParamDB->DECREASE_REFTOL_FACTOR
                              double increase_coarsetol_factor,                // TDatabase::ParamDB->INCREASE_COARSETOL_FACTOR
                              double fraction_of_error,                        // TDatabase::ParamDB->FRACTION_OF_ERROR
                              int max_cell_level,                              // TDatabase::ParamDB->MAX_CELL_LEVEL
                              int current_estimator);                          // TDatabase::ParamDB->ADAPTIVE_REFINEMENT_CRITERION

    bool shouldRefineCell(size_t cellIndex);

    void applyEstimator(CDErrorEstimator2D &d_estimator, NSEErrorEstimator2D &nse_estimator, TCollection *finestCollection);
};


#endif //PARMOON_STODA2DREFINEMENTSTRATEGY_H

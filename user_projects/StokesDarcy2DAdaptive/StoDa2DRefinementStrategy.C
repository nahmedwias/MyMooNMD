// =============================================================================
// 
// Purpose: 
// 
// Author: Moritz Hoffmann, 15/12/15.
//
// =============================================================================

#include "StoDa2DRefinementStrategy.h"

StoDa2DRefinementStrategy::StoDa2DRefinementStrategy(int refine_strategy, double reftol, double coarsetol,
                                                     double min_fraction_to_change, double decrease_reftol_factor,
                                                     double increase_coarsetol_factor, double fraction_of_error,
                                                     int max_cell_level, int current_estimator) :
        RefinementStrategy(refine_strategy, reftol, coarsetol, min_fraction_to_change, decrease_reftol_factor,
                           increase_coarsetol_factor, fraction_of_error, max_cell_level, current_estimator) {
}

StoDa2DRefinementStrategy::StoDa2DRefinementStrategy() : RefinementStrategy() {
}

void StoDa2DRefinementStrategy::applyEstimator(ErrorEstimator2D &estimator2D) {
    std::cerr << "Called apply estimator(estimator2d) in stoda estimator, not supported!" << std::endl;
    exit(EXIT_FAILURE);
}

bool StoDa2DRefinementStrategy::shouldRefineCell(size_t cellIndex) {
    if (!finest_collection) {
        std::cerr << "no finest collection set, cannot check if cell should be refined or not." << std::endl;
        exit(EXIT_FAILURE);
    }
    auto cell = finest_collection->GetCell((int) cellIndex);
    if (cell->GetReference_ID() == this->reference_id_darcy) {
        auto darcy_strategy = this->d_strategy.get();
        return darcy_strategy->shouldRefineCell(darcy_collection->GetIndex(cell));
    } else if (cell->GetReference_ID() == this->reference_id_stokes) {
        auto stokes_strategy = this->nse_strategy.get();
        return stokes_strategy->shouldRefineCell(stokes_collection->GetIndex(cell));
    } else {
        std::cerr << "cell with index=" << cellIndex << " had reference_id=" << cell->GetReference_ID()
        << ", which is neither the darcy-reference-id=" << reference_id_darcy
        << " nor the stokes-reference-id=" << reference_id_stokes << std::endl;
        exit(EXIT_FAILURE);
    }
}

void StoDa2DRefinementStrategy::applyEstimator(CDErrorEstimator2D &d_estimator, NSEErrorEstimator2D &nse_estimator, TCollection *finestCollection) {
    Output::print<1>("---- refinement strategies start ----");
    d_strategy.reset(new RefinementStrategy(int(refine_strategy), reftol, coarsetol, min_fraction_to_change, decrease_reftol_factor, increase_coarsetol_factor, fraction_of_error, max_cell_level, current_estimator));
    nse_strategy.reset(new RefinementStrategy(int(refine_strategy), reftol, coarsetol, min_fraction_to_change, decrease_reftol_factor, increase_coarsetol_factor, fraction_of_error, max_cell_level, current_estimator));
    this->finest_collection = finestCollection;
    this->stokes_collection = nse_estimator.GetCollection();
    this->darcy_collection = d_estimator.GetCollection();

    {
        Output::print<1>("Applying Darcy estimator:");
        auto first_darcy_cell = d_estimator.GetCollection()->GetCell(0);
        reference_id_darcy = first_darcy_cell->GetReference_ID();
        d_strategy.get()->applyEstimator(d_estimator);
    }

    {
        Output::print<1>("Applying Stokes estimator:");
        auto first_cell_stokes = nse_estimator.GetCollection()->GetCell(0);
        reference_id_stokes = first_cell_stokes->GetReference_ID();
        nse_strategy.get()->applyEstimator(nse_estimator);
    }
    Output::print<1>("---- refinement strategies end ----");
}

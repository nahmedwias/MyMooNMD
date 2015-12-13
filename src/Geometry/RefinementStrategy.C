// =============================================================================
//
// Purpose: Implementation of RefinementStrategy.h.
//
// Author: Moritz Hoffmann, 06/07/15.
//
// =============================================================================

#include "RefinementStrategy.h"

RefinementStrategy::RefinementStrategy() : RefinementStrategy::RefinementStrategy(
        TDatabase::ParamDB->REFINE_STRATEGY,
        TDatabase::ParamDB->REFTOL,
        TDatabase::ParamDB->COARSETOL,
        TDatabase::ParamDB->MIN_FRACTION_TO_CHANGE,
        TDatabase::ParamDB->DECREASE_REFTOL_FACTOR,
        TDatabase::ParamDB->INCREASE_COARSETOL_FACTOR,
        TDatabase::ParamDB->FRACTION_OF_ERROR,
        TDatabase::ParamDB->MAX_CELL_LEVEL,
        TDatabase::ParamDB->ADAPTIVE_REFINEMENT_CRITERION) {}

bool RefinementStrategy::shouldRefineCell(size_t cellIndex) {
    return refinements[cellIndex];
}

void RefinementStrategy::applyEstimator(ErrorEstimator2D &estimator2D) {
    // get the collection of the finest available mesh
    auto collection = estimator2D.GetCollection();
    // number of cells in that finest mesh
    auto n_cells = collection->GetN_Cells();
    // resize indicator vector accordingly
    refinements.resize((unsigned long) n_cells);
    // initialize it with false, i.e., do not refine anything
    std::fill(refinements.begin(), refinements.end(), false);

    // number of marked cells, current iteration and maximal number of iterations
    int changed = -1, it = 0, max_it = 100;
    // maximal local error
    double eta_max = estimator2D.GetMaximalLocalError();
    // reference tolerance as percentage value of the maximal local error
    double reftol = this->reftol * eta_max;
    // since local estimates are stored as square
    reftol *= reftol;
    // reference coarsening tolerance as percentage value of the maximal local error
    double coarsetol = this->coarsetol * eta_max;
    // again, stored as square
    coarsetol = coarsetol * coarsetol;
    // minimal number of cells to change
    double min_changed = n_cells * min_fraction_to_change;
    // current cell
    TBaseCell *cell;

    OutPut("Refine strategy: " << refine_strategy << ", min cells to change: " << ((int) min_changed));
    switch (refine_strategy) {
        case RefinementStrategyType::MaximalEstimatedLocalError: {
            // loop while the number of changed cells is below the minimal number and we have not reached max_it
            while ((changed < min_changed) && (it < max_it)) {
                // reset number of changed cells
                // if this is not the first pass, the conditions have been relaxed and the cells that have been marked
                // will now be marked in particular
                changed = 0;
                // loop over all cells
                #pragma omp parallel for
                for (size_t m = 0; m < n_cells; m++) {
                    cell = collection->GetCell((int) m);
                    // mark cell for refinement
                    if ((estimator2D.GetEta_K()[m] >= reftol) && (cell->GetGeoLevel() < max_cell_level)) {
                        refinements[m] = true;
                        #pragma omp atomic update
                        changed++;
                    }
                    if (estimator2D.GetEta_K()[m] <= coarsetol)
                        ;   //TODO: mark cell for coarsening
                }
                // not enough cells marked, relax conditions
                if (changed < min_changed) {
                    reftol *= decrease_reftol_factor;
                    coarsetol *= increase_coarsetol_factor;
                }
                if (TDatabase::ParamDB->SC_VERBOSE > 1) OutPut("total " << n_cells << " changed " << changed << endl);
                it++;
            }
            break;
        }
        case RefinementStrategyType::PortionOfEstimatedGlobalError: {
            // mark smallest set of cells whose sum of local errors is a prescribed fraction of the global
            double sum_local_errors = 0.0;
            // sum of local errors of marked cells
            double eta = estimator2D.GetEstimatedGlobalError()[current_estimator];
            // because we compare to the sum over the eta_K^2
            eta *= eta;
            // if number of changed cells is smaller than min_changed and the sum of local errors is smaller than a
            // fraction of the global error, refine more cells
            while ((changed < min_changed) && (sum_local_errors <= fraction_of_error * eta) && (it < max_it)) {
                // reset number of changed cells
                // if this is not the first pass, the conditions have been relaxed and the cells that have been marked
                // will now be marked in particular
                changed = 0;
                // loop over all cells
                #pragma omp parallel for
                for (size_t m = 0; m < n_cells; m++) {
                    // consider m-th cell
                    cell = collection->GetCell((int) m);
                    // mark cell for refinement
                    if ((estimator2D.GetEta_K()[m] >= reftol) && (cell->GetGeoLevel() < max_cell_level)) {
                        refinements[m] = true;
                        // accumulate sum_local_errors
                        #pragma omp atomic update
                        sum_local_errors += estimator2D.GetEta_K()[m];
                        // number of changed cells increased
                        #pragma omp atomic update
                        changed++;
                    }
                    if (estimator2D.GetEta_K()[m] <= coarsetol)
                        ;   //TODO: mark cell for derefinement
                }
                if (TDatabase::ParamDB->SC_VERBOSE > 1) OutPut("total " << n_cells << " changed " << changed << " global error "<< eta << " error of refined cells " << sum_local_errors << endl);
                // if the number of marked cells is still smaller than min_changed or the accumulated errors are still
                // smaller than a presribed fraction of the global error, relax conditions
                if ((changed < min_changed) || (sum_local_errors <= fraction_of_error * eta)) {
                    reftol *= decrease_reftol_factor;
                    coarsetol *= increase_coarsetol_factor;
                    sum_local_errors = 0.0;
                }
                it++;
            }
            break;
        }
        case RefinementStrategyType::EquilibrationStrategy: {
            double sum_local_errors = 0.0;
            double eta = estimator2D.GetEstimatedGlobalError()[current_estimator];
            eta *= eta;
            while(sum_local_errors < (1-getReftol())*eta) {
                // find maximal eta_K
                double eta_max_value = -1;
#ifdef _OMP
                #pragma omp parallel
                {
                    double eta_max_private = -1;
                    #pragma omp for
                    for(size_t m = 0; m < n_cells; ++m) {
                        if(!refinements[m] && estimator2D.GetEta_K()[m] > eta_max_value && (collection->GetCell((int) m)->GetGeoLevel() < max_cell_level)) {
                            eta_max_private = estimator2D.GetEta_K()[m];
                        }
                    }
                    #pragma omp flush (eta_max_value)
                    if(eta_max_private > eta_max_value) {
                        #pragma omp critical
                        {
                            if(eta_max_private > eta_max_value) eta_max_value = eta_max_private;
                        }
                    }
                }
#else
                {
                    for (size_t m = 0; m < n_cells; ++m) {
                        if(!refinements[m] && estimator2D.GetEta_K()[m] > eta_max_value && (collection->GetCell((int) m)->GetGeoLevel() < max_cell_level)) {
                            eta_max_value = estimator2D.GetEta_K()[m];
                        }
                    }
                }
#endif
                // prevent infinite loop
                if(eta_max_value == -1) break;
                // mark all cells with value of eta_K
                #pragma omp parallel for
                for (size_t m = 0; m < n_cells; ++m) {
                    if(!refinements[m] && estimator2D.GetEta_K()[m] == eta_max_value && (collection->GetCell((int) m)->GetGeoLevel() < max_cell_level)) {
                        #pragma omp atomic update
                        sum_local_errors += eta_max_value;
                        #pragma omp atomic update
                        changed++;
                        refinements[m] = true;
                    }
                }
            }
            if (TDatabase::ParamDB->SC_VERBOSE > 1) OutPut("total " << n_cells << " changed " << changed << " global error " << eta << " error of refined cells " << sum_local_errors << endl);
            break;
        }
    }
    Output::print<1>("Cells marked for refinement: ", changed);
    if(it >= max_it && changed < min_changed) {
        Output::print<1>("Failed to mark ", int(min_changed), " cells, since parameters were relaxed it=", it," >= max_it times. Consider decreasing DECREASE_REFTOL_FACTOR.");
    }
}

bool shouldRefineCell(RefinementStrategy *strategy, size_t cell) {
    return strategy->shouldRefineCell(cell);
}

RefinementStrategyType  RefinementStrategy::getRefineStrategy() const {
    return refine_strategy;
}

void  RefinementStrategy::setRefineStrategy(RefinementStrategyType refine_strategy) {
    RefinementStrategy::refine_strategy = refine_strategy;
}

double  RefinementStrategy::getReftol() const {
    return reftol;
}

void  RefinementStrategy::setReftol(double reftol) {
    RefinementStrategy::reftol = reftol;
}

double  RefinementStrategy::getCoarsetol() const {
    return coarsetol;
}

void  RefinementStrategy::setCoarsetol(double coarsetol) {
    RefinementStrategy::coarsetol = coarsetol;
}

double  RefinementStrategy::getMinFractionToChange() const {
    return min_fraction_to_change;
}

void  RefinementStrategy::setMinFractionToChange(double min_fraction_to_change) {
    RefinementStrategy::min_fraction_to_change = min_fraction_to_change;
}

double  RefinementStrategy::getDecreaseReftolFactor() const {
    return decrease_reftol_factor;
}

void  RefinementStrategy::setDecreaseReftolFactor(double decrease_reftol_factor) {
    RefinementStrategy::decrease_reftol_factor = decrease_reftol_factor;
}

double  RefinementStrategy::getIncreaseCoarsetolFactor() const {
    return increase_coarsetol_factor;
}

void  RefinementStrategy::setIncreaseCoarsetolFactor(double increase_coarsetol_factor) {
    RefinementStrategy::increase_coarsetol_factor = increase_coarsetol_factor;
}

double  RefinementStrategy::getFractionOfError() const {
    return fraction_of_error;
}

void  RefinementStrategy::setFractionOfError(double fraction_of_error) {
    RefinementStrategy::fraction_of_error = fraction_of_error;
}

int  RefinementStrategy::getMaxCellLevel() const {
    return max_cell_level;
}

void  RefinementStrategy::setMaxCellLevel(int max_cell_level) {
    RefinementStrategy::max_cell_level = max_cell_level;
}

int  RefinementStrategy::getCurrentEstimator() const {
    return current_estimator;
}

void  RefinementStrategy::setCurrentEstimator(int current_estimator) {
    RefinementStrategy::current_estimator = current_estimator;
}

std::ostream &operator<<(std::ostream &os, RefinementStrategyType &type) {
    switch (type) {
        case RefinementStrategyType::MaximalEstimatedLocalError: {
            os << "Maximal estimated local error";
            break;
        }
        case RefinementStrategyType::PortionOfEstimatedGlobalError: {
            os << "Prescribed portion of estimated global error";
            break;
        }
        case RefinementStrategyType::EquilibrationStrategy: {
            os << "Equilibration strategy";
            break;
        }
    }
    int typeIndex = int(type);
    os << " (database value = " << typeIndex << ")";
    return os;
}
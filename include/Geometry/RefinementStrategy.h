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


#include <Database.h>
#include <MooNMD_Io.h>
#include <vector>
#include <ErrorEstimator2D.h>

enum class RefinementStrategyType {
    // this refinement strategy marks cells if the error is above some tolerance and relaxes the tolerance if
    // not enough cells have been marked yet
            MaximalEstimatedLocalError = 0,
    // this refinement strategy marks cells if the error is above some tolerance and relaxes the tolerance if
    // not enough cells have been marked yet or the accumulated error of the cells is not some prescribed
    // portion of the global error
            PortionOfEstimatedGlobalError,
    // this refinement strategy is computationally more heavy than the other strategies and aims at refining
    // as few cells as possible while sustaining uniform convergence of the adaptive process
            EquilibrationStrategy
};

// enable output of refinement strategy type
std::ostream &operator<<(std::ostream &os, RefinementStrategyType &type);

class RefinementStrategy {
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
    // this constructor calls the other constructor with the paramDB values
    RefinementStrategy();
    RefinementStrategy(int refine_strategy,                             // TDatabase::ParamDB->REFINE_STRATEGY
                       double reftol,                                   // TDatabase::ParamDB->REFTOL
                       double coarsetol,                                // TDatabase::ParamDB->COARSETOL
                       double min_fraction_to_change,                   // TDatabase::ParamDB->MIN_FRACTION_TO_CHANGE
                       double decrease_reftol_factor,                   // TDatabase::ParamDB->DECREASE_REFTOL_FACTOR
                       double increase_coarsetol_factor,                // TDatabase::ParamDB->INCREASE_COARSETOL_FACTOR
                       double fraction_of_error,                        // TDatabase::ParamDB->FRACTION_OF_ERROR
                       int max_cell_level,                              // TDatabase::ParamDB->MAX_CELL_LEVEL
                       int current_estimator)                           // TDatabase::ParamDB->ADAPTIVE_REFINEMENT_CRITERION
            : refine_strategy(RefinementStrategyType(refine_strategy)), reftol(reftol), coarsetol(coarsetol),
              min_fraction_to_change(min_fraction_to_change), decrease_reftol_factor(decrease_reftol_factor),
              increase_coarsetol_factor(increase_coarsetol_factor), fraction_of_error(fraction_of_error),
              max_cell_level(max_cell_level), current_estimator(current_estimator) {
        if (coarsetol >= reftol) this->coarsetol = 0.001 * reftol;
    }

    ~RefinementStrategy() {
    }


    virtual void applyEstimator(ErrorEstimator2D &estimator2D);

    virtual bool shouldRefineCell(size_t cellIndex);

    RefinementStrategyType getRefineStrategy() const;

    void setRefineStrategy(RefinementStrategyType refine_strategy);

    double getReftol() const;

    void setReftol(double reftol);

    double getCoarsetol() const;

    void setCoarsetol(double coarsetol);

    double getMinFractionToChange() const;

    void setMinFractionToChange(double min_fraction_to_change);

    double getDecreaseReftolFactor() const;

    void setDecreaseReftolFactor(double decrease_reftol_factor);

    double getIncreaseCoarsetolFactor() const;

    void setIncreaseCoarsetolFactor(double increase_coarsetol_factor);

    double getFractionOfError() const;

    void setFractionOfError(double fraction_of_error);

    int getMaxCellLevel() const;

    void setMaxCellLevel(int max_cell_level);

    int getCurrentEstimator() const;

    void setCurrentEstimator(int current_estimator);
};


#endif //PARMOON_REFINEMENTSTRATEGY_H

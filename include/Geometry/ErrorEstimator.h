// =============================================================================
//
// Purpose: Base class of all 2D error estimators.
//
// Author: Moritz Hoffmann, 05/07/15.
//
// =============================================================================

#ifndef INCLUDE_GEOMETRY_ERRORESTIMATOR_H
#define INCLUDE_GEOMETRY_ERRORESTIMATOR_H


class TCollection;
class TFEFunction2D;
#include "Enumerations.h"
#include "templateNames.h"
#include <vector>

template <int d>
class ErrorEstimator 
{
  public:
    using FEFunction = typename Template_names<d>::FEFunction;
    using MultiIndex_vector = typename Template_names<d>::MultiIndex_vector;
    
  protected:
  
    const TCollection* currentCollection;
    std::vector<double> eta_K;
    double maximal_local_error;
    std::vector<double> estimated_global_error;
    double estimated_afc_error;

  public:
    ErrorEstimator() {}
    virtual ~ErrorEstimator() = default;
    
    virtual void info() = 0;

    const std::vector<double> get_eta_K() const
    { return eta_K; }
    
    double get_eta_K(int cell_index) const
    { return eta_K.at(cell_index); }

    double get_maximal_local_error() const 
    { return maximal_local_error; }
    
    double get_estimated_afc_error() const 
    { return estimated_afc_error; }

    const std::vector<double> get_estimated_global_error() const
    { return estimated_global_error; }

    const TCollection* get_collection() const
    { return currentCollection; }
};


#endif // INCLUDE_GEOMETRY_ERRORESTIMATOR_H

// =============================================================================
//
// Purpose: Base class of all 2D error estimators.
//
// Author: Moritz Hoffmann, 05/07/15.
//
// =============================================================================

#ifndef ERRORESTIMATOR_2D_H
#define ERRORESTIMATOR_2D_H


#include <Example2D.h>
#include <Enumerations.h>

class TDomain;
class TDatabase;

class ErrorEstimator2D {
protected:
    Example2D& example2D;
    TDomain& domain;
    double *eta_K = nullptr, maximal_local_error;
    std::vector<double> estimated_global_error;

public:
    ErrorEstimator2D(Example2D &ex, TDomain& domain) : example2D(ex), domain(domain) {}
    ~ErrorEstimator2D() {
        if(eta_K && eta_K != nullptr) delete[] eta_K;
    }

    virtual void estimate(const std::vector<MultiIndex2D> &derivatives, const TFEFunction2D &fe_function2D) = 0;

    TDomain &GetDomain() const {
        return this->domain;
    }

    double *GetEta_K() const {
        return eta_K;
    }

    double GetMaximalLocalError() const {
        return maximal_local_error;
    }


    const std::vector<double> GetEstimatedGlobalError() const {
        return estimated_global_error;
    }
};


#endif //ERRORESTIMATOR_2D_H

#ifndef LOCALASSEMBLINGPR_H
#define LOCALASSEMBLINGPR_H

#include "../include/AssembleRoutines/LocalAssembling.h"

/**
 * @todo write docs
 */

class LocalAssemblingPR : public LocalAssembling<2>
{
public:
    /**
     * Default constructor
     */
    LocalAssemblingPR(ParameterDatabase db, CoeffFct2D coeffs, std::string name);

};

#endif // LOCALASSEMBLINGPR_H

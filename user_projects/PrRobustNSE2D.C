#include <PrRobustNSE2D.h>
#include <NSE2D.h>
#include <Domain.h>

PrRobustNSE2D::SystemPerGrid::SystemPerGrid(const Example_NSE2D& example, 
TCollection& coll, unsigned int order)
: projection_space(&coll, (char*)"u", (char*)"Darcy velocity", example.get_bc(0),
                  order, nullptr),
 matrix(&projection_space, &projection_space)
{
  
}

PrRobustNSE2D::PrRobustNSE2D(const TDomain& domain, const Example_NSE2D& _example,
                             unsigned int reference_id)
 : NSE2D(domain, _example, reference_id), Systems()
{
}




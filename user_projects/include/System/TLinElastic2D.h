/** ************************************************************************
 * 
 * @name         TLinElastic2D
 * @brief        store everything needed to solve a 2D linear elastic
 *               problem
 *               
 * @author       Najib Alia
 * @History      16.03.2017
**************************************************************************/

#ifndef __TLinElastic2D__
#define __TLinElastic2D__

#include <BlockFEMatrix.h>
#include <BlockVector.h>
#include <FESpace2D.h>
#include <FEFunction2D.h>

#include <MainUtilities.h>
#include <ParameterDatabase.h>
#include <vector>
#include <deque>
#include <utility>
#include <array>

class TLinElastic2D
{
  public:
    int test = 1;

    /** @brief constructor*/
    TLinElastic2D();

    /*************************************************************/
   /**
    * Special member functions mostly deleted
    * ...needs to be optimized
    */
   //! Delete copy constructor.
    TLinElastic2D(const TLinElastic2D&) = delete;

   //! Delete move constructor.
    TLinElastic2D(TLinElastic2D&&) = delete;

   //! Delete copy assignment operator.
    TLinElastic2D& operator=(const TLinElastic2D&) = delete;

   //! Delete move assignment operator.
    TLinElastic2D& operator=(TLinElastic2D&&) = delete;

   //! Default destructor. Most likely causes memory leaks.
   ~TLinElastic2D() = default;

};

#endif // __TLinElastic2D__

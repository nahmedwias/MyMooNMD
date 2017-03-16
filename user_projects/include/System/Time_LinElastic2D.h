/** ************************************************************************
 * 
 * @name         Time_LinElastic2D
 * @brief        store everything needed to solve a 2D linear elastic
 *               problem
 *               
 * @author       Najib Alia
 * @History      16.03.2017
**************************************************************************/

#ifndef __Time_LinElastic2D__
#define __Time_LinElastic2D__

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

class Time_LinElastic2D
{
  public:
    int test = 1;

    /** @brief constructor*/
    Time_LinElastic2D();

    /*************************************************************/
   /**
    * Special member functions mostly deleted
    * ...needs to be optimized
    */
   //! Delete copy constructor.
    Time_LinElastic2D(const Time_LinElastic2D&) = delete;

   //! Delete move constructor.
    Time_LinElastic2D(Time_LinElastic2D&&) = delete;

   //! Delete copy assignment operator.
    Time_LinElastic2D& operator=(const Time_LinElastic2D&) = delete;

   //! Delete move assignment operator.
    Time_LinElastic2D& operator=(Time_LinElastic2D&&) = delete;

   //! Default destructor. Most likely causes memory leaks.
   ~Time_LinElastic2D() = default;

};

#endif // __Time_LinElastic2D__

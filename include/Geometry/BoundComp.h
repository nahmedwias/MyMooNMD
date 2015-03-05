// =======================================================================
// @(#)BoundComp.h        1.5 08/12/99
//
// Class:       TBoundComp
// Purpose:     components of boundary faces
//
// Author:      Volker Behns  18.06.97
//
// =======================================================================

#ifndef __BOUNDCOMP__
#define __BOUNDCOMP__

#include <Constants.h>
#include <fstream>

#define TOL_SECANT_BOUND 1.01
#define MAXVERTS 20
#define RECURS_DEPTH 8

enum BoundTypes {Line = 1, Circle, Spline2D, Polygon, NonUniformSpline2D,
                 Plane = 10, Sphere, Cylinder, 
                 Wall = 100,
                 NoPRM = 4711 };

/** components of boundary faces */
class TBoundComp
{
  protected:
    /** component identifier */
    int ID;
    /** type of component */
    BoundTypes Type;
    /** TRUE if component is on free boundary */
    bool FreeBoundaryStatus;
    
  public:
    // Constructor
    TBoundComp(int id);

    // Methods
    /** read parameter from input file */
    virtual int ReadIn(std::ifstream &dat) = 0;

    /** return ID */
    int GetID()
    { return ID; }

    /** get type of component */
    BoundTypes GetType()
    { return Type; }

    /** get free boundary status */
    bool IsFreeBoundary()
    { return FreeBoundaryStatus; }

    /** set free boundary status */
    void SetFreeBoundaryStatus(bool status)
    { FreeBoundaryStatus = status; }
    
    void ChangeType(BoundTypes New_Type)
     { Type = New_Type; }    
};

#endif

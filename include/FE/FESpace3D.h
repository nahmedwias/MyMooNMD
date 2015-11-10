// =======================================================================
// %W% %G%
// 
// Class:       TFESpace3D
// Purpose:     class for all 3D finite element spaces
//
// Author:      Gunar Matthies (22.11.97)
//
// History:     start of implementation 22.11.97 (Gunar Matthies)
//
// =======================================================================

#ifndef __FESPACE3D__
#define __FESPACE3D__

#include <FESpace.h>
#include <FE3D.h>

class THangingNode;

/** class for all 3D finite element spaces */
class TFESpace3D : public TFESpace
{
  protected:
    /** number of active degrees of freedom */
    int N_ActiveDegrees;

    /** number of slave degrees of freedom */
    int N_SlaveDegrees;

    /** NeumannBound <= i < HangingBound for all hanging nodes i */
    /** => HangingBound <= i < DirichletBound for all Dirichlet node i */
    int HangingBound;

    /** array containing the used elements */
    FE3D *UsedElements; 

    /** array with an element for each shape */
    FE3D *ElementForShape;

    /** array of hanging nodes */
    THangingNode **HangingNodeArray;

    /** array storing the fe for each element, if necessary */
    FE3D *AllElements;

    /**
     *  @brief Boundary condition used to create this space.
     *
     *  In order to reproduce the information, which boundary condition was used
     *  in the creation of the space, store a function pointer to it.
     */
   BoundCondFunct3D* boundCondition_;

  public:
    /** constructor */
    TFESpace3D(TCollection *coll, char *name, char *description);

    /** constructor for building a space with elements of order k */
    TFESpace3D(TCollection *coll, char *name, char *description, 
               BoundCondFunct3D *BoundaryCondition, int k);

    /** constructor for building a space with the given elements */
    TFESpace3D(TCollection *coll, char *name, char *description,
               BoundCondFunct3D *BoundaryCondition,
               FE3D *fes);

    TFESpace3D(TCollection *coll, char *name, char *description, 
               BoundCondFunct3D *BoundaryCondition, SpaceType type,
               int ord);

    /** destructor */
    ~TFESpace3D();

    /** find used elements */
    void FindUsedElements();

    /** construct space */
    void ConstructSpace(BoundCondFunct3D *BoundaryCondition);

   /** @return The boundary condition function pointer. */
   BoundCondFunct3D* getBoundCondition() const
   { return boundCondition_; }

    /** return number of active degrees of freedom */
    int GetN_ActiveDegrees() const
    { return N_ActiveDegrees; }

    /** return number of slave degrees of freedom */
    int GetN_SlaveDegrees() const
    { return N_SlaveDegrees; }

    /** return HangingBound */
    int GetHangingBound() const
    { return HangingBound; }

    /** return N_Hanging=N_SlaveDegrees */
    int GetN_Hanging() const
    { return N_SlaveDegrees; }

    /** return identifiers of used elements */
    FE3D *GetUsedElements() const
    { return UsedElements; }

    /** return array with all hanging nodes */
    THangingNode **GetHangingNodes() const
    { return HangingNodeArray; }

    /** return the FE Id for element i, corresponding to cell */
    FE3D GetFE3D(int i, TBaseCell *cell) const;

    /** return position of one given DOF */
    void GetDOFPosition(int dof, double &x, double &y, double &z);

    /** return position of all dofs */
    void GetDOFPosition(double *x, double *y, double *z);


};

#endif

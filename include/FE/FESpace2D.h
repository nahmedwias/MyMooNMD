// =======================================================================
// @(#)FESpace2D.h        1.5 06/13/00
// 
// Class:       TFESpace2D
// Purpose:     class for all 2D finite element spaces
//
// Author:      Gunar Matthies (03.11.97)
//
// History:     start of implementation 03.11.97 (Gunar Matthies)
//
//              split FESpace into TFESpacexD
//              15.04.1998 (Volker Behns)
//
//              start of reimplementation
//              30.07.98 (Gunar Matthies)
//
// =======================================================================

#ifndef __FESPACE2D__
#define __FESPACE2D__

#include <FESpace.h>
#include <FE2D.h>

class THangingNode;

/** class for all 2D finite element spaces */
class TFESpace2D : public TFESpace
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
    FE2D *UsedElements; 

    /** array with an element for each shape */
    FE2D *ElementForShape;

    /** array of hanging nodes */
    THangingNode **HangingNodeArray;

    /** array storing the fe for each element, if necessary */
    FE2D *AllElements;
    
    /** boundary condition used to create this space */
   BoundCondFunct2D *BoundCondition;

   /** indices for mapping between Nodalfunctional/Nodal-interpolation point 
      in operator-splitting methods ---  Sashikumaar Ganesan */
   int *IntlPtIndexOfPts;


  public:
    /** constructor */
    TFESpace2D(TCollection *coll, std::string name, std::string description);

    /** constructor for building a space with elements of order k */
    TFESpace2D(TCollection *coll, std::string name, std::string description,
               BoundCondFunct2D *BoundaryCondition, int k);

    TFESpace2D(TCollection *coll, std::string name, std::string description,
               BoundCondFunct2D *BoundaryCondition, SpaceType type, int k);

    /** constructor for building a space with the given elements */
    TFESpace2D(TCollection *coll, std::string name, std::string description,
               BoundCondFunct2D *BoundaryCondition, FE2D *fes);

    TFESpace2D(const  TFESpace2D&)=delete;
    TFESpace2D& operator=(TFESpace2D) = delete;
    /** destructor */
    ~TFESpace2D();

    /** find used elements */
    void FindUsedElements();

    /** construct space */
    void ConstructSpace(BoundCondFunct2D *BoundaryCondition);
    
    /** @brief get dimension of the vector basis function */
    virtual int GetBaseVectDim() const;  

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
    FE2D *GetUsedElements() const
    { return UsedElements; }

    /** return array with all hanging nodes */
    THangingNode **GetHangingNodes() const
    { return HangingNodeArray; }

    /** return the FE Id for element i, corresponding to cell */
    FE2D GetFE2D( int i, const TBaseCell* cell ) const;
    
    /** @brief return the Finite Element on a given cell */
    const TFE2D& get_fe(unsigned int cell_number) const;

    /** return position of one given DOF */
    void GetDOFPosition(int dof, double &x, double &y) const;

    /** return position of all dofs */
    void GetDOFPosition(double *x, double *y);

    void SetIntlPtIndexOfPts(int *intlPtIndexOfPts)
     { IntlPtIndexOfPts = intlPtIndexOfPts; }

    int *GetIntlPtIndexOfPts() const
     { return IntlPtIndexOfPts; }
     
     FE2D *GetAllElements() const
     { return AllElements; }
     
     /**
      * @brief print some information on this fe space
      * 
      * Depending on the verbosity level more or less is printed.
      */
     void info() const;
     /** return boundary condition */
    BoundCondFunct2D *get_boundary_condition() const
    { return BoundCondition; }
    
    friend  bool operator== (const TFESpace2D &lhs, const TFESpace2D &rhs);
    friend  bool operator!= (const TFESpace2D &lhs, const TFESpace2D &rhs);
};


#endif

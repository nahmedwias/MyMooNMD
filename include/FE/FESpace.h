// =======================================================================
// @(#)FESpace.h        1.1 10/30/98
// 
// Class:       TFESpace
// Purpose:     general super class for all finite element spaces
//              special spaces are implemented in subclasses
//
// Author:      Gunar Matthies (03.11.97)
//
// History:     start of implementation 03.11.97 (Gunar Matthies)
//
//              add data for storing information on the used elements
//              29.11.1997 (Gunar Matthies)
//
//              split FESpace into TFESpacexD
//              15.04.1998 (Volker Behns)
//
//              start of reimplementation
//              30.07.1998 (Gunar Matthies)
//
// =======================================================================

#ifndef __FESPACE__
#define __FESPACE__

#include <Collection.h>
#include <string.h>

/** general super class for all finite element spaces, special spaces are
    implemented in subclasses */
class TFESpace
{
  protected:
// =======================================================================
// administration
// =======================================================================
    /** name of the space */
    std::string Name;

    /** some more words describing the space */
    std::string Description;

// =======================================================================
// information of cell collection
// =======================================================================
    /** Collection containing the cells used for building this space */
    TCollection *Collection;

    /** number of cells in the triangulation used for building this space */
    int N_Cells;

// =======================================================================
// general information on degrees of freedom and their numbers
// =======================================================================
    /** number of all degrees of freedom */
    int N_DegreesOfFreedom;

    /** array containing all global numbers of local degrees of freedom
        for all elements */
    int *GlobalNumbers;

    /** array containing the begin index in GlobalNumbers array for each
        element */
    int *BeginIndex;

    /** number of used elements */
    int N_UsedElements;

// =======================================================================
// counts and bounds for different types of degrees of freedom
// =======================================================================
    /** number of different boundary node type, except Dirichlet */
    int N_DiffBoundNodeTypes;

    /** type for each of the different types */
    BoundCond *BoundaryNodeTypes;

    /** number of Dirichlet nodes */
    int N_Dirichlet;

    /** number of nodes for each boundary node type */
    int *N_BoundaryNodes;

    /** number of inner nodes */
    int N_Inner;

    /** 0 <= i < InnerBound for all inner degrees of freedom i */
    int InnerBound;

    /** InnerBound <= i < NeumannBound for all Neumann nodes i */
    int *BoundaryNodesBound;

    /** NeumannBound <= i < DirichletBound for all Dirichlet node i */
    int DirichletBound;

    /** number of inner and non-Dirichlet boundary nodes are less than */
    int ActiveBound;

    /// True if this space contains discontinuous elements.
    bool is_discontinuous_galerkin_space;

  public:
    /**
     * Constructor, setting name, description, cellgrid, and a lot of
     * dummy variables. Note that FESpace is more like an interface class,
     * only its daughter classes (with fixed dimension) are usable.
     *
     * @param[in] coll The cell grid (i.e. the finite element mesh).
     * @param[in] name The name of the space, used in printout etc.
     * @param[in] description A description of the space, used in printout etc.
     */
    TFESpace(TCollection *coll, std::string name, std::string description);

    /** destrcutor */
    virtual ~TFESpace();

    /** return name */
    std::string GetName() const
    { return Name; }

    /** return description */
    std::string GetDescription() const
    { return Description; }

    /** return number of cells in the triangulation used for building 
        this space */
    int GetN_Cells() const
    { return N_Cells; }

    /** return the collection of this space */
    TCollection *GetCollection() const
    { return Collection; }
    
    /** @brief get dimension of the vector basis function */
    virtual int GetBaseVectDim() const = 0;

    /** return global numbers of local degrees of freedom */
    int *GetGlobalNumbers() const
    { return GlobalNumbers; }
    
    void SetGlobalNumbers(int* NewGN)
    { GlobalNumbers=NewGN; }

    /** return begin index for each element */
    int *GetBeginIndex() const
    { return BeginIndex; }
    
    /** @brief return correspondence map from local to global degrees of freedom
     * 
     * set int * DOF=feSpace->GetGlobalDOF(i); then DOF[j]-th global degree of 
     * freedom corresponds to the j-th local degree of freedom.
    */
    int* GetGlobalDOF(int i) const
    { return GlobalNumbers+BeginIndex[i];}

    /** return number of used elements */
    int GetN_UsedElements() const
    { return N_UsedElements; }

    /** return number of all degrees of freedom */
    int GetN_DegreesOfFreedom() const
    { return N_DegreesOfFreedom; }

// =======================================================================
// counts and bounds for different types of degrees of freedom
// =======================================================================
    /** get number of different boundary node types, except Dirichlet */
    int GetN_DiffBoundaryNodeTypes() const
    { return N_DiffBoundNodeTypes; }

    /** return type for each of the different types */
    BoundCond *GetBoundaryNodeTypes() const
    { return BoundaryNodeTypes; }

    /** return number of nodes for each boundary node type */
    int *GetN_BoundaryNodes() const
    { return N_BoundaryNodes; }

    /** return N_Dirichlet */
    int GetN_Dirichlet() const
    { return N_Dirichlet; }

    /** @return The number of Neumann boundary d.o.f. */
    int get_n_neumann_dof() const
    {
      return N_BoundaryNodes[0];
    }

    /** @return The number of Neumann boundary d.o.f. */
    int get_n_robin_dof() const
    { return N_BoundaryNodes[1]; }

    /** return N_Inner */
    int GetN_Inner() const
    { return N_Inner; }

    /** return InnerBound */
    int GetInnerBound() const
    { return InnerBound; }

    /** return BoundaryNodesBound */
    int *GetBoundaryNodesBound() const
    { return BoundaryNodesBound; }

    /** return DirichletBound */
    int GetDirichletBound() const
    { return DirichletBound; }

    /** return ActiveBound */
    int GetActiveBound() const
    { return ActiveBound; }

    /** write info on fespace into file */
    int Write(const std::string filename);

    void SetAsDGSpace() { is_discontinuous_galerkin_space = true; }

    bool IsDGSpace() const { return is_discontinuous_galerkin_space; }


};

#endif

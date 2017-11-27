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

#include <memory>
#include <FESpace.h>
#include <FE3D.h>

class THangingNode;

#ifdef _MPI
class TParFECommunicator3D;
class TParFEMapper3D;
#endif

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

#ifdef  _MPI
     /** Maximum number of subdomains associated with any dof */
    int MaxSubDomainPerDof;

   /// There belongs a ParFECommunicator to this space. Store it!
   std::shared_ptr<TParFECommunicator3D> comm_;

   ///The communicator needs a mapper TODO Move mapper into comm as member!
   std::shared_ptr<TParFEMapper3D> mapper_;

#endif

  public:
    /** constructor */
    TFESpace3D(TCollection *coll, std::string name, std::string description);

    /** constructor for building a space with elements of order k */
    TFESpace3D(TCollection *coll, std::string name, std::string description,
               BoundCondFunct3D *BoundaryCondition, int k);

    /** constructor for building a space with the given elements */
    TFESpace3D(TCollection *coll, std::string name, std::string description,
               BoundCondFunct3D *BoundaryCondition,
               FE3D *fes);

    TFESpace3D(TCollection *coll, std::string name, std::string description,
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
    
    /** @brief return the Finite Element on a given cell */
    const TFE3D& get_fe(unsigned int cell_number) const;


    /** return position of one given DOF */
    void GetDOFPosition(int dof, double &x, double &y, double &z) const;

    /** return position of all dofs */
    void GetDOFPosition(double *x, double *y, double *z) const;

    /**
     * velocity space, if there is a element that only has dirichlet dof's??
     */
    bool CheckMesh() const;

#ifdef  _MPI
    /**
     * As soon as the maxSubDomainPerDof is known, the parallel infrastructure
     * can be initialized.
     */
    void initialize_parallel(int maxSubDomainPerDof);

    /** return  MaxSubDomainPerDof */
    int GetMaxSubDomainPerDof()
    { return MaxSubDomainPerDof; }

    const TParFECommunicator3D& get_communicator() const
    { return *comm_; }

    TParFECommunicator3D& get_communicator()
    { return *comm_; }

    const TParFEMapper3D& get_mapper() const
    { return *mapper_; }

    TParFEMapper3D& get_mapper()
    { return *mapper_; }
#endif

};

#endif

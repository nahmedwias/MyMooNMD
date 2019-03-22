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
   /// There belongs a ParFECommunicator to this space. Store it!
   std::shared_ptr<TParFECommunicator3D> comm_;

   ///The communicator needs a mapper TODO Move mapper into comm as member!
   std::shared_ptr<TParFEMapper3D> mapper_;

#endif

  public:
    /** constructor */
    TFESpace3D(TCollection *coll, const std::string& name,
               const std::string& description);

    /** constructor for building a space with elements of order k */
    TFESpace3D(TCollection *coll, const std::string& name,
               const std::string& description,
               BoundCondFunct3D *BoundaryCondition, int k);

    /** constructor for building a space with the given elements */
    TFESpace3D(TCollection *coll, const std::string& name,
               const std::string& description,
               BoundCondFunct3D *BoundaryCondition, FE3D *fes);

    TFESpace3D(TCollection *coll, const std::string& name,
               const std::string& description,
               BoundCondFunct3D *BoundaryCondition, SpaceType type, int ord);

    /** destructor */
    ~TFESpace3D();

    /** find used elements */
    void FindUsedElements();

    /** construct space */
    void ConstructSpace(BoundCondFunct3D *BoundaryCondition);
        
    /** @brief get dimension of the vector basis function */
    virtual int GetBaseVectDim() const; 

   /** @return The boundary condition function pointer. */
   BoundCondFunct3D* get_boundary_condition() const
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
    FE3D GetFE3D(int i, const TBaseCell *cell) const;
    
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

    /**
     * @brief get degree of FE space on a given cell
     **/
    int getFEDegree(TBaseCell *cell) const;
      
    /**
     * @brief get quadrature data on the m-th face of a given cell
    **/
    void getFaceQuadratureData(TBaseCell *cell, int m,
			       std::vector<double>& qWeights,std::vector<double>& qPointsT,
			       std::vector<double>& qPointsS,
			       std::vector< std::vector<double> >& basisFunctionsValues,
			       int _deg = -1) const;

    /**
     * @brief get the right quadrature formula for a given space on m-th face of a given cell
    **/
    QuadFormula2D getFaceQuadratureFormula(TBaseCell *cell, int m, int d=-1) const;
      
    /**
     * @brief get quadrature value on the m-th face of a given cell for a given quad formula
    **/
    void getFaceQuadratureValue(TBaseCell *cell, int m, QuadFormula2D FaceQuadFormula,
				std::vector< std::vector<double> >& basisFunctionsValues) const;

#ifdef  _MPI
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

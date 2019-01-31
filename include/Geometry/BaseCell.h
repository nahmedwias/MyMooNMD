/** ************************************************************************ 
*
* @class      TBaseCell
* @brief      prototype of a cell
*

* @author     Volker Behns (97) & Gunar Matthies & Sashikumaar Ganesan
* @date       09.07.97
* @History    Making some methods inline (Gunar Matthies,  10.10.97)
              Adding ClipBoard (Gunar Matthies, 10.10.97)
              Twophase flows (Sashikumaar Ganesan. 20.09.09)
              Added edge info in 3D (Sashikumaar Ganesan, 04.09.10)
              MPI methods (Sashikumaar Ganesan, 08.08.14)
************************************************************************  */

#ifndef __BASECELL__
#define __BASECELL__

class TBaseCell;
class TDomain;

#include <vector>
#include <Edge.h>
#include <Joint.h>
#include <RefDesc.h>

 /**  @brief information for finite element data structure */
class TBaseCell
{
  protected:
    /**  @brief current property of refinement (including shape descriptor) */
    TRefDesc *RefDesc;

    /**  @brief array of all joints */
    TJoint **Joints;

    /**  @brief an integer for storing clipboard information. Bad style!*/
    mutable int ClipBoard;

    /** @brief an integer for storing physical reference of cell (e.g. material properties) **/
    int Reference_ID;
    
    /**  @brief an integer for storing which phase(multiphase domains) this cell belongs*/
    int Phase_ID;

    /**  @brief an integer for storing boundary part (surface meshes) */
    int Bd_Part;

    /**  @brief cell index value in the collection */
    int CellIndex;

  /**  @brief an integer for storing the global cell number **/
    int GlobalCellNo;
    
    /** @brief normal orientation of joints (for each joint this is 1 or -1) */
    int *normalOrientation;
    
  /**  @brief an integer for storing the region of this cell number **/
    int region;

  /**  @brief an integer for indicating layer cells **/
    int LayerCell;
    
#ifdef __3D__
    /**  @brief array of all Edges in 3D */
    TEdge **Edges;
#endif

#ifdef  _MPI
  /**  @brief an integer for storing clipboard information in parallel FEspace 
   * mapping.
   * Bad style!
   */
  mutable int ClipBoard_Par;
    
  /**  @brief an integer for storing which subdomain contains this cell **/
    int SubDomainNumber;

  /**  @brief an integer for storing the global cell number **/
    int SubDomainLocalCellNo;

  /**  @brief a bool to check this cell is own cell of the SubDomain or not **/
    bool OwnCell;

  /**  @brief a bool to check this cell is a halo cell for the SubDomain or not **/
    bool HaloCell;

  /**  @brief a bool to check this cell is neibs' halo cell or not **/
    bool DependentCell;
    
  /**  @brief a bool to check this cell contains SubDomainInterface(s) or not **/
    bool SubDomainInterfaceCell;

  /**  @brief a bool to check this cell contains  cross edges or not **/
    bool CrossEdgeCell;
    
  /**  @brief a bool to check this cell contains cross vertices or not **/
    bool CrossVertexCell;

  /**  @brief Number of neib processes for this cell **/
    int N_NeibProcesses;

  /**  @brief if(N_NeibProcesses), the rank ID of neib processes **/
    int *NeibProcessesIds;
 
#endif

  /**  @brief an integer for storing the local cell number (2Phase flows) **/
    int LocalCellNo;   
    
    virtual TBaseCell *GetChild(int C_i) = 0;
    virtual TBaseCell *GetParent() = 0;
    
    friend class TDomain;
    
  public:
    // Constructor
    explicit TBaseCell(TRefDesc *refdesc);

    // Destructor
    virtual ~TBaseCell();

    // Methods
    /**  @brief set refinement descriptor to newrefdesc */
    int SetRefDesc(TRefDesc *newrefdesc)
    {
      if (newrefdesc->GetShapeDesc() == GetShapeDesc())
      {
        RefDesc = newrefdesc;
        return 0;
      }
      else
        if (newrefdesc->GetShapeDesc()->GetN_Vertices() == 
              GetShapeDesc()->GetN_Vertices())
        {
          RefDesc = newrefdesc;
          return 0;
        }
        else
          return -1;
    }

    /**  @brief return refinement descriptor */
    const TRefDesc *GetRefDesc() const
    { return RefDesc; }
    /**  @brief return shape descriptor of refinement descriptor */
    const TShapeDesc *GetShapeDesc() const
    { return RefDesc->GetShapeDesc(); }
    /**  @brief return shape type of refinement descriptor */
    Shapes GetType() const
    { return RefDesc->GetShapeDesc()->GetType(); }

    /**  @brief return refinement descriptor of edge i */
    Refinements GetEdgeRef(int i) const
    { return RefDesc->GetEdgeRef(i); }

#ifdef __3D__
    /**  @brief return refinement descriptor of face i */
    Refinements GetFaceRef(int i) const
    { return RefDesc->GetFaceRef(i); }

    /**  @brief set the pointer to edge E_i to E */
    int SetEdge(int E_i, TEdge *E)
    {
      Edges[E_i] = E;
      return 0;
    }

    /**  @brief return the pointer to edge with number E_i */
    TEdge *GetEdge(int E_i)
    { return Edges[E_i]; }
    const TEdge *GetEdge(int E_i) const
    { return Edges[E_i]; }

#endif

    /**  @brief set the pointer of vertex Vert\_i to Vert */
    virtual int SetVertex(int Vert_i, TVertex *Vert) = 0;
    /**  @brief return the pointer to vertex with number i */
    virtual TVertex *GetVertex(int Vert_i) = 0;
    virtual const TVertex *GetVertex(int Vert_i) const = 0;

    /**  @brief set the pointer to face J\_i to J */
    int SetJoint(int J_i, TJoint *J)
    {
      Joints[J_i] = J;
      return 0;
    }

    /**  @brief return the pointer to face with number i */
    TJoint *GetJoint(int J_i)
    { return Joints[J_i]; }
    const TJoint *GetJoint(int J_i) const
    { return Joints[J_i]; }

    /**  @brief return the number of vertices of the cell */
    int GetN_Vertices() const
    { return RefDesc->GetN_OrigVertices(); }
    /**  @brief return the number of edges of the cell */
    int GetN_Edges() const
    {  return RefDesc->GetN_OrigEdges(); }
    /**  @brief return the number of joints */
    int GetN_Joints() const
    {  return RefDesc->GetShapeDesc()->GetN_Joints(); }

    #ifdef __3D__
      /**  @brief return the number of faces of the cell */
      int GetN_Faces() const
      { return RefDesc->GetN_OrigFaces(); }
    #endif

    /**  @brief return the number of children of the cell */
    virtual int GetN_Children() const = 0;
    /**  @brief return the number of parents of the cell */
    virtual int GetN_Parents() const = 0;

    /**  @brief return the child with the number C\_i */
    virtual const TBaseCell *GetChild(int C_i) const = 0;
    /**  @brief return the parent cell */
    virtual const TBaseCell *GetParent() const = 0;
    /**  @brief set the parent to parent */
    virtual int SetParent(TBaseCell *parent) = 0;
    /**  @brief return the child number of cell Me */
    virtual int GetChildNumber(TBaseCell *Me) const =  0;

    /**  @brief write the postscript cell data to stream dat.
     * A positive cell_index is written to the postscript file, negative values
     * are ignored.
     */
    virtual void PS(std::ofstream &dat, double scale, double StartX,
                    double StartY, int cell_index = -1) const = 0;
    
    /**  @brief refine the current cell on level RefLevel according actual
        refinement descriptor */
    virtual int Refine(int RefLevel) = 0;

    /**  @brief derefine the current cell, remove the children */
    virtual int Derefine() = 0;

    /**  @brief make refinement or derefinement according to cell's clipboard */
    virtual int RefDeref() = 0;
    /**  @brief set marks in neighbour cell in order to maintain 1 regularity */
    virtual int Gen1RegMarks() = 0;

    /**  @brief generate such refinement information in the neighbouring cells
        so that the traingulation will become 1-regular */
    virtual int Gen1RegGrid() = 0;
    /**  @brief regular refinement of joint LocJointNum */
    virtual int Ref1Reg(int LocJointNum, TBaseCell *&RefCell) = 0;
    /**  @brief check whether the surroundings of cell is 1-regular */
    virtual int Check1Reg() = 0;
    /**  @brief set refinement descriptor to no refinement */
    virtual int SetNoRefinement() = 0;
    /**  @brief set refinement descriptor to regular refinement */
    virtual int SetRegRefine() = 0;
    /**  @brief set refinement descriptor to adaptive refinement */
    virtual int Set1Refine(int i)= 0;    
    /**  @brief is the cell to refine */
    virtual int IsToRefine() const = 0;
    /**  @brief are there any children of this cell */
    virtual int ExistChildren() const = 0;
    /**  @brief generate conforming closures */
    virtual int MakeConfClosure() = 0;

    #ifdef __2D__
      /**  @brief return (x,y) coordinates */
      virtual int LineMidXY(int J_i, int P_j, double &X, double &Y) = 0;
      /**  @brief return parameter values */
      virtual int LineMidT(int J_i, int SJ_j, double &T_0, double &T_1) = 0;
    #else
    #endif

    /**  @brief set value in ClipBoard. Don't rely on this! */
    void SetClipBoard(int value) const
    { ClipBoard=value; }
    /**  @brief get value from ClipBoard */
    int GetClipBoard() const
    { return ClipBoard; }

    /**  @brief get diameter of a cell */
    virtual double GetDiameter() const = 0;

    /**  @brief return shortest edge of a cell */
    virtual double GetShortestEdge() const = 0;

    /**  @brief return the length of the cell defined with the reference map */
    virtual double GetLengthWithReferenceMap() const = 0;

     /**  @brief get measure of a cell */
    virtual double GetMeasure() const = 0;
    
    /** @brief get the value of hK 
     * 
     * This function calls either this->GetDiameter(), thid->GetShortestEdge, or
     * this->GetMeasure(), depending on cell_measure.
     * 
     * Typically you should set cell_measure = TDatabase::ParamDB->CELL_MEASURE.
     */
    double Get_hK(int cell_measure) const;

#ifdef __3D__
    /**  @brief check if the line define by the position P and the direction is
     * intersecting this cell.
     * If so, lmin and lmax are the intersecting points */     
    virtual bool IsLineCutingCell(int                  direction,
                                  std::array<double,3> P,
                                  double&              lmin,
                                  double&              lmax) const = 0;
#endif /** #ifdef __3D__ */

    /**  @brief return whether a point is inside a cell */
    virtual bool PointInCell(double X, double Y) const = 0;

#ifdef __3D__
    virtual bool PointInCell(double X, double Y, double Z) const = 0;

     // added 25.04.2010 for fixing refinement problem
     void CorrectBoundaryVertices(TVertex **NewVertices, TJoint **NewJoints);
#endif

    /**  @brief get geometry level */
    virtual int GetGeoLevel() = 0;

    /**  @brief get subgrid ID */
    virtual int GetSubGridID() const = 0;

    /** @brief set reference number to this cell   */
    void SetReference_ID(int val)
    { Reference_ID = val;}
    
    /** @brief get reference number of this cell   */
    int GetReference_ID() const
    {return Reference_ID;}
    
    /**  @brief set phase number to this cell   */
    void SetPhase_ID(int val)
       {Phase_ID = val;}

    /**  @brief get phase number to this cell   */
    int GetPhase_ID() const
       {return Phase_ID;}

    /**  @brief set phase number to this cell   */
    void SetBd_Part(int val)
       {Bd_Part = val;}

    /**  @brief get phase number to this cell   */
    int GetBd_Part() const
       {return Bd_Part;}

    /**  @brief set subdomain number to this cell   */
    void SetCellIndex(int val)
       {CellIndex = val;}

    /**  @brief set subdomain number to this cell   */
    int GetCellIndex() const
       {return CellIndex;}

    /**
     * @brief Set a global cell number to this cell.
     *
     * This must be consistent over all processors, and no excuses.
     */
    void SetGlobalCellNo(int val)
       {GlobalCellNo = val;}

    /**
     * @brief Get the global cell number of this cell.
     *
     *  The GlobalCellNo being consistent on all processors is crucial
     *  for a functioning domain decomposition method.
     */
    int GetGlobalCellNo() const
       {return GlobalCellNo;}

    /**  @brief set subdomain number to this cell   */
    void SetLocalCellNo(int val)
       {LocalCellNo = val;}

    /**  @brief set LocalCellNo number to this cell   */
    int GetLocalCellNo() const 
       {return LocalCellNo;}

    /**  @brief set region number to this cell   */
    void SetRegionID(int val)
       {region = val;}

    /**  @brief set region number to this cell   */
    int GetRegionID() const
       {return region;}       

     /**  @brief set as LayerCell cell   */
    void SetAsLayerCell(int val)
       {LayerCell = val;}

    /**  @brief get LayerCell info  */
    int IsLayerCell() const
       {return LayerCell;}          

    /** @brief compute normal orientation w.r.t cell */
    void SetNormalOrientation();
    
    /** @brief get normal orientation w.r.t cell at i-th joint */
    int GetNormalOrientation(int i) const
    { return normalOrientation[i]; }
    
    //LB ====================================================
    /** @brief test if current cell is a boundary cell on BoundComp id (optional) */
    bool IsBoundaryCell( int BoundComp_id = -4711) const;

    
//    bool Get_IsBoundaryCellOnBoundComp()
//    {
//        return IsBoundaryCellOnBoundComp;
//    }
    //LB ====================================================

#ifdef __3D__

    /** @brief 
        compute normal to m-th face and transformation (for face integral) 
        @attention The normal to boundary faces might not be directed outward. For this, use
        the function get_normal_vector for TBoundFace objects.  
     */
    void computeNormalAndTransformationData(int m,
                                            std::vector<double>& normal,
                                            double &transformationDeterminant) const;
    
    /** @brief get number of vertices on m-th face */
    int getNumberOfFaceVertices(int m) const
    {
      const int *faceVertexMap, *faceVertexMapLength;
      int maxNVerticesPerFace;
      GetShapeDesc()->GetFaceVertex(faceVertexMap,faceVertexMapLength,maxNVerticesPerFace);
      // simplify: number of vertices on face m
      return faceVertexMapLength[ m ];
    }
#endif
    
#ifdef  _MPI

    /**  @brief set value in ClipBoard */
    void SetClipBoard_Par(int value) const
    { ClipBoard_Par=value; }
    
    /**  @brief get value from ClipBoard */
    int GetClipBoard_Par() const
    { return ClipBoard_Par; }
    
    /**  @brief set subdomain number to this cell   */
    void SetSubDomainNo(int val)
       {SubDomainNumber = val;}

    /**  @brief set subdomain number to this cell   */
    int GetSubDomainNo() const
       {return SubDomainNumber;}

    void SetAsOwnCell()
     { OwnCell=true; }

    void SetAsDependentCell()
     { DependentCell=true; }

    void SetAsHaloCell()
     { 
      OwnCell=false;
      HaloCell=true;
     }

    bool IsHaloCell() const
     {  return HaloCell; }

    void SetAsSubDomainInterfaceCell()
     { SubDomainInterfaceCell=true; }

    bool IsSubDomainInterfaceCell() const
     {  return SubDomainInterfaceCell; }

    void SetAsCrossEdgeCell()
     { CrossEdgeCell=true; }

    bool IsCrossEdgeCell() const
     {  return CrossEdgeCell; }


    void SetAsCrossVertexCell()
     { CrossVertexCell=true; }

    bool IsCrossVertexCell() const
     {  return CrossVertexCell; }

    bool IsDependentCell() const
     {  return DependentCell; }

    /// \todo TODO This method is never called - that is a problem, since
    /// information on N_NeibProcesses might be of interest in other parts of the program
    /// (I need it for Vanka Smoothers in multgrid! CB)
    void SetN_NeibProcesses(int n)
     { N_NeibProcesses = n; }

    int GetN_NeibProcesses() const
     { return N_NeibProcesses; }

    void SetNeibProcessesIds(int *Neiblist);

    int *GetNeibProcessesIds() const
     { return NeibProcessesIds; } 
     
#endif

    virtual void check() const = 0;
    
#ifdef __3D__
   /** @brief the vertices of a 3D cell have a specific order
    * (for tets the right hand rule holds; for hexas it holds using vert (0,1,2,4) )
    * check the sign of the triple product of three vectors - for the right hand rule it has to be positive*/
   bool check_orientation() const;

   /** @brief check whether the ShapeDesc actually describes the shape of the cell
    * (does the cell have a volume?)*/
   bool check_shape() const;
#endif

};

#endif

/** ************************************************************************ 
*
* @class     TGridCell
* @brief     represent geometric information of the cell
* @author    Volker Behns  09.07.97
* @History 
 ************************************************************************  */

#ifndef __GRIDCELL__
#define __GRIDCELL__

#include <BaseCell.h>

/**  @brief represent geometric information of the cell */
class TGridCell : public TBaseCell
{
  protected:
    /**  @brief field of pointer to children */
    TBaseCell **Children;
    /**  @brief pointer to parent cell */
    TBaseCell *Parent;

    /**  @brief field of all vertices */
    TVertex  **Vertices;

    /**  @brief grid level on with this cell was generated */
    int RefLevel;

  public:
    // Constructor
    TGridCell(TRefDesc *refdesc, int reflevel);

    // Destructor
    ~TGridCell();

    // Methods
    /**  @brief set the pointer to vertex with number i */
    virtual int SetVertex(int Vert_i, TVertex *Vert) override;
    /**  @brief return the pointer to vertex with number i */
    virtual TVertex *GetVertex(int Vert_i) override;
    virtual const TVertex *GetVertex(int Vert_i) const override;
    /**  @brief return field of pointers to all vertices */
    TVertex **GetVertices() const
    { return Vertices; }

    /**  @brief return number of children */
    virtual int GetN_Children() override;
    /**  @brief return number of parents */
    virtual int GetN_Parents() override;

    /**  @brief return pointer to child cell with number C\_i */
    virtual TBaseCell *GetChild(int C_i) override;
    /**  @brief return pointer to parent cell */
    virtual TBaseCell *GetParent() override;
    /**  @brief set parent */
    virtual int SetParent(TBaseCell *parent) override;
    /**  @brief return local number of child Me */
    virtual int GetChildNumber(TBaseCell *Me) override;

    /**  @brief put out postscript data to a file */
    virtual int PS(std::ofstream &dat, double scale, double StartX,
                   double StartY, bool gridcell_with_numbers) const override;
    
    /**  @brief derefine the cell */
    virtual int Derefine() override;
    /**  @brief refine or derefine the cell according to cell's clipboard */
    virtual int RefDeref() override;
    /**  @brief set marks in neighbour cells in order to maintain 1-regularity */
    virtual int Gen1RegMarks() override;
    /**  @brief generate conforming closures */
    virtual int MakeConfClosure() override;
    /**  @brief refine a cell */
    virtual int Refine(int RefLevel) override;

    #ifdef __MORTAR__
      /**  @brief refine a mortar cell */
      virtual int RefineMortar(int RefLevel);
    #endif

    /**  @brief generate a 2-regular grid */
    virtual int Gen1RegGrid() override;
    /**  @brief set refinement for the neighbour of your parent on joint LocJointNum */
    virtual int Ref1Reg(int LocJointNum, TBaseCell *&RefCell) override;
    /**  @brief check whether the surroundings of cell is 1-regular */
    virtual int Check1Reg() override;
    /**  @brief set RefDesc to no refinement */
    virtual int SetNoRefinement() override;
    /**  @brief set RefDesc to regular refinement */
    virtual int SetRegRefine() override;
    /**  @brief set RefDesc to adaptive refinement */
    virtual int Set1Refine(int i) override;      
    /**  @brief check whether a cell should be refined */
    virtual int IsToRefine() const override;
    /**  @brief check whether exist some children */
    virtual int ExistChildren() override
    { return Children == nullptr ? false : true; }

#ifdef __2D__
    /**  @brief return coordinates of mid point P\_j on edge J\_i */
    virtual int LineMidXY(int J_i, int P_j, double &X, double &Y) override;
    /**  @brief return parameters on boundary of subedge SJ\_j on edge J\_i */
    virtual int LineMidT(int J_i, int SJ_j, double &T_0, double &T_1) override;
#else
    /**  @brief return whether a point is inside a cell */
    virtual bool PointInCell(double X, double Y, double Z) const override;
#endif

    /**  @brief return whether a point is inside a cell */
    virtual bool PointInCell(double X, double Y) const override;

    /**  @brief get diameter of a cell */
    virtual double GetDiameter()  override
    { return RefDesc->GetShapeDesc()->GetDiameter(Vertices); }

    /**  @brief get shortest edge of a cell */
    virtual double GetShortestEdge() override
    { return RefDesc->GetShapeDesc()->GetShortestEdge(Vertices); }

    /**  @brief return the length of the cell defined with the reference map */
    virtual double GetLengthWithReferenceMap() override
    { return RefDesc->GetShapeDesc()->GetLengthWithReferenceMap(Vertices); }

     /**  @brief get measure of a cell */
    virtual double GetMeasure() override
    { return RefDesc->GetShapeDesc()->GetMeasure(Vertices); }

    /**  @brief get geometry level */
    virtual int GetGeoLevel() override;

    /**  @brief return subgrid ID */
    virtual int GetSubGridID() override;

    /**  @brief compute number of edges at the boundary */
    virtual int GetN_BoundaryEdges();
#ifdef __3D__
    /**  @brief compute number of faces at the boundary */
    virtual int GetN_BoundaryFaces();    
#endif 
    /**  @brief compute number of vertices at the boundary */
    virtual int GetN_BoundaryVertices();
    
    virtual void check() const override;
};

#endif

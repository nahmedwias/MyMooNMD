/** ************************************************************************ 
*
* @class TCollection 
* @date  14.10.97
* @brief store cells in an array
* @author Gunar Matthies & Sashikumaar Ganesan
* @History: MPI methods (Sashikumaar Ganesan, 08.08.14)
   
****************************************************************************/

#ifndef __COLLECTION__
#define __COLLECTION__

#include <vector>
#include <BaseCell.h>
#include <JointCollection.h>
#include <BoundEdge.h>

/** @brief store cells in an array, used by cell iterators */
class TCollection
{
  protected:
    /** @brief number of cells stored */
    int N_Cells;

    /** @brief array containing the pointers to the cells */
    TBaseCell **Cells;

    /** @brief array with all cells sorted by pointer */
    TBaseCell **SortedCells;

    /** @brief array with index of SortedCells in Cells */
    int *Index;

    
    
#ifdef  _MPI
    /** @brief Number of own cells (excluding Halo cells) */
    int N_OwnCells;
#endif
    
    /** @brief array for Globalcell number in Cells */
    int *GlobalIndex;
    
  public:
    /** @brief constructor */
    TCollection(int n_cells, TBaseCell **cells);

    /** @brief return number of cells */
    int GetN_Cells() const
    { return N_Cells; }

    /** @brief return Cell with index i in Cells-array */
    TBaseCell *GetCell(int i) const
    { return Cells[i]; }

    /** @brief destructor: delete arrays */
    ~TCollection();

    /** @brief get maximal and minimal diameter */
    int GetHminHmax(double *hmin, double *hmax);

    /** @brief return Index of cell in Cells-array */
    int GetIndex(TBaseCell *cell);

    /** @brief mark the vertices that are on the boundary */
    int MarkBoundaryVertices();

    /** @brief return Index of joints in Cells-array */
    TJointCollection  *GetJointCollection();

    /** @brief Generate   Vertex Neibs for all cells in the collection */
    void GenerateCellVertNeibs();

    /** @brief return the Index of the vertex in the sorted array */
    int GetIndex(TVertex **Array, int Length, TVertex *Element);
    
#ifdef  _MPI
    void SetN_OwnCells(int n_OwnCells)
     { N_OwnCells = n_OwnCells; }

    int GetN_OwnCells()
     { return N_OwnCells; }

    int GetN_HaloCells()
     { return (N_Cells - N_OwnCells); }

    int *GetGlobalIndex()
    {
      return GlobalIndex;
    }

    /**
     * Find the lowest-number process which contains the given point (x,y,z)
     * in an OwnCell.
     *
     * Throws an error if the point was not found anywhere.
     *
     * @param x x value
     * @param y y value
     * @param z z value
     * @return The lowest number process of those processes which contain
     * the given point in an "OwnCell".
     */
    int find_process_of_point(double x, double y, double z) const;
#endif

    void Replace_Coll(int n_cells, TBaseCell **cells)
     {
      N_Cells = n_cells;
      Cells = cells;
     }
   
   // ------------------------------------------------
   ///@brief create a list of nodes, vertices, elements
   int createElementLists();
   std::vector<double> NodesCoords;
   std::vector<int> NodesReferences;
   std::vector< std::vector<int> > ElementNodes;
   std::vector<int> ElementReferences;
   std::vector<int> BdFacesNodes;
   std::vector<int> BdFacesReferences;
   std::vector<int> DomainVertexNumbers;
   int getIndexInCollection(TBaseCell *cell);

   ///@brief write the geometry in .mesh format
   /// @param dimension If the file should be written in 2D or 3D. 3D is necessary
   /// for opening the file with gmsh. 2D must be used if the file is supposed to
   /// be used by, e.g., Brush.
   int writeMesh(const char *meshFileName, int dimension = 3);
    
   /**@brief Write a list of boundary edges
    */
   void get_edge_list_on_component(int i,std::vector<TBoundEdge*> &edges);
   ///@todo it is better to return the vector?
   // std::vector<TBoundEdge*> get_edge_list_on_component(int i);
   // ------------------------------------------------

 private:
    /** @brief provide additional arrays */
    void GenerateSortedArrays();

    /** @brief return Index of cell in SortedCells-array */
    int GetSortedIndex(TBaseCell *cell);

};

#endif

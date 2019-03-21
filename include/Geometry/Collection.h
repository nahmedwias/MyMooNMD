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
#include <map>
class TBaseCell;
class TBoundEdge;
class TBoundFace;

/** @brief store cells in an array, used by cell iterators */
class TCollection
{
  protected:
    /** @brief number of cells stored */
    int N_Cells;

    /** @brief array containing the pointers to the cells */
    TBaseCell **Cells;

    /** @brief map each cell to its index within this collection.
     * This enables the method get_cell_index. It serves as a cache is therefore
     * mutable.
     */
    mutable std::map<const TBaseCell*, int> cell_to_index_map;
    
#ifdef  _MPI
    /** @brief Number of own cells (excluding Halo cells) */
    int N_OwnCells;
#endif
    
    struct ElementLists
    {
      //----- used in DataWriter ----
      std::vector<double> NodesCoords;
      std::vector<unsigned int> NodesReferences;
      std::vector< std::vector<unsigned int> > ElementNodes;
      std::vector<unsigned int> ElementReferences;
      std::vector<unsigned int> BdFacesNodes;
      std::vector<unsigned int> BdFacesReferences;
      std::vector<unsigned int> DomainVertexNumbers;
      unsigned int NLocVertices;
      bool empty() const { return NodesCoords.empty(); }
    };
    /// @brief store redundant data for faster access, mostly needed for 
    /// DataWriter. It serves as a cache is therefore mutable.
    mutable ElementLists element_lists;
    /// @brief create a list of nodes, vertices, elements.
    /// This method is 'const' because it only alters the mutable member 
    /// 'element_lists'.
    void createElementLists() const;
   
    
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
    int GetHminHmax(double *hmin, double *hmax) const;

    /** @brief return index of cell in Cells-array */
    int get_cell_index(const TBaseCell *cell) const;
    
#ifdef  _MPI
    void SetN_OwnCells(int n_OwnCells)
     { N_OwnCells = n_OwnCells; }

    int GetN_OwnCells() const
     { return N_OwnCells; }

    int GetN_HaloCells() const
     { return (N_Cells - N_OwnCells); }

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

   ///@brief write the geometry in .mesh format
   int writeMesh(const char *meshFileName);
    
   /**@brief Write a list of boundary edges
    */
   void get_edge_list_on_component(int i,std::vector<TBoundEdge*> &edges) const;
   ///@todo it is better to return the vector?
   // std::vector<TBoundEdge*> get_edge_list_on_component(int i);
   // ------------------------------------------------
   
   void get_boundary_edge_list(std::vector<TBoundEdge*> &edges) const;
   
   //New LB 11.10.18
#ifdef __3D__
   void get_face_list_on_component(int boundary_component_id,
                                   std::vector<TBoundFace*> &faces) const;
#endif

   //####################################################
   //--------------- used in DataWriter -----------------
   //####################################################
   
   unsigned int GetN_Vertices() const;
   
   unsigned int GetN_BdFaces() const;
   
   double GetCoord(unsigned int vert) const;
   
   unsigned int GetBdFacesNode(unsigned int face) const;
   
   unsigned int GetGlobalVerNo(unsigned int cell, unsigned int locvert) const;
   
   unsigned int GetNLocVertices() const;
};

#endif

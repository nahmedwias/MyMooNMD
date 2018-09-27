// =======================================================================
// @(#)Vertex.h        1.1 10/30/98
// 
// Class:       TVertex
// Purpose:     a vertex in a grid
//
// Author:      Volker Behns  09.07.97
//              Sashikumaar Ganesan 05.11.09 (added parallel methods)
//              Sashikumaar Ganesan 08.09.2010 (added 3D parallel methods)
// =======================================================================

#ifndef __VERTEX__
#define __VERTEX__

#include <MooNMD_Io.h>
class TBaseCell;

/** a vertex in a grid */
class TVertex
{
  protected:
    /** first coordinate */
    double X;
    /** second coordinate */
    double Y;
    /** third coordinate (3D) */
    double Z = 0;

    /** an integer for storing clipboard information*/
    int ClipBoard;

#ifdef _MPI

    /** Number of 3D cells containing this cells **/
    /** Note !this info only set for dependent cells !!!!!!! */
    int N_Cells;

    /** cells */
    /** Note ! this info only set for dependent cells !!!!!!!!!!*/
    TBaseCell **Cells;

    /** marking this vertex as subdomain vertex */
    bool SubDomainVert;

    /** marking this vertex as cross vertex */
    bool CrossVert;

    /** an integer which stores the number of ranks (SubDomains) contain this vertex */
    int N_SubDomains;

    /** an integer which stores the rank of SubDomains, which contain this vertex */
    int *SubDomain_Ranks;

    /** list of neib cell Globalnumbers, which incident this vertex */
    int *SubDomainGlobalCellNo;

    /** list of neib cell local vert no */
    int *SubDomainLocVertNo;

    /** an integer which stores the number of Cross neib cells, which incident this vertex */
    int N_CrossNeibCells;    
#endif

    
    /** marking this vertex as Bound vertex */
    bool BoundVert;   
    
  public:
    // Constructors

#ifdef __3D__
      /** 3D vertex */
      TVertex(double initX, double initY, double initZ);
#else
      /** 2D vertex */ 
      TVertex(double initX, double initY);
#endif

    // Destructor
    ~TVertex();
    TVertex& operator=(const TVertex&) = delete;
    TVertex(const TVertex&) = delete;

    // Methods

    // set coordinates
#ifdef __3D__
      /** set the coordinates in 3D */
      void SetCoords(double initX, double initY, double initZ);
#else
      /** set the coordinate in 2D */
      void SetCoords(double initX, double initY);
#endif

    /** return the x coordinate */
    double GetX() const
    { return X; }
    /** return the y coordinate */
    double GetY() const
    { return Y; }

    /** return the z coordinate (3D) */
    double GetZ() const
    { return Z; }
    /** return all three coordinates */
    void GetCoords(double& x, double& y, double& z) const
    {
      x = X;
      y = Y;
      z = Z; // 0 in 2D
    }
#ifdef __2D__
    /** return all two coordinates */
    void GetCoords(double& x, double& y) const
    {
      x = X;
      y = Y;
    }
#endif

    /** write some information of the vertex in stream s */
    friend std::ostream& operator << (std::ostream& s, const TVertex *v);
    /**
     * This operator introduces an alphanumeric order on the vertices. It will
     * compare first the x component, then the y and finally the z component, if the
     * other ones are equal.
     * The vertices are regarded as equal with a tolerance of 1e-8.
     */
    friend bool operator < (const TVertex& V, const TVertex& W);

    /** set value in ClipBoard */
    void SetClipBoard(int value)
    { ClipBoard=value; }
    /** get value from ClipBoard */
    int GetClipBoard() const
    { return ClipBoard; }

     void SetAsBoundVert()
      { BoundVert = true; }
      
     bool IsBoundVert() const
     { return BoundVert; }

     
#ifdef _MPI

    /** Note ! this info only set for dependent cells !!!!!! */
    void SetVertexCells(int n_Cells, TBaseCell **cells);

    void SetSubDomainInfo(int n_SubDomains, int *subDomain_Ranks, int *subDomainGlobalCellNo, 
                          int *subDomainLocVertNo);

    void AddCrossNeib(int Neib_ID);

    void SetAsSubDomainVert()
     { SubDomainVert = true; }

     bool IsSubDomainVert()
     { return SubDomainVert; }

     void SetAsCrossVert()
     {  CrossVert = true; }
     
     bool IsCrossVert()
     { return CrossVert; }

     void GetCrossNeibs(int &n_VertCrossNeibs, int *&vertCrossNeibs) 
      {
       n_VertCrossNeibs = N_CrossNeibCells;
       vertCrossNeibs = SubDomain_Ranks;
      }

     void GetCrossNeibsInfo(int &N_NeibCells, int *&NeibCellRank, 
                            int *&GlobalNo, int *&LocVertNo)
      {
       N_NeibCells = N_CrossNeibCells;
       NeibCellRank = SubDomain_Ranks;
       GlobalNo = SubDomainGlobalCellNo;
       LocVertNo = SubDomainLocVertNo;
      }

     void GetNeibs(int &n_Neibs, TBaseCell **&neighbs)
      {
       n_Neibs = N_Cells;
       neighbs = Cells;
      }
      
      int GetNNeibs()
      {  return N_Cells; }     
#endif
};

#endif

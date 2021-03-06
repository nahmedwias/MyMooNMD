// =======================================================================
// @(#)Joint.h        1.7 11/15/99
// 
// Class:       TJoint
// Purpose:     superclass for edges and faces
//
// Author:      Volker Behns  23.07.97
//
// History:     Add method for getting a certain neighbour
//              (Gunar Matthies 17.10.97)
//
// =======================================================================

#ifndef __JOINT__
#define __JOINT__

#include <Constants.h>

#ifndef __3D__
  #define MAXN_nVpoJ  3
  #define MAXN_nJpoJ  2
#else
  #define MAXN_nVpoJ  9
  #define MAXN_nJpoJ  4
#endif

class TJoint;

struct StoreGeom
{
  const TVertex *Vertices[MAXN_nVpoJ];
  const TJoint *Joints[MAXN_nJpoJ];
  bool Filled;
};

/** supercall for edges and faces */
class TJoint
{
  protected:
    JointType ID;

    /** first neighbour */
    TBaseCell *Neighb0;
    /** second neighbour */
    TBaseCell *Neighb1;

    /** value in ClipBoard (dG method)*/
    int ClipBoard;
    
    /** */
    int NeibSubDomainLocalJointNo;
    
    /** The index of this joint in the two neighbors */
    ///@todo set the size of the array as a function of the Joint type (e.g. 1 for boundaries)
    int IndexInNeighbor[2];
    
#ifdef __3D__
    int MapType;
#endif

  public:
    // Constructors
    TJoint();

    // Methods
    /** return type */
    JointType GetType() const
    { return ID; }
    
    void ChangeType(JointType New_ID)
     {ID = New_ID;}

#ifdef __3D__
    /** set mapper type automatically */
    void SetMapType();

    /** set mapper type */
    void SetMapType(int maptype)
    { MapType = maptype; }

    /** return mapper type */
    int GetMapType() const
    { return MapType; }
    
    /** Function is used to get local edge index on neighboured element */
    int GetNeighbourEdgeIndex(const TBaseCell*, int) const;
#endif

    /** check the refinement pattern on both sides for matching,
        return already existing object on the joint in Tmp */
    virtual int CheckMatchingRef(TBaseCell *Me, int J_i,
                  struct StoreGeom &Tmp) = 0;

    /** create a new instance of the same class */
    virtual TJoint *NewInst(double T_0, double T_1, TBaseCell *Me) = 0;
    virtual TJoint *NewInst() = 0;

    /** set the neighbour to Neighb */
    int SetNeighbour(TBaseCell *Neighb);
    /** return the neighbour of this joint which is not equal to Me */
    TBaseCell *GetNeighbour(const TBaseCell *Me) const;

    /** set neighbour i to Neighb */
    int SetNeighbour(int i, TBaseCell *Neighb);
    /** return neighbour with number i */
    TBaseCell *GetNeighbour(int i) const;

    /** remove a neighbour */
    void Delete(TBaseCell *Neighb);

    /** function for debug purpose only */
    TBaseCell *GetNeighb(int i) const
    { return(i ? Neighb1 : Neighb0); }

    /** return whether this is an interior joint */
    virtual bool InnerJoint() const = 0;

    #ifdef __3D__
      /** return mapper of refined vertices and faces */
      void GetMapperRef(const int *&MapVerts, const int *&MapFaces) const;

      /** return mapper of original vertices and edges */
      void GetMapperOrig(const int *&MapVerts, const int *&MapEdges) const;
    #endif
            
    /** set value in ClipBoard */
    void SetClipBoard(int value)
    { ClipBoard=value; }
    /** get value from ClipBoard */
    int GetClipBoard()
    { return ClipBoard; }
    
    void SetNeibSubDomainLocalJointNo(int value)
    { NeibSubDomainLocalJointNo = value; }
    
    int GetNeibSubDomainLocalJointNo()
    { return NeibSubDomainLocalJointNo; }  
    
    // Destructor
    virtual ~TJoint();

};

#endif

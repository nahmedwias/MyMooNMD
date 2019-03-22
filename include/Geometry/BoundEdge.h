// =======================================================================
// @(#)BoundEdge.h        1.2 08/27/99
// 
// Class:       TBoundEdge
// Purpose:     edge on a boundary component
//
// Author:      Volker Behns  02.08.97
//
// =======================================================================

#ifndef __BOUNDEDGE__
#define __BOUNDEDGE__

#include <Joint.h>
#include <BoundComp2D.h>
#include <Vertex.h>

/** edge on a boundary component */
class TBoundEdge : public TJoint
{
  protected:
    /** boundary component to which this edge belongs to */
    TBoundComp2D *BoundComp;

    /** paramter of starting point */
    double T_0;
    /** parameter of end point */
    double T_1;

  public:
    // Constructors
    /** initialize the edge with the boundary component bdcomp and the
        paramter of starting and end point t\_0, t\_1 */
    TBoundEdge(TBoundComp2D *bdcomp, double t_0, double t_1);

    // Methods
    /** check whether the refinement pattern on both side patch,
        dummy here: there is no neighbour */
    virtual int CheckMatchingRef(TBaseCell *Me, int J_i,
                  struct StoreGeom &Tmp);

    /** create a new instance of this class */
    virtual TJoint *NewInst(double T_0, double T_1, TBaseCell *Me);
    virtual TJoint *NewInst();

    /** return start parameter T0 */
    double GetStartParameter() const
    { return T_0; }

    /** return end paramter T1 */
    double GetEndParameter() const
    { return T_1; }

    /** return parameters */
    void GetParameters(double& t0, double& t1) const
    {
      t0 = T_0;
      t1 = T_1;
    }
    
    void get_vertices(double &x0, double &y0, double &x1, double &y1)
    {
        GetXYofT( T_0, x0, y0);
        GetXYofT( T_1, x1, y1);
    }

    double get_length()
    {
      double x0, x1, y0, y1;
      GetXYofT( this->T_0, x0, y0);
      GetXYofT( this->T_1, x1, y1);
      return sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0));
    }

    void get_normal(double &nx, double &ny)
    {
      double x0, x1, y0, y1;
      GetXYofT( this->T_0, x0, y0);
      GetXYofT( this->T_1, x1, y1);
      double length =  sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0));
      nx =  (y1-y0)/length;
      ny = (x0-x1)/length; //(x1-x0)/length; //
    }

    void get_tangent(double &tx, double &ty)
    {
      double x0, x1, y0, y1;
      GetXYofT( this->T_0, x0, y0);
      GetXYofT( this->T_1, x1, y1);
      double length = sqrt( (x1-x0) * (x1-x0) + (y1-y0) * (y1-y0) );
      tx = (x1-x0)/length;
      ty = (y1-y0)/length;
    }

    void set_index_in_neighbour(TBaseCell *neigh, int index);
    int get_index_in_neighbour(const TBaseCell*const neigh) const;


#ifdef __2D__
    /** update parameters according to the new vertex positions */
    void UpdateParameters(const TVertex *Begin, const TVertex *End);
#endif

    /** @brief return the coordinates {X,Y} of parameter value T */
    void GetXYofT(double T, double &X, double &Y)
    { BoundComp->GetXYofT(T, X, Y); }
    /** @brief return parameter value T of the coordinates {X,Y} */
    void GetTofXY(double X, double Y, double& T)
    { BoundComp->GetTofXY(X, Y, T); }

    /** return boundary component */
    const TBoundComp2D *GetBoundComp() const
    { return BoundComp; }

    /** return whether this is an interior joint */
    virtual bool InnerJoint() const
    { return false; }
    
    /** change the boundary component */
    void ChangeBoundComp(TBoundComp2D *New_BoundComp)
     {BoundComp = New_BoundComp;}    
};

#endif

#ifndef __REFTRANS2D__
#define __REFTRANS2D__

// forward declaration
class TBaseCell;

/** @brief reference transformations for 2D geometric objects */
class TRefTrans2D
{
  protected:
    const TBaseCell *Cell;

  public:
    /** @brief constuctor */
    TRefTrans2D() {};

    /** @brief transfer form reference element to original element */
    void GetOrigFromRef(double eta, double xi, double &x, double &y);

    /** @brief transfer form reference element to original element */
    void GetOrigFromRef(double *ref, double *orig);

    /** @brief transfer from original element to reference element */
    void GetRefFromOrig(double x, double y, double &eta, double &xi);

    /** @brief transfer from original element to reference element */
    void GetRefFromOrig(double *orig, double *ref);

    /** @brief calculate functions and derivatives from reference element
        to original element */
    void GetOrigValues(const TBaseCell *cell);

    /** @brief set original element to cell */
    virtual void SetCell(const TBaseCell *cell)
    {  Cell = cell; }

    /** @brief return outer normal vector */
    virtual void GetOuterNormal(int j, double zeta, double &n1, double &n2) = 0;

    /** @brief return tangent */
    virtual void GetTangent(int j, double zeta, double &t1, double &t2) = 0;

    /** @brief return volume of cell according to reference transformation */
    double GetVolume();
    
    /** @brief Piola map, needed for vector values basis functions such as 
     *         Raviart-Thomas (RT) or Brezzi-Douglas-Marini (BDM).
     */
    virtual void PiolaMapOrigFromRef(int N_Functs, double *refD00,
                                     double *origD00);
    
    /// @brief Piola map for derivatives
    virtual void PiolaMapOrigFromRef(int N_Functs, double *refD10, 
                                     double *refD01, double *origD10, 
                                     double *origD01);
};

#endif

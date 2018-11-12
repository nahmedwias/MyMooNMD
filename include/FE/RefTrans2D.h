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
    TRefTrans2D() = default;

    /** @brief transfer form reference element to original element */
    virtual void GetOrigFromRef(double eta, double xi, double &x, double &y)
      = 0;
    
    /** @brief transfer form reference element to original element */
    virtual void GetOrigFromRef(int N_Points, const double *eta,
                                const double *xi, double *x, double *y,
                                double *absdetjk) = 0;

    /** @brief transfer form reference element to original element */
    virtual void GetOrigFromRef(const double *ref, double *orig) = 0;

    /** @brief transfer from original element to reference element */
    virtual void GetRefFromOrig(double x, double y, double &eta, double &xi)
      = 0;

    /** @brief transfer from original element to reference element */
    virtual void GetRefFromOrig(const double *orig, double *ref) = 0;

    /** @brief set original element to cell */
    virtual void SetCell(const TBaseCell *cell)
    {  Cell = cell; }

    /** @brief return outer normal vector */
    virtual void GetOuterNormal(int j, double zeta, double &n1, double &n2)
      const = 0;

    /** @brief return tangent */
    virtual void GetTangent(int j, double zeta, double &t1, double &t2) const
      = 0;

    /** @brief return volume of cell according to reference transformation */
    virtual double GetVolume() const = 0;
    
    /** @brief Piola map, needed for vector values basis functions such as 
     *         Raviart-Thomas (RT) or Brezzi-Douglas-Marini (BDM).
     */
    virtual void PiolaMapOrigFromRef(double xi, double eta, int N_Functs,
                                     const double *refD00, double *origD00) = 0;
    
    /// @brief Piola map for derivatives
    virtual void PiolaMapOrigFromRef(double xi, double eta, int N_Functs,
                                     const double *refD00, const double *refD10,
                                     const double *refD01, double *origD10,
                                     double *origD01) = 0;
};

#endif

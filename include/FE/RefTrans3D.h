#ifndef __REFTRANS3D__
#define __REFTRANS3D__

// forward declaration
class TBaseCell;

/** @brief reference transformations for 3D geometric objects */
class TRefTrans3D
{
  protected:
    const TBaseCell *Cell;

  public:
    /** @brief constuctor */
    TRefTrans3D() = default;

    /** @brief transfer form reference element to original element */
    virtual void GetOrigFromRef(double xi, double eta, double zeta,
                        double &x, double &y, double &z) = 0;

    /** @brief transfer form reference element to original element */
    virtual void GetOrigFromRef(const double *ref, double *orig) = 0;
    
    virtual void GetOrigFromRef(int N_Points, const double *eta,
                                const double *xi, const double *zeta,
                                double *x, double *y, double *z,
                                double *absdetjk) = 0;

    /** @brief transfer from original element to reference element */
    virtual void GetRefFromOrig(double x, double y, double z,
                                double &xi, double &eta, double &zeta) = 0;

    /** @brief transfer from original element to reference element */
    virtual void GetRefFromOrig(const double *orig, double *ref) = 0;

    /** @brief set original element to cell */
    virtual void SetCell(const TBaseCell *cell)
    {  Cell = cell; }

    /** @brief return outer normal unit vector */
    virtual void GetOuterNormal(int j, double s, double t,
                                double &n1, double &n2, double &n3) const = 0;

    /** @brief return two tangent vectors */
    virtual void GetTangentVectors(int j, double p1, double p2,
                                   double &t11, double &t12, double &t13,
                                   double &t21, double &t22, double &t23)
      const = 0;
    
    /** @brief Piola map, needed for vector values basis functions such as 
     *         Raviart-Thomas (RT) or Brezzi-Douglas-Marini (BDM).
     */
    virtual void PiolaMapOrigFromRef(double xi, double eta, double zeta,
                                     int N_Functs, const double *refD00, 
                                     double *origD00) = 0;
};

#endif

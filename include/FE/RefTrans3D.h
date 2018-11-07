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
    TRefTrans3D() {};

    /** @brief transfer form reference element to original element */
    void GetOrigFromRef(double xi, double eta, double zeta,
                        double &x, double &y, double &z);

    /** @brief transfer form reference element to original element */
    void GetOrigFromRef(double *ref, double *orig);

    /** @brief transfer from original element to reference element */
    void GetRefFromOrig(double x, double y, double z,
                        double &xi, double &eta, double &zeta);

    /** @brief transfer from original element to reference element */
    void GetRefFromOrig(double *orig, double *ref);

    /** @brief calculate functions and derivatives from reference element
        to original element */
    void GetOrigValues(const TBaseCell *cell);

    /** @brief set original element to cell */
    virtual void SetCell(const TBaseCell *cell)
    {  Cell = cell; }

    /** @brief return outer normal unit vector */
    void GetOuterNormal(int j, double s, double t,
                        double &n1, double &n2, double &n3);

    /** @brief return two tangent vectors */
    void GetTangentVectors(int j, double p1, double p2,
        double &t11, double &t12, double &t13,
        double &t21, double &t22, double &t23);
    
    /** @brief Piola map, needed for vector values basis functions such as 
     *         Raviart-Thomas (RT) or Brezzi-Douglas-Marini (BDM).
     */
    virtual void PiolaMapOrigFromRef(int N_Functs, double *refD00, 
                                     double *origD00);
};

#endif

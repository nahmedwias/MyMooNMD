#ifndef __REFTRANS1D__
#define __REFTRANS1D__

// forward declaration
class TBaseCell;

/** @brief reference transformations for 1D geometric objects */
class TRefTrans1D
{
  protected:
    TBaseCell *Cell;

  public:
    /** @brief constuctor */
    TRefTrans1D() {};

    /** @brief transfer form reference element to original element */
    void GetOrigFromRef(double xi, double &x);


    /** @brief transfer form reference element to original element */
    void GetOrigFromRef(double *ref, double *orig);

    /** @brief transfer from original element to reference element */
    void GetRefFromOrig(double x, double &eta);

    /** @brief transfer from original element to reference element */
    void GetRefFromOrig(double *orig, double *ref);

    /** @brief calculate functions and derivatives from reference element
        to original element */
    void GetOrigValues(TBaseCell *cell);

    /** @brief set original element to cell */
    virtual void SetCell(TBaseCell *cell)
    {  Cell = cell; }

    /** @brief return volume of cell according to reference transformation */
    double GetVolume();
};

#endif


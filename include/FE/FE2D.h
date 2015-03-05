// =======================================================================
// @(#)FE2D.h        1.2 05/04/99
//
// Class:       TFE2D
// Purpose:     store all information for one finite element class
//
// Author:      Gunar Matthies  09.07.98
//
// =======================================================================

#ifndef __FE2D__
#define __FE2D__

#include <AllClasses.h>
#include <Constants.h>
#include <Enumerations.h>
#include <BaseFunct2D.h>
#include <FEDesc2D.h>

/** store all information for one finite element class */
class TFE2D
{
  protected:
   /** ID for set of basis functions */
   BaseFunct2D BaseFunct_ID;

   /** set of basis function */
   TBaseFunct2D *BaseFunct;

   /** ID for set of nodal functional */
   NodalFunctional2D NodalFunctional_ID;

   /** set of nodal functional */
   TNodalFunctional2D *NodalFunctional;

   /** ID for reference transformation */
   RefTrans2D RefTransID;

   /** ID for element description */
   FEDesc2D FEDesc_ID;

   /** element description */
   TFEDesc2D *FEDesc;

   /** number of needed integer entries (numbers + infos) */
   int Size;

   /** number of degrees of freedom */
   int N_DOF;

   /** number of info blocks */
   int N_Info;

  public:
    /** constructor */
    TFE2D();

    /** constructor with data */
    TFE2D(BaseFunct2D basefunct_id, NodalFunctional2D nodalfunctional_id,
        RefTrans2D reftransid, FEDesc2D fedesc_id,
        int info);

    /** return BaseFunct2D_ID */
    BaseFunct2D GetBaseFunct2D_ID()
      { return BaseFunct_ID; };

    /** return BaseFunct2D */
    TBaseFunct2D *GetBaseFunct2D()
      { return BaseFunct; };

    /** return BaseFunct2D_ID and BaseFunct2D */
    void GetBaseFunct2D(BaseFunct2D &ID,
                        TBaseFunct2D* &Obj)
      { ID = BaseFunct_ID; Obj = BaseFunct; };

    /** return NodalFunctional2D_ID */
    NodalFunctional2D GetNodalFunctional2D_ID()
      { return NodalFunctional_ID; };

    /** return NodalFunctional2D */
    TNodalFunctional2D *GetNodalFunctional2D()
      { return NodalFunctional; };

    /** return NodalFunctional2D_ID and NodalFunctional2D */
    void GetNodalFunctional2D(NodalFunctional2D &ID,
                              TNodalFunctional2D* &Obj)
      { ID = NodalFunctional_ID; Obj = NodalFunctional; };

    /** return RefTransID */
    RefTrans2D GetRefTransID()
      { return RefTransID; };

    /** return FEDesc2D_ID */
    FEDesc2D GetFEDesc2D_ID()
      { return FEDesc_ID; };

    /** return FEDesc2D */
    TFEDesc2D *GetFEDesc2D()
      { return FEDesc; };

    /** return FEDesc2D_ID and FEDesc2D */
    void GetFEDesc2D(FEDesc2D &ID, TFEDesc2D* &Obj)
      { ID = FEDesc_ID; Obj = FEDesc; };

    /** return size */
    int GetSize()
      { return Size; };

    /** return number of degrees of freedom */
    int GetN_DOF()
      { return N_DOF; };

    /** return number of info blocks */
    int GetN_Info()
      { return N_Info; };

    /** check N[i](b[j]) = delta[ij] */
    void CheckNFandBF();
};

#endif

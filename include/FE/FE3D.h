// =======================================================================
// %W% %G%
//
// Class:       TFE3D
// Purpose:     store all information for one finite element class
//
// Author:      Gunar Matthies  19.11.99
//
// =======================================================================

#ifndef __FE3D__
#define __FE3D__

class TBaseFunct3D;
class TNodalFunctional3D;
class TFEDesc3D;
#include <Constants.h>
#include <Enumerations.h>
#include <BaseFunct3D.h>
#include <FEDesc3D.h>
#include <NodalFunctional3D.h>

/** store all information for one finite element class */
class TFE3D
{
  protected:
   /** ID for set of basis functions */
   BaseFunct3D BaseFunct_ID;

   /** set of basis function */
   TBaseFunct3D *BaseFunct;

   /** ID for set of nodal functional */
   NodalFunctional3D NodalFunctional_ID;

   /** set of nodal functional */
   TNodalFunctional3D *NodalFunctional;

   /** ID for reference transformation */
   RefTrans3D RefTransID;

   /** ID for element description */
   FEDesc3D FEDesc_ID;

   /** element description */
   TFEDesc3D *FEDesc;

   /** number of needed integer entries (numbers + infos) */
   int Size;

   /** number of degrees of freedom */
   int N_DOF;

   /** number of info blocks */
   int N_Info;

  public:
    /** constructor */
    TFE3D();

    /** constructor with data */
    TFE3D(BaseFunct3D basefunct_id, NodalFunctional3D nodalfunctional_id,
        RefTrans3D reftransid, FEDesc3D fedesc_id,
        int info);

    /** return BaseFunct3D_ID */
    BaseFunct3D GetBaseFunct3D_ID() const
      { return BaseFunct_ID; };
    
    BaseFunct3D get_id() const
      { return BaseFunct_ID; };

    /** return BaseFunct3D */
    TBaseFunct3D *GetBaseFunct3D() const
      { return BaseFunct; };

    /** return BaseFunct3D_ID and BaseFunct3D */
    void GetBaseFunct3D(BaseFunct3D &ID,
                        TBaseFunct3D* &Obj) const
      { ID = BaseFunct_ID; Obj = BaseFunct; };

    /** return NodalFunctional3D_ID */
    NodalFunctional3D GetNodalFunctional3D_ID() const
      { return NodalFunctional_ID; };

    /** return NodalFunctional3D */
    TNodalFunctional3D *GetNodalFunctional3D() const
      { return NodalFunctional; };

    /** return NodalFunctional3D_ID and NodalFunctional3D */
    void GetNodalFunctional3D(NodalFunctional3D &ID,
                              TNodalFunctional3D* &Obj) const
      { ID = NodalFunctional_ID; Obj = NodalFunctional; };

    /** return RefTransID */
    RefTrans3D GetRefTransID() const
      { return RefTransID; };

    /** return FEDesc3D_ID */
    FEDesc3D GetFEDesc3D_ID() const
      { return FEDesc_ID; };

    /** return FEDesc3D */
    TFEDesc3D *GetFEDesc3D() const
      { return FEDesc; };

    /** return FEDesc3D_ID and FEDesc3D */
    void GetFEDesc3D(FEDesc3D &ID, TFEDesc3D* &Obj) const
      { ID = FEDesc_ID; Obj = FEDesc; };

    /** return size */
    int GetSize() const
      { return Size; };

    /** return number of degrees of freedom */
    int GetN_DOF() const
      { return N_DOF; };

    /** return number of info blocks */
    int GetN_Info() const
      { return N_Info; };

    /** check N[i](b[j]) = delta[ij] */
    void CheckNFandBF()  const;
};

#endif

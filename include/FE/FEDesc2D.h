// =======================================================================
// @(#)FEDesc2D.h        1.1 10/30/98
//
// Class:       TFEDesc2D
// Purpose:     store a finite element descriptor for a 2D element
//
// Author:      Gunar Matthies  23.07.98
//
// =======================================================================

#ifndef __FEDESC2D__
#define __FEDESC2D__

#include <Enumerations.h>

/** store a finite element descriptor for a 2D element */
class TFEDesc2D
{
  protected:
    /** description for the object */
    char *Description;

    /** number of degrees of freedom */
    int N_DOF;

    /** number of degrees of freedom on closure of each joint */
    int N_JointDOF;

    /** local numbers of all degrees of freedom on the joints */
    int **JointDOF;

    /** number of inner degrees of freedom */
    int N_InnerDOF;

    /** local numbers of all inner degrees of freedom */
    int *InnerDOF;

    /** number of degrees of freedom on cell boundary */
    int N_OuterDOF;

    /** local numbers of all degrees of freedom on cell boundary */
    int *OuterDOF;

  public:
    /** constructor, setting all data without dof on cell boundary */
    TFEDesc2D(char *description, int n_dof, int n_jointdof,
              int **jointdof, int n_innerdof, int *innerdof);

    /** constructor, setting all data with dof on cell boundary */
    TFEDesc2D(char *description, int n_dof, int n_jointdof,
              int **jointdof, int n_innerdof, int *innerdof,
              int n_outerdof, int *outerdof);

    /** return description */
    char *GetDescription()
      { return Description; }

    /** return number of degrees of freedom */
    int GetN_DOF()
      { return N_DOF; }

    /** return number of degrees of freedom per closure of each joint */
    int GetN_JointDOF()
      { return N_JointDOF; }

    /** return number of inner degrees of freedom */
    int GetN_InnerDOF()
      { return N_InnerDOF; }

    /** return local numbers of inner degrees of freedom */
    int *GetInnerDOF()
      { return InnerDOF; }

    /** return number of degrees of freedom on cell boundary */
    int GetN_OuterDOF()
      { return N_OuterDOF; }

    /** return local numbers of degrees of freedom on cell boundary */
    int *GetOuterDOF()
      { return OuterDOF; }

    /** return total number and local numbers of degrees of freedom
        on cell boundary */ 
    void GetOuterDOF(int &n_outerdof, int* &outerdof)
      { n_outerdof = N_OuterDOF; outerdof = OuterDOF; }

    /** return local numbers of degrees of freedom on each joint */
    int **GetJointDOF()
      { return JointDOF; }

    /** return local numbers of degrees of freedom on joint i */
    int *GetJointDOF(int i)
      { return JointDOF[i]; }

};

#endif

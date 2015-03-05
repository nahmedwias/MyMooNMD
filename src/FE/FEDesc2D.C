// =======================================================================
// @(#)FEDesc2D.C        1.1 10/30/98
//
// Class:       TFEDesc2D
// Purpose:     store a finite element descriptor for a 2D element
//
// Author:      Gunar Matthies  23.07.98
//
// =======================================================================

#include <FEDesc2D.h>
#include <string.h>

/** constructor, setting all data with dof on cell boundary */
TFEDesc2D::TFEDesc2D(char *description, int n_dof, int n_jointdof,
                     int **jointdof, int n_innerdof, int *innerdof,
                     int n_outerdof, int *outerdof)
{
  Description = strdup(description);
  N_DOF = n_dof;
  N_JointDOF = n_jointdof;
  JointDOF = jointdof;
  N_InnerDOF = n_innerdof;
  InnerDOF = innerdof;
  N_OuterDOF = n_outerdof;
  OuterDOF = outerdof;
}

/** constructor, setting all data with dof on cell boundary */
TFEDesc2D::TFEDesc2D(char *description, int n_dof, int n_jointdof,
                     int **jointdof, int n_innerdof, int *innerdof)
{
  Description = strdup(description);
  N_DOF = n_dof;
  N_JointDOF = n_jointdof;
  JointDOF = jointdof;
  N_InnerDOF = n_innerdof;
  InnerDOF = innerdof;
  N_OuterDOF = 0;
  OuterDOF = NULL;
}

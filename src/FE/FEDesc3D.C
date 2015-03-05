// =======================================================================
// %W% %G%
//
// Class:       TFEDesc3D
// Purpose:     store a finite element descriptor for a 3D element
//
// Author:      Gunar Matthies  23.07.98
//
// =======================================================================

#include <FEDesc3D.h>
#include <string.h>
#include <iostream>

/** constructor, setting all data */
TFEDesc3D::TFEDesc3D(char *description, int n_dof, int n_jointdof,
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
  
#ifdef _MPI    
  N_EdgeDOF = 0;
  EdgeDOF = NULL;
  N_VertDOF = 0;
  VertDOF = NULL;
  EdgeVertData_Filled = 0;  
#endif

}

/** constructor, setting all data */
TFEDesc3D::TFEDesc3D(char *description, int n_dof, int n_jointdof,
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
  
#ifdef _MPI    
  N_EdgeDOF = 0;
  EdgeDOF = NULL;
  N_VertDOF = 0;
  VertDOF = NULL;
  EdgeVertData_Filled = 0;
#endif
  
}


#ifdef _MPI 
TFEDesc3D::TFEDesc3D(char *description, int n_dof, int n_jointdof,
                     int **jointdof, int n_innerdof, int *innerdof,
                     int n_edgeDOF, int **edgeDOF, int n_vertDOF, int *vertDOF )
{
  Description = strdup(description);
  N_DOF = n_dof;
  N_JointDOF = n_jointdof;
  JointDOF = jointdof;
  N_InnerDOF = n_innerdof;
  InnerDOF = innerdof;
  N_OuterDOF = 0;
  OuterDOF = NULL;

  N_EdgeDOF = n_edgeDOF;
  EdgeDOF = edgeDOF;
  N_VertDOF = n_vertDOF;
  VertDOF = vertDOF; 
  EdgeVertData_Filled = 1;
}


TFEDesc3D::TFEDesc3D(char *description, int n_dof, int n_jointdof,
                     int **jointdof, int n_innerdof, int *innerdof,
                     int n_outerdof, int *outerdof,
                     int n_edgeDOF, int **edgeDOF, int n_vertDOF, int *vertDOF )
{
  Description = strdup(description);
  N_DOF = n_dof;
  N_JointDOF = n_jointdof;
  JointDOF = jointdof;
  N_InnerDOF = n_innerdof;
  InnerDOF = innerdof;
  N_OuterDOF = n_outerdof;
  OuterDOF = outerdof;

  N_EdgeDOF = n_edgeDOF;
  EdgeDOF = edgeDOF;
  N_VertDOF = n_vertDOF;
  VertDOF = vertDOF; 
  EdgeVertData_Filled = 1;  
}

#endif  